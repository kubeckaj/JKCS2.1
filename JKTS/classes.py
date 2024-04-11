from numpy import array, dot, cos, sin, cross, zeros_like, arccos, degrees, exp, pi, radians, all
from numpy.linalg import norm
import re
import shutil
import os
import pickle
import random
from rdkit import Chem
from rdkit.Chem import AllChem
##################################################WORKFLOWS#################################################################
# Change as you feel like
OH_H2O_workflow = ['optimization', 'DLPNO']
reactant_product_workflow = ['crest_sampling', 'optimization', 'DLPNO']
TS_workflow = ['opt_constrain', 'TS_opt', 'crest_sampling', 'opt_constrain_conf', 'TS_opt_conf', 'DLPNO']
# TS_workflow = ['crest_sampling', 'opt_constrain_conf', 'TS_opt_conf', 'DLPNO', 'Done']
############################################################################################################################

class Vector:
    @staticmethod
    def calculate_vector(coord1, coord2):
        return array(coord1) - array(coord2)

    @staticmethod
    def vector_length(vector):
        return norm(vector)

    @staticmethod
    def atom_distance(atom1, atom2):
        return norm(array(atom2) - array(atom1))

    @staticmethod
    def normalize_vector(vector):
        norm_ = norm(vector)
        if norm_ < 1e-8:
            return zeros_like(vector)
        return vector / norm_

    @staticmethod
    def rotate_vector(vector, axis, angle):
        cos_theta = cos(angle)
        sin_theta = sin(angle)
        cross_product = cross(axis, vector)
        return (vector * cos_theta +
                cross_product * sin_theta +
                axis * dot(axis, vector) * (1 - cos_theta))

    @staticmethod
    def calculate_angle(coord1, coord2, coord3):
        vector1 = Vector.calculate_vector(coord2, coord1)
        vector2 = Vector.calculate_vector(coord2, coord3)

        dot_product = dot(vector1, vector2)
        magnitude1 = Vector.vector_length(vector1)
        magnitude2 = Vector.vector_length(vector2)

        angle_rad = arccos(dot_product / (magnitude1 * magnitude2))
        angle_deg = degrees(angle_rad)

        return angle_deg

    @staticmethod
    def get_upwards_perpendicular_axis(self, direction):
        # direction is the vector lying between 2 carbons in C=C
        pass


class Molecule(Vector):
    def __init__(self, file_path=None, log_file_path="", name="", directory="", atoms=None, coordinates=None, electronic_energy=None, reactant=False, product=False, program=None, indexes=None, smiles=None, method=''):
        self.smiles = smiles
        self.name = name
        self.method = method
        self.directory = directory
        self.file_path = file_path
        self.log_file_path = log_file_path
        self.reactant = True if 'reactant' in self.name.split("_") or 'OH' in self.name else reactant
        self.product = True if 'product' in self.name.split("_") or 'H2O' in self.name else product
        self.job_id = ""
        self.atoms = atoms if atoms is not None else []
        self.coordinates = coordinates if coordinates is not None else []
        self.constrained_indexes = {'C': indexes[0], 'H': indexes[1], 'O': indexes[2], 'OH': indexes[2]+1} if indexes else {}
        self.converged = False
        self.mult = 1 if self.reactant else 2 
        self.charge = 0
        self.vibrational_frequencies = []
        self.zero_point = None
        self.electronic_energy = electronic_energy
        self.free_energy = None
        self.Q = None
        self.dipole_moment = None
        self._program = None
        self.input = None
        self.output = None
        self.status = None
        self.error_termination_count = 0
        self.reactant_pair = None
        self.product_pair = None
        self.program = program if program is not None else 'g16'
        self.init_molecule(indexes)
                        

    def init_molecule(self, indexes):
        # If XYZ file is given
        if self.file_path and self.file_path.split(".")[-1] == 'xyz':
            self.read_xyz_file()
            self.workflow = self.determine_workflow()
            self.set_current_step(self.workflow[0])

        if self.smiles:
            self.name, self.atoms, self.coordinates = self.smiles_to_atoms_coordinates(self.smiles)
            self.workflow = self.determine_workflow()
            self.set_current_step(self.workflow[0])

        elif self.file_path and self.file_path.split(".")[-1] in ['log', 'out', 'com', 'inp', 'pkl']:
            self.name = f"{self.file_path.split('/')[-1].split('.')[0]}"
            self.directory = os.path.dirname(self.file_path)
            if self.file_path.split(".")[-1] in ['log', 'out']: # If a log file is given
                self.log_file_path = self.file_path
                self.program = self.log2program()
                self.atoms, self.coordinates = self.log2xyz(atoms=True)
                self.update_energy()
            elif self.file_path.split(".")[-1] in ['com', 'inp']:
                self.program = 'g16' if self.file_path.split(".")[-1] == 'com' else 'orca'
                self.atoms, self.coordinates = self.inp2xyz()

            if self.reactant or 'reactant' in self.name or 'OH' in self.name:
                self.reactant = True
                if 'OH' in self.name:
                    self.mult = 2
                else:
                    self.mult = 1
                
            elif self.product or 'product' in self.name or 'H2O' in self.name:
                self.product = True
                if 'H2O' in self.name:
                    self.mult = 1
                else:
                    self.mult = 2
                        
            else:
                self.find_active_site(indexes)

            self.workflow = self.determine_workflow()
            self.set_current_step()
        if not self.method:
            self.method = self.log2method()


    def set_current_step(self, step=None):
        if step:
            # If the step is provided and it's in the workflow, update the current and next steps accordingly
            if step in self.workflow:
                self.current_step_index = self.workflow.index(step)
                self.current_step = step
            else:
                print(f"Step '{step}' not found in the workflow.")
                return  # Optionally handle the error case differently
        else:
            self.current_step = self.determine_current_step()
            if self.current_step in self.workflow:
                self.current_step_index =  self.workflow.index(self.current_step)
            # If no step is provided, reset to the start of the workflow
            else:
                self.current_step_index = 0
                self.current_step = self.workflow[0] if self.workflow else None

        # Update the next step based on the current step index
        self.next_step_index = self.current_step_index + 1 if self.current_step_index + 1 < len(self.workflow) else None
        self.next_step = self.workflow[self.next_step_index] if self.next_step_index is not None and self.next_step_index < len(self.workflow) else None


    def determine_workflow(self):
        if self.reactant or self.product or 'OH' in self.name or 'H2O' in self.name:
            if 'OH' in self.name or 'H2O' in self.name:
                return OH_H2O_workflow
            else:
                return reactant_product_workflow
        else:
            return TS_workflow


#################################################__init__#############################################################

    @classmethod
    def create_OH(cls, directory=os.getcwd()):
        # Creates an OH radical molecule
        molecule = cls(name='OH', directory=directory, atoms=['O', 'H'], coordinates=[[0.0, 0.0, 0.0], [0.97, 0.0, 0.0]], reactant=True)
        molecule.mult = 2
        molecule.reactant = True
        molecule.workflow = molecule.determine_workflow()
        molecule.converged = False
        molecule.set_current_step()
        molecule.method = molecule.log2method()
        return molecule

    @classmethod
    def create_H2O(cls, directory=os.getcwd()):
        # Creates a water (H2O) molecule
        molecule = cls(name='H2O', directory=directory, atoms=['O', 'H', 'H'], coordinates=[[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0], [-0.239664, 0.926711, 0.0]], product=True)
        molecule.mult = 1
        molecule.product = True
        molecule.workflow = molecule.determine_workflow()
        molecule.converged = False
        molecule.set_current_step()
        molecule.method = molecule.log2method()
        return molecule


    def read_xyz_file(self):
        try:
            with open(self.file_path, 'r') as file:
                lines = file.readlines()[2:]  
                for line in lines:
                    parts = line.split()
                    if len(parts) == 4:
                        self.atoms.append(parts[0])
                        self.coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except FileNotFoundError:
            print(f"File not found: {self.file_path}")

    def save_to_pickle(self, file_path):
        with open(file_path, 'wb') as file:
            pickle.dump(self, file)

    def move_files(self, ignore_file=''):
        dest_directory = os.path.join(self.directory, "log_files")
        if not os.path.exists(dest_directory):
            os.makedirs(dest_directory, exist_ok=True)

        for filename in os.listdir(self.directory):
            if filename.endswith(self.output):
                file_path = os.path.join(self.directory, filename)
                dest_file_path = os.path.join(dest_directory, filename)
                if os.path.exists(dest_file_path):
                    os.remove(dest_file_path)
                shutil.move(file_path, dest_file_path)

    def move_converged(self):
        if not os.path.exists(os.path.join(self.directory, "log_files")):
            os.makedirs(os.path.join(self.directory, "log_files"), exist_ok=True)
        destination = os.path.join(self.directory, "log_files", f"{self.name}{self.output}")
        if os.path.exists(destination):
            os.remove(destination)
        shutil.move(os.path.join(self.directory, f"{self.name}{self.output}"), destination)
        self.log_file_path = destination

    def move_inputfile(self):
        input_files_dir = os.path.join(self.directory, "input_files")
        if not os.path.exists(input_files_dir):
            os.makedirs(input_files_dir, exist_ok=True)
        source = os.path.join(self.directory, f"{self.name}{self.input}")
        destination = os.path.join(input_files_dir, f"{self.name}{self.input}")
        if os.path.exists(destination):
            os.remove(destination)
        if os.path.exists(source):
            shutil.move(source, destination)

    def move_failed(self):
        if not os.path.exists(os.path.join(self.directory, "failed_logs")):
            os.makedirs(os.path.join(self.directory, "failed_logs"), exist_ok=True)
        destination = os.path.join(self.directory, "failed_logs", f"{self.name}_{self.error_termination_count}{self.output}")
        if os.path.exists(destination):
            os.remove(destination)
        shutil.move(os.path.join(self.directory, f"{self.name}{self.output}"), destination)
            
    def write_xyz_file(self, output_file_path):
        with open(output_file_path, 'w') as file:
            file.write(str(len(self.atoms)) + '\n\n')
            for atom, coord in zip(self.atoms, self.coordinates):
                coord_str = ' '.join(['{:.6f}'.format(c) for c in coord])  
                file.write(f'{atom} {coord_str}\n')


    def update_current_and_next_step(self):
        # Update the current and next step based on the current step index
        self.current_step = self.workflow[self.current_step_index]
        self.next_step_index = self.current_step_index + 1
        if self.next_step_index < len(self.workflow):
            self.next_step = self.workflow[self.next_step_index]
        else:
            self.next_step = None  # Indicates the end of the workflow

    def update_step(self):
        if self.next_step_index < len(self.workflow):
            self.current_step_index += 1
            self.update_current_and_next_step()
        else:
            print("Reached the end of the workflow.")


    def find_active_site(self, indexes=None):
        if self.constrained_indexes:
            return
        if indexes:
            C_index, H_index, O_index = indexes
            H_OH_index = O_index + 1
            self.constrained_indexes = {'C': C_index, 'H': H_index, 'O': O_index, 'OH': H_OH_index}
            return
        else:
            constrain_file_path = os.path.join(self.directory, ".constrain")
            if os.path.exists(constrain_file_path):
                try:
                    with open(constrain_file_path, "r") as file:
                        constraints = {}
                        for line in file:
                            atom, index = line.split(":")
                            constraints[atom.strip()] = int(index.strip())
                        C_index = constraints.get('C', None)
                        H_index = constraints.get('H', None)
                        O_index = constraints.get('O', None)
                        H_OH_index = O_index + 1
                        self.constrained_indexes = {'C': C_index, 'H': H_index, 'O': O_index, 'OH': H_OH_index}
                        return
                except Exception as e:
                    print(f"Error reading .constrain file: {e}")
                    # If reading .constrain fails, proceed to loop through the molecule
                    indexes = None

            if not indexes:  # Fallback to original logic if indexes not set
                # Directly access the indexes of the OH radical
                O_index = len(self.atoms) - 2  # Second-last atom is the oxygen of the OH radical
                H_OH_index = len(self.atoms) - 1  # Last atom is the hydrogen of the OH radical

                for C, atom in enumerate(self.atoms):
                    if atom == 'C':
                        for H, neighbor in enumerate(self.atoms):
                            if neighbor == 'H' and self.is_nearby(C, H):
                                # Since we already know the indexes of O and H in the OH radical,
                                # we directly use them instead of iterating through all atoms again.
                                if self.calculate_angle(self.coordinates[C], self.coordinates[H], self.coordinates[O_index]) > 120:
                                    # Check if the H from the molecule is close enough to the O from the OH radical
                                    # indicating a potential active site for hydrogen abstraction.
                                    if self.is_nearby(H, O_index, threshold_distance=2):
                                        if self.is_nearby(C, O_index, threshold_distance=2.8):
                                            self.constrained_indexes = {'C': C+1, 'H': H+1, 'O': O_index+1, 'OH': H_OH_index+1}
                                        return

        if not self.constrained_indexes:
            raise ValueError('Indexes for atoms involved in transition state site could not be determined\n Consider using "-CHO c_index h_index o_index" in the input argument. Open visualizer to manually get indexes')

    def is_nearby(self, atom_index1, atoms_index2, threshold_distance=1.7):
        distance = norm(array(self.coordinates[atom_index1]) - array(self.coordinates[atoms_index2]))
        return distance < threshold_distance

    @property
    def zero_point_corrected(self):
        if self.zero_point is not None and self.electronic_energy is not None:
            return float(self.zero_point) + float(self.electronic_energy)
        return None
    
    @staticmethod
    def load_from_pickle(file_path):
        with open(file_path, 'rb') as file:
            return pickle.load(file)

    @staticmethod
    def molecules_to_pickle(molecules, file_path):
        '''molecules = [molecule1, molecule2, ...]
        usage: Molecule.molecules_to_pickle(molecules, 'path_to_temp_pkl_file')'''
        with open(file_path, 'wb') as file:
            pickle.dump(molecules, file)

    @staticmethod
    def load_molecules_from_pickle(file_path):
        '''usage: loaded_molecules = Molecule.load_molecules_from_pickle('path_to_temp_pkl_file')'''
        with open(file_path, 'rb') as file:
            return pickle.load(file)

    @property
    def program(self):
        return self._program

    @program.setter
    def program(self, value):
        self._program = value.lower()
        if self._program == 'orca':
            self.input = '.inp'
            self.output = '.out' #.out
        elif self._program == 'g16':
            self.input = '.com'
            self.output = '.log'
        elif self._program == 'crest':
            self.input = '.xyz'
            self.output = '.output'
        else:
            raise ValueError(f"Unsupported program: {value}")


    def update_energy(self, logger=None, log_file_path=None, DLPNO=False, program=None):
        file_path = log_file_path if log_file_path else self.log_file_path
        program = program if program else self.program

        with open(file_path, 'r') as f:
            log_content = f.read()

            if program.lower() == 'g16':
                zero_point_correction = re.findall(r'Zero-point correction=\s+([-.\d]+)', log_content)
                free_energy = re.findall(r'Sum of electronic and thermal Free Energies=\s+([-.\d]+)', log_content)
                dipole_moment = re.findall(r"Tot=\s+([-\d.]+)", log_content)
                if zero_point_correction:
                    self.zero_point = float(zero_point_correction[-1])
                    self.free_energy = float(free_energy[-1])
                    self.dipole_moment = float(dipole_moment[-1])
                electronic_energy = re.findall(r'(SCF Done:  E\(\S+\) =)\s+([-.\d]+)', log_content)
                if electronic_energy:
                    self.electronic_energy = float(electronic_energy[-1][-1])
                    freq_matches = re.findall(r"Frequencies --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)?\s+(-?\d+\.\d+)?", log_content)
                    if freq_matches:
                        self.vibrational_frequencies = [float(freq) for match in freq_matches for freq in match if freq]
                        rot_temp = re.findall(r"Rotational temperatures? \(Kelvin\)\s+(-?\d+\.\d+)(?:\s+(-?\d+\.\d+))?(?:\s+(-?\d+\.\d+))?", log_content)
                        self.rot_temps = [float(rot) for rot in rot_temp[-1] if rot != '']
                        symmetry_num = re.search(r"Rotational symmetry number\s*(\d+)", log_content)
                        if symmetry_num:
                            self.symmetry_num = int(symmetry_num.group(1))
                        else: 
                            if logger:
                                logger.log(f"No symmetry number found in {self.name}. assuming 1")
                            self.symmetry_num = 1
                        mol_mass = re.search(r"Molecular mass:\s+(-?\d+\.\d+)", log_content)
                        self.mol_mass = float(mol_mass.group(1))
                        mult = re.search(r"Multiplicity =\s*(\d+)", log_content)
                        if mult:    
                            self.mult = int(mult.group(1))
                        else: self.mult = 2

                        self.partition_function()
                    elif 'TS' in self.name:
                        if logger:
                            logger.log(f"No frequencies found in {self.name}")

                # partition_function = re.search() # implement reading Q from log file

            elif program.lower() == 'orca' or self.current_step == 'DLPNO' or DLPNO:
                zero_point_correction = re.findall(r"Zero point energy\s+...\s+([-+]?\d*\.\d+|\d+)", log_content)
                free_energy = re.findall(r"Final Gibbs free energy\s+...\s+([-+]?\d*\.\d+|\d+)", log_content)
                electronic_energy = re.findall(r'(FINAL SINGLE POINT ENERGY)\s+([-.\d]+)', log_content)
                dipole_moment = re.findall(r"Magnitude \(Debye\)\s*:\s*([\d.]+)", log_content)
                # find ORCA zero point zorrected
                if electronic_energy:
                    self.electronic_energy = float(electronic_energy[-1][-1])
                    if zero_point_correction and free_energy:
                        self.zero_point = float(zero_point_correction[-1])
                        self.free_energy = float(free_energy[-1])
                        self.dipole_moment = float(dipole_moment[-1])
                    freq_matches = re.findall(r'([-+]?\d*\.\d+)\s*cm\*\*-1', log_content)
                    if freq_matches:
                        n = 3*len(self.atoms)-6 # Utilizing the fact that non-linear molecules has 3N-6 degrees of freedom
                        self.vibrational_frequencies = [float(freq) for freq in freq_matches][-n:]
                        rot_temp = re.findall(r"Rotational constants in cm-1: \s*[-+]?(\d*\.\d*)  \s*[-+]?(\d*\.\d*) \s*[-+]?(\d*\.\d*)", log_content)
                        self.rot_temps = [float(rot) for rot in rot_temp[-1]]
                        symmetry_num = re.search(r'Symmetry Number:\s*(\d*)', log_content)
                        if symmetry_num:
                            self.symmetry_num = int(symmetry_num.group(1))
                        else: 
                            self.symmetry_num = 1

                        mol_mass = re.search(r'Total Mass\s*...\s*\d*\.\d+', log_content)
                        if mol_mass:
                            self.mol_mass = float(mol_mass.group().split()[-1])
                        multiplicity = re.search(r'Mult\s* ....\s*(\d*)', log_content)
                        if multiplicity:
                            self.mult = int(multiplicity.group().split()[-1])
                        else:
                            self.mult = 2

                        self.partition_function()
                
                # partition_function = re.search() # implement reading Q from log file


    def partition_function(self,  T=298.15):
        h = 6.62607015e-34  # Planck constant in J.s
        k_b = 1.380649e-23  # Boltzmann constant in J/K
        c = 299792458       # Speed of light in m/s
        R = 8.314462618153  # m³ Pa K⁻¹ mol⁻¹
        P = 100000          # Pa

        # Vibrational partition function
        qvib = 1 
        for freq in self.vibrational_frequencies:
            if freq > 0:
                f = freq * 100 * c
                qvib *= 1 / (1 - exp(-(h * f) / (k_b * T))) 

        # Rotational partition function
        if self.program.lower() == 'g16':
            rot_temps = self.rot_temps
            if len(rot_temps) == 1:
                qrot = ((1/self.symmetry_num) * (T/rot_temps[0]))
            else:
                qrot = (pi**0.5/self.symmetry_num) * ((T**3 / (rot_temps[0] * rot_temps[1] * rot_temps[2]))**0.5)
        elif self.program.lower() == 'orca':
            rot_constants = self.rot_temps
            if 0.0 in rot_constants or len(rot_constants) == 1:
                rot_constant = [e for e in rot_constants if e != 0.0]  
                qrot = ((k_b*T)/(self.symmetry_num*h*c*100*rot_constant[0]))
            else:
                qrot = (1/self.symmetry_num) * ((k_b*T)/(h*c*100))**1.5 * (pi / (rot_constants[0] * rot_constants[1] * rot_constants[2]))**0.5


        # Translational partition function
        mol_mass = self.mol_mass * 1.66053906660e-27
        V = k_b*T/P # R*T/P
        qtrans = ((2 * pi * mol_mass * k_b * T) / h**2)**(3/2) * V

        if 'OH' in self.name:
            qelec = 3.019 # OH radical with 2 low lying near degenerate energy level
        else:
            qelec = self.mult

        self.Q = qvib*qrot*qtrans*qelec


    def perturb_active_site(self, indexes=None, scaling_factor=0.1):
        perturbed_coords = []
        if not self.constrained_indexes:
            self.find_active_site(indexes)
        C_index = self.constrained_indexes['C']
        H_index = self.constrained_indexes['H']
        O_index = self.constrained_indexes['O']
        OH_index = O_index + 1
        C_index, H_index, O_index, OH_index = C_index-1, H_index-1, O_index-1, OH_index-1
        original_coordinates = self.coordinates
        perturbed_coords = original_coordinates.copy()
        for index in [C_index, H_index, O_index]:  # Adjusted for 0 based indexing in ORCA
            random_perturbation = [random.uniform(-scaling_factor, scaling_factor) for _ in range(3)]
            perturbed_coords[index] = [coord + delta for coord, delta in zip(original_coordinates[index], random_perturbation)]
        self.coordinates = perturbed_coords


    def set_active_site(self, indexes=None, distance_CH=None, distance_HO=None, distance_OH_H=None, reaction_angle=None):
        if not self.constrained_indexes:
            self.find_active_site(indexes)
        C_index = self.constrained_indexes['C']
        H_index = self.constrained_indexes['H']
        O_index = self.constrained_indexes['O']
        OH_index = O_index + 1
        C_index, H_index, O_index, OH_index = C_index-1, H_index-1, O_index-1, OH_index-1

        # If distances and angle are not given, use predefined values
        if distance_CH is None:
            distance_CH = 1.35
        if distance_HO is None:
            distance_HO = 1.38
        if distance_OH_H is None:
            distance_OH_H = 0.97
        if reaction_angle is None:
            reaction_angle = 173.0

        # Set the C-H distance
        vector_CH = self.calculate_vector(self.coordinates[C_index], self.coordinates[H_index])
        norm_vector_CH = self.normalize_vector(vector_CH)
        new_H_position = self.coordinates[C_index] - (norm_vector_CH * distance_CH)
        if dot(self.calculate_vector(self.coordinates[C_index], new_H_position), vector_CH) < 0:
            new_H_position = self.coordinates[C_index] + (norm_vector_CH * distance_CH)

        # Set the H-C-O angle
        complement_angle = 180.0 - reaction_angle
        rotation_angle = radians(complement_angle)
        rotation_axis = self.normalize_vector(cross(norm_vector_CH, [1, 0, 0]))
        if all(rotation_axis == 0):
            rotation_axis = self.normalize_vector(cross(norm_vector_CH, [0, 1, 0]))
        rotated_vector = self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle)

        # Set the H-O distance
        new_O_position = new_H_position - (rotated_vector * distance_HO)
        if dot(self.calculate_vector(self.coordinates[C_index], new_O_position), vector_CH) < 0:
            new_O_position = new_H_position + rotated_vector * distance_HO 
        rotation_axis = self.normalize_vector(rotation_axis)

        rotation_angle_H = radians(104.5)
        new_OH_H_position = new_O_position - self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H
        HOH_angle = self.calculate_angle(new_OH_H_position, new_O_position, new_H_position)
        if HOH_angle < 95:
            new_OH_H_position = new_O_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H

        # Update positions
        self.coordinates[H_index] = new_H_position
        self.coordinates[O_index] = new_O_position
        self.coordinates[OH_index] = new_OH_H_position


    def H_abstraction(self, Cl=False, products=False, num_molecules=None):
        original_coords = self.coordinates.copy()
        atoms = self.atoms
        abstraction_molecules = []
        product_molecules = []
        methyl_C_indexes, aldehyde_groups, ketone_methyl_groups = self.identify_functional_groups()
        aldehyde_H_indexes = {group['H'] for group in aldehyde_groups}
        carbon_iteration_counter = {index: 0 for index in range(len(atoms)) if atoms[index] == 'C'}

        for i, atom in enumerate(atoms):
            if atom != 'H':
                continue
            # Reset parameters for each hydrogen iteration
            distance_CH = 1.35
            distance_OH = 1.38
            distance_OH_H = 0.97
            water_angle = 104.5
            perp_axis = [1, 0, 0]
            perp_axis2 = [0, 1, 0]
            reaction_angle = 170

            is_aldehyde_H = i in aldehyde_H_indexes
            if is_aldehyde_H:
                # Find the specific aldehyde group this H belongs to
                for group in aldehyde_groups:
                    if group['H'] == i:
                        aldehyde_O = group['O']
                        aldehyde_C = group['C']
                        # Adjust settings specifically for aldehyde H
                        # distance_CH = 1.19481
                        distance_CH = 1.41
                        distance_OH = 1.40971
                        reaction_angle = 153.569
                        water_angle = 99.589
                        perp_axis = self.normalize_vector(self.calculate_vector(original_coords[aldehyde_O], original_coords[aldehyde_C]))
                        break
            for j, other_atom in enumerate(atoms):
                if other_atom == "C" and self.atom_distance(self.coordinates[i], self.coordinates[j]) < 1.3:
                    if j in methyl_C_indexes and carbon_iteration_counter[j] >= 1:
                        continue
                    carbon_iteration_counter[j] += 1

                    vector_CH = self.calculate_vector(original_coords[j], original_coords[i])
                    dist_CH = self.vector_length(vector_CH)
                    norm_vector_CH = self.normalize_vector(vector_CH)

                    new_coords = original_coords.copy()
                    new_atoms = atoms.copy()

                    if products:
                        product_coords = original_coords.copy()
                        product_atoms = atoms.copy()
                        product_coords.pop(i)
                        product_atoms.pop(i)
                        # Add the product molecule to the list
                        product_molecule = Molecule(self.file_path, product=True, smiles=self.smiles)
                        product_molecule.atoms = product_atoms
                        product_molecule.coordinates = product_coords
                        product_molecules.append(product_molecule)

                    new_H_position = array(original_coords[i]) + norm_vector_CH * (dist_CH - distance_CH)
                    new_distance_CH = self.vector_length(self.calculate_vector(new_H_position, new_coords[j]))
                    if new_distance_CH < dist_CH:
                        new_H_position = array(original_coords[i]) - norm_vector_CH * (dist_CH - new_distance_CH)

                    rotation_axis = self.normalize_vector(cross(perp_axis, norm_vector_CH))
                    if all(rotation_axis == 0):
                        rotation_axis = self.normalize_vector(cross(perp_axis2, norm_vector_CH))

                    rotation_angle = radians(180 - reaction_angle)
                    rotated_vector = self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle)

                    new_oxygen_position = new_H_position - (rotated_vector * distance_OH)
                    if self.atom_distance(new_oxygen_position, original_coords[j]) < distance_OH:
                        new_oxygen_position = new_H_position + (rotated_vector * distance_OH)

                    rotation_angle_H = radians(180-water_angle)
                    norm_vector_OH = self.normalize_vector(self.calculate_vector(new_H_position, new_oxygen_position))
                    new_OH_H_position = new_oxygen_position - self.rotate_vector(norm_vector_OH, rotation_axis, rotation_angle_H) * distance_OH_H
                    if self.atom_distance(new_OH_H_position, new_H_position) < 1.5:
                        new_OH_H_position = new_oxygen_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H

                    new_coords[i] = (new_H_position).tolist()

                    if Cl:
                        new_atoms.append('Cl')
                        new_coords.append(new_oxygen_position.tolist())
                        constrained_indexes = {'C': j+1, 'H': i+1, 'O': len(new_coords)}
                    else:
                        new_atoms.append('O')
                        new_atoms.append('H')
                        new_coords.append(new_oxygen_position.tolist())
                        new_coords.append(new_OH_H_position.tolist())
                        constrained_indexes = {'C': j+1, 'H': i+1, 'O': len(new_coords)-1, 'OH': len(new_coords)}

                    new_molecule = Molecule(self.file_path, smiles=self.smiles)
                    new_molecule.atoms = new_atoms
                    new_molecule.coordinates = new_coords
                    new_molecule.constrained_indexes = constrained_indexes
                    abstraction_molecules.append(new_molecule)

        if num_molecules is not None and 0 <= num_molecules <= len(abstraction_molecules):
            return abstraction_molecules[:num_molecules], product_molecules[:num_molecules]
        else:
            return abstraction_molecules, product_molecules


    def OH_addition(self, distance=1.45, double_bond_distance=1.36, dist_OH=0.97):
        atoms = self.atoms
        coordinates = self.coordinates
        num_atoms = len(atoms)
        modified_coords = []

        for i in range(num_atoms):
            if atoms[i] == "C":
                for j in range(num_atoms):
                    if atoms[j] == "C" and i != j:
                        vector_cc = self.calculate_vector(coordinates[i], coordinates[j])
                        dist_cc = self.vector_length(vector_cc)
                        if dist_cc <= double_bond_distance:
                            norm_vector_cc = self.normalize_vector(vector_cc)
                            new_coords = coordinates.copy()

                            perpendicular_axis = cross(norm_vector_cc, [0,1,0])
                            perpendicular_axis = self.normalize_vector(perpendicular_axis)
                            # shift carbon
                            new_coords[i][1:] = array(new_coords[i][1:]) + norm_vector_cc * 0.1
                            # update oxygen in OH coordiantes
                            oxygen_coords = array(new_coords[j][1:]) + perpendicular_axis * distance

                            rotation_axis = norm_vector_cc
                            rotation_angle_h = radians(45) 

                            rotated_vector_h = self.rotate_vector(perpendicular_axis, rotation_axis, rotation_angle_h)

                            hydrogen_coords = oxygen_coords + rotated_vector_h * dist_OH
                            
                            new_coords.append(['O', *oxygen_coords])
                            new_coords.append(['H', *hydrogen_coords])
                            
                            C_index = j+1
                            H_index = i+1
                            O_index = len(new_coords)-1

                            modified_coords.append((new_coords, (C_index, H_index, O_index)))

        coordinates = modified_coords

    
    def addition(self, other_molecule, distance=1.55, double_bond_distance=1.36, elongation_factor=0.1):
        pass 

         
    def identify_functional_groups(self, distance=1.5, angle_tolerance=5):
        methyl_C_indexes = []
        aldehyde_groups = []
        ketone_methyl_groups = []

        for i, atom in enumerate(self.atoms):
            if atom == 'C':
                H_neighbors = [j for j, other_atom in enumerate(self.atoms) if other_atom == 'H' and self.atom_distance(self.coordinates[i], self.coordinates[j]) < distance]
                if len(H_neighbors) >= 3:
                    # Check for neighboring carbons that could be part of a ketone
                    C_neighbors = [j for j, other_atom in enumerate(self.atoms) if other_atom == 'C' and self.atom_distance(self.coordinates[i], self.coordinates[j]) < distance]
                    for C_neighbor in C_neighbors: # Look for an oxygen atom double-bonded to the neighboring carbon
                        O_neighbors = [k for k, other_atom in enumerate(self.atoms) if other_atom == 'O' and self.atom_distance(self.coordinates[C_neighbor], self.coordinates[k]) < distance]
                        for O in O_neighbors: # Check the angle to see if it's roughly 109.5 degrees, indicative of a ketone
                            angle = self.calculate_angle(self.coordinates[i], self.coordinates[C_neighbor], self.coordinates[O])
                            if abs(angle - 109.5) <= angle_tolerance:
                                ketone_methyl_groups.append({'methyl_C': i, 'ketone_C': C_neighbor, 'O': O})
                    methyl_C_indexes.append(i)
                elif len(H_neighbors) <= 2:  # Potential for being part of an aldehyde group
                    # Look for an oxygen atom double-bonded to this carbon
                    O_neighbors = [j for j, other_atom in enumerate(self.atoms) if other_atom == 'O' and self.atom_distance(self.coordinates[i], self.coordinates[j]) < distance]
                    for O in O_neighbors:
                        aldehyde_groups.append({'C': i, 'H': H_neighbors[0], 'O': O})

        return methyl_C_indexes, aldehyde_groups, ketone_methyl_groups


    def smiles_to_atoms_coordinates(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        conf = mol.GetConformer()
        coordinates = [(conf.GetAtomPosition(atom.GetIdx()).x, conf.GetAtomPosition(atom.GetIdx()).y, conf.GetAtomPosition(atom.GetIdx()).z) for atom in mol.GetAtoms()]
        
        molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        return molecular_formula, atoms, coordinates

    def get_terminal_O(self, distance=1.5):
        O_index = []
        for i, atom in enumerate(self.coordinates):
            if atom[0] == "O":
                neighbors = sum(1 for j, other_atom in enumerate(self.coordinates) if i != j and self.atom_distance(atom[1:], other_atom[1:]) < distance)
                if neighbors == 1:
                    O_index.append(i)
        return O_index[0]

    def log2program(self):
        try:
            with open(self.log_file_path, 'r') as file:
                for _ in range(10):
                    line = file.readline()
                    if "Entering Gaussian System" in line or "Gaussian(R)" in line:
                        return 'g16'
                    elif '* O   R   C   A *' in line:
                        return 'orca'
                    elif '|                 C R E S T                  |' in line:
                        return 'crest'
        except FileNotFoundError:
            print(f"File not found: {self.log_file_path}")
            return None
        except Exception as e:
            print(f"An error occurred: {e}")
            return None
        return None


    def log2method(self):
        methods = [
            "wb97xd", "wb97x-d3", "wb97x-d3bj", "wb97x-d4", "b97-3c", "r2scan-3c", "pm3", "am1",
            "pm6", "pm7", 'g3mp2', 'g3', "b3lyp", "m062x", "m06-2x", "m08", "dlpno-ccsd(t)"
        ]
        # if self.program.lower() == 'crest' or 'collection' in self.name or 'crest' in self.name.lower():
        #     return 'crest'
        
        file_to_read = self.log_file_path if self.log_file_path else self.file_path if self.file_path else None
        if not file_to_read:
            method_file = os.path.join(self.directory, ".method")
            if os.path.exists(method_file):
                with open(method_file, 'r') as file:
                    method = file.read().strip()
                    if method in methods:
                        return method
            return 'wb97xd' if self.program.lower() == 'g16' else 'wb97x-d3bj'
        
        try:
            with open(file_to_read, 'r') as file:
                for line in file:
                    line_content = ''
                    if self.program.lower() == 'orca':
                        if '|  1> !' in line:
                            line_content = line.split('! ', 1)[-1]
                    elif self.program.lower() == 'g16':
                        if line.strip().startswith('#'):
                            line_content = line
                    
                    for method in methods:
                        if method.lower() in line_content.lower():
                            return method
        except FileNotFoundError:
            method_file = os.path.join(self.directory, ".method")
            if os.path.exists(method_file):
                with open(method_file, 'r') as file:
                    method = file.read().strip()
                    if method in methods:
                        return method
            return 'wb97xd' if self.program.lower() == 'g16' else 'wb97x-d3bj'
            
        return 'wb97x-d3bj' if self.program.lower() == 'orca' else 'wb97xd'


    # def log2method(self):
    #     methods = ["wb97xd", "wb97x-d3", "wb97x-d3bj", "wb97x-d4", "b97-3c", "r2scan-3c", "pm3", "am1", "pm6", "pm7", 'g3mp2', 'g3', "b3lyp", "m062x", "m06-2x", "m08", "dlpno-ccsd(t)"]
    #     if self.program.lower() == 'crest' or 'collection' in self.name or 'crest' in self.name.lower():
    #         return 'crest'
    #     
    #     file_to_read = self.log_file_path if self.log_file_path else self.file_path if self.file_path else None
    #     if not file_to_read:
    #         return 'wb97xd' if self.program.lower() == 'g16' else 'wb97x-d3bj'
    #     
    #     with open(file_to_read, 'r') as file:
    #         for line in file:
    #             line_content = ''
    #             if self.program.lower() == 'orca':
    #                 if '|  1> !' in line:
    #                     line_content = line.split('! ', 1)[-1]
    #             elif self.program.lower() == 'g16':
    #                 if line.strip().startswith('#'):
    #                     line_content = line
    #             
    #             for method in methods:
    #                 if method.lower() in line_content.lower():
    #                     return method
    #
    #     return 'wb97x-d3bj' if self.program.lower() == 'orca' else 'wb97xd'

    def determine_current_step(self):
        # Handle cases not dependent on file contents first
        if self.program.lower() == 'crest' or 'collection' in self.name or 'crest' in self.name.lower():
            return 'crest_sampling'
        if self.reactant or self.product:
            return 'DLPNO' if 'DLPNO' in self.name else 'optimization'

        file_to_read = self.log_file_path if self.log_file_path else self.file_path if self.file_path else None

        if not file_to_read:
            # If there's no file, attempt to determine the current step from the molecule's name
            if 'DLPNO' in self.name:
                return 'DLPNO'
            elif 'TS' in self.name:
                return 'TS_opt_conf' if 'conf' in self.name else 'TS_opt'
            elif 'conf' in self.name:
                return 'opt_constrain_conf'
            return 'opt_constrain'  # Fallback step if no other conditions are met

        with open(file_to_read, 'r') as file:
            for line in file:
                if self.program.lower() == 'orca':
                    line_content = line.split('! ', 1)[-1] if '|  1> !' in line else ''
                elif self.program.lower() == 'g16':
                    line_content = line if line.strip().startswith('#') else ''
                else:
                    continue

                lower_line_content = line_content.lower()
                if 'dlpno' in lower_line_content or 'f12' in lower_line_content:
                    return 'DLPNO'
                elif re.search(r'\bts\b', lower_line_content) or re.search(r'\boptts\b', lower_line_content):
                    return 'TS_opt_conf' if 'conf' in self.name else 'TS_opt'
                elif 'opt' in lower_line_content:
                    return 'opt_constrain_conf' if 'conf' in self.name else 'optimization' if self.product or self.reactant else 'opt_constrain'


    def inp2xyz(self):
        with open(self.file_path, 'r') as file:
            coordinates = []
            element = []
            for line in file:
                parts = line.split()
                if len(parts) >= 4 and parts[0] in ['H', 'C', 'O', 'N', 'S']:
                    element.append(parts[0])
                    coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])

        return (element, coordinates)


    def log2xyz(self, atoms=False):
        if self.program.lower() == "g16":
            match_string='Standard orientation'
            atomic_number_to_symbol = {
                1: 'H',
                6: 'C',
                7: 'N',
                8: 'O',
                16: 'S',
                17: 'Cl'
            }

            with open(self.log_file_path, 'r') as file:
                coordinates = []
                element = []
                start_reading = False
                for line in file:
                    if start_reading:
                        parts = line.split()
                        if len(parts) >= 6 and parts[1].isdigit() and all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in parts[-3:]):
                            element_symbol = atomic_number_to_symbol.get(int(parts[1]), 'Unknown')
                            element.append(element_symbol)
                            coords = [float(parts[3]), float(parts[4]), float(parts[5])]
                            coordinates.append(coords)
                    if match_string in line:
                        start_reading = True
                        coordinates = []
                        element = []
                    if "Rotational" in line:
                        start_reading = False
                if not coordinates:
                    raise ValueError("Coordinates could not be found in the G16 log file.")
                if atoms:
                    return (element, coordinates)
                else: return coordinates

        elif self.program.lower() == "orca":
            match_string='CARTESIAN COORDINATES (ANGSTROEM)'
            with open(self.log_file_path, 'r') as file:
                coordinates = []
                element = []
                start_reading = False

                for line in file:
                    if start_reading:
                        parts = line.split()
                        if len(parts) == 4 and parts[0].isalpha():  # Check if line starts with element symbol
                                element_symbol = parts[0]
                                coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                                element.append(element_symbol)
                                coordinates.append(coords)
                    if match_string in line:
                        start_reading = True
                        coordinates = []
                        element = []
                        parts = line.split()
                        if len(parts) == 4 and parts[0].isalpha():  # Check if line starts with element symbol
                                element_symbol = parts[0]
                                coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                                element.append(element_symbol)
                                coordinates.append(coords)
                    if match_string in line:
                        start_reading = True
                        coordinates = []
                        element = []
                    if "------" in line:
                        continue
                    if "CARTESIAN COORDINATES (A.U.)" in line:
                        start_reading = False
                if not coordinates:
                    raise ValueError("Coordinates could not be found in the ORCA log file.")

                if atoms:
                    return (element, coordinates)
                else: return coordinates


    def print_items(self, logger=None):
        def output(message):
            if logger:
                logger.log(message)
            else:
                print(message)

        output(f"Molecule: {self.name}")
        output(f"Method: {self.method.upper()}")
        output(f"File Path: {self.file_path}")
        output(f"Directory: {self.directory}")
        output(f"Log File Path: {self.log_file_path}")
        output(f"Program: {self.program.upper()}")
        output(f"Reactant: {self.reactant}")
        output(f"Product: {self.product}")
        output(f"Multiplicity: {self.mult}")
        output(f"Charge: {self.charge}")
        output(f"Dipole Moment: {self.dipole_moment}")
        output(f"Workflow: {self.workflow}")
        if 'Cl' in self.atoms[-2:]:
            output(f"Constrained Indexes: [C: {self.constrained_indexes['C']}, H: {self.constrained_indexes['H']}, Cl: {self.constrained_indexes['O']}]")
        else:
            output(f"Constrained Indexes: {self.constrained_indexes}")
        output(f"Electronic Energy: {self.electronic_energy}")
        output(f"Zero Point Correction: {self.zero_point}")
        output(f"Gibbs free energy: {self.free_energy}")
        output(f"Partition Function: {self.Q}")
        output(f"Vibrational Frequencies: {self.vibrational_frequencies}")
        output(f"Current Step: {self.current_step}")
        output(f"Next Step: {self.next_step}")
        output("Atoms and Coordinates:")
        for atom, coord in zip(self.atoms, self.coordinates):
            coord_str = f"{coord[0]:>9.6f} {coord[1]:>9.6f} {coord[2]:>9.6f}"
            output(f"  {atom:2} {coord_str}")
        output("-----------------------------------------------------------------------")


class Logger:
    def __init__(self, log_file):
        self.log_file = log_file

    def log(self, message):
        with open(self.log_file, 'a') as file:
            file.write(message + '\n')

    def log_with_stars(self, message):
        wrapped_message = self.wrap_in_stars(message)
        self.log(wrapped_message)

    def log_results(self, message):
        pass # TODO

    @staticmethod
    def wrap_in_stars(s):
        star_length = len(s) + 14  # 6 stars and 2 spaces padding on each side

        top_bottom_line = '*' * star_length
        middle_line = f"****** {s} ******"
        wrapped_string = f"\n{top_bottom_line}\n{middle_line}\n{top_bottom_line}\n"

        return wrapped_string
