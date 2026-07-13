from enum import Enum
from numpy import array, dot, cos, sin, cross, zeros_like, arccos, degrees, exp, pi, radians, all
from numpy.linalg import norm
import glob
import re
import shutil
import os
import pickle
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from output import console

_G16_SCRATCH_EXTENSIONS = ('.int', '.d2e', '.rwf', '.skr', '.chk')

def _remove_g16_scratch(directory):
    for path in glob.glob(os.path.join(directory, 'Gau-*')):
        if path.endswith(_G16_SCRATCH_EXTENSIONS):
            try:
                os.remove(path)
            except OSError:
                pass

##################################################WORKFLOWS#################################################################
class Step(str, Enum):
    OPT_CONSTRAIN      = 'opt_constrain'
    TS_OPT             = 'TS_opt'
    CREST_SAMPLING     = 'crest_sampling'
    OPT_CONSTRAIN_CONF = 'opt_constrain_conf'
    TS_OPT_CONF        = 'TS_opt_conf'
    OPTIMIZATION       = 'optimization'
    DLPNO              = 'DLPNO'

# Change as you feel like
OH_H2O_workflow          = (Step.OPTIMIZATION, Step.DLPNO)
reactant_product_workflow = (Step.CREST_SAMPLING, Step.OPTIMIZATION, Step.DLPNO)
TS_workflow               = (Step.OPT_CONSTRAIN, Step.TS_OPT, Step.CREST_SAMPLING,
                             Step.OPT_CONSTRAIN_CONF, Step.TS_OPT_CONF, Step.DLPNO)
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
        self.reactant = True if 'reactant' in self.name.split("_") else reactant
        self.product = True if 'product' in self.name.split("_") else product
        self.job_id = ""
        self.atoms = atoms if atoms is not None else []
        self.coordinates = coordinates if coordinates is not None else []
        self.constrained_indexes = {'C': indexes[0], 'H': indexes[1], 'X': indexes[2], 'XH': indexes[2]+1} if indexes else {}
        self.converged = False
        self.mult = self._determine_mult()
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
        self.reaction_path_degeneracy = 1  # σᵢ: equivalent H's this TS represents
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

            if self.reactant or 'reactant' in self.name or self.name in ('OH', 'OH_DLPNO', 'Cl', 'Cl_DLPNO', 'NO3', 'NO3_DLPNO'):
                self.reactant = True
            elif self.product or 'product' in self.name or self.name in ('H2O', 'H2O_DLPNO', 'HCl', 'HCl_DLPNO', 'HNO3', 'HNO3_DLPNO'):
                self.product = True
            else:
                self.find_active_site(indexes)
            self.mult = self._determine_mult()

            self.workflow = self.determine_workflow()
            self.set_current_step()
        if not self.method:
            self.method = self.log2method()


    def set_current_step(self, step=None):
        if step is not None:
            if step in self.workflow:
                self.current_step = step
            else:
                console.warning(f"Step '{step}' not found in the workflow of {self.name}.")
                return
        else:
            detected = self.determine_current_step()
            self.current_step = detected if detected in self.workflow else (self.workflow[0] if self.workflow else None)

        idx = list(self.workflow).index(self.current_step) if self.current_step in self.workflow else 0
        self.next_step = self.workflow[idx + 1] if idx + 1 < len(self.workflow) else None


    @property
    def small_molecule(self):
        return self.name in ('OH', 'H2O', 'H2O_DLPNO', 'OH_DLPNO',
                             'Cl', 'HCl', 'Cl_DLPNO', 'HCl_DLPNO',
                             'NO3', 'HNO3', 'NO3_DLPNO', 'HNO3_DLPNO')

    def determine_workflow(self):
        if self.small_molecule:
            return OH_H2O_workflow
        elif self.reactant or self.product:
            return reactant_product_workflow
        else:
            return TS_workflow

    def _determine_mult(self):
        if 'H2O' in self.name or self.name in ('HCl', 'HCl_DLPNO', 'HNO3', 'HNO3_DLPNO'):
            return 1   # water/HCl/HNO3 — singlet
        if 'OH' in self.name and not self.product:
            return 2   # OH radical — doublet
        if self.name in ('Cl', 'Cl_DLPNO', 'NO3', 'NO3_DLPNO'):
            return 2   # Cl/NO3 radical — doublet
        if self.reactant:
            return 1   # closed-shell organic reactant
        if self.product:
            return 2   # alkyl radical product — doublet
        return 2       # TS — open-shell doublet


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

    @classmethod
    def create_Cl(cls, directory=os.getcwd()):
        # Creates a Cl radical molecule (reactant for Cl-radical H abstraction)
        molecule = cls(name='Cl', directory=directory, atoms=['Cl'], coordinates=[[0.0, 0.0, 0.0]], reactant=True)
        molecule.mult = 2
        molecule.reactant = True
        molecule.workflow = molecule.determine_workflow()
        molecule.converged = False
        molecule.set_current_step()
        molecule.method = molecule.log2method()
        return molecule

    @classmethod
    def create_HCl(cls, directory=os.getcwd()):
        # Creates a HCl molecule (product of Cl-radical H abstraction)
        molecule = cls(name='HCl', directory=directory, atoms=['H', 'Cl'], coordinates=[[0.0, 0.0, 0.0], [1.275, 0.0, 0.0]], product=True)
        molecule.mult = 1
        molecule.product = True
        molecule.workflow = molecule.determine_workflow()
        molecule.converged = False
        molecule.set_current_step()
        molecule.method = molecule.log2method()
        return molecule

    @classmethod
    def create_NO3(cls, directory=os.getcwd()):
        # Creates a NO3 radical molecule (reactant for NO3-radical H abstraction)
        # Planar D3h: N at center, 3 O atoms at 120°, N-O bond = 1.24 Å
        N_O = 1.24
        molecule = cls(name='NO3', directory=directory,
                       atoms=['N', 'O', 'O', 'O'],
                       coordinates=[[0.0, 0.0, 0.0],
                                    [N_O, 0.0, 0.0],
                                    [-N_O * 0.5,  N_O * 0.866, 0.0],
                                    [-N_O * 0.5, -N_O * 0.866, 0.0]],
                       reactant=True)
        molecule.mult = 2
        molecule.reactant = True
        molecule.workflow = molecule.determine_workflow()
        molecule.converged = False
        molecule.set_current_step()
        molecule.method = molecule.log2method()
        return molecule

    @classmethod
    def create_HNO3(cls, directory=os.getcwd()):
        # Creates a HNO3 molecule (product of NO3-radical H abstraction)
        # Planar structure (Cs symmetry), N-O(H) along +x axis.
        # Exp. geometry: N-OH=1.406, N=O=1.206/1.211 Å, O-H=0.964 Å
        # Angles: ∠O=N=O=130.3°, ∠O(H)-N=O=115.9°/113.8°, ∠H-O-N=102.2°
        # Cis conformer (H cis to the shorter N=O oxygen).
        molecule = cls(name='HNO3', directory=directory,
                       atoms=['H', 'O', 'N', 'O', 'O'],
                       coordinates=[[ 1.609,  0.943, 0.0],   # H  (cis to O-nitro-1)
                                    [ 1.406,  0.000, 0.0],   # O  (hydroxyl, N-O-H)
                                    [ 0.000,  0.000, 0.0],   # N
                                    [-0.524,  1.086, 0.0],   # O  (nitro-1, N=O 1.206 Å, cis to H)
                                    [-0.490, -1.107, 0.0]],  # O  (nitro-2, N=O 1.211 Å, trans to H)
                       product=True)
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
            console.error(f"File not found: {self.file_path}")

    def save_to_pickle(self, file_path):
        with open(file_path, 'wb') as file:
            pickle.dump(self, file)

    def move_files(self, ignore_file=''):
        try:
            dest_directory = os.path.join(self.directory, "log_files")
            if not os.path.exists(dest_directory):
                os.makedirs(dest_directory, exist_ok=True)

            for filename in os.listdir(self.directory):
                if filename.endswith(self.output) and "JKTS" not in filename:
                    file_path = os.path.join(self.directory, filename)
                    dest_file_path = os.path.join(dest_directory, filename)
                    if os.path.exists(dest_file_path):
                        os.remove(dest_file_path)
                    shutil.move(file_path, dest_file_path)
        except (OSError, shutil.Error):
            pass
        _remove_g16_scratch(self.directory)

    def move_converged(self):
        try:
            if not os.path.exists(os.path.join(self.directory, "log_files")):
                os.makedirs(os.path.join(self.directory, "log_files"), exist_ok=True)
            destination = os.path.join(self.directory, "log_files", f"{self.name}{self.output}")
            if os.path.exists(destination):
                os.remove(destination)
            shutil.move(os.path.join(self.directory, f"{self.name}{self.output}"), destination)
            self.log_file_path = destination
        except Exception:
            pass

    def move_inputfile(self):
        try:
            input_files_dir = os.path.join(self.directory, "input_files")
            if not os.path.exists(input_files_dir):
                os.makedirs(input_files_dir, exist_ok=True)
            source = os.path.join(self.directory, f"{self.name}{self.input}")
            destination = os.path.join(input_files_dir, f"{self.name}{self.input}")
            if os.path.exists(destination):
                os.remove(destination)
            if os.path.exists(source):
                shutil.move(source, destination)
        except Exception:
            pass

    def move_failed(self):
        try:
            if not os.path.exists(os.path.join(self.directory, "failed_logs")):
                os.makedirs(os.path.join(self.directory, "failed_logs"), exist_ok=True)
            destination = os.path.join(self.directory, "failed_logs", f"{self.name}_{self.error_termination_count}{self.output}")
            if os.path.exists(destination):
                os.remove(destination)
            shutil.move(os.path.join(self.directory, f"{self.name}{self.output}"), destination)
        except Exception:
            pass
        _remove_g16_scratch(self.directory)

    def write_xyz_file(self, output_file_path):
        with open(output_file_path, 'w') as file:
            file.write(str(len(self.atoms)) + '\n\n')
            for atom, coord in zip(self.atoms, self.coordinates):
                coord_str = ' '.join(['{:.6f}'.format(c) for c in coord])  
                file.write(f'{atom} {coord_str}\n')


    def update_step(self):
        if self.next_step is not None:
            self.set_current_step(self.next_step)
        else:
            console.info(f"{self.name}: reached the end of the workflow.")


    def find_active_site(self, indexes=None):
        # Migrate old 'O'/'OH' keys to new 'X'/'XH' naming if needed
        if self.constrained_indexes and 'O' in self.constrained_indexes and 'X' not in self.constrained_indexes:
            self.constrained_indexes['X'] = self.constrained_indexes.pop('O')
            if 'OH' in self.constrained_indexes:
                self.constrained_indexes['XH'] = self.constrained_indexes.pop('OH')
        if self.constrained_indexes:
            return
        if indexes:
            C_index, H_index, abstractor_index = indexes
            abstractor = self.atoms[abstractor_index - 1] if self.atoms else 'O'
            if abstractor == 'O' and len(self.atoms) >= abstractor_index and self.atoms[abstractor_index] == 'H':
                # OH TS: the atom right after the abstracting O is the OH-H
                self.constrained_indexes = {'C': C_index, 'H': H_index, 'X': abstractor_index, 'XH': abstractor_index + 1}
            else:
                # Cl or NO3 TS (or OH TS where H is not immediately after O)
                self.constrained_indexes = {'C': C_index, 'H': H_index, 'X': abstractor_index}
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
                        # Support both new ('X'/'XH') and old ('O'/'OH') key names
                        X_index = constraints.get('X', constraints.get('O', None))
                        XH_index = constraints.get('XH', constraints.get('OH', None))
                        if X_index is not None and self.atoms:
                            abstractor = self.atoms[X_index - 1]
                            if abstractor == 'Cl' or XH_index is None:
                                # Cl TS or NO3 TS: no secondary radical-H atom
                                self.constrained_indexes = {'C': C_index, 'H': H_index, 'X': X_index}
                            else:
                                # OH TS: secondary H on the radical
                                self.constrained_indexes = {'C': C_index, 'H': H_index, 'X': X_index, 'XH': XH_index}
                        return
                except Exception as e:
                    console.warning(f"Error reading .constrain file: {e}")
                    # If reading .constrain fails, proceed to loop through the molecule
                    indexes = None

            if not indexes:  # Fallback to original logic if indexes not set
                if self.atoms and self.atoms[-1] == 'Cl':
                    # Cl radical: single abstractor atom (Cl) appended at the end
                    Cl_index = len(self.atoms) - 1  # 0-based
                    for C, atom in enumerate(self.atoms[:-1]):
                        if atom == 'C':
                            for H, neighbor in enumerate(self.atoms[:-1]):
                                if neighbor == 'H' and self.is_nearby(C, H):
                                    if self.calculate_angle(self.coordinates[C], self.coordinates[H], self.coordinates[Cl_index]) > 120:
                                        if self.is_nearby(H, Cl_index, threshold_distance=3.0):
                                            self.constrained_indexes = {'C': C+1, 'H': H+1, 'X': Cl_index+1}
                                            return
                elif len(self.atoms) >= 4 and self.atoms[-3] == 'N' and self.atoms[-4] == 'O':
                    # NO3 radical: last 4 atoms are O(abstracting), N, O, O
                    X_0based = len(self.atoms) - 4  # 0-based index of abstracting O
                    for C, atom in enumerate(self.atoms[:-4]):
                        if atom == 'C':
                            for H, neighbor in enumerate(self.atoms[:-4]):
                                if neighbor == 'H' and self.is_nearby(C, H):
                                    if self.calculate_angle(self.coordinates[C], self.coordinates[H], self.coordinates[X_0based]) > 120:
                                        if self.is_nearby(H, X_0based, threshold_distance=2.5):
                                            self.constrained_indexes = {'C': C+1, 'H': H+1, 'X': X_0based+1}
                                            return
                else:
                    # OH radical: last two atoms are O and H of the OH radical
                    abstractor_0based = len(self.atoms) - 2  # Second-last atom is the abstracting O
                    H_OH_index = len(self.atoms) - 1  # Last atom is the H of the OH radical

                    for C, atom in enumerate(self.atoms):
                        if atom == 'C':
                            for H, neighbor in enumerate(self.atoms):
                                if neighbor == 'H' and self.is_nearby(C, H):
                                    if self.calculate_angle(self.coordinates[C], self.coordinates[H], self.coordinates[abstractor_0based]) > 120:
                                        if self.is_nearby(H, abstractor_0based, threshold_distance=2):
                                            if self.is_nearby(C, abstractor_0based, threshold_distance=2.8):
                                                self.constrained_indexes = {'C': C+1, 'H': H+1, 'X': abstractor_0based+1, 'XH': H_OH_index+1}
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
        with open(file_path, 'wb') as file:
            pickle.dump(molecules, file)

    @staticmethod
    def load_molecules_from_pickle(file_path):
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
                partition_function = re.findall(r"Total V=0\s+([-\d.]+D[+-]\d+)", log_content)
                if zero_point_correction:
                    self.zero_point = float(zero_point_correction[-1])
                    self.free_energy = float(free_energy[-1])
                    self.dipole_moment = float(dipole_moment[-1])
                    self.Q = float(partition_function[-1].replace('D', 'E'))
                electronic_energy = re.findall(r'(SCF Done:  E\(\S+\) =)\s+([-.\d]+)', log_content)
                if electronic_energy:
                    self.electronic_energy = float(electronic_energy[-1][-1])
                    freq_matches = re.findall(r"Frequencies --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)?\s+(-?\d+\.\d+)?", log_content)
                    if freq_matches:
                        flat = [float(freq) for match in freq_matches for freq in match if freq]
                        rot_temp = re.findall(r"Rotational temperatures? \(Kelvin\)\s+(-?\d+\.\d+)(?:\s+(-?\d+\.\d+))?(?:\s+(-?\d+\.\d+))?", log_content)
                        self.rot_temps = [float(rot) for rot in rot_temp[-1] if rot != '']
                        # Retry-path logs contain the spectrum twice; keep only the last complete 3N-6 (3N-5 linear) block
                        n = 3*len(self.atoms) - (5 if len(self.rot_temps) == 1 else 6)
                        if len(flat) != n:
                            (logger or console).warning(f"{self.name}: parsed {len(flat)} frequencies, expected {n}; keeping last {n}")
                        self.vibrational_frequencies = flat[-n:]
                        symmetry_num = re.search(r"Rotational symmetry number\s*(\d+)", log_content)
                        if symmetry_num:
                            self.symmetry_num = int(symmetry_num.group(1))
                        else:
                            (logger or console).warning(f"No symmetry number found in {self.name}; assuming 1")
                            self.symmetry_num = 1
                        mol_mass = re.search(r"Molecular mass:\s+(-?\d+\.\d+)", log_content)
                        self.mol_mass = float(mol_mass.group(1))
                        mult = re.search(r"Multiplicity =\s*(\d+)", log_content)
                        if mult:
                            self.mult = int(mult.group(1))
                        else: self.mult = 2

                        if self.Q is not None:
                            # Gaussian's Q assumes q_elec = mult; correct it for OH/Cl
                            correct_qelec = self._electronic_partition_function()
                            if correct_qelec != self.mult:
                                self.Q *= correct_qelec / self.mult
                    elif 'TS' in self.name:
                        (logger or console).warning(f"No frequencies found in {self.name}")


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
                        rot_temp = re.findall(r"Rotational constants in cm-1: \s*[-+]?(\d*\.\d*)  \s*[-+]?(\d*\.\d*) \s*[-+]?(\d*\.\d*)", log_content)
                        self.rot_temps = [float(rot) for rot in rot_temp[-1]]
                        linear = len(self.rot_temps) == 1 or 0.0 in self.rot_temps
                        n = 3*len(self.atoms) - (5 if linear else 6)
                        self.vibrational_frequencies = [float(freq) for freq in freq_matches][-n:]
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

        qelec = self._electronic_partition_function(T)

        self.Q = qvib*qrot*qtrans*qelec


    def _electronic_partition_function(self, T=298.15):
        h = 6.62607015e-34  # Planck constant in J.s
        k_b = 1.380649e-23  # Boltzmann constant in J/K
        c = 299792458       # Speed of light in m/s
        if self.name in ('OH', 'OH_DLPNO'):
            return 3.019  # ²Π spin-orbit splitting (~140 cm⁻¹)
        elif self.name in ('Cl', 'Cl_DLPNO'):
            return 4 + 2 * exp(-(882 * 100 * h * c) / (k_b * T))  # ²P₃/₂ + ²P₁/₂ at 882 cm⁻¹
        else:
            return self.mult


    def perturb_active_site(self, indexes=None, scaling_factor=0.1):
        if not self.constrained_indexes:
            self.find_active_site(indexes)
        C_index = self.constrained_indexes['C'] - 1
        H_index = self.constrained_indexes['H'] - 1
        abstractor_index = self.constrained_indexes['X'] - 1
        original_coordinates = self.coordinates
        perturbed_coords = original_coordinates.copy()
        for index in [C_index, H_index, abstractor_index]:
            random_perturbation = [random.uniform(-scaling_factor, scaling_factor) for _ in range(3)]
            perturbed_coords[index] = [coord + delta for coord, delta in zip(original_coordinates[index], random_perturbation)]
        self.coordinates = perturbed_coords


    def set_active_site(self, indexes=None, distance_CH=None, distance_HX=None, distance_XH=None, reaction_angle=None):
        if not self.constrained_indexes:
            self.find_active_site(indexes)
        C_index = self.constrained_indexes['C'] - 1
        H_index = self.constrained_indexes['H'] - 1
        abstractor_index = self.constrained_indexes['X'] - 1
        abstractor = self.atoms[abstractor_index]

        # If distances and angle are not given, use predefined values
        if distance_CH is None:
            distance_CH = 1.35
        if distance_HX is None:
            distance_HX = 1.84 if abstractor == 'Cl' else 1.38
        if distance_XH is None:
            distance_XH = 0.97
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

        # Set the H-X distance (X = abstracting atom)
        new_X_position = new_H_position - (rotated_vector * distance_HX)
        if dot(self.calculate_vector(self.coordinates[C_index], new_X_position), vector_CH) < 0:
            new_X_position = new_H_position + rotated_vector * distance_HX
        rotation_axis = self.normalize_vector(rotation_axis)

        # Set the XH position (only for OH TS: the H that stays on the radical's O)
        if 'XH' in self.constrained_indexes:
            radical_H_index = self.constrained_indexes['XH'] - 1
            rotation_angle_XH = radians(104.5)
            new_XH_position = new_X_position - self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_XH) * distance_XH
            XH_angle = self.calculate_angle(new_XH_position, new_X_position, new_H_position)
            if XH_angle < 95:
                new_XH_position = new_X_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_XH) * distance_XH
            self.coordinates[radical_H_index] = new_XH_position

        # Update positions
        self.coordinates[H_index] = new_H_position
        self.coordinates[abstractor_index] = new_X_position


    def H_abstraction(self, Cl=False, NO3=False, products=False, num_molecules=None):
        original_coords = self.coordinates.copy()
        atoms = self.atoms
        method = self.method
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
            distance_HX = 1.38    # H to abstracting atom (X) distance in TS
            distance_XH = 0.97    # X-H bond in radical (only used for OH)
            radical_angle = 104.5 # H-X-XH angle in radical product (only used for OH)
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
                        distance_CH = 1.12092
                        distance_HX = 1.84524
                        reaction_angle = 145.391
                        radical_angle = 94.051
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
                        product_molecule = Molecule(self.file_path, product=True, smiles=self.smiles, method=method)
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

                    # Place the abstracting atom X (O for OH/NO3, Cl for Cl) along the reaction vector
                    new_X_position = new_H_position - (rotated_vector * distance_HX)
                    if self.atom_distance(new_X_position.tolist(), original_coords[j]) < distance_HX:
                        new_X_position = new_H_position + (rotated_vector * distance_HX)

                    # Compute XH position for OH TS (the H that stays bonded to the abstracting O)
                    rotation_angle_XH = radians(180 - radical_angle)
                    norm_vector_XH = self.normalize_vector(self.calculate_vector(new_H_position, new_X_position))
                    new_XH_position = new_X_position - self.rotate_vector(norm_vector_XH, rotation_axis, rotation_angle_XH) * distance_XH
                    if self.atom_distance(new_XH_position.tolist(), new_H_position.tolist()) < 1.5:
                        new_XH_position = new_X_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_XH) * distance_XH

                    new_coords[i] = new_H_position.tolist()

                    if Cl:
                        # Cl TS: place Cl at H-Cl TS distance (~1.84 Å)
                        distance_HCl = 1.84
                        new_Cl_position = new_H_position - (rotated_vector * distance_HCl)
                        if self.atom_distance(new_Cl_position.tolist(), original_coords[j]) < distance_HCl:
                            new_Cl_position = new_H_position + (rotated_vector * distance_HCl)
                        new_atoms.append('Cl')
                        new_coords.append(new_Cl_position.tolist())
                        constrained_indexes = {'C': j+1, 'H': i+1, 'X': len(new_coords)}
                    elif NO3:
                        # NO3 TS: place abstracting O, then build N and two terminal O atoms
                        N_O_bond = 1.24
                        new_atoms.append('O')   # abstracting O of NO3
                        new_coords.append(new_X_position.tolist())
                        X_coord_index = len(new_coords)  # 1-based index of abstracting O
                        # N placed along rotation_axis (perpendicular to C-H-X reaction plane)
                        new_N_position = new_X_position + rotation_axis * N_O_bond
                        new_atoms.append('N')
                        new_coords.append(new_N_position.tolist())
                        # Two terminal O atoms at ±120° from (N→abstracting O) around the C-H axis
                        N_to_X = new_X_position - new_N_position
                        new_O2_position = new_N_position + self.rotate_vector(N_to_X, norm_vector_CH, radians(120))
                        new_O3_position = new_N_position + self.rotate_vector(N_to_X, norm_vector_CH, radians(-120))
                        new_atoms.append('O')
                        new_coords.append(new_O2_position.tolist())
                        new_atoms.append('O')
                        new_coords.append(new_O3_position.tolist())
                        constrained_indexes = {'C': j+1, 'H': i+1, 'X': X_coord_index}
                    else:
                        # OH TS: abstracting O + H of OH radical
                        new_atoms.append('O')
                        new_atoms.append('H')
                        new_coords.append(new_X_position.tolist())
                        new_coords.append(new_XH_position.tolist())
                        constrained_indexes = {'C': j+1, 'H': i+1, 'X': len(new_coords)-1, 'XH': len(new_coords)}

                    new_molecule = Molecule(self.file_path, smiles=self.smiles, method=method)
                    new_molecule.atoms = new_atoms
                    new_molecule.coordinates = new_coords
                    new_molecule.constrained_indexes = constrained_indexes
                    # σᵢ: a methyl carbon collapses its equivalent H's into one TS
                    if j in methyl_C_indexes:
                        new_molecule.reaction_path_degeneracy = sum(
                            1 for k, a in enumerate(atoms)
                            if a == 'H' and self.atom_distance(original_coords[j], original_coords[k]) < 1.3)
                    else:
                        new_molecule.reaction_path_degeneracy = 1
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
                elif len(H_neighbors) in {1, 2}:
                    # Carbonyl C=O only (~1.20 Å), not a carbinol C-OH
                    O_neighbors = [j for j, other_atom in enumerate(self.atoms) if other_atom == 'O' and self.atom_distance(self.coordinates[i], self.coordinates[j]) < 1.3]
                    for O in O_neighbors:
                        O_has_H = any(other_atom == 'H' and self.atom_distance(self.coordinates[O], self.coordinates[k]) < 1.1
                                      for k, other_atom in enumerate(self.atoms))
                        if O_has_H:
                            continue
                        for H in H_neighbors:  # every H (formaldehyde has two)
                            aldehyde_groups.append({'C': i, 'H': H, 'O': O})

        return methyl_C_indexes, aldehyde_groups, ketone_methyl_groups


    def smiles_to_atoms_coordinates(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            # If it fails, try interpreting the input as a common name
            mol = Chem.MolFromName(smiles)
            if not mol:
                return "Invalid input: neither a valid SMILES nor a recognized chemical name", [], []

        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            return "Could not embed molecule", [], []

        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        conf = mol.GetConformer()
        coordinates = [(conf.GetAtomPosition(atom.GetIdx()).x,
                        conf.GetAtomPosition(atom.GetIdx()).y,
                        conf.GetAtomPosition(atom.GetIdx()).z) for atom in mol.GetAtoms()]

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
            console.error(f"File not found: {self.log_file_path}")
            return None
        except Exception as e:
            console.error(f"Error detecting QC program from {self.log_file_path}: {e}")
            return None
        return None


    def log2method(self):
        methods = [
            "wb97xd", "wb97x-d3", "wb97x-d3bj", "wb97x-d4", "b97-3c", "r2scan-3c", "pm3", "am1",
            "pm6", "pm7", 'g3mp2', 'g3', "b3lyp", "m062x", "m06-2x", "m08", "dlpno-ccsd(t)", "mp2", "bhandhlyp"
        ]

        if self.name in ('OH', 'H2O', 'OH_DLPNO', 'H2O_DLPNO'):
            method_file = os.path.join(self.directory, ".method")
            if os.path.exists(method_file):
                with open(method_file, 'r') as file:
                    method = file.read().strip().lower()
                    if method in methods:
                        return method
            return 'method could not be determined'

        
        file_to_read = self.log_file_path if self.log_file_path else self.file_path if self.file_path else None
        if not file_to_read:
            method_file = os.path.join(self.directory, ".method")
            if os.path.exists(method_file):
                with open(method_file, 'r') as file:
                    method = file.read().strip().lower()
                    if method in methods:
                        return method
            return 'method could not be determined'
        
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
            return 'method could not be determined'
            
        return 'method could not be determined'


    def compare_structures(self, mol2):
        if len(self.atoms) != len(mol2.atoms):
            raise ValueError("Molecules have different number of atoms, comparison not possible.")

        for atom_self, atom_other in zip(self.atoms, mol2.atoms):
            if atom_self != atom_other:
                raise ValueError("Order or atoms between molecules is not the same, make sure order of atoms is identical")
        
        bond_length_errors = []
        bond_angle_errors = []

        # Calculate bond length errors
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                dist1 = Vector.atom_distance(self.coordinates[i], self.coordinates[j])
                dist2 = Vector.atom_distance(mol2.coordinates[i], mol2.coordinates[j])
                bond_length_errors.append(abs(dist1 - dist2))

        # Calculate bond angle errors
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                for k in range(j + 1, len(self.atoms)):
                    if i != j and i != k:
                        angle1 = Vector.calculate_angle(self.coordinates[i], self.coordinates[j], self.coordinates[k])
                        angle2 = Vector.calculate_angle(mol2.coordinates[i], mol2.coordinates[j], mol2.coordinates[k])
                        bond_angle_errors.append(abs(angle1 - angle2))
        
        total_weights = len(bond_length_errors)
        wmae_bond_lengths = sum(bond_length_errors) / total_weights if bond_length_errors else 0
        maxe_bond_lengths = max(bond_length_errors) if bond_length_errors else 0

        total_weights = len(bond_angle_errors)
        wmae_bond_angles = sum(bond_angle_errors) / total_weights if bond_angle_errors else 0
        maxe_bond_angles = max(bond_angle_errors) if bond_angle_errors else 0

        return {
            'WMAE_bond_length': wmae_bond_lengths,
            'MaxE_bond_length': maxe_bond_lengths,
            'WMAE_bond_angle': wmae_bond_angles,
            'MaxE_bond_angle': maxe_bond_angles
        }


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
                        # Standard orientation lines always have exactly 6 fields:
                        # Center# AtomicNum Type X Y Z
                        # parts[0] and parts[1] are both integers (center# and atomic#)
                        if len(parts) == 6 and parts[0].isdigit() and parts[1].isdigit():
                            try:
                                element_symbol = atomic_number_to_symbol.get(int(parts[1]), 'Unknown')
                                element.append(element_symbol)
                                coords = [float(parts[3]), float(parts[4]), float(parts[5])]
                                coordinates.append(coords)
                            except ValueError:
                                pass
                    if match_string in line:
                        start_reading = True
                        coordinates = []
                        element = []
                    if "Rotational" in line:
                        start_reading = False
                if not coordinates:
                    return False
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
                    return False

                if atoms:
                    return (element, coordinates)
                else: return coordinates


    def print_items(self):
        print(f"Molecule: {self.name}")
        print(f"Method: {self.method.upper()}")
        print(f"File Path: {self.file_path}")
        print(f"Directory: {self.directory}")
        print(f"Log File Path: {self.log_file_path}")
        print(f"Program: {self.program.upper()}")
        print(f"Reactant: {self.reactant}")
        print(f"Product: {self.product}")
        print(f"Multiplicity: {self.mult}")
        print(f"Charge: {self.charge}")
        print(f"Dipole Moment: {self.dipole_moment}")
        print(f"Workflow: {self.workflow}")
        if self.constrained_indexes and 'X' in self.constrained_indexes:
            X_idx = self.constrained_indexes['X']
            abstractor_sym = self.atoms[X_idx - 1] if self.atoms and X_idx <= len(self.atoms) else '?'
            label = 'Cl' if abstractor_sym == 'Cl' else ('NO3-O' if abstractor_sym == 'O' and 'XH' not in self.constrained_indexes and len(self.atoms) >= 4 and self.atoms[X_idx] == 'N' else 'O')
            print(f"Constrained Indexes: [C: {self.constrained_indexes['C']}, H: {self.constrained_indexes['H']}, {label}: {X_idx}]")
        else:
            print(f"Constrained Indexes: {self.constrained_indexes}")
        print(f"Electronic Energy: {self.electronic_energy}")
        print(f"Zero Point Correction: {self.zero_point}")
        print(f"Gibbs free energy: {self.free_energy}")
        print(f"Partition Function: {self.Q}")
        print(f"Vibrational Frequencies: {self.vibrational_frequencies}")
        print(f"Current Step: {self.current_step}")
        print(f"Next Step: {self.next_step}")
        print("Atoms and Coordinates:")
        for atom, coord in zip(self.atoms, self.coordinates):
            coord_str = f"{coord[0]:>9.6f} {coord[1]:>9.6f} {coord[2]:>9.6f}"
            print(f"  {atom:2} {coord_str}")
        print("-----------------------------------------------------------------------")


