#!/usr/bin/env python3
'''Dynamic Approach for Transition State'''
###############################LIBRARIES#####################################
import numpy as np
import pandas as pd
import argparse
import os
import re
import subprocess
import time
import threading

###############################CONSTANS#####################################
T = 298.15
h = 6.62607015e-34
k_b = 1.380649e-23
c = 29979245800
R = 8.314 # J/mol*K 

np.set_printoptions(suppress=True, precision=6)
###########################VECTOR MANIPULATION################################
class VectorManipulation:
    @staticmethod
    def calculate_vector(coord1, coord2):
        return np.array(coord1) - np.array(coord2)

    @staticmethod
    def vector_length(vector):
        return np.linalg.norm(vector)

    @staticmethod
    def atom_distance(atom1, atom2):
        return np.linalg.norm(np.array(atom2) - np.array(atom1))

    @staticmethod
    def normalize_vector(vector):
        norm = np.linalg.norm(vector)
        if norm < 1e-8:
            return np.zeros_like(vector)
        return vector / norm

    @staticmethod
    def rotate_vector(vector, axis, angle):
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        cross_product = np.cross(axis, vector)
        return (vector * cos_theta +
                cross_product * sin_theta +
                axis * np.dot(axis, vector) * (1 - cos_theta))

    @staticmethod
    def get_upwards_perpendicular_axis(self, direction):
        # direction is the vector lying between 2 carbons in C=C
        pass


class Logger:
    def __init__(self, log_file):
        self.log_file = log_file

    def log(self, message):
        with open(self.log_file, 'a') as file:
            file.write(message + '\n')

    def log_with_stars(self, message):
        wrapped_message = self.wrap_in_stars(message)
        self.log(wrapped_message)

    @staticmethod
    def wrap_in_stars(s):
        star_length = len(s) + 14  # 6 stars and 2 spaces padding on each side

        top_bottom_line = '*' * star_length
        middle_line = f"****** {s} ******"
        wrapped_string = f"\n{top_bottom_line}\n{middle_line}\n{top_bottom_line}\n"

        return wrapped_string


   #########################################SUBMIT JOBS############################
def submit_job(molecule, input_file, ncpus, mem, partition, time, nnodes=1):
    job_name = molecule.name
    file_extension = molecule.input
    pwd = molecule.directory    
    path_submit_script = os.path.join(pwd, f"{job_name}_submit.sh")
    submit_file = os.path.join(pwd, "qsub.tmp")

    if molecule.program.lower() == "orca" or molecule.program.lower() == "g16":
        program_mem = mem + 2000
    else:
        program_mem = mem

    with open(path_submit_script, 'w') as file:
        file.write("#!/bin/bash\n")
        file.write(f"SUBMIT={submit_file}\n\n") 
        file.write("cat > $SUBMIT <<!EOF\n")
        file.write("#!/bin/sh\n")
        file.write(f"#SBATCH --job-name={job_name}\n")
        if molecule.program.lower() == 'g16':
            file.write(f"#SBATCH --nodes={nnodes}\n")
            file.write(f"#SBATCH --cpus-per-task={ncpus}\n")
            file.write(f"#SBATCH --ntasks={nnodes}\n")
        else:
            file.write(f"#SBATCH --nodes={nnodes}\n")
            file.write(f"#SBATCH --cpus-per-task=1\n")
            file.write(f"#SBATCH --ntasks={ncpus}\n")
        file.write(f"#SBATCH --error={pwd}/{job_name}_%j.err\n")
        file.write(f"#SBATCH --output={pwd}/{job_name}_%j.out\n")
        # file.write(f"#SBATCH --error=error.log\n")
        # file.write(f"#SBATCH --output=output.log\n")
        file.write(f"#SBATCH --time={time}\n")
        file.write(f"#SBATCH --partition={partition}\n")
        file.write(f"#SBATCH --no-requeue\n")
        file.write(f"#SBATCH --mem={program_mem}\n\n")

        if molecule.program.lower() == "g16" or file_extension == '.com':
            file.write("mkdir /scratch/\$SLURM_JOB_ID\n\n")
            file.write(f"cd {pwd}\n")
            file.write("export GAUSS_SCRDIR=/scratch/\$SLURM_JOB_ID\n\n")
            file.write(f"srun \$(which g16) {input_file} > {job_name}.log\n")
        elif molecule.program.lower() == "orca" or file_extension == '.inp':
            file.write("source /comm/groupstacks/chemistry/bin/modules.sh\n")
            file.write("ml orca/5.0.4\n")
            file.write("SCRATCH=/scratch/\$SLURM_JOB_ID\n")
            file.write("mkdir -p \$SCRATCH || exit $?\n")
            file.write("cd \$SCRATCH\n")
            file.write(f"cp {pwd}/{job_name}.inp .\n")
            file.write(f"\$(which orca) {job_name}.inp > {pwd}/{job_name}.log\n")
        elif molecule.program.lower() == "crest" or file_extension == '.xyz':
            file.write("source /comm/groupstacks/chemistry/bin/modules.sh\n")
            file.write("ml xtb/6.3.3\n\n")
            file.write("SCRATCH=/scratch/\$SLURM_JOB_ID\n\n")
            file.write("mkdir -p \$SCRATCH || exit $?\n")
            file.write("cd \$SCRATCH\n")
            file.write(f"cp {pwd}/{job_name}.xyz .\n")
            file.write(f"/home/kubeckaj/Applications/crest/crest {job_name}.xyz -gfn2 -ewin 2 -noreftopo -cinp constrain.inp > {job_name}.log\n")
            file.write(f"cp crest_conformers.xyz {pwd}/.\n")
            file.write(f"cp *log {pwd}/.\n")

        file.write("rm -rf /scratch/\$SLURM_JOB_ID\n")
        file.write("!EOF\n\n")
        file.write("sbatch $SUBMIT\n")

    subprocess.run(['sh', path_submit_script])


def submit_array_job(dir, job_files, input_array_list_name, job_name, job_program, partition, time, ncpus, mem, nnodes=1):
    path_submit_script = os.path.join(dir, f"{job_name}_submit.sh")
    array_path = os.path.join(dir, input_array_list_name)
    
    if job_program.lower() == "orca" or job_program.lower() == "g16":
        program_mem = mem + 2000
    else:
        program_mem = mem

    with open(array_path, 'w') as f:
        for file in job_files:
            f.write(f"{file}\n")

    with open(path_submit_script, 'w') as file:
        file.write("#!/bin/bash\n\n")
        file.write("submit=sbatch\n\n")
        file.write(f"IN=$1\n")
        file.write("[ `cut -c1 <<< $IN` == '-' ] && { submit=cat; IN=`cut -c2- <<< $IN`; }\n\n")

        file.write('[ -f "${IN:-##}" ] || { echo "File not found. Good Bye"; exit 1; }\n')
        file.write("PAR=${2:-10}\n")
        file.write("egrep -q '^X[0-9]*$' <<< \"X${PAR}\" || { echo 'Illegal number, $PAR. Good Bye'; exit 1; }\n")

        file.write("MAX=$(wc -l < $IN)\n")
        file.write("[ $PAR -gt $MAX ] && PAR=$MAX\n")
        file.write('ARRAY="1-${MAX}%${PAR}"\n\n')

        file.write('JOB=${IN%.*}\n\n')

        file.write('SUBMIT=qsub.tmp\n')

        file.write("REMDIR=`pwd`\n\n")

        file.write(f"cat <<!EOF  | $submit\n")
        file.write("#!/bin/sh\n")
        file.write(f"#SBATCH --job-name={job_name}\n")
        if job_program.lower() == 'g16':
            file.write(f"#SBATCH --nodes={nnodes}\n")
            file.write(f"#SBATCH --cpus-per-task={ncpus}\n")
            file.write(f"#SBATCH --ntasks={nnodes}\n")
        else:
            file.write(f"#SBATCH --nodes={nnodes}\n")
            file.write(f"#SBATCH --cpus-per-task=1\n")
            file.write(f"#SBATCH --ntasks={ncpus}\n")
        file.write(f"#SBATCH --error={dir}/{job_name}_%j.err\n")
        file.write(f"#SBATCH --output={dir}/{job_name}_%j.out\n")
        # file.write(f"#SBATCH --error=error.log\n")
        # file.write(f"#SBATCH --output=output.log\n")
        file.write(f"#SBATCH --time={time}\n")
        file.write(f"#SBATCH --partition={partition}\n")
        file.write(f"#SBATCH --no-requeue\n")
        file.write(f"#SBATCH --mem={program_mem}\n")
        file.write(f"#SBATCH --array=$ARRAY\n\n")

        file.write("# Create scratch folder\n")
        file.write("SCRATCH=/scratch/\${SLURM_JOB_ID}/\${SLURM_ARRAY_TASK_ID}\n")
        file.write("mkdir -p \$SCRATCH\n\n")

        if job_program.lower() == 'g16':
            file.write(f"cd {dir}\n")
            file.write("export GAUSS_SCRDIR=\$SCRATCH\n\n")

            file.write("source /comm/groupstacks/gaussian/bin/modules.sh\n")
            file.write("ml gaussian16/Rev.B.01\n")
            file.write("ml gcc/9.2.0\n")
            file.write("ml openmpi/4.0.1\n\n")

            file.write('GJ=\$(awk "NR == \$SLURM_ARRAY_TASK_ID" $IN)\n')
            file.write('LOG=\${GJ%.*}.log\n\n')
            file.write("srun $(which g16) \$GJ > \$LOG\n")
            file.write("#\n")
            file.write("!EOF")

        elif job_program.lower() == "orca":
            file.write("  ulimit -c 0\n")
            file.write("source /comm/groupstacks/gaussian/bin/modules.sh\n")
            file.write("ml orca/5.0.4\n\n")

            file.write("cd \$SCRATCH\n")
            file.write(f"cp {dir}/{input_array_list_name} .\n")
            file.write(f"cp {dir}/*.inp .\n\n")

            file.write('GJ=\$(awk "NR == \$SLURM_ARRAY_TASK_ID" $IN)\n')
            file.write('LOG=\${GJ%.*}.log\n\n')

            file.write(f"\$(which orca) \$GJ > \$LOG\n\n")
            file.write("!EOF")

    subprocess.run(['sh', path_submit_script, array_path])


def log2xyz(molecule, atoms=False):
    if molecule.program.lower() == "g16":
        match_string='Standard orientation'
        atomic_number_to_symbol = {
            1: 'H',
            6: 'C',
            7: 'N',
            8: 'O',
            16: 'S'
        }

        with open(molecule.log_file_path, 'r') as file:
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
            if atoms:
                return (element, coordinates)
            else: return coordinates

    elif molecule.program.lower() == "orca":
        match_string='CARTESIAN COORDINATES (ANGSTROEM)'
        with open(molecule.log_file_path, 'r') as file:
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
                if "------" in line:
                    continue
                if "CARTESIAN COORDINATES (A.U.)" in line:
                    start_reading = False
            return coordinates

    elif molecule.program.lower() == "crest": 
        crest_conformer_path = os.path.join(molecule.directory, "crest_conformers.xyz")
        conformers_list = []
        with open(crest_conformer_path, 'r') as file:
            conformer = []
            for line in file:
                stripped_line = line.strip()
                if stripped_line:
                    parts = stripped_line.split()
                    if parts[0].isdigit():
                        if conformer: 
                            conformers_list.append(conformer)
                            conformer = []
                    elif len(parts) == 4: 
                        element, x, y, z = parts
                        conformer.append([element, float(x), float(y), float(z)])
            if conformer:  
                conformers_list.append(conformer)
        return conformers_list


#########################################FILES MANIPULATION############################
class Molecule(VectorManipulation):
    def __init__(self, file_path=None, log_file_path=None, name="", directory="", atoms=None, coordinates=None, reactant=False, program=None, current_step=None, next_step=None):
        self.name = name
        self.directory = directory
        self.file_path = file_path
        self.log_file_path = log_file_path
        self.reactant = reactant
        self.atoms = atoms if atoms is not None else []
        self.coordinates = coordinates if coordinates is not None else []
        self.constrained_indexes = {}
        self.current_step = current_step
        self.next_step = next_step
        self.converged = False
        self.mult = 2
        self.charge = 0
        self.vibrational_frequencies = []
        self.single_point = None
        self.thermal_constribution = None
        self.Q = None
        self._program = None
        self.input = None
        self.output = None
        self.program = program if program is not None else global_program
        self._current_step = None
        self._next_step = None

        
        if self.file_path:
            self.read_xyz_file()

        if self.name == 'OH':
            self.atoms = ['O','H']
            self.coordinates = self.coordinates = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0]]

    def set_name(self, name):
        self.name = name

    def set_directory(self, directory):
        self.directory = directory

    def update_file_path(self, new_path):
        self.file_path = new_path

    def read_xyz_file(self):
        with open(self.file_path, 'r') as file:
            lines = file.readlines()[2:]  
            for line in lines:
                parts = line.split()
                if len(parts) == 4:
                    self.atoms.append(parts[0])
                    self.coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])

    @property
    def program(self):
        return self._program

    @program.setter
    def program(self, value):
        self._program = (value or global_program).lower()
        if self._program == 'orca':
            self.input = '.inp'
            self.output = '.log'
        elif self._program == 'g16':
            self.input = '.com'
            self.output = '.log'
        elif self._program == 'crest':
            self.input = '.xyz'
            self.output = '.log'
        else:
            raise ValueError(f"Unsupported program: {value}")

    @property
    def current_step(self):
        return self._current_step

    @current_step.setter
    def current_step(self, value):
        self._current_step = value
        self._update_next_step()

    @property
    def next_step(self):
        return self._next_step

    @next_step.setter
    def next_step(self, value):
        self._next_step = value
        self._update_current_step()

    def _get_workflow(self):
        if self.reactant:
            if self.name == 'OH':
                return ['optimization', 'DLPNO_SP']
            else:
                return ['crest_sampling', 'opt_conf', 'DLPNO_conf']
        else:
            return ['preopt', 'TS_opt', 'crest_sampling', 'TS_opt_conf', 'DLPNO_conf']

    def _update_next_step(self):
        workflow = self._get_workflow()
        if self._current_step in workflow:
            current_index = workflow.index(self._current_step)
            if current_index + 1 < len(workflow):
                self._next_step = workflow[current_index + 1]
            else:
                self._next_step = None  # End of workflow

    def _update_current_step(self):
        workflow = self._get_workflow()
        if self._next_step in workflow:
            next_index = workflow.index(self._next_step)
            if next_index - 1 >= 0:
                self._current_step = workflow[next_index - 1]
            else:
                self._current_step = None  # Beginning of workflow

    def update_step(self):
        if self._current_step is None and self._next_step is None:
            self._current_step = 'preopt' if not self.reactant else 'crest_sampling'
            self._next_step = self._get_workflow()[1] if len(self._get_workflow()) > 1 else None
        elif self._next_step is not None:
            self.current_step = self._next_step


    def _get_initial_next_step(self):
            workflow = self._get_workflow()
            if self._current_step in workflow:
                current_index = workflow.index(self._current_step)
                if current_index + 1 < len(workflow):
                    return workflow[current_index + 1]
            return None


    def update_energy(self, logger):
        with open(self.log_file_path, 'r') as f:
            log_content = f.read()

            if self.program.lower() == 'g16':
                energy_matches = re.findall(r'(SCF Done:  E\(\S+\) =)\s+([-.\d]+)', log_content)
                if energy_matches:
                    self.single_point = float(energy_matches[-1][-1])
                    freq_matches = re.findall(r"Frequencies --\s+(-?\d+\.\d+)", log_content)
                    if freq_matches:
                        self.vibrational_frequencies = [float(freq) for freq in freq_matches]
                        rot_constant = re.findall(r"Rotational constants \(GHZ\):\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", log_content)
                        self.rot_constant = [float(rot) for rot in rot_constant[-1]]
                        symmetry_num = re.search(r"Rotational symmetry number\s*(\d+)", log_content)
                        if symmetry_num:
                            self.symmetry_num = int(symmetry_num.group(1))
                        else: 
                            logger.log(f"No symmetry number found in {self.name}. assuming 1")
                            self.symmetry_num = 1
                        mol_mass = re.search(r"Molecular mass:\s+(-?\d+\.\d+)", log_content)
                        self.mol_mass = float(mol_mass.group(1))
                        mult = re.search(r"Multiplicity =\s*(\d+)", log_content)
                        if mult:    
                            self.mult = int(mult.group(1))
                        else: self.mult = 2
                        self.partition_function()
                    else:
                        logger.log(f"No frequencies found in {self.name}")

            elif self.program.lower() == 'orca':
                energy_matches = re.findall(r'(FINAL SINGLE POINT ENERGY)\s+([-.\d]+)', log_content)
                if energy_matches:
                    self.single_point = float(energy_matches[-1][-1])

        self.log_items(logger) # DELETE


    def partition_function(self,  T=293.15):
        h = 6.62607015e-34  # Planck constant in J.s
        k_b = 1.380649e-23  # Boltzmann constant in J/K
        c = 299792458       # Speed of light in m/s

        # Vibrational partition function
        qvib = 1 
        for freq in self.vibrational_frequencies:
            if freq > 0:
                f = freq * 100 * c
                qvib *= 1 / (1 - np.exp(-(h * f) / (k_b * T))) 

        # Rotational partition function
        ghztowave = 29.978869
        rot_constants = [rot * ghztowave for rot in self.rot_constant]
        if 0.0 in rot_constants:
            rot_constant = [e for e in rot_constants if e != 0.0]  
            qrot = ((k_b*T)/(self.symmetry_num*h*c*rot_constant[0]))
        else:
            qrot = (1/self.symmetry_num) * ((k_b*T)/(h*c))**1.5 * (np.pi / (rot_constants[0] * rot_constants[1] * rot_constants[2]))**0.5


        # Translational partition function
        mol_mass = self.mol_mass * 1.6605*(10**(-27))
        lmbda = h/(np.sqrt(2*np.pi*mol_mass*k_b*T))
        qtrans = 0.02479/(lmbda**3)

        qelec = self.mult

        self.Q = qvib*qrot*qtrans*qelec

            
    def write_xyz_file(self, output_file_path):
        with open(output_file_path, 'w') as file:
            file.write(str(len(self.atoms)) + '\n\n')
            for atom, coord in zip(self.atoms, self.coordinates):
                coord_str = ' '.join(['{:.6f}'.format(c) for c in coord])  
                file.write(f'{atom} {coord_str}\n')


    def print_items(self):
        print(f"Molecule: {self.name}")
        print(f"File Path: {self.file_path}")
        print(f"Directory: {self.directory}")
        print(f"Log File Path: {self.log_file_path}")
        print(f"Program: {self.program}")
        print(f"Reactant: {self.reactant}")
        print(f"Constrained Indexes: {self.constrained_indexes}")
        print(f"Single Point Energy: {self.single_point}")
        print(f"Partition Function: {self.Q}")
        print(f"Vibrational Frequencies: {self.vibrational_frequencies}")
        print(f"Current Step: {self.current_step}")
        print(f"Next Step: {self.next_step}")
        print("Atoms and Coordinates:")
        for atom, coord in zip(self.atoms, self.coordinates):
            print(f"  {atom}: {coord}")
   
    def log_items(self, logger):
        logger.log(f"Molecule: {self.name}")
        logger.log(f"File Path: {self.file_path}")
        logger.log(f"Directory: {self.directory}")
        logger.log(f"Log File Path: {self.log_file_path}")
        logger.log(f"Reactant: {self.reactant}")
        logger.log(f"Constrained Indexes: {self.constrained_indexes}")
        logger.log(f"Single Point Energy: {self.single_point}")
        logger.log(f"Partition Function: {self.Q}")
        logger.log(f"Vibrational Frequencies: {self.vibrational_frequencies}")
        logger.log(f"Current Step: {self.current_step}")
        logger.log(f"Next Step: {self.next_step}")
        logger.log("Atoms and Coordinates:")
        for atom, coord in zip(self.atoms, self.coordinates):
            logger.log(f"  {atom}: {coord}")

    def H_abstraction(self, distance=1.35, dist_OH=0.97, NEB=False):
        original_coords = self.coordinates.copy()
        atoms = self.atoms
        num_atoms = len(atoms)
        abstraction_molecules = []
        methyl_C_indexes = self.methyl_C_index()  
        carbon_iteration_counter = {index: 0 for index in range(num_atoms) if atoms[index] == 'C'}

        for i in range(num_atoms):
            if atoms[i] == "H":
                for j in range(num_atoms):
                    if atoms[j] == "C":
                        if j in methyl_C_indexes and carbon_iteration_counter[j] >= 1:
                            continue
                        vector_CH = self.calculate_vector(original_coords[j], original_coords[i])
                        dist_CH = self.vector_length(vector_CH)

                        if dist_CH < distance:
                            new_dist_CH = 1.215  # Adjusted C-H distance
                            new_dist_OH = 1.325
                            CHO_angle = 174.845
                            carbon_iteration_counter[j] += 1
                            norm_vector_CH = self.normalize_vector(vector_CH)
                            H_perturb_axis = np.cross(norm_vector_CH, [0, 1, 0])
                            H_perturb_axis = self.normalize_vector(H_perturb_axis)

                            new_coords = original_coords.copy()
                            new_atoms = atoms.copy()

                            new_H_position = np.array(original_coords[i]) + norm_vector_CH * (distance - dist_CH)
                            new_dist_CH = self.vector_length(self.calculate_vector(new_H_position, new_coords[j]))
                            if new_dist_CH < dist_CH:
                                norm_vector_CH = -norm_vector_CH
                                new_H_position = np.array(original_coords[i]) + norm_vector_CH * (distance - dist_CH)
                            oxygen_position = new_H_position + norm_vector_CH * dist_OH
                            norm_vector_OH = self.normalize_vector(oxygen_position - new_H_position)
                            rotation_axis = np.cross(norm_vector_OH, H_perturb_axis)
                            rotation_axis = self.normalize_vector(rotation_axis)

                            rotation_angle_H = np.radians(75.5)
                            new_OH_H_position = oxygen_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * dist_OH
                            # Perturb initial hydrogen to avoid linear C--H--O configuration
                            new_coords[i] = (np.array(new_coords[i]) + H_perturb_axis * 0.3).tolist()

                            new_atoms.append('O')
                            new_atoms.append('H')
                            new_coords.append(oxygen_position.tolist())
                            new_coords.append(new_OH_H_position.tolist())

                            if NEB:
                                # NEB specific code goes here
                                pass

                            new_molecule = Molecule(self.file_path, reactant=self.reactant, program=self.program)
                            new_molecule.atoms = new_atoms
                            new_molecule.coordinates = new_coords
                            new_molecule.constrained_indexes = {'C': j+1, 'H': i+1, 'O': len(new_coords)-1}
                            abstraction_molecules.append(new_molecule)

        return abstraction_molecules 


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

                            perpendicular_axis = np.cross(norm_vector_cc, [0,1,0])
                            perpendicular_axis = self.normalize_vector(perpendicular_axis)
                            # shift carbon
                            new_coords[i][1:] = np.array(new_coords[i][1:]) + norm_vector_cc * 0.1
                            # update oxygen in OH coordiantes
                            oxygen_coords = np.array(new_coords[j][1:]) + perpendicular_axis * distance

                            rotation_axis = norm_vector_cc
                            rotation_angle_h = np.radians(45) 

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


    def methyl_C_index(self, distance=1.5):
            indexes = []
            for i, atom in enumerate(self.atoms):
                if atom == 'C':
                    H_neighbors = sum(1 for j, other_atom in enumerate(self.atoms) if other_atom == 'H' and self.vector_length(self.calculate_vector(self.coordinates[i], self.coordinates[j])) < distance)
                    if H_neighbors >= 3:
                        indexes.append(i)
            return indexes
     
    def get_terminal_O(self, distance=1.5):
        O_index = []
        for i, atom in enumerate(self.coordinates):
            if atom[0] == "O":
                neighbors = sum(1 for j, other_atom in enumerate(self.coordinates) if i != j and self.atom_distance(atom[1:], other_atom[1:]) < distance)
                if neighbors == 1:
                    O_index.append(i)
        return O_index[0]


def get_time(molecule):
    atoms = molecule.atoms
    heavy_count = 0
    H_count = 0
    for atom in atoms:
        if atom == 'H':
            H_count += 1
        else:
            heavy_count += 1
    time_seconds = 0.4 * ((heavy_count**3) + (H_count**2)) + 100
    time_seconds = max(60, min(time_seconds, 3600))
    return time_seconds


def mkdir(molecule, CREST=True):
    if not os.path.exists(molecule.directory):
        os.makedirs(molecule.directory)
    if CREST:
        crest_constrain(molecule)


def pkl_to_xyz(file):
    all_conformers = []
    with open(file, 'rb') as f:
        df = pd.read_pickle(f) 
        for xyz in df[('xyz', 'structure')]:
            coordinates_list = []
            np.set_printoptions(precision=6, suppress=True)
            coords = xyz.get_positions()
            atom = xyz.get_chemical_symbols()
            for atom, coord in zip(atom, coords):
                coordinates = [atom, coord[0], coord[1], coord[2]]
                coordinates_list.append(coordinates)
            all_conformers.append(coordinates_list)
    return all_conformers


#########################################GENERATE INPUT FILES#############################
def crest_constrain(molecule, force_constant=0.95):
    '''Force constant: How tight to constrain the atoms i.e. the magnitude of the oscillation between the constrained atoms'''
    C_index = molecule.constrained_indexes['C']
    H_index = molecule.constrained_indexes['H']
    O_index = molecule.constrained_indexes['O']
    with open (molecule.directory + "/constrain.inp","w") as f:
        f.write("$constrain\n")
        f.write(f"  force constant={force_constant}\n") 
        f.write(f"  distance: {C_index}, {H_index}, auto\n")
        f.write(f"  distance: {H_index}, {O_index}, auto\n")
        f.write(f"  angle: {C_index}, {H_index}, {O_index}, auto\n")
        f.write("$end\n")


def QC_input(molecule, constrain,  method, basis_set, TS, C_index=None, H_index=None, O_index=None):
    file_name = f"{molecule.name}{molecule.input}"
    file_path = os.path.join(molecule.directory, file_name)
    atoms = molecule.atoms
    coords = molecule.coordinates

    if constrain:
        C_index = molecule.constrained_indexes['C']
        H_index = molecule.constrained_indexes['H']
        O_index = molecule.constrained_indexes['O']

    if molecule.program.lower() == "orca":
        with open(file_path, "w") as f:
            if method == "DLPNO" or basis_set == "DLPNO":
                f.write(f"! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF RI-JK aug-cc-pVTZ/JK\n")
                args.cpu = 1 # TEMP solution
            elif TS:
                f.write(f"! {method} {basis_set} OptTS freq\n")
            else:
                f.write(f"! {method} {basis_set} OPT FREQ\n")
            f.write(f"%pal nprocs {args.cpu} end\n")
            f.write(f"%maxcore {args.mem}\n")
            if constrain:
                f.write("%geom\n")
                f.write("Constraints\n")
                # ORCA start atom indexing from 0
                f.write(f"{{B {C_index-1} {H_index-1} C}}\n") # Bond being broken
                f.write(f"{{B {H_index-1} {O_index-1} C}}\n") # Bond being formed
                f.write(f"{{A {C_index-1} {H_index-1} {O_index-1} C}}\n")
                f.write("end\nend\n")
            f.write("\n")
            f.write(f"* xyz {molecule.charge} {molecule.mult}\n") 
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("*")
    elif molecule.program.lower() == "g16":
        with open(file_path, "w") as f:
            f.write(f"%nprocshared={args.cpu}\n")
            f.write(f"%mem={args.mem}mb\n")
            if TS and constrain:
                f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,modredundant) freq\n\n') # freq may be redundant
            elif TS and constrain is False:
                f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen) freq\n\n')
            elif TS is False and constrain:
                f.write(f'# {method} {basis_set} opt=modredundant\n\n')
            else:
                f.write(f"# {method} {basis_set} opt freq\n\n")
            f.write("Title\n\n")
            f.write(f"{molecule.charge} {molecule.mult}\n")
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("\n")
            if constrain:
                f.write(f"B {C_index} {H_index} F\n")
                f.write(f"B {H_index} {O_index} F\n")
                f.write(f"A {C_index} {H_index} {O_index} F\n")
                f.write("\n")
    else:
        print(f"QC_input was called but no program was specified for {molecule.name}")


def NEP_input(file_path, file_name):
    if args.NEB:
        with open(file_path + "/NEB_TS.inp", "w") as f:
            f.write(f'! B3LYP 6-31+g(d,p)  NEB-TS FREQ\n')
            f.write(f'%NEB PREOPT_ENDS TRUE NEB_END_XYZFILE "{file_path + "/" + file_name}_product.xyz" END\n')
            f.write(f'* XYZfile 0 2 {file_path + "/" + file_name}_reactant.xyz\n')


def read_last_lines(filename, num_lines):
    with open(filename, 'rb') as f:
        f.seek(0, os.SEEK_END)
        buffer_size = 8192
        content = ''
        while len(content.splitlines()) < num_lines + 1:
            byte_offset = f.tell() - buffer_size
            if byte_offset < 0: byte_offset = 0
            f.seek(byte_offset)
            content = f.read().decode('utf-8')
            if f.tell() == 0: break  # Reached the beginning of the file
        return content.splitlines()[-num_lines:]


def resubmit_job(molecule, logger, job_type):
    if job_type == 'TS_opt_conf':
        QC_input(molecule, constrain=False,  method=high_method, basis_set=high_basis, TS=True)
        input_file_name = f"{molecule.name}{molecule.input}"
        input_file_path = os.path.join(molecule.directory, input_file_name)
        logger.log(f"Submitting file {input_file_name} for calculation in path {molecule.directory}")
        submit_job(molecule, input_file_path, args.cpu, args.mem, args.par, args.time)
    elif job_type == 'DLPNO_conf':
        molecule.program = 'ORCA'
        QC_input(molecule, constrain=False, method="DLPNO", basis_set="DLPNO", TS=False)
        input_file_name = f"{molecule.name}.inp"
        input_file_path = os.path.join(molecule.directory, input_file_name)
        logger.log(f"Submitting file {input_file_name} for calculation in path {molecule.directory}")
        submit_job(molecule, input_file_path, args.cpu, args.mem, args.par, args.time)


def check_convergence_for_conformers(molecules, logger, threads):
    time_seconds = get_time(molecules[0])
    if args.attempts: max_attempts = args.attempts
    else: max_attempts = 100
    if args.initial_delay: initial_delay = args.initial_delay
    else: initial_delay = int(time_seconds)
    if args.interval: interval = args.interval
    else: interval = initial_delay/2
    attempts = 0
    pending = 0
    freq_cutoff = -250
    job_type = molecules[0].current_step
    
    termination = termination_strings.get(molecules[0].program.lower(), None)
    error_termination = error_strings.get(molecules[0].program.lower(), None)

    all_converged = False
    for m in molecules: # Initiate with the all molecules not being converged
        m.converged = False

    while attempts < max_attempts:
        for molecule in molecules:
            log_file_name = f"{molecule.name}.log"
            log_file_path = os.path.join(molecule.directory, log_file_name)
            molecule.log_file_path = log_file_path
            if molecule.converged:
                continue

            if not os.path.exists(log_file_path):
                all_converged = False
                pending += 1
                time.sleep(interval)
                continue

            last_lines = read_last_lines(log_file_path, 10)
            termination_detected = any(termination in line for line in last_lines)
            error_termination_detected = any(error_termination in line for line in last_lines)

            if termination_detected:
                molecule.converged = True
                logger.log(f"{molecule.name} converged.")
                xyz_coordinates = log2xyz(molecule)
                if xyz_coordinates: 
                    molecule.coordinates = xyz_coordinates
                    molecule.update_energy(logger)
                    if job_type == 'TS_opt_conf':
                        negative_freq = any(freq < freq_cutoff for freq in molecule.vibrational_frequencies)
                        if not negative_freq:
                            logger.log(f"However, no imaginary frequency found for {molecule.name}. Resubmitting TS calculation")
                            molecule.converged = False
                            resubmit_job(molecule, logger, job_type)
                            time.sleep(interval)
                            continue
            elif error_termination_detected:
                molecule.converged = False
                all_converged = False
                logger.log(f"Error termination found in {log_file_name}.")
                xyz_coordinates = log2xyz(molecule)
                if xyz_coordinates:
                    molecule.coordinates = xyz_coordinates
                    logger.log(f"Resubmitting job for {molecule.name} due to error termination")
                    resubmit_job(molecule, logger, job_type)
                    time.sleep(interval)
                    continue
            else:
                print(f"No termination detected in {molecule.name}")
                molecule.converged = False
                all_converged = False

        if all({m.converged for m in molecules}):
            all_converged = True
            break

        if pending >= len(molecules)/3:
            logger.log(f"Some log files have not been generated. Jobs may be pending")
            pending = 0
            time.sleep(initial_delay)
        else:
            attempts += 1
            if attempts % 1 == 0:
                logger.log(f"Conformers not converged. Checking every {interval} seconds. Attempt {attempts}/{max_attempts}")
            time.sleep(interval)

    if all_converged:
        logger.log("All conformer jobs have converged.")
        if job_type in ['TS_opt_conf', 'opt_conf']:
            logger.log("Doing DLPNO single point calculations on top of DFT geometry")
            converged_job(molecules, logger, threads)
        else:
            [global_molecules.append(molecule) for molecule in molecules]


def check_convergence(molecule, logger, threads):
    time_seconds = get_time(molecule)
    log_file_name = f"{molecule.name}.log"
    directory = molecule.directory
    log_file_path = os.path.join(directory, log_file_name)
    molecule.log_file_path = log_file_path
    job_type = molecule.current_step
    molecule.print_items()

    if args.attempts: max_attempts = args.attempts
    else: max_attempts = 100
    if args.initial_delay: initial_delay = args.initial_delay
    else: initial_delay = time_seconds
    if args.interval: interval = args.interval
    else: interval = initial_delay/2

    if job_type in ['preopt', 'crest_sampling'] or molecule.name == 'OH':
        initial_delay = time_seconds//2
        interval = initial_delay//2

    termination = termination_strings.get(molecule.program.lower(), None)
    error_termination = error_strings.get(molecule.program.lower(), None)
    if not termination:
        logger.log("Invalid program specified. Unable to check convergence of log file")
        return False

    initial_mod_time = os.path.getmtime(log_file_path) if os.path.exists(log_file_path) else None
    logger.log(f"Waiting {initial_delay} seconds before first check.")
    time.sleep(initial_delay)
    attempts = 0
    find_log_attempts = 0
    checks = 0

    if args.OH:
        freq_cutoff = -250
    elif args.CC:
        freq_cutoff = -50 # TEST
    else:
        freq_cutoff = -50 # TEST

    while attempts < max_attempts and find_log_attempts < max_attempts:
        try:
            if os.path.exists(log_file_path):
                current_mod_time = os.path.getmtime(log_file_path)
                file_changed = current_mod_time != initial_mod_time
                checks += 1
                if file_changed or checks % 10 == 0:
                    initial_mod_time = current_mod_time
                    last_lines = read_last_lines(log_file_path, 15)
                    termination_detected = any(termination in line for line in last_lines)
                    error_termination_detected = any(error_termination in line for line in last_lines)
                    link1_detected = any("Link1:  Proceeding to internal job step number  2" in line for line in last_lines)
                    if link1_detected:
                        time.sleep(interval)
                        continue
                    if termination_detected:
                        if molecule.name == 'OH' and job_type == 'DLPNO_SP':
                            global_molecules.append(molecule)
                            return True
                        next_step = determine_next_step(molecule)
                        if job_type == 'crest_sampling' or next_step == 'TS_opt_conf':
                            xyz_conformers = log2xyz(molecule)
                            if xyz_conformers:
                                conformer_molecules = []
                                for n, conformer_coords in enumerate(xyz_conformers, start=1):
                                    conformer_molecule = Molecule(name=molecule.name.replace('_CREST', f'_conf{n}'),
                                    directory=molecule.directory,
                                    atoms=[atom[0] for atom in conformer_coords],
                                    coordinates=[atom[1:] for atom in conformer_coords])
                                    conformer_molecule.constrained_indexes = molecule.constrained_indexes
                                    conformer_molecule.reactant = molecule.reactant
                                    conformer_molecule.program = global_program
                                    conformer_molecules.append(conformer_molecule)
                                logger.log_with_stars(f"Yay CREST {log_file_name} converged!")
                                if molecule.reactant:
                                    for m in conformer_molecules: m.mult = 1; m.next_step = 'optimization'
                                    converged_job(conformer_molecules, logger, threads)
                                else:
                                    for m in conformer_molecules: m.next_step = 'TS_opt_conf'
                                    converged_job(conformer_molecules, logger, threads)
                                return True
                            else:
                                logger.log(f"Normal CREST termination in {log_file_name}, but failed to extract xyz coordinates of confomers")
                                non_converged_job(molecule, logger, threads)
                                return False
                        else:
                            xyz_coordinates = log2xyz(molecule)
                            if xyz_coordinates:
                                logger.log(f"Extracting XYZ coordinates and energies from {log_file_name}")
                                molecule.coordinates = xyz_coordinates
                                molecule.update_energy(logger)
                                if job_type == 'TS_opt': #Future: For other types of reactions than H abstraction, the -250 cutoff may not be accurate 
                                    if any(freq < freq_cutoff for freq in molecule.vibrational_frequencies):
                                        negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if freq < freq_cutoff)
                                        logger.log_with_stars(f"Yay {log_file_name} converged with imaginary frequency: {negative_freqs}")
                                        molecule.converged = True
                                        converged_job(molecule, logger, threads)
                                        return True
                                    elif any(freq_cutoff <= freq < 0 for freq in molecule.vibrational_frequencies):
                                        negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if -freq_cutoff <= freq < 0)
                                        logger.log(f"Small negative frequency found in the range {freq_cutoff} to 0: {negative_freqs}. This may be indicate wrong transition state. Resubmitting job")
                                        # Add perturbation of atoms here
                                        non_converged_job(molecule, logger, threads)
                                        return False
                                    else:
                                        logger.log(f"No negative frequency found for transition state: {molecule.vibrational_frequencies}. Resubmitting job")
                                        non_converged_job(molecule, logger, threads)
                                        return False
                                else:
                                    logger.log_with_stars(f"Yay {log_file_name} converged")
                                    molecule.converged = True
                                    converged_job(molecule, logger, threads) 
                                    return True
                            else:
                                logger.log(f"Normal termination of {molecule.program}, but no XYZ coordinates found. Check log file.")
                                return False
                    elif error_termination_detected:
                        if file_changed is False:
                            logger.log("Error termiantion from last log file. Meaning log file has not been updatet, maybe because new job is pending")
                            time.sleep(interval)
                            continue
                        logger.log(f"Error termination in {log_file_name}. Gathering last XYZ coordinates")
                        molecule.log_file_path = log_file_path
                        molecule.update_energy(logger)
                        xyz_coordinates = log2xyz(molecule)
                        if xyz_coordinates:
                            logger.log(f"XYZ coordinates found in failed log file {log_file_name}. Trying to resubmit job")
                            molecule.coordinates = xyz_coordinates
                            non_converged_job(molecule, logger, threads)
                            return False
                        else:
                            logger.log(f"No XYZ coordinates found in {log_file_name}. Check log file for errortype.")
                            # Maybe non_coverged_job here
                            return False
                    else:
                        attempts += 1
                        if attempts % 1 == 0:
                            logger.log(f"No termination yet in {log_file_name}. Checking every {interval} seconds. Attempt: {attempts}/{max_attempts}")
                        time.sleep(interval)
                    if file_changed: 
                        checks = 0
                else:
                    time.sleep(interval)
                    continue
        except FileNotFoundError:
            if find_log_attempts < 10:
                time.sleep(interval)
                logger.log(f"Log file {log_file_name} not found. Perhaps job is pending. Checking again in {interval} seconds.")
                find_log_attempts += 1
            else:
                logger.log(f"Log file {log_file_name} still not found. Check whether calculation is pending or terminated.")
                find_log_attempts += 1
                interval *= 2 

    logger.log("Max attempts reached. Calculation may be stuck. Check for convergence.")
    return False


def submit_and_monitor(molecule, logger, threads):
    if isinstance(molecule, list):
        dir = molecule[0].directory  # Assuming all molecules share the same directory
        job_program = molecule[0].program
        job_files = [os.path.join(dir, f"{conf.name}{conf.input}") for conf in molecule]
        input_array_list_name = "array.txt"
        if molecule[0].current_step == 'TS_opt_conf':
            job_name = f"{molecule[0].name}_conformers"
        elif molecule[0].current_step == 'opt_conf':
            job_name = f'{molecule[0].name}_conformers' 
        else:
            job_name = f"{molecule[0].name}_conformers"

        submit_array_job(dir, job_files, input_array_list_name, job_name, job_program, args.par, args.time, args.cpu, args.mem)
        logger.log(f"Submitted SLURM array job for TS optimization of conformers in {dir}")

        thread = threading.Thread(target=check_convergence_for_conformers, args=(molecule, logger, threads))
        threads.append(thread)
        thread.start()

    elif molecule.current_step not in ['TS_opt_conf', 'DLPNO_conf']:
        input_file_name = f"{molecule.name}{molecule.input}"
        input_file_path = os.path.join(molecule.directory, input_file_name)
        logger.log(f"Submitting file {input_file_name} for calculation in path {molecule.directory}")
        submit_job(molecule, input_file_path, args.cpu, args.mem, args.par, args.time)
        thread = threading.Thread(target=check_convergence, args=(molecule, logger, threads))
        threads.append(thread)
        thread.start()

    else:
        logger.log("Error in molecules list")


def converged_job(molecule, logger, threads):
    if not isinstance(molecule, list):
        molecule.update_step()
        job_type = molecule.current_step

        if job_type == 'TS_opt':
            molecule.name += '_TS'
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.output}")
            molecule.update_file_path(output_file_path)
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=True)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'crest_sampling':
            molecule.name = molecule.name.replace("_TS", "_CREST")
            molecule.program = 'CREST'
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}.xyz")
            molecule.write_xyz_file(output_file_path)
            submit_and_monitor(molecule, logger, threads)
    
        elif job_type == 'DLPNO_SP':
            molecule.program = 'ORCA'
            QC_input(molecule, constrain=False, method='DLPNO', basis_set='DLPNO', TS=False)
            submit_and_monitor(molecule, logger, threads)

    elif isinstance(molecule, list):
        conformer_molecules = molecule
        job_type = conformer_molecules[0].next_step
        for conf in conformer_molecules:
            conf.program = global_program
            conf.update_step()
            if job_type == 'TS_opt_conf':
                conf.name += '_TS'
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=True)

            elif job_type == 'DLPNO_conf':
                conf.program = 'ORCA'
                conf.name = conf.name.replace("_TS","_DLPNO")
                QC_input(conf, constrain=False, method='DLPNO', basis_set='DLPNO', TS=False)

            elif job_type == 'opt_conf':
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=False)
        submit_and_monitor(conformer_molecules, logger, threads)

    else:
        logger.log(f"Error: next step could not be determined for molecule {molecule.name}")


def non_converged_job(molecule, logger, threads):
    job_type = molecule.current_step
    logger.log(f"Handling non-converged job {job_type} for molecule: {molecule.name}")

    if job_type == "preopt":
        logger.log(f"Failed preoptimization for molecule: {molecule.name}. Redoing calculation")
        QC_input(molecule, constrain=True, method=low_method, basis_set=low_basis, TS=False)
        submit_and_monitor(molecule, logger, threads)

    elif job_type == "TS_opt":
        logger.log(f"Failed optimization towards first order saddle point for molecule: {molecule.name}. Redoing calculation")
        QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=True)
        submit_and_monitor(molecule, logger, threads)

    elif job_type == "crest_sampling":
        logger.log(f"Failed CREST sampling for molecule: {molecule.name}. Redoing calculation")
        output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.input}")
        molecule.write_xyz_file(output_file_path)
        submit_and_monitor(molecule, logger, threads)

    elif job_type == "TS_opt_conf":
        logger.log(f"Failed CREST sampling for molecule: {molecule.name}. Redoing calculation")
        if not isinstance(molecule, list):
            conformer_molecules = [molecule]
        else: conformer_molecules = molecule
        for conf in conformer_molecules:
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=True)
        submit_and_monitor(conformer_molecules, logger, threads)

    elif job_type == 'DLPNO_conf':
        if not isinstance(molecule, list):
            conformer_molecules = [molecule]
        else: conformer_molecules = molecule
        for conf in conformer_molecules:
            conf.program = 'ORCA'
            QC_input(conf, constrain=False, method='DLPNO', basis_set='DLPNO', TS=False)
        submit_and_monitor(conformer_molecules, logger, threads)

    else:
        logger.log(f"Error: next step could not be determined for non-converged job of molecule: {molecule.name}")


def determine_next_step(molecule):
    # For reactants
    if molecule.name == "OH": # Should also add that OH needs freqs
        return 'DLPNO_SP'
    if 'reactant_CREST' in molecule.name:
        return 'opt_conf'
    elif 'reactant_opt' in molecule.name:
        return 'DLPNO_conf'
    # For TS molecules
    try:
        with open(molecule.log_file_path, 'r') as file:
            for _ in range(30, 200):  # Check for method section between line 30 and 200. 
                try:
                    line = next(file).strip()
                    if line.startswith('#'):
                        components = line.split(",")
                        if 'ts' not in components:
                            return 'TS_opt'
                        elif 'ts' in components: 
                            if 'conf' in molecule.name:
                                vibrations = log2vib(molecule)
                                if vibrations:
                                    return 'DLPNO_conf'
                            else: 
                                return 'crest_sampling' 
                    elif line.startswith('xTB'):
                        return 'TS_opt_conf'
                except StopIteration:
                    break  
        return 'Method section not found'
    except FileNotFoundError:
        print(f"Error: Log file for molecule {molecule.name} not found.")
        return 'error'



def log2vib(molecule):
    with open(molecule.log_file_path, 'r') as file:
        content = file.read()
        if molecule.program.lower() == "g16":
            vibrations = re.findall(r"Frequencies --\s+(-?\d+\.\d+)", content)
        elif molecule.program.lower() == "orca":
            vibrations = []
            vib = re.search(r'[-+]?\d*\.\d+\s*cm\*\*-1', content)
            if vib:
                vibration = float(vib.group().split()[0])
                vibrations.append(vibration)
        else:
            return 'No vibrations found'
    return vibrations
        
def log2program(log_file_path):
    try:
        with open(log_file_path, 'r') as file:
            for _ in range(10):
                line = file.readline()
                if "Entering Gaussian System" in line or "Gaussian(R)" in line:
                    return 'g16'
                elif '* O   R   C   A *' in line:
                    return 'orca'
                elif '|                 C R E S T                  |' in line:
                    return 'crest'
    except FileNotFoundError:
        print(f"File not found: {log_file_path}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    return None
    

def main():
    parser = argparse.ArgumentParser(description='''    'Dynamic Approach for Transition State'
    Automated tool for generating input files, primarily 
    for transition state geometry optimization. 
    Calculation of tunneling corrected multi-configurational 
    rate constants can also be calculated from log files.''',
                                     prog="JKTS",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
    Examples of use:
                JKTS CH4.xyz -H -auto --low_method "am1 3-21g" --high_method "B3LYP 6-31+g(d,p)"
                JKTS reactant.log product.log -k
                JKTS pinonaldehyde.xyz -H -auto
                                     ''')


    parser.add_argument('input_files', metavar='reactant.xyz', nargs='+', help='XYZ files (e.g., pinonaldehyde.xyz or even *.xyz)')

    reaction_options = parser.add_argument_group("Types of reactions")

    reaction_options.add_argument('-CC', action='store_true', help='Perform addition to C=C bonds')
    reaction_options.add_argument('-H', action='store_true', help='Perform H abstraction with OH radical')
    reaction_options.add_argument('-OH', action='store_true', help='Perform OH addition to C=C bonds')

    parser.add_argument('-G16', action='store_true', help='Gaussian16 is used for QC calculations')
    parser.add_argument('-ORCA', action='store_true', help='ORCA is used for QC calculations')
    parser.add_argument('-crest', action='store_true', help='If CREST sampling should be performed on input file(.xyz or .log)')
    parser.add_argument('-constrain', action='store_true', help='Constrain is integrated into relevant input file')
    parser.add_argument('-reactants', action='store_true', help='Prepare folder for reactants')
    parser.add_argument('-NEB', action='store_true', help='Prepare input file for Nudged Elsatic Band')
    parser.add_argument('-auto', action='store_true', help='Automated process with the workflow: \n-> Preoptimization of geometry (Low Level of Theory) \n-> Optimization towards transition state (High Level of Theory) \n-> CREST TS conformer sampling (GFN2-xTB -ewin=2kcal/mol) \n-> Reoptimization of CREST conformers to TS (High Level of Theory)\n-> DLPNO-CCSD(T) SP energy calculations on top of TS conformers \n-> calculate rate constants and branching rations for reaction type \n* Resubmission of failed calculations is automatically done until convergence or wall-time reached')

    additional_options = parser.add_argument_group("Additional arguments")

    additional_options.add_argument('--no-TS', action='store_true', help='Input files for geometry relaxation are generated')
    additional_options.add_argument('--no-xyz', action='store_true', help='No XYZ files generated')
    additional_options.add_argument('-k', action='store_true', help='Calculate Multiconformer Transition State rate constant')
    additional_options.add_argument('--high_level', nargs='+', metavar='', help='Specify the high level of theory for QC method TS optimization [def = wB97X-D aug-cc-pVTZ]')
    additional_options.add_argument('--low_level', nargs='+', metavar='', help='Specify the low level of theory for preoptimization [def = B3LYP 6-31+G(d,p)]')
    additional_options.add_argument('-cpu', metavar="int", nargs='?', const=1, type=int, default=4, help='Number of CPUs [def = 4]')
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=4000, help='Amount of memory allocated for job [def = 400mb]')
    additional_options.add_argument('-par', metavar="partition", nargs='?', const=1, default="qany", help='Partition to use [def = qany]')
    additional_options.add_argument('-time', metavar="hours:minutes:seconds", nargs='?', const=1, default="72:00:00", help='Specify total time for calculation [def = 72 Hours]')
    additional_options.add_argument('-interval', metavar="int", nargs='?', const=1, type=int, help='Set time interval between checks of log files')
    additional_options.add_argument('-initial_delay', metavar="int", nargs='?', const=1, type=int, help='Set an initial delay before checking log files')
    additional_options.add_argument('-attempts', metavar="int", nargs='?', const=1, type=int, help='Set how mnay times a log files should be checked')

    global args, start_dir
    args = parser.parse_args()
    start_dir = os.getcwd()

    ################################ARGUMENT SPECIFICATIONS############################
    if args.auto:
        args.constrain=True; args.no_TS=True; args.reactants=True; args.no_xyz=True;
    # Program and method 
    global low_method, low_basis, high_method, high_basis

    methods_no_basis = {"B97-3c", "r2scan-3c", "pm3", "am1", "pm6", "pm7"}

    def extract_method_basis(input_args, default):
        if not input_args:
            return default

        if len(input_args) == 1:
            if input_args[0] in methods_no_basis:
                return input_args[0], ""
            else:
                raise ValueError(f"Basis set required for method {input_args[0]}")
        
        return input_args[0], input_args[1]

    high_method, high_basis = extract_method_basis(args.high_level, ["wb97xd", "aug-cc-pVTZ"])
    low_method, low_basis = extract_method_basis(args.low_level, ["B3LYP", "6-31+G(d,p)"])

    global global_program
    if args.G16:
        global_program = "G16"
    elif args.ORCA:
        global_program = "ORCA"
        if low_method in methods_no_basis:
            low_basis = ""
        if high_method in methods_no_basis:
            high_basis = ""
        if high_method.lower() == "wb97xd":
            high_method = "WB97X-D3"
        if low_method.lower() == "wb97xd":
            low_method = "WB97X-D3" 
    else: 
        global_program = None

    global termination_strings, error_strings
    termination_strings = {
        "g16": "Normal termination",
        "orca": "****ORCA TERMINATED NORMALLY****",
        "crest": "CREST terminated normally"
    }
    error_strings = {
        "g16": "Error termination",
        "orca": "The optimization did not converge but reached the maximum",
        "crest": "Find my error message"
    }

    ####################################################################################################
    global global_molecules
    global_molecules = []
    threads = []

    for n, input_file in enumerate(args.input_files):
        file_name, file_type = os.path.splitext(input_file)
        input_file_path = os.path.join(start_dir, input_file)

        if file_type == ".xyz":
            input_molecule = Molecule(input_file_path)

            if args.reactants: # and args.OH:
                OH = Molecule(name="OH", reactant=True)
                OH.set_directory(os.path.join(start_dir, "reactants"))
                OH.mult = 2
                OH.next_step = 'optimization'
                reactant = Molecule(input_file_path, name=f"{file_name}_reactant_CREST", reactant=True, program='CREST')
                reactant.set_directory(os.path.join(start_dir, "reactants"))
                reactant.mult = 1
                reactant.current_step = 'crest_sampling' 
                reactant_logger = Logger(os.path.join(reactant.directory, "log"))
                mkdir(reactant, CREST=False)

                QC_input(OH, constrain=False, method=high_method, basis_set=high_basis, TS=False)
                reactant.write_xyz_file(os.path.join(reactant.directory, f"{reactant.name}.xyz"))
                submit_and_monitor(reactant, reactant_logger, threads)
                submit_and_monitor(OH, reactant_logger, threads)

            if args.OH:
                # Hydrogen abstraction process
                abstraction_molecules = input_molecule.H_abstraction(NEB=args.NEB)

                for count, molecule in enumerate(abstraction_molecules, start=1):
                    molecule.set_name(f"{file_name}_H{count}")
                    molecule.set_directory(os.path.join(start_dir, molecule.name))
                    molecule.program = global_program
                    molecule.current_step = 'preopt'
                    molecule.reactant = False
                    mkdir(molecule)
                    logger = Logger(os.path.join(molecule.directory, "log"))

                    if args.no_xyz is False:
                        molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz"))

                    if molecule.program is not None:
                        QC_input(molecule=molecule, constrain=args.constrain, method=low_method, basis_set=low_basis, TS=False)
                        submit_and_monitor(molecule, logger, threads)


            # elif args.crest:
            #     executing_path = os.getcwd()
            #     logger = Logger(os.path.join(executing_path, "log"))
            #     input_file_path = os.path.join(executing_path, input_file)
            #     submit_and_monitor(executing_path, input_file_path, logger, convergence_info, threads, job_type='crest', job_program=program)
            #
            # elif args.CC and len(args.input_files) == 2:
            #     # C=C addition
            #     addition(args.input_files[n], args.input_files[n+1], constrain=args.constrain, program=program, method=low_method, TS=False, no_xyz=args.no_xyz, basis_set=low_basis)
            #
            # elif not args.H and not args.CC:
            #     # Default case for XYZ files without specific reactions
            #     coords = read_xyz_file(input_file)
            #     QC_input(file_name=input_file.split(".")[0], destination=os.getcwd(), coords=coords, constrain=args.constrain, program=program, method=args.method, basis_set=args.basis, TS=TS)
            #
            else:
                print("Unsupported reaction type or combination of arguments.")

        elif file_type in [".log", ".out"]:
            executing_path = os.getcwd()
            log_file_path = os.path.join(executing_path, input_file)
            logger = Logger(os.path.join(executing_path, "log"))
            log_program = log2program(log_file_path)
            molecule = Molecule(name=file_name, directory=executing_path, program=log_program)
            if 'reactant' in molecule.name: molecule.reactant = True; molecule.mult = 1
            molecule.log_file_path = log_file_path
            molecule.program = log_program
            next_step = determine_next_step(molecule)
            last_lines = read_last_lines(input_file_path, 10) 
            termination = termination_strings.get(log_program.lower(), "")
            error_termination = error_strings.get(log_program.lower(), "")

            if any(termination in line for line in last_lines):
                logger.log(f"Detected log file with normal termination. Next step is: {next_step}")
                if log_program in ['orca', 'g16']:
                    atoms, xyz_coordinates = log2xyz(molecule, True)
                    molecule.atoms, molecule.coordinates = atoms, xyz_coordinates
                    molecule.next_step = next_step
                    molecule.update_energy(logger)
                    converged_job(molecule, logger, threads)
                elif log_program == 'crest':
                    xyz_conformers = log2xyz(molecule)
                    if xyz_conformers:
                        conformer_molecules = []
                        for n, conformer_coords in enumerate(xyz_conformers, start=1):
                            conformer_molecule = Molecule(name=molecule.name.replace('_CREST', f'_conf{n}'),
                            directory=molecule.directory,
                            atoms=[atom[0] for atom in conformer_coords],
                            coordinates=[atom[1:] for atom in conformer_coords])
                            conformer_molecule.constrained_indexes = molecule.constrained_indexes
                            conformer_molecule.reactant = molecule.reactant
                            conformer_molecule.next_step = 'TS_opt_conf' 
                            conformer_molecule.program = log_program
                            conformer_molecules.append(conformer_molecule)
                        if molecule.reactant:
                            for m in conformer_molecules: m.mult = 1
                            converged_job(conformer_molecules, logger, threads)
                        else:
                            converged_job(conformer_molecules, logger, threads)
                
            elif any(error_termination in line for line in last_lines):
                logger.log("Detected log file with error termination.")
                # non_converged_job(molecule, logger, next_step, threads, next_step)

        else:
            parser.error(f"Unsupported file type for input: {input_file}")

    # Monitor and handle convergence of submitted jobs
    while threads:
        for thread in list(threads):
            thread.join(timeout=0.1)
            if not thread.is_alive():
                threads.remove(thread)

    print("Final molecules:")
    for molecule in global_molecules:
        molecule.print_items()

    reactant_molecules = [molecule for molecule in global_molecules if molecule.reactant]
    TS_molecules = [molecule for molecule in global_molecules if molecule.reactant is False]

    if reactant_molecules and TS_molecules:
        lowest_SP_reactant = min(reactant_molecules, key=lambda molecule: molecule.single_point)
        lowest_SP_TS = min(TS_molecules, key=lambda molecule: molecule.single_point)

        HtoJ = 4.3597447222*(10**(-18))
        sum_TS = sum([np.exp((lowest_SP_TS*HtoJ - molecule.single_point*HtoJ) / (k_b*T) * molecule.Q) for molecule in TS_molecules])
        sum_reactant = sum([np.exp((lowest_SP_reactant*HtoJ - molecule.single_point*HtoJ) / (k_b*T) * molecule.Q) for molecule in reactant_molecules])

        print(f"sum_TS: {sum_TS}")
        print(f"sum_reactant: {sum_reactant}")

        k = [(k_b*T)/h * (sum_TS * np.exp(-(molecule_i.single_point - lowest_SP_TS)/k_b*T)*molecule_i.Q) / (sum_reactant*(np.exp(-(molecule_j.single_point - lowest_SP_reactant)/k_b*T)*molecule_j.Q)) * np.exp(-(lowest_SP_TS - lowest_SP_reactant)/k_b*T) for molecule_i, molecule_j in zip(TS_molecules, reactant_molecules)]
        result_logger = Logger(start_dir)
        result_logger.log(k)


if __name__ == "__main__":
    main()
