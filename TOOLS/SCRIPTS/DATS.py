#!/usr/bin/env python3
'''Dynamic Approach for Transition State'''
###############################LIBRARIES#####################################
from os.path import exists
import numpy as np
import pandas as pd
import argparse
import os
import re
import subprocess
import time
import threading
import copy
import shutil
import glob

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
def extract_job_id(output):
    match = re.search(r'Submitted batch job (\d+)', output)
    if match:
        return int(match.group(1))
    else:
        print("Could not extract job ID.")
        return None

def check_job_status(job_id):
    try:
        result = subprocess.run(['squeue', '-j', str(job_id)], capture_output=True, text=True, check=True)
        status = parse_job_status(result.stdout)
        return status
    except subprocess.CalledProcessError as e:
        print(f"Error checking job status: {e}")
        return "Unknown"

def parse_job_status(output):
    lines = output.strip().split("\n")
    # Check if there's more than just the header line
    for line in lines[1:]:  # Skip the header line
        parts = line.split()
        if len(parts) > 4:
            status_code = parts[4]  # Status code is expected to be the fifth element
            if status_code == "PD":
                return "Pending"
            elif status_code == "R":
                return "Running"
            elif status_code == 'CP':
                return 'Completed'
    return "Completed or Not Found"


def submit_job(molecule, ncpus, mem, partition, time, nnodes=1):
    input_file_name = f"{molecule.name}{molecule.input}"
    input_file_path = os.path.join(molecule.directory, input_file_name)
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
        file.write(f"#SBATCH --error={pwd}/out_and_err_files/{job_name}_%j.err\n")
        file.write(f"#SBATCH --output={pwd}/out_and_err_files/{job_name}_%j.out\n")
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
            file.write(f"srun \$(which g16) {input_file_path} > {job_name}.log\n")

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
            file.write(f"cp {pwd}/constrain.inp .\n")
            if molecule.reactant and 'OH' in molecule.name:
                file.write(f"/home/kubeckaj/Applications/crest/crest {job_name}.xyz -gfn2 -ewin 2 -noreftopo -uhf 1 > {job_name}.log\n")
            elif molecule.reactant and 'OH' not in molecule.name:
                file.write(f"/home/kubeckaj/Applications/crest/crest {job_name}.xyz -gfn2 -ewin 2 -noreftopo > {job_name}.log\n")
            else:
                file.write(f"/home/kubeckaj/Applications/crest/crest {job_name}.xyz -gfn2 -ewin 2 -noreftopo -uhf 1 -cinp constrain.inp > {job_name}.log\n")
            file.write(f"cp crest_conformers.xyz {pwd}/.\n")
            file.write(f"cp *log {pwd}/.\n")

        file.write("rm -rf /scratch/\$SLURM_JOB_ID\n")
        file.write("!EOF\n\n")
        file.write("sbatch $SUBMIT\n")

    submit_command = ['sh', path_submit_script]
    try:
        result = subprocess.run(submit_command, capture_output=True, text=True, check=True)
        job_id = extract_job_id(result.stdout)
        molecule.job_id = job_id
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"Error in job submission: {e}")
        return None



def submit_array_job(molecules, partition, time, ncpus, mem, nnodes=1):
    dir = molecules[0].directory
    job_files = [os.path.join(dir, f"{conf.name}{conf.input}") for conf in molecules]
    job_name = molecules[0].current_step
    job_program = molecules[0].program
    path_submit_script = os.path.join(dir, f"{job_name}_submit.sh")
    input_array_list_name = "array.txt"
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
        file.write("PAR=${2:-50}\n")
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
        file.write(f"#SBATCH --error={dir}/out_and_err_files/{job_name}_%j.err\n")
        file.write(f"#SBATCH --output={dir}/out_and_err_files/{job_name}_%j.out\n")
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

    submit_command = ['sh', path_submit_script, array_path]
    try:
        result = subprocess.run(submit_command, capture_output=True, text=True, check=True)
        job_id = extract_job_id(result.stdout)
        for n, molecule in enumerate(molecules, start=1): 
            molecule.job_id = f"{job_id}_{n}"
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"Error in job submission: {e}")
        return None


def log2xyz(molecule, atoms=False, max_conformers=50):
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
            if atoms:
                return (element, coordinates)
            else: return coordinates

    if molecule.program.lower() == "crest": 
        if args.max_conformers:
            max_conformers = args.max_conformers
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
        return conformers_list[:max_conformers]

#########################################FILES MANIPULATION############################
class Molecule(VectorManipulation):
    def __init__(self, file_path=None, log_file_path=None, name="", directory="", atoms=None, coordinates=None, reactant=False, program=None, current_step=None, next_step=None):
        self.name = name
        self.directory = directory
        self.file_path = file_path
        self.log_file_path = log_file_path
        self.reactant = reactant
        self.job_id = None
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
                return ['crest_sampling', 'opt_conf_reactants', 'DLPNO_conf']
        else:
            return ['crest_sampling', 'opt_conf', 'TS_opt_conf', 'DLPNO_conf']

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


    def update_energy(self, logger, log_file_path=None, DLPNO=False, program=None):
        file_path = log_file_path if log_file_path else self.log_file_path
        program = program if program else self.program

        with open(file_path, 'r') as f:
            log_content = f.read()

            if program.lower() == 'g16':
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
                    elif 'TS' in self.name:
                        logger.log(f"No frequencies found in {self.name}")

            elif program.lower() == 'orca' or DLPNO:
                energy_matches = re.findall(r'(FINAL SINGLE POINT ENERGY)\s+([-.\d]+)', log_content)
                if energy_matches:
                    self.single_point = float(energy_matches[-1][-1])
                    freq_matches = re.findall(r'[-+]?\d*\.\d+\s*cm\*\*-1', log_content)
                    if freq_matches and DLPNO is False:
                        self.vibrational_frequencies = [float(freq) for freq in freq_matches]
                        rot_constant = re.findall(r"Rotational constants in cm-1: \s*[-+]?\d*\.\d*  \s*[-+]?\d*\.\d* \s*[-+]?\d*\.\d*", log_content)
                        self.rot_constant = [float(rot) for rot in rot_constant[-1]]
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


        # self.log_items(logger) # DELETE


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
        print("-----------------------------------------------")
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
        print("-----------------------------------------------")
   
    def log_items(self, logger):
        logger.log("-----------------------------------------------")
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
        logger.log("-----------------------------------------------")

    def H_abstraction(self, NEB=False):
        original_coords = self.coordinates.copy()
        atoms = self.atoms
        num_atoms = len(atoms)
        abstraction_molecules = []
        methyl_C_indexes = self.methyl_C_index()  
        carbon_iteration_counter = {index: 0 for index in range(num_atoms) if atoms[index] == 'C'}
        distance_CH = 1.215  # Adjusted C-H distance
        distance_OH = 1.325
        distance_OH_H = 0.97
        CHO_angle = 174.845

        for i in range(num_atoms):
            if atoms[i] == "H":
                for j in range(num_atoms):
                    if atoms[j] == "C":
                        if j in methyl_C_indexes and carbon_iteration_counter[j] >= 1:
                            continue
                        vector_CH = self.calculate_vector(original_coords[j], original_coords[i])
                        dist_CH = self.vector_length(vector_CH)

                        if dist_CH <= distance_OH:
                            carbon_iteration_counter[j] += 1
                            norm_vector_CH = self.normalize_vector(vector_CH)
                            H_perturb_axis = np.cross(norm_vector_CH, [0, 1, 0])
                            H_perturb_axis = self.normalize_vector(H_perturb_axis)

                            new_coords = original_coords.copy()
                            new_atoms = atoms.copy()

                            new_H_position = np.array(original_coords[i]) + norm_vector_CH * (dist_CH - distance_CH)
                            distance_CH = self.vector_length(self.calculate_vector(new_H_position, new_coords[j]))
                            if distance_CH < dist_CH:
                                norm_vector_CH = -norm_vector_CH
                                new_H_position = np.array(original_coords[i]) + norm_vector_CH * (dist_CH - distance_CH)
                            oxygen_position = new_H_position - norm_vector_CH * distance_OH 
                            norm_vector_OH = self.normalize_vector(oxygen_position - new_H_position)
                            rotation_axis = np.cross(norm_vector_OH, H_perturb_axis)
                            rotation_axis = self.normalize_vector(rotation_axis)

                            rotation_angle_H = np.radians(104.5)
                            new_OH_H_position = oxygen_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H
                            # Perturb initial hydrogen to avoid linear C--H--O configuration
                            new_coords[i] = (new_H_position + H_perturb_axis * 0.07).tolist()

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
        if not os.path.exists(os.path.join(molecule.directory, "input_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "input_files")))
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

def extract_normal_coordinates(molecule):
    '''Return format: [[0.02, -0.02, 0.06], [0.07, -0.11, -0.02], ...]''' 
    with open(molecule.log_file_path, 'r') as file:
        lines = file.readlines()
    negative_freq_found = False
    xyz_coordinates = []
    for line in lines:
        if "Frequencies --" in line and '-' in line:
            negative_freq_found = True
            continue
        if any(substring in line for substring in ['Red. masses', 'Frc consts', 'IR Inten', 'Atom  AN']):
            continue
        if negative_freq_found:
            if re.match(r'\s*\d+\s+\d+\s+(-?\d+\.\d+\s+)*', line):
                xyz = re.findall(r'-?\d+\.\d+', line)[:3]  # Extract first three floats
                xyz = [float(coord) for coord in xyz]
                xyz_coordinates.append(xyz)
            else:
                break

    return xyz_coordinates


def perturb_atoms(molecule):
    C_index = molecule.constrained_indexes['C']
    H_index = molecule.constrained_indexes['H']
    O_index = molecule.constrained_indexes['O']

    original_coordinates = molecule.coordinates
    perturbation_directions = extract_normal_coordinates(molecule)
    perturbed_coords = []
    for i,j in zip(original_coordinates, perturbation_directions):
        coords = list(map(lambda x,y: x+y, i, j))
        perturbed_coords.append(coords)

    molecule.coordinates = perturbed_coords

#########################################GENERATE INPUT FILES#############################
def crest_constrain(molecule, force_constant=1.00):
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
            if method == "DLPNO":
                f.write(f"%maxcore {args.mem+3000}\n") # Find some way to make sure enough memory
            else:
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


def move_files(directory, files_type):
    if files_type in ['.com', '.inp']:
        folder = "input_files"
        file_type = files_type 
    elif files_type in ['.out', '.log']:
        folder = "log_files"
        file_type = files_type
    else:
        folder = "temp" 
        file_type = files_type
    destination = os.path.join(directory, folder)
    os.makedirs(destination, exist_ok=True)

    for file in glob.glob(os.path.join(directory, f"*{file_type}")):
        shutil.move(file, destination)


def read_last_lines(filename, logger, num_lines):
    attempts = 0
    max_attempts = 10

    while attempts < max_attempts:
        try:
            with open(filename, 'rb') as f:
                f.seek(0, os.SEEK_END)
                buffer_size = 8192
                content = ''
                while len(content.splitlines()) < num_lines + 1:
                    byte_offset = f.tell() - buffer_size
                    if byte_offset < 0: 
                        byte_offset = 0
                    f.seek(byte_offset)
                    content = f.read().decode('utf-8')
                    if f.tell() == 0: 
                        break  # Reached the beginning of the file
                return content.splitlines()[-num_lines:]
        except FileNotFoundError:
            attempts += 1
            # logger.log(f"Attempt {attempts}: File {filename} not found. Sleeping for 120 seconds before retrying.")
            time.sleep(120)

    logger.log(f"Failed to find the file {filename} after {max_attempts} attempts.")
    return ""

def resubmit_job(molecule, logger, job_type):
    if job_type == 'opt_conf':
        QC_input(molecule, constrain=True, method=low_method, basis_set=low_basis, TS=False)

    elif job_type == 'TS_opt_conf':
        QC_input(molecule, constrain=False,  method=high_method, basis_set=high_basis, TS=True)

    elif job_type == 'DLPNO_conf':
        molecule.program = 'ORCA'
        QC_input(molecule, constrain=False, method="DLPNO", basis_set="DLPNO", TS=False)

    elif job_type == 'opt_conf_reactants':
        QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=False)

    job_id = submit_job(molecule, args.cpu, args.mem, args.par, args.time)
    molecule.job_id = job_id
    logger.log(f"submitted file {molecule.name}{molecule.input} with job type: {job_type} and job id {molecule.job_id} in path {molecule.directory}")


def check_convergence_for_conformers(molecules, logger, threads):
    time_seconds = get_time(molecules[0])
    if args.attempts: max_attempts = args.attempts
    else: max_attempts = 100
    if args.initial_delay: initial_delay = args.initial_delay
    else: initial_delay = int(time_seconds)
    if args.interval: interval = args.interval
    else: interval = initial_delay/2
    attempts = 0
    count = 0
    pending = []
    running = [] 
    job_type = molecules[0].current_step

    freq_cutoff = -250
    
    termination = termination_strings.get(molecules[0].program.lower(), "")
    error_termination = error_strings.get(molecules[0].program.lower(), "")

    all_converged = False
    for m in molecules: # Initiate with the all molecules not being converged
        m.converged = False

    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)
    while attempts < max_attempts and all_converged is False:
        count += 1
        for molecule in molecules:
            if molecule.converged and count % 10 != 0:
                continue
            status = check_job_status(molecule.job_id) 
            if status in ['Running', 'Completed or Not Found'] or molecule.converged is False:
                if molecule.job_id not in running:
                    logger.log(f"Job {job_type} for {molecule.name} with job id {molecule.job_id} is running.")
                    input_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.input}")
                    destination = os.path.join(molecule.directory, "input_files")
                    if os.path.exists(destination): os.remove(destination)
                    shutil.move(input_file_path, destination)
                    running.append(molecule.job_id)

                log_file_name = f"{molecule.name}.log"
                log_file_path = os.path.join(molecule.directory, log_file_name)
                molecule.log_file_path = log_file_path
                last_lines = read_last_lines(log_file_path, logger, 10)
                termination_detected = any(termination in line for line in last_lines)
                error_termination_detected = any(error_termination in line for line in last_lines)

                if termination_detected:
                    xyz_coordinates = log2xyz(molecule)
                    molecule.coordinates = xyz_coordinates
                    molecule.update_energy(logger)
                    if job_type == 'TS_opt_conf':
                        if any(freq < freq_cutoff for freq in molecule.vibrational_frequencies):
                            negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if freq < freq_cutoff)
                            if molecule.converged == False:
                                logger.log_with_stars(f"Yay {log_file_name} converged with imaginary frequency: {negative_freqs}")
                            molecule.converged = True
                            converged_job(molecule, logger, threads)
                        elif any(freq_cutoff <= freq < 0 for freq in molecule.vibrational_frequencies):
                            negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if -freq_cutoff <= freq < 0)
                            logger.log(f"Small negative frequency found in the range {freq_cutoff} to 0: {negative_freqs}. This may be indicate wrong transition state. Resubmitting job")
                            molecule.converged = False
                            # Add perturbation of atoms here
                            # Perhaps constraining methyl groups may work
                            resubmit_job(molecule, logger, job_type)
                            time.sleep(interval)
                            continue
                        else:
                            logger.log(f"However, no imaginary frequency found for {molecule.name}. Resubmitting TS calculation")
                            molecule.converged = False
                            # perturb atoms
                            resubmit_job(molecule, logger, job_type)
                            time.sleep(interval)
                            continue
                    else:
                        molecule.converged = True
                        logger.log(f"{molecule.name} converged.")
                        xyz_coordinates = log2xyz(molecule)
                        if xyz_coordinates: 
                            molecule.coordinates = xyz_coordinates
                        else:
                            logger.log(f"Error gathering XYZ cordinates from {log_file_name}. Check log file.")
                elif error_termination_detected:
                    molecule.converged = False
                    all_converged = False
                    logger.log(f"Error termination found in {log_file_name}.")
                    xyz_coordinates = log2xyz(molecule)
                    if xyz_coordinates:
                        molecule.coordinates = xyz_coordinates
                        logger.log(f"Resubmitting job for {molecule.name} due to error termination")
                        os.remove(log_file_path)
                        resubmit_job(molecule, logger, job_type)
                        time.sleep(interval)
                        continue
                else:
                    molecule.converged = False
                    all_converged = False
                    time.sleep(interval)
                    continue
            elif status == "Pending":
                if molecule.job_id not in pending:
                    logger.log(f"Job {molecule.job_id} is pending in the queue.")
                    pending.append(molecule.job_id)
                time.sleep(interval)
                continue

            elif status == "Unknown":
                logger.log(f"Status of job {molecule.job_id} is unknown. Checking log file.")
                log_file_name = f"{molecule.name}.log"
                log_file_path = os.path.join(molecule.directory, log_file_name)
                last_lines = read_last_lines(log_file_path, logger, 10)
                termination_detected = any(termination in line for line in last_lines)
                error_termination_detected = any(error_termination in line for line in last_lines)
                logger.log(f"Termination status for {molecule.name}: {termination_detected}")

            if all({m.converged for m in molecules}):
                all_converged = True
                break

        time.sleep(interval)
        attempts += 1
        logger.log(f"Log files of the {len(molecules)} conformers has been checked. Checking again in {interval} seconds. Attempt: {attempts}/{max_attempts}")

    if all_converged:
        logger.log_with_stars("Yay! All conformer jobs have converged.")
        move_files(molecules[0].directory, molecules[0].input)
        if job_type == "DLPNO_conf":
            [global_molecules.append(molecule) for molecule in molecules]
            return True
        else:
            logger.log(f"Job type {job_type} for conformers finished. Next step is: {molecules[0].next_step}")
            converged_job(molecules, logger, threads)
            return True


def check_convergence(molecule, logger, threads):
    time_seconds = get_time(molecule)
    log_file_name = f"{molecule.name}.log"
    directory = molecule.directory
    log_file_path = os.path.join(directory, log_file_name)
    molecule.log_file_path = log_file_path
    job_type = molecule.current_step

    if args.attempts: max_attempts = args.attempts
    else: max_attempts = 100
    if args.initial_delay: initial_delay = args.initial_delay
    else: initial_delay = time_seconds
    if args.interval: interval = args.interval
    else: interval = initial_delay/2
    attempts = 0
    running = 1
    pending = 1
    crest_running = 1

    if job_type in ['preopt', 'crest_sampling'] or molecule.name == 'OH':
        initial_delay = time_seconds//2
        interval = initial_delay//2

    if args.OH:
        freq_cutoff = -250
    elif args.CC:
        freq_cutoff = -50 # TEST
    else:
        freq_cutoff = -50 # TEST

    while attempts < max_attempts:
        status = check_job_status(molecule.job_id)

        if status in ["Running", "Completed or Not Found"]:
            if running:
                logger.log(f"Job {job_type} for {molecule.name} with job id {molecule.job_id} is {status}.")
                running = 0
            if status == "Completed or Not Found":
                running = 1
            if not os.path.exists(log_file_path):
                if crest_running:
                    logger.log("CREST log file not generated yet")
                    crest_running = 0
                time.sleep(initial_delay) # crest log file and conformer file is first generated when the CREST calculation finishes
                continue

            last_lines = read_last_lines(log_file_path, logger, 10) 
            termination = termination_strings.get(molecule.program.lower(), "")
            error_termination = error_strings.get(molecule.program.lower(), "")
            termination_detected = any(termination in line for line in last_lines)
            error_termination_detected = any(error_termination in line for line in last_lines)
            link1_detected = any("Link1:  Proceeding to internal job step number  2" in line for line in last_lines)
            if link1_detected:
                time.sleep(interval)
                continue
            if termination_detected:
                if 'OH' in molecule.name and job_type == 'DLPNO_SP':
                    molecule.update_energy(logger, DLPNO=True)
                    global_molecules.append(molecule)
                    logger.log_with_stars(f"Yay! DLPNO single point for OH molecule converged")
                    return True
                if job_type == 'crest_sampling':
                    xyz_conformers = log2xyz(molecule)
                    if xyz_conformers:
                        conformer_molecules = []
                        for n, conformer_coords in enumerate(xyz_conformers, start=1):
                            conformer_molecule = copy.deepcopy(molecule)
                            conformer_molecule.name = f"{molecule.name}_conf{n}"
                            conformer_molecule.atoms = [atom[0] for atom in conformer_coords]
                            conformer_molecule.coordinates = [atom[1:] for atom in conformer_coords]
                            conformer_molecules.append(conformer_molecule)
                        logger.log_with_stars(f"Yay CREST job {log_file_name} converged")
                        converged_job(conformer_molecules, logger, threads)
                        return True
                    else:
                        logger.log(f"ERROR: could not extract coordinates from crest_conformers.xyz")
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
                                # Perhaps constraining methyl groups may work
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
                logger.log(f"Error termination in {log_file_name}. Gathering last XYZ coordinates")
                molecule.update_energy(logger)
                xyz_coordinates = log2xyz(molecule)
                if xyz_coordinates:
                    logger.log(f"XYZ coordinates found in failed log file {log_file_name}. Trying to resubmit job")
                    molecule.coordinates = xyz_coordinates
                    os.remove(log_file_path)
                    non_converged_job(molecule, logger, threads)
                    return False
                else:
                    logger.log(f"No XYZ coordinates found in {log_file_name}. Check log file for errortype.")
                    return False
            else:
                crest_conformer_file = os.path.join(molecule.directory, "crest_conformers.xyz")
                if job_type == 'crest_sampling' and os.path.exists(log_file_path) and os.path.exists(crest_conformer_file):
                    logger.log("Log file found, but no sort of termination yet.")
                    logger.log("Checking for crest_conformers.xyz file")
                    crest_xyz_file = os.path.exists(os.path.join(molecule.directory, "crest_conformers.xyz"))
                    if crest_xyz_file:
                        logger.log("crest_conformers.xyz found! Extracting conformers")
                        xyz_conformers = log2xyz(molecule)
                        if xyz_conformers:
                            logger.log("Conformers extracted. Trying to submit next step")
                            conformer_molecules = []
                            for n, conformer_coords in enumerate(xyz_conformers, start=1):
                                conformer_molecule = copy.deepcopy(molecule)
                                conformer_molecule.name = f"{molecule.name}_conf{n}"
                                conformer_molecule.atoms = [atom[0] for atom in conformer_coords]
                                conformer_molecule.coordinates = [atom[1:] for atom in conformer_coords]
                                conformer_molecules.append(conformer_molecule)
                            logger.log(f"No termination found in {log_file_name}. However, conformers could still be extracted from crest_conformers.xyz")
                            converged_job(conformer_molecules, logger, threads)
                            return True
                        else:
                            logger.log(f"ERROR: could not extract coordinates from crest_conformers.xyz")
                            return False

                # In some cases the log file is crashed and a termination will never occur.
                # Implement time stamp for log file and check if it changes

        elif status == "Pending":
            if pending:
                logger.log(f"Job {molecule.job_id} is pending in the queue.")
                pending = 0
        time.sleep(interval)
        attempts += 1

    logger.log("Max attempts reached or job status is not running. Calculation may be stuck or completed.")
    return False



def submit_and_monitor(molecule, logger, threads):
    if isinstance(molecule, list):
        logger.log(f"Submitted SLURM array job for conformers in {molecule[0].directory}")
        job_id = submit_array_job(molecule, args.par, args.time, args.cpu, args.mem)

        if job_id is not None:
            thread = threading.Thread(target=check_convergence_for_conformers, args=(molecule, logger, threads))
            threads.append(thread)
            thread.start()

    else:
        logger.log(f"Submitting file {molecule.name}{molecule.input} for calculation in path {molecule.directory}")
        job_id = submit_job(molecule, args.cpu, args.mem, args.par, args.time)

        if job_id is not None:
            thread = threading.Thread(target=check_convergence, args=(molecule, logger, threads))
            threads.append(thread)
            thread.start()


def converged_job(molecule, logger, threads):
    if not isinstance(molecule, list):
        molecule.update_step()
        job_type = molecule.current_step

        if job_type == 'crest_sampling':
            molecule.name += "_CREST"
            molecule.program = 'CREST'
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}.xyz")
            molecule.write_xyz_file(output_file_path)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'TS_opt':
            molecule.name.replace("_CREST","")
            molecule.name += '_TS'
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.output}")
            molecule.update_file_path(output_file_path)
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=True)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'DLPNO_SP':
            molecule.name.replace("_TS", "")
            molecule.name += "_DLPNO"
            molecule.program = 'ORCA'
            QC_input(molecule, constrain=False, method='DLPNO', basis_set='DLPNO', TS=False)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'optimization': #Usually for OH radical
            molecule.name.replace("_CREST", "")
            molecule.program = global_program
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.output}")
            molecule.update_file_path(output_file_path)
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=False)
            submit_and_monitor(molecule, logger, threads)


    elif isinstance(molecule, list):
        conformer_molecules = molecule
        for conf in conformer_molecules: conf.update_step()
        job_type = conformer_molecules[0].current_step
        for conf in conformer_molecules:
            conf.program = global_program

            if job_type == 'opt_conf': # Changed
                QC_input(conf, constrain=True, method=low_method, basis_set=low_basis, TS=False)

            elif job_type == 'opt_conf_reactants':
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=False)
            
            elif job_type == 'TS_opt_conf':
                conf.name += '_TS'
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=True)

            elif job_type == 'DLPNO_conf':
                conf.program = 'ORCA'
                conf.name = conf.name.replace("_TS","")
                conf.name += '_DLPNO'
                QC_input(conf, constrain=False, method='DLPNO', basis_set='DLPNO', TS=False)

            elif job_type  == None:
                logger.log("DLPNO calculation has converged")
                break

            else:
                logger.log("Job type could not be determined for conformer list")
                break
        if job_type: submit_and_monitor(conformer_molecules, logger, threads)

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


def determine_current_step(molecule):
    # For reactants
    if molecule.program.lower == 'crest':
        return 'crest_sampling'
    else:
        if molecule.reactant:
            if 'OH' in molecule.name and 'DLPNO' in molecule.name: # Should also add that OH needs freqs
                return 'DLPNO_SP'
            elif 'OH' in molecule.name:
                return 'optimization'
            elif 'reactant_CREST' in molecule.name:
                return 'opt_conf_reactants'
            elif 'reactant' in molecule.name and 'DLPNO' in molecule.name:
                return 'DLPNO_conf'
        # For TS molecules
        else:
            if 'TS' in molecule.name:
                return 'TS_opt_conf'
            elif 'DLPNO' in molecule.name:
                return 'DLPNO_conf'
            else:
                return 'opt_conf'

        # try:
        #     with open(molecule.log_file_path, 'r') as file:
        #         for _ in range(30, 200):  # Check for method section between line 30 and 200. 
        #             try:
        #                 line = next(file).strip()
        #                 if molecule.program.lower == 'g16':
        #                     if line.startswith('#'):
        #                         components = line.split(",")
        #                         if 'ts' not in components:
        #                             return 'TS_opt'
        #                         elif 'ts' in components: 
        #                             if 'conf' in molecule.name:
        #                                 vibrations = log2vib(molecule)
        #                                 if vibrations:
        #                                     return 'DLPNO_conf'
        #                             else: 
        #                                 return 'crest_sampling' 
        #                     elif line.startswith('xTB'):
        #                         return 'TS_opt_conf'
        #                 elif molecule.program.lower == 'orca':
        #                     if line.startswith('|  1>'):
        #                         components = line.split()
        #                         if 'Opt' in components:
        #                             if 'conf' in molecule.name:
        #                                 vibrations = log2vib(molecule)
        #                                 if vibrations:
        #                                     return 'DLPNO_conf'
        #                                 else:
        #                                     return 'crest_sampling'
        #                             else:
        #                                 return 'TS_opt'
        #                 elif molecule.program.lower == 'crest':
        #                     return 'opt_conf'
        #             except StopIteration:
        #                 break  
        #     return 'Method section not found'
        # except FileNotFoundError:
        #     print(f"Error: Log file for molecule {molecule.name} not found.")
        #     return 'error'


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
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=4000, help='Amount of memory allocated for job [def = 4000mb]')
    additional_options.add_argument('-par', metavar="partition", nargs='?', const=1, default="q24,q28,q36,q40,q48,q64", help='Partition to use [def = q24,q28,q36,q40,q48,q64]')
    additional_options.add_argument('-time', metavar="hours:minutes:seconds", nargs='?', const=1, default="72:00:00", help='Specify total time for calculation [def = 72 Hours]')
    additional_options.add_argument('-interval', metavar="int", nargs='?', const=1, type=int, help='Set time interval between checks of log files [def = based on molecule size]')
    additional_options.add_argument('-initial_delay', metavar="int", nargs='?', const=1, type=int, help='Set an initial delay before checking log files [def = based on molecule size]')
    additional_options.add_argument('-attempts', metavar="int", nargs='?', const=1, type=int, help='Set how many times a log files should be checked [def = 100]')
    additional_options.add_argument('-max_conformers', metavar="int", nargs='?', const=1, type=int, help='Set max number of conformers from CREST [def = 50]')

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
    if args.ORCA:
        global_program = "ORCA"
        if low_method in methods_no_basis:
            low_basis = ""
        if high_method in methods_no_basis:
            high_basis = ""
        if high_method.lower() == "wb97xd":
            high_method = "WB97X-D3"
        if low_method.lower() == "wb97xd":
            low_method = "WB97X-D3" 
    elif args.G16:
        global_program = "G16"
    else: 
        global_program = "G16" # Default case

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
    converged_molecules = []
    non_converged_molecules = []

    for n, input_file in enumerate(args.input_files, start=1):
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
                    molecule.set_name(f"{file_name}_H{count}") # Changed
                    molecule.set_directory(os.path.join(start_dir, molecule.name))
                    molecule.program = 'CREST' # Changed
                    molecule.current_step = 'crest_sampling' #Changed
                    molecule.reactant = False
                    mkdir(molecule)
                    logger = Logger(os.path.join(molecule.directory, "log"))

                    # if args.no_xyz is False: #Changed
                    molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz"))
                    submit_and_monitor(molecule, logger, threads)

                    # if molecule.program is not None: # Changed
                    #     QC_input(molecule=molecule, constrain=args.constrain, method=low_method, basis_set=low_basis, TS=False)
                    #     submit_and_monitor(molecule, logger, threads)


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

        elif file_type in [".log", ".out"]:
            # Initialize molecule from log file
            log_logger = Logger(os.path.join(start_dir, "log"))
            log_program = log2program(os.path.join(start_dir, input_file))
            molecule = Molecule(name=f"{input_file.split('.')[0]}", directory=start_dir, program=log_program)
            if 'reactant' in molecule.name or 'OH' in molecule.name: 
                molecule.reactant = True
                if 'OH' not in molecule.name: molecule.mult = 1
            molecule.log_file_path = os.path.join(start_dir, input_file)
            molecule.program = log_program
            current_step = determine_current_step(molecule)
            molecule.current_step = current_step
            last_lines = read_last_lines(input_file_path, log_logger, 10) 
            atoms, xyz_coordinates = log2xyz(molecule, True)
            molecule.atoms, molecule.coordinates = atoms, xyz_coordinates
            if molecule.current_step in ['DLPNO_conf', 'DLPNO_SP']:
                n = re.search(r'conf(\d+)', input_file)
                if n or 'OH' in molecule.name:
                    if molecule.reactant:
                        if 'OH' in molecule.name: thermo_data_log = "OH.log"
                        else: thermo_data_log = re.sub(r"_conf\d+_DLPNO", f"_conf{int(n.group(1))}", input_file) 
                    else:
                        thermo_data_log = re.sub(r"_conf\d+_DLPNO", f"_conf{int(n.group(1))}_TS", input_file)
                    thermo_data_path = os.path.join(start_dir, thermo_data_log)            
                    thermo_program = log2program(thermo_data_path)
                    molecule.update_energy(log_logger, log_file_path=thermo_data_path, program=thermo_program) # first update thermochemistry
                    molecule.update_energy(log_logger, DLPNO=True) # Then update single point energy for DLPNO
                    global_molecules.append(molecule)
                continue

            else:
                termination = termination_strings.get(log_program.lower(), "")
                error_termination = error_strings.get(log_program.lower(), "")

                if any(termination in line for line in last_lines):
                    log_logger.log(f"Detected log file with normal termination. Next step is: {molecule.next_step}")
                    if log_program in ['orca', 'g16']:
                        molecule.converged = True
                        converged_molecules.append(molecule)
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
                                converged_job(conformer_molecules, log_logger, threads)
                            else:
                                converged_job(conformer_molecules, log_logger, threads)
                    
                elif any(error_termination in line for line in last_lines):
                    log_logger.log("Detected log file with error termination.")
                    molecule.converged = False
                    non_converged_molecules.append(molecule)
                    non_converged_job(molecule, log_logger, threads)

                if converged_molecules:
                    converged_job(converged_molecules, log_logger, threads)

                if non_converged_molecules:
                    non_converged_job(non_converged_molecules, log_logger, threads)

        else:
            parser.error(f"Unsupported file type for input: {input_file}")

    # Monitor and handle convergence of submitted jobs
    while threads:
        for thread in list(threads):
            thread.join(timeout=0.1)
            if not thread.is_alive():
                threads.remove(thread)

    for molecule in global_molecules:
        logger = Logger(os.path.join(molecule.directory, "molecules"))
        molecule.log_items(logger)

    k_b = 1.380649e-23  # Boltzmann constant in J/K
    T = 298.15  # Temperature in K
    h = 6.62607015e-34  # Planck constant in J*s
    HtoJ = 4.3597447222e-18  # Conversion factor from Hartree to Joules

    # Molecule lists
    reactant_molecules = [molecule for molecule in global_molecules if molecule.reactant]
    OH_reactant = [molecule for molecule in global_molecules if 'OH' in molecule.name and molecule.reactant]
    TS_molecules = [molecule for molecule in global_molecules if not molecule.reactant]

    # Calculation of rate constant
    if reactant_molecules and TS_molecules:
        # Find the lowest single point energies for reactant (excluding OH) and TS
        lowest_SP_reactant = min([mol for mol in reactant_molecules if 'OH' not in mol.name], key=lambda molecule: molecule.single_point).single_point
        lowest_SP_TS = min(TS_molecules, key=lambda molecule: molecule.single_point).single_point
        OH = OH_reactant[0]

        # Convert energies from Hartree to Joules
        lowest_SP_reactant_J = np.float64(lowest_SP_reactant) * HtoJ
        lowest_SP_TS_J = np.float64(lowest_SP_TS) * HtoJ

        # Sum terms for reactants and transition states
        sum_TS = np.sum([np.exp((lowest_SP_TS_J - np.float64(molecule.single_point * HtoJ)) / np.float64((k_b * T))) * np.float64(molecule.Q) for molecule in TS_molecules])
        sum_reactant = np.sum([np.exp((lowest_SP_reactant_J - np.float64(molecule.single_point * HtoJ)) / np.float64((k_b * T))) * np.float64(molecule.Q) for molecule in reactant_molecules if 'OH' not in molecule.name])        

        # Pre-exponential factor
        pre_exponential = (k_b * T) / h

        # Exponential term for the reference energies
        energy_diff = lowest_SP_TS_J - (lowest_SP_reactant_J + (OH.single_point*HtoJ))

        # Calculate the rate constant, k
        k = pre_exponential * (sum_TS / (sum_reactant*OH.Q)) * np.exp(-energy_diff / (k_b * T))

        # Log the calculated rate constant
        result_logger = Logger(os.path.join(start_dir, "results.log"))
        result_logger.log(f"Calculated rate constant k: {k}")



if __name__ == "__main__":
    main()
