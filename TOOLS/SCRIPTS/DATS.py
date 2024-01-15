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
import copy
import shutil
import pickle
import random

global_molecules = []
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
    def calculate_angle(coord1, coord2, coord3):
        vector1 = VectorManipulation.calculate_vector(coord2, coord1)
        vector2 = VectorManipulation.calculate_vector(coord2, coord3)

        dot_product = np.dot(vector1, vector2)
        magnitude1 = VectorManipulation.vector_length(vector1)
        magnitude2 = VectorManipulation.vector_length(vector2)

        angle_rad = np.arccos(dot_product / (magnitude1 * magnitude2))
        angle_deg = np.degrees(angle_rad)

        return angle_deg


    @staticmethod
    def get_upwards_perpendicular_axis(self, direction):
        # direction is the vector lying between 2 carbons in C=C
        pass


class Logger:
    def __init__(self, log_file):
        self.log_file = log_file
        self.lock = threading.Lock()

    def log(self, message):
        with self.lock:
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

class SCFAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, f"SCF={option_string.strip('-')}")

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
        result = subprocess.run(['squeue', '-j', f'{job_id}'], capture_output=True, text=True, check=True)
        status = parse_job_status(result.stdout)
        return status
    except subprocess.CalledProcessError as e:
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


def submit_array_job(molecules, partition, time, ncpus, mem, nnodes=1):
    dir = molecules[0].directory
    job_files = [f"{conf.name}{conf.input}" for conf in molecules] 
    if molecules[0].product:
        job_name = f"{molecules[0].current_step.split('_')[0]}_prod_{os.path.basename(dir)}"
    elif molecules[0].reactant:
        job_name = f"{molecules[0].current_step.split('_')[0]}_reac_{os.path.basename(dir)}"
    else:
        job_name = f"{molecules[0].current_step.split('_')[0]}_{os.path.basename(dir)}"

    job_program = molecules[0].program
    submit_name = f"{molecules[0].current_step}_submit.sh"
    path_submit_script = os.path.join(dir, submit_name)
    array_txt = "array.txt"
    array_path = os.path.join(dir, array_txt)
    
    if job_program.lower() == "orca" or job_program.lower() == "g16":
        program_mem = mem + 2000 # Headspace for program
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

        file.write("PWD=`pwd`\n\n")

        if job_program.lower() == 'crest':
            file.write(f"cat > $SUBMIT <<!EOF\n")
        else:
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
        # file.write(f"#SBATCH --error=./slurm_output/{job_name}_%A_%a.err\n")
        file.write(f"#SBATCH --output=./slurm_output/{job_name}_%A_%a.out\n")
        file.write(f"#SBATCH --time={time}\n")
        file.write(f"#SBATCH --partition={partition}\n")
        file.write(f"#SBATCH --no-requeue\n")
        file.write(f"#SBATCH --mem={program_mem}\n")
        file.write(f"#SBATCH --array=$ARRAY\n\n")

        file.write("# Create scratch folder\n")
        file.write("SCRATCH=/scratch/\\${SLURM_JOB_ID}/\\${SLURM_ARRAY_TASK_ID}\n")
        file.write("mkdir -p \\$SCRATCH\n\n")

        if job_program.lower() == 'g16':
            file.write(f"cd $PWD\n")
            file.write("export GAUSS_SCRDIR=\\$SCRATCH\n\n")

            file.write("source /comm/groupstacks/gaussian/bin/modules.sh\n")
            file.write("ml gaussian16/Rev.B.01\n")
            file.write("ml gcc/9.2.0\n")
            file.write("ml openmpi/4.0.1\n\n")

            file.write('GJ=\\$(awk "NR == \\$SLURM_ARRAY_TASK_ID" $IN)\n')
            file.write('LOG=\\${GJ%.*}.log\n\n')
            file.write("srun $(which g16) \\$GJ > \\$LOG\n")
            file.write("#\n")
            file.write("!EOF")

        elif job_program.lower() == "orca":
            file.write("  ulimit -c 0\n")
            file.write("source /comm/groupstacks/gaussian/bin/modules.sh\n")
            file.write("ml openmpi\n")
            file.write("ml orca/5.0.4\n\n")

            file.write("cd \\$SCRATCH\n")
            file.write(f"cp \\$SLURM_SUBMIT_DIR/{array_txt} .\n")
            file.write(f"cp \\$SLURM_SUBMIT_DIR/*.inp .\n\n")

            file.write('GJ=\\$(awk "NR == \\$SLURM_ARRAY_TASK_ID" $IN)\n')
            file.write('LOG=\\${GJ%.*}.log\n\n')

            file.write(f"\\$(which orca) \\$GJ > \\$SLURM_SUBMIT_DIR/\\$LOG\n\n")
            file.write("!EOF")

        else:
            file.write('source ~/.JKCSusersetup.txt\n')
            file.write('MOLECULE_NAME=\\$(awk "NR == \\${SLURM_ARRAY_TASK_ID}" "$IN")\n')
            file.write(f'cp \\$SLURM_SUBMIT_DIR/\\$MOLECULE_NAME \\$SCRATCH/\n')
            file.write(f'cd \\$SCRATCH\n')
            file.write(f"program_CREST \\$MOLECULE_NAME -gfn2 -ewin {args.ewin} -noreftopo > \\${{MOLECULE_NAME%.xyz}}.log\n")
            file.write(f"cp *pkl \\$SLURM_SUBMIT_DIR/.\n")
            file.write(f"cp *log \\$SLURM_SUBMIT_DIR/.\n")
            file.write(f"cp *output {dir}/.\n")
            file.write(f"!EOF\n\n")
            file.write(f"sbatch $SUBMIT")

    submit_command = ['sh', path_submit_script, array_txt]
    try:
        result = subprocess.run(submit_command, capture_output=True, text=True, check=True)
        job_id = extract_job_id(result.stdout)
        for n, molecule in enumerate(molecules, start=1): 
            molecule.job_id = f"{job_id}_{n}"
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"Error in job submission: {e}")
        return None


def submit_job(molecule, ncpus, mem, partition, time, nnodes=1):
    input_file_name = f"{molecule.name}{molecule.input}"
    job_name = f"{molecule.name}"
    dir = molecule.directory    
    submit_name = f"{molecule.program}_submit.sh"
    path_submit_script = os.path.join(dir, submit_name)

    if molecule.reactant:
        crest_input = f"{molecule.name}.xyz --gfn2 --ewin {args.ewin} --noreftopo > {molecule.name}.log"
    elif molecule.product:
        crest_input = f"{molecule.name}.xyz --gfn2 --ewin {args.ewin} --noreftopo --uhf 1 > {molecule.name}.log"
    else:
        crest_input = f"{molecule.name}.xyz --gfn2 --ewin {args.ewin} --noreftopo --uhf 1 --cinp {dir}/constrain.inp > {molecule.name}.log"

    if molecule.program.lower() == "orca" or molecule.program.lower() == "g16":
        program_mem = mem + 2000 # Headspace for program
    else:
        program_mem = mem


    script_content_crest = f"""#!/bin/bash
IN=$1
JOB=${{IN%.*}}

SUBMIT=qsub.tmp
PWD=`pwd`

cat > $SUBMIT <<!EOF
#!/bin/sh
#SBATCH --job-name=${{JOB}}_crest
#SBATCH --nodes={nnodes}
#SBATCH --cpus-per-task=1
#SBATCH --ntasks={ncpus}
#SBATCH --output=./slurm_output/{job_name}_%A_%a.out
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --no-requeue
#SBATCH --mem={mem}

source ~/.JKCSusersetup.txt
SCRATCH=/scratch/\\$SLURM_JOB_ID

mkdir -p \\$SCRATCH || exit $?
cd \\$SCRATCH
cp \\$SLURM_SUBMIT_DIR/{molecule.name}.xyz .
program_CREST {crest_input}
cp *pkl \\$SLURM_SUBMIT_DIR/.
cp *output \\$SLURM_SUBMIT_DIR/.
rm -rf /scratch/\\$SLURM_JOB_ID
!EOF

sbatch $SUBMIT
    """
######################################################################################################################

    script_content_G16 = f"""#!/bin/bash
IN=$1
JOB=${{IN%.*}}

SUBMIT=qsub.tmp
PWD=`pwd`

cat > $SUBMIT <<!EOF
#!/bin/sh
#SBATCH --job-name=$JOB
#SBATCH --nodes={nnodes}
#SBATCH --cpus-per-task={ncpus}
#SBATCH --ntasks={nnodes}
#SBATCH --output=./slurm_output/{job_name}_%A_%a.out
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --no-requeue
#SBATCH --mem={program_mem}  # requested total memory in MB

# Create scratch folder
mkdir /scratch/\\$SLURM_JOB_ID
cd $PWD
export GAUSS_SCRDIR=/scratch/\\$SLURM_JOB_ID

srun /comm/groupstacks/gaussian/gaussian/gaussian16/Rev.B.01/g16/g16 $JOB.com > $JOB.log

# Remove scratch folder
rm -rf /scratch/\\$SLURM_JOB_ID

!EOF

sbatch $SUBMIT
"""
######################################################################################################################
    script_content_orca = f"""#!/bin/bash
IN=$1
JOB=${{IN%.*}}

SUBMIT=qsub.tmp
PWD=`pwd`

cat > $SUBMIT <<!EOF
#!/bin/bash
#SBATCH --job-name=$JOB
#SBATCH --nodes={nnodes}
#SBATCH --cpus-per-task=1
#SBATCH --ntasks={ncpus}
#SBATCH --output=./slurm_output/{job_name}_%A_%a.out
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --no-requeue
#SBATCH --mem={mem}  # requested memory per process in MB

source /comm/groupstacks/chemistry/bin/modules.sh
ml openmpi
ml orca/5.0.4

# create and set scratch dir
SCRATCH=/scratch/\\$SLURM_JOB_ID
mkdir -p \\$SCRATCH || exit $?

cd \\$SCRATCH
cp \\$SLURM_SUBMIT_DIR/$JOB.inp .
\\$(which orca) $JOB.inp > \\$SLURM_SUBMIT_DIR/$JOB.log

# Copy back results
cp *log \\$SLURM_SUBMIT_DIR/.
!EOF

sbatch $SUBMIT
    """

    with open(path_submit_script, 'w') as file:
        if molecule.program.lower() == 'g16':
            file.write(script_content_G16)
        elif molecule.program.lower() == 'orca':
            file.write(script_content_orca)
        elif molecule.program.lower() == 'crest':
            file.write(script_content_crest)

    submit_command = ['sh', path_submit_script, input_file_name]
    try:
        result = subprocess.run(submit_command, capture_output=True, text=True, check=True)
        job_id = extract_job_id(result.stdout)
        molecule.job_id = f"{job_id}"
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
    def __init__(self, file_path=None, log_file_path=None, name="", directory="", atoms=None, coordinates=None, reactant=False, product=False, program=None, current_step=None, next_step=None):
        self.name = name
        self.directory = directory
        self.file_path = file_path
        self.log_file_path = log_file_path
        self.reactant = reactant
        self.product = product
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
        self.Q = None
        self._program = None
        self.input = None
        self.output = None
        self.program = program if program is not None else global_program
        self._current_step = None
        self._next_step = None
        self.error_termination_count = 0
        self.wrong_TS = 0

        
        if self.file_path and self.file_path.split(".")[-1] == 'xyz':
            self.read_xyz_file()

        if self.file_path and self.file_path.split(".")[-1] in ['log', 'out']:
            self.name = f"{self.file_path.split('/')[-1].split('.')[0]}"
            self.directory = os.path.dirname(self.file_path)
            self.program = log2program(file_path)
            self.log_file_path = file_path
            self.atoms, self.coordinates = log2xyz(self, True)
            self.update_energy()

            if self.reactant or 'reactant' in self.name:
                self.reactant = True
                self.mult = 1
                self.current_step = determine_current_step(self)

            elif self.product or 'product' in self.name:
                self.product = True
                self.current_step = determine_current_step(self)
        
            else:
                self.current_step = determine_current_step(self)
                self.find_active_site()

        if 'OH' in self.name:
            self.atoms = ['O','H']
            self.reactant = True
            self.directory = start_dir
            self.current_step = 'optimization'
            self.log_file_path = file_path
            if args.high_level is None: # Use optimized geometry from wB97X-D and single point from DLPNO-CCSD(T)
                self.coordinates = [[0.000000, 0.000000, 0.107778], [0.000000, 0.000000, -0.862222]]
                self.single_point = -75.645450513535 
                self.vibrational_frequencies = [3778.5332]
                self.rot_temps = [27.18932]
                self.symmetry_num = 1
                self.mol_mass = 17.00274
                self.partition_function()
            else:
                self.coordinates = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0]]
    
        elif 'H2O' in self.name:
            self.atoms = ['O', 'H', 'H']
            self.mult = 1
            self.product = True
            self.directory = start_dir
            self.current_step = 'optimization'
            if args.high_level is None: # Use optimized geometry from wB97X-D and single point from DLPNO-CCSD(T)
                O_coords = [0.000000, -0.000000, 0.116473]
                H1_coords = [0.000000, -0.760262, -0.465892]
                H2_coords = [0.000000, 0.760262, -0.465892]
                self.single_point = -76.342050858187
                self.vibrational_frequencies = [1638.1838, 3874.6394, 3981.4927]
                self.rot_temps = [39.95124, 20.81842, 13.68646]
                self.symmetry_num = 1
                self.mol_mass = 18.01056
                self.partition_function()
            else:
                O_coords = [0.0, 0.0, 0.0]
                H1_coords = [0.9572, 0.0, 0.0]
                H2_coords = [-0.239664, 0.926711, 0.000000]
            self.coordinates = [O_coords, H1_coords, H2_coords]


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

    def move_converged(self):
        if not os.path.exists(os.path.join(self.directory, "log_files")):
            os.makedirs(os.path.join(self.directory, "log_files"), exist_ok=True)
        destination = os.path.join(self.directory, "log_files", f"{self.name}{self.output}")
        if os.path.exists(destination):
            os.remove(destination)
        shutil.move(os.path.join(self.directory, f"{self.name}{self.output}"), destination)
        self.log_file_path = destination

    def move_failed(self):
        if not os.path.exists(os.path.join(self.directory, "failed_logs")):
            os.makedirs(os.path.join(self.directory, "failed_logs"), exist_ok=True)
        destination = os.path.join(self.directory, "failed_logs", f"{self.name}{self.output}")
        if os.path.exists(destination):
            os.remove(destination)
        shutil.move(os.path.join(self.directory, f"{self.name}{self.output}"), destination)
            
    def write_xyz_file(self, output_file_path):
        with open(output_file_path, 'w') as file:
            file.write(str(len(self.atoms)) + '\n\n')
            for atom, coord in zip(self.atoms, self.coordinates):
                coord_str = ' '.join(['{:.6f}'.format(c) for c in coord])  
                file.write(f'{atom} {coord_str}\n')

    def read_constrained_indexes(self):
        try:
            with open(os.path.join(self.directory, ".constrain"), "r") as f:
                self.constrained_indexes = {}
                for line in f:
                    parts = line.strip().split(":")
                    if len(parts) == 2:
                        atom, index = parts[0].strip(), int(parts[1].strip())
                        self.constrained_indexes[atom] = index
        except FileNotFoundError:
            print(f"No .constrain file found in {self.directory}")
        except Exception as e:
            print(f"Error reading .constrain file: {e}")

    def find_active_site(self):
        for C, atom in enumerate(self.atoms):
            if atom == 'C':
                for H, neighbor in enumerate(self.atoms):
                    if neighbor == 'H' and self.is_nearby(C, H):
                        for O, radical in enumerate(self.atoms):
                            if radical == 'O' and self.is_nearby(H, O) and self.calculate_angle(self.coordinates[C], self.coordinates[H], self.coordinates[O]) > args.reaction_angle - 30:
                                self.constrained_indexes = {'C': C+1, 'H': H+1, 'O': O+1, 'OH': O+2}


    def is_nearby(self, atom_index1, atoms_index2, threshold_distance=1.5):
        distance = np.linalg.norm(np.array(self.coordinates[atom_index1]) - np.array(self.coordinates[atoms_index2]))
        return distance < threshold_distance

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
        self._program = (value or global_program).lower()
        if self._program == 'orca':
            self.input = '.inp'
            self.output = '.log' #.out
        elif self._program == 'g16':
            self.input = '.com'
            self.output = '.log'
        elif self._program == 'crest':
            self.input = '.xyz'
            self.output = '.output'
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

    def _update_next_step(self):
        workflow = self._get_workflow()
        if self._current_step in workflow:
            current_index = workflow.index(self._current_step)
            self._next_step = workflow[current_index + 1] if current_index + 1 < len(workflow) else None

    def update_step(self):
        if self._current_step is None and self._next_step is None:
            self._current_step = 'optimization' if self.name in ['OH', 'H2O'] else 'crest_sampling'
            self._next_step = self._get_initial_next_step()
        elif self._next_step is not None:
            self.current_step = self._next_step

    def _get_initial_next_step(self):
        workflow = self._get_workflow()
        if self._current_step in workflow:
            current_index = workflow.index(self._current_step)
            return workflow[current_index + 1] if current_index + 1 < len(workflow) else None
        return None

    def _get_workflow(self):
        if self.reactant or self.product:
            if self.name == 'OH' or self.name == 'H2O':
                return ['optimization', 'DLPNO', 'Done']
            else:
                return ['crest_sampling', 'opt', 'DLPNO', 'Done']
        else:
            return ['crest_sampling', 'opt_constrain', 'TS_opt', 'DLPNO', 'Done']


    def update_energy(self, logger=None, log_file_path=None, DLPNO=False, program=None):
        file_path = log_file_path if log_file_path else self.log_file_path
        program = program if program else self.program

        with open(file_path, 'r') as f:
            log_content = f.read()

            if program.lower() == 'g16':
                energy_matches = re.findall(r'(SCF Done:  E\(\S+\) =)\s+([-.\d]+)', log_content)
                if energy_matches:
                    self.single_point = float(energy_matches[-1][-1])
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
                energy_matches = re.findall(r'(FINAL SINGLE POINT ENERGY)\s+([-.\d]+)', log_content)
                if energy_matches:
                    self.single_point = float(energy_matches[-1][-1])
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
                qvib *= 1 / (1 - np.exp(-(h * f) / (k_b * T))) 

        # Rotational partition function
        if self.program.lower() == 'g16':
            rot_temps = self.rot_temps
            if len(rot_temps) == 1:
                qrot = ((1/self.symmetry_num) * (T/rot_temps[0]))
            else:
                qrot = (np.pi**0.5/self.symmetry_num) * ((T**3 / (rot_temps[0] * rot_temps[1] * rot_temps[2]))**0.5)
        elif self.program.lower() == 'orca':
            rot_constants = self.rot_temps
            if 0.0 in rot_constants:
                rot_constant = [e for e in rot_constants if e != 0.0]  
                qrot = ((k_b*T)/(self.symmetry_num*h*c*100*rot_constant[0]))
            else:
                qrot = (1/self.symmetry_num) * ((k_b*T)/(h*c*100))**1.5 * (np.pi / (rot_constants[0] * rot_constants[1] * rot_constants[2]))**0.5


        # Translational partition function
        mol_mass = self.mol_mass * 1.66053906660e-27
        V = k_b*T/P # R*T/P
        qtrans = ((2 * np.pi * mol_mass * k_b * T) / h**2)**(3/2) * V

        if 'OH' in self.name:
            qelec = 3 # OH radical with 2 low lying near degenerate energy level
        else:
            qelec = self.mult

        self.Q = qvib*qrot*qtrans*qelec


    def set_active_site(self, perturb=False, distance_CH=None, distance_HO=None, distance_OH_H=None, reaction_angle=None):
        # python indexing starting from 0
        C_index = self.constrained_indexes['C']-1
        H_index = self.constrained_indexes['H']-1
        O_index = self.constrained_indexes['O']-1
        OH_index = self.constrained_indexes['OH']-1

        if perturb:
            scaling_factor = 0.1
            original_coordinates = self.coordinates
            perturbed_coords = original_coordinates.copy()
            for index in [C_index, H_index, O_index]:  # Adjusted for 0 based indexing in ORCA
                random_perturbation = [random.uniform(-scaling_factor, scaling_factor) for _ in range(3)]
                perturbed_coords[index] = [coord + delta for coord, delta in zip(original_coordinates[index], random_perturbation)]
            self.coordinates = perturbed_coords

        else:
            # If distances and angle are not given, use molecules current values
            if distance_CH is None:
                distance_CH = self.atom_distance(self.coordinates[C_index], self.coordinates[H_index])
            if distance_HO is None:
                distance_HO = self.atom_distance(self.coordinates[H_index], self.coordinates[O_index])
            if distance_OH_H is None:
                distance_OH_H = self.atom_distance(self.coordinates[O_index], self.coordinates[OH_index])
            if reaction_angle is None:
                reaction_angle = self.calculate_angle(self.coordinates[C_index], self.coordinates[H_index], self.coordinates[O_index])

            # Set the C-H distance
            vector_CH = self.calculate_vector(self.coordinates[C_index], self.coordinates[H_index])
            norm_vector_CH = self.normalize_vector(vector_CH)
            new_H_position = self.coordinates[C_index] - (norm_vector_CH * distance_CH)
            if np.dot(self.calculate_vector(self.coordinates[C_index], new_H_position), vector_CH) < 0:
                new_H_position = self.coordinates[C_index] + (norm_vector_CH * distance_CH)

            # Set the H-C-O angle
            complement_angle = 180.0 - reaction_angle
            rotation_angle = np.radians(complement_angle)
            rotation_axis = self.normalize_vector(np.cross(norm_vector_CH, [1, 0, 0]))
            if np.all(rotation_axis == 0):
                rotation_axis = self.normalize_vector(np.cross(norm_vector_CH, [0, 1, 0]))
            rotated_vector = self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle)

            # Set the H-O distance
            new_O_position = new_H_position - (rotated_vector * distance_HO)
            if np.dot(self.calculate_vector(self.coordinates[C_index], new_O_position), vector_CH) < 0:
                new_O_position = new_H_position + rotated_vector * distance_HO 
            rotation_axis = self.normalize_vector(rotation_axis)

            rotation_angle_H = np.radians(104.5)
            new_OH_H_position = new_O_position - self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H
            HOH_angle = self.calculate_angle(new_OH_H_position, new_O_position, new_H_position)
            if HOH_angle < 95:
                new_OH_H_position = new_O_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H

            # Update positions
            self.coordinates[H_index] = new_H_position
            self.coordinates[O_index] = new_O_position
            self.coordinates[OH_index] = new_OH_H_position

    def H_abstraction(self, NEB=False, num_molecules=None):
        original_coords = self.coordinates.copy()
        atoms = self.atoms
        num_atoms = len(atoms)
        abstraction_molecules = []
        product_molecules = []
        methyl_C_indexes = self.methyl_C_index()  
        carbon_iteration_counter = {index: 0 for index in range(num_atoms) if atoms[index] == 'C'}
        distance_CH = 1.3  # Adjusted C-H distance
        distance_OH = 1.25
        distance_OH_H = 0.97
        CHO_angle = args.reaction_angle
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

                            new_coords = original_coords.copy()
                            new_atoms = atoms.copy()

                            if args.products:
                                product_coords = original_coords.copy()
                                product_atoms = atoms.copy()
                                product_coords.pop(i)
                                product_atoms.pop(i)
                                # Add the product molecule to the list
                                product_molecule = Molecule(self.file_path, reactant=self.reactant, program=self.program)
                                product_molecule.atoms = product_atoms
                                product_molecule.coordinates = product_coords
                                product_molecules.append(product_molecule)

                            new_H_position = np.array(original_coords[i]) + norm_vector_CH * (dist_CH - distance_CH)
                            distance_CH = self.vector_length(self.calculate_vector(new_H_position, new_coords[j]))
                            if distance_CH < dist_CH:
                                norm_vector_CH = -norm_vector_CH
                                new_H_position = np.array(original_coords[i]) + norm_vector_CH * (dist_CH - distance_CH)


                            complement_angle = 180 - CHO_angle
                            rotation_angle = np.radians(complement_angle)
                            rotation_axis = self.normalize_vector(np.cross(norm_vector_CH, [1, 0, 0]))
                            if np.all(rotation_axis == 0):
                                rotation_axis = self.normalize_vector(np.cross(norm_vector_CH, [0, 1, 0]))
                            rotated_vector = self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle)

                            oxygen_position = new_H_position - rotated_vector * distance_OH 
                            if self.atom_distance(oxygen_position, original_coords[j]) < distance_OH:
                                oxygen_position = new_H_position + rotated_vector * distance_OH 

                            # norm_vector_OH = self.normalize_vector(oxygen_position - new_H_position)
                            rotation_axis = self.normalize_vector(rotation_axis)

                            rotation_angle_H = np.radians(104.5)
                            new_OH_H_position = oxygen_position + self.rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H) * distance_OH_H

                            # Perturb initial hydrogen to avoid linear C--H--O configuration
                            new_coords[i] = (new_H_position).tolist()

                            new_atoms.append('O')
                            new_atoms.append('H')
                            new_coords.append(oxygen_position.tolist())
                            new_coords.append(new_OH_H_position.tolist())

                            if NEB:
                                # NEB specific code goes here
                                pass

                            new_molecule = Molecule(self.file_path, reactant=self.reactant, product=self.product, program=self.program)
                            new_molecule.atoms = new_atoms
                            new_molecule.coordinates = new_coords
                            new_molecule.constrained_indexes = {'C': j+1, 'H': i+1, 'O': len(new_coords)-1, 'OH': len(new_coords)}
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


    def print_items(self):
        print(f"Molecule: {self.name}")
        print(f"File Path: {self.file_path}")
        print(f"Directory: {self.directory}")
        print(f"Log File Path: {self.log_file_path}")
        print(f"Program: {self.program.upper()}")
        print(f"Reactant: {self.reactant}")
        print(f"Product: {self.product}")
        print(f"Multiplicity: {self.mult}")
        print(f"Charge: {self.charge}")
        print(f"Constrained Indexes: {self.constrained_indexes}")
        print(f"Single Point Energy: {self.single_point}")
        print(f"Partition Function: {self.Q}")
        print(f"Vibrational Frequencies: {self.vibrational_frequencies}")
        print(f"Current Step: {self.current_step}")
        print(f"Next Step: {self.next_step}")
        print("Atoms and Coordinates:")
        for atom, coord in zip(self.atoms, self.coordinates):
            print(f"  {atom}: {coord}")
        print("-----------------------------------------------------------------------")
   
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
        logger.log("-----------------------------------------------------------------------")


def get_time(molecule):
    job_type = molecule.current_step
    atoms = molecule.atoms
    if molecule.reactant:
        heavy_count = 0
    else:
        heavy_count = -1
    H_count = 0
    for atom in atoms:
        if atom == 'H':
            H_count += 1
        else:
            heavy_count += 1
    
    if 'TS' in job_type:
        a = 42.1; b = -120.4; c = 732.5
        interval = (a*heavy_count**2) + (b*heavy_count) + c
    elif job_type == 'crest_sampling':
        a = 27; b = 10
        interval = (a*heavy_count) + b
    elif job_type == 'optimization':
        a=0.4; b=2.5; c=13.4; d=4.8
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    else:
        a=2.26; b=-30.8; c=249.8; d=-114
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d

    return max(interval, 10)

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def mkdir(molecule, CREST=True):
    if not os.path.exists(molecule.directory):
        os.makedirs(molecule.directory, exist_ok=True)
        if not os.path.exists(os.path.join(molecule.directory, "input_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "input_files")))
        if not os.path.exists(os.path.join(molecule.directory, "log_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "log_files")))
    if CREST:
        crest_constrain(molecule)


def move_files(directory, file_type, ignore_file='constrain.inp'):
    """
    Moves all files of a given type in the specified directory to a folder named 'input_files'.
    If 'input_files' folder does not exist, it is created.

    Parameters:
    directory (str): The path of the directory to search for files.
    file_type (str): The file extension to search for, e.g., '.log'.
    """

    if file_type in ['.com', '.inp']:
        input_folder = os.path.join(directory, 'input_files')
    else:
        input_folder = os.path.join(directory, 'log_files')

    if not os.path.exists(input_folder):
        os.makedirs(input_folder)

    for filename in os.listdir(directory):
        if filename.endswith(file_type) and filename != ignore_file:
            file_path = os.path.join(directory, filename)
            dest_file_path = os.path.join(input_folder, filename)
            if os.path.exists(dest_file_path):
                os.remove(dest_file_path)
            shutil.move(file_path, dest_file_path)


def pkl_to_xyz(pkl_file_path, max_conformers=50):
    if args.max_conformers:
        max_conformers = args.max_conformers
    all_conformers = []
    with open(pkl_file_path, 'rb') as f:
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
    return all_conformers[:max_conformers]

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


def handle_incorrect_TS(molecule):
    '''Approach: analyse the normal coordinates associated with the imag freq and determine whether it is a methyl rotation, wiggling molecule or something else.
       If methyl rotation perhaps fix involving atoms; small wiggling: move C--H--O atoms into a better resemples of transition state.
       Ensure other parts of the molecule are not positioned in a way unwanted interactions can occur'''
    C_index = molecule.constrained_indexes['C']
    H_index = molecule.constrained_indexes['H']
    O_index = molecule.constrained_indexes['O']
    if molecule.program.lower() == 'orca':
        C_index -= 1
        H_index -= 1
        O_index -= 1


def check_transition_state(molecule):
    pass
                    


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

    with open (molecule.directory + "/.constrain","w") as f:
        f.write(f"C: {C_index}\n")
        f.write(f"H: {H_index}\n")
        f.write(f"O: {O_index}\n")

def QC_input(molecule, constrain,  method, basis_set, TS, C_index=None, H_index=None, O_index=None):
    file_name = f"{molecule.name}{molecule.input}"
    file_path = os.path.join(molecule.directory, file_name)
    atoms = molecule.atoms
    coords = molecule.coordinates

    # Ensure correct usage of unretricted method for Gaussian16
    if molecule.program.lower() == 'g16':
        if molecule.mult == 2 and method[0].lower() != "u":
            method = f"u{method}"
        elif molecule.mult == 1 and method[0].lower() == "u":
            method = method[1:]

    if constrain or TS:
        if not molecule.constrained_indexes:
            molecule.find_active_site()
        if molecule.program.lower() == 'orca': # ORCA indexes from 0
            C_index = molecule.constrained_indexes['C']-1
            H_index = molecule.constrained_indexes['H']-1
            O_index = molecule.constrained_indexes['O']-1
        else:
            C_index = molecule.constrained_indexes['C']
            H_index = molecule.constrained_indexes['H']
            O_index = molecule.constrained_indexes['O']

    if 'H2O' in molecule.name or 'OH' in molecule.name or molecule.reactant: # dont need for H2O?
        freq = "freq"
    else:
        freq = ""

    if molecule.program.lower() == "orca":
        with open(file_path, "w") as f:
            if method == "DLPNO":
                f.write(f"! {basis_set} {basis_set}/C DLPNO-CCSD(T) TightSCF RI-JK {basis_set}/JK\n")
            elif TS:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OptTS freq\n")
            else:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OPT {freq}\n")
            if method == 'DLPNO':
                f.write(f"%pal nprocs 1 end\n")
                f.write(f"%maxcore {args.mem+14000}\n") # Find some scaling with the molecule to ensure enough memory
            else:
                f.write(f"%pal nprocs {args.cpu} end\n")
                f.write(f"%maxcore {args.mem}\n")
            if constrain:
                f.write("%geom\n")
                f.write("Constraints\n")
                f.write(f"{{B {C_index} {H_index} C}}\n") # Bond being broken
                f.write(f"{{B {H_index} {O_index} C}}\n") # Bond being formed
                # f.write(f"{{A {C_index-1} {H_index-1} {O_index-1} C}}\n") # CHO angle
                f.write("end\nend\n")
            elif TS:
                f.write("%geom\n")
                f.write("maxiter 80\n")
                f.write("Calc_Hess true\n")
                f.write("Recalc_Hess 5\n")
                # f.write(f"TS_Mode {{ B {C_index} {H_index} }} end\n")
                f.write(f"TS_Active_Atoms {{ {C_index} {H_index} {O_index} }} end\n")
                f.write("TS_Active_Atoms_Factor 3\n")
                f.write("end\n")
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
                f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,modredundant,MaxCycles=60) freq {args.SCF}\n\n')
            elif TS and constrain is False:
                if molecule.wrong_TS == 1:
                    f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,CalcAll,MaxCycles=60) freq {args.SCF}\n\n')
                elif molecule.wrong_TS == 2:
                    f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,CalcAll,MaxCycles=60, MaxStep=10) freq {args.SCF}\n\n')
                elif molecule.wrong_TS == 3:
                    molecule.set_active_site(perturb=True)
                    f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,CalcAll,MaxCycles=60) freq {args.SCF}\n\n')
                else:
                    f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,MaxCycles=60) freq {args.SCF}\n\n')
            elif TS is False and constrain:
                f.write(f'# {method} {basis_set} opt=modredundant {args.SCF}\n\n')
            else:
                f.write(f"# {method} {basis_set} opt {freq} {args.SCF}\n\n")
            f.write("Title\n\n")
            f.write(f"{molecule.charge} {molecule.mult}\n")
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("\n")
            if constrain:
                f.write(f"B {C_index} {H_index} F\n")
                f.write(f"B {H_index} {O_index} F\n")
                # f.write(f"A {C_index} {H_index} {O_index} {args.reaction_angle} F\n") # CHO angle
                f.write("\n")
    else:
        print(f"QC_input was called but no program was specified for {molecule.name}")


def NEP_input(file_path, file_name):
    if args.NEB:
        with open(file_path + "/NEB_TS.inp", "w") as f:
            f.write(f'! B3LYP 6-31+g(d,p)  NEB-TS FREQ\n')
            f.write(f'%NEB PREOPT_ENDS TRUE NEB_END_XYZFILE "{file_path + "/" + file_name}_product.xyz" END\n')
            f.write(f'* XYZfile 0 2 {file_path + "/" + file_name}_reactant.xyz\n')

def energy_cutoff(molecules, logger, initial_cutoff=args.energy_cutoff, max_cutoff_increment=30.0, increment_step=5.0):
    """
    Filter molecules to those within a certain energy range of the lowest energy conformer.
    Adjusts cutoff to avoid removing more than 50% of molecules.

    molecules: List of molecule objects.
    initial_cutoff: Initial energy cutoff in kcal/mol (default is 5 kcal/mol).
    max_cutoff_increment: Maximum additional cutoff to add (default is 10 kcal/mol).
    increment_step: Step size for increasing the cutoff (default is 1 kcal/mol).
    """
    hartree_to_kcalmol = 627.509

    # Convert to kcal/mol from Hartree and sort by energy
    energies = [molecule.single_point * hartree_to_kcalmol for molecule in molecules]
    min_energy = min(energies)

    cutoff = initial_cutoff
    while cutoff <= initial_cutoff + max_cutoff_increment:
        # Filter molecules within the current cutoff range
        filtered_molecules = [molecule for molecule, energy in zip(molecules, energies)
                              if (energy - min_energy) <= cutoff]

        # Check if more than 50% of molecules are retained
        if len(filtered_molecules) >= 0.5 * len(molecules):
            return filtered_molecules

        # Increase cutoff
        cutoff += increment_step

    # Return filtered list if cutoff reaches the maximum limit
    return filtered_molecules


def read_last_lines(filename, logger, num_lines, interval=120):
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
            # logger.log(f"Attempt {attempts}: File {filename} not found. Sleeping for {interval} seconds before retrying.")
            time.sleep(interval)

    logger.log(f"Failed to find the file {filename} after {max_attempts} attempts.")
    return False

def resubmit_job(molecule, logger, job_type):
    molecule.move_failed()
    if job_type == 'opt_constrain':
        QC_input(molecule, constrain=True, method=low_method, basis_set=low_basis, TS=False)

    elif job_type == 'TS_opt':
        QC_input(molecule, constrain=False,  method=high_method, basis_set=high_basis, TS=True)

    elif job_type == 'DLPNO':
        molecule.program = 'ORCA'
        QC_input(molecule, constrain=False, method="DLPNO", basis_set=high_basis, TS=False)

    elif job_type == 'opt':
        QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=False)

    if args.par == 'qtest':
        job_time = '1:00:00'
        time_seconds = 3600 / args.attempts
    else:
        time_seconds = get_time(molecule)
        total_seconds = time_seconds * (args.attempts + 5) # Taking into account initial_delay
        job_time = seconds_to_hours(total_seconds)

    job_id = submit_job(molecule, args.cpu, args.mem, args.par, job_time)
    molecule.job_id = f"{job_id}"
    logger.log(f"submitted file {molecule.name}{molecule.input} with job type: {job_type} and job id {molecule.job_id} in path {molecule.directory}")


def check_convergence_for_conformers(molecules, logger, threads, time_seconds, max_attempts):
    if args.initial_delay:
        initial_delay = args.initial_delay
    else:
        initial_delay = int(time_seconds*3)
    if args.interval:
        interval = args.interval
    else:
        interval = int(time_seconds) 
    attempts = 0
    pending = []
    running = []
    job_type = molecules[0].current_step

    if args.freq_cutoff:
        freq_cutoff = args.freq_cutoff
    elif args.OH:
        freq_cutoff = -150
    elif args.CC:
        freq_cutoff = -100
    else:
        freq_cutoff = -150

    termination = termination_strings.get(molecules[0].program.lower(), "")
    error_termination = error_strings.get(molecules[0].program.lower(), "")

    all_converged = False
    for m in molecules:  # Initialize with all molecules not being converged and no terminations counted
        m.error_termination_count = 0
        m.wrong_TS = 0

    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)

    while attempts < max_attempts and not all_converged:
        i = 0
        while i < len(molecules):
            molecule = molecules[i]
            if molecule.error_termination_count >= 3 or molecule.wrong_TS >= 4:
                logger.log(f"Molecule {molecule.name} is being dropped due to repeated error terminations.")
                molecules.pop(i)
                # Check if all remaining molecules are converged
                if all(m.converged for m in molecules):
                    all_converged = True
                    break
            else:
                i += 1

        if all_converged:
            break

        for molecule in molecules:
            if molecule.converged:
                continue
            status = check_job_status(molecule.job_id)
            if status in ['Running', 'Completed or Not Found'] or not molecule.converged:
                if molecule.job_id not in running:
                    logger.log(f"Job {job_type} for {molecule.name} with job id {molecule.job_id} is running.")
                    running.append(molecule.job_id)

                log_file_name = f"{molecule.name}{molecule.output}"
                log_file_path = os.path.join(molecule.directory, log_file_name)
                molecule.log_file_path = log_file_path
                last_lines = read_last_lines(log_file_path, logger, 30)
                if not last_lines:
                    molecule.error_termination_count += 1
                    continue
                termination_detected = any(termination in line for line in last_lines)
                error_termination_detected = any(error_termination in line for line in last_lines)

                if termination_detected:
                    logger.log("termination_detected")
                    xyz_coordinates = log2xyz(molecule)
                    molecule.coordinates = xyz_coordinates
                    molecule.update_energy(logger)
                    if job_type == 'TS_opt':
                        if any(freq < freq_cutoff for freq in molecule.vibrational_frequencies):
                            negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if freq < freq_cutoff)
                            logger.log_with_stars(f"Yay {log_file_name} converged with imaginary frequency: {negative_freqs}")
                            molecule.converged = True
                            molecule.wrong_TS = 0
                            molecule.move_converged()
                        elif any(freq < freq_cutoff for freq in molecule.vibrational_frequencies):
                            negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if -freq_cutoff <= freq < 0)
                            logger.log(f"Small negative frequency {negative_freqs} found. This may be indicate wrong transition state. Resubmitting job")
                            molecule.converged = False
                            molecule.wrong_TS += 1
                            if molecule.wrong_TS >= 4:
                                continue
                            resubmit_job(molecule, logger, job_type)
                            time.sleep(interval / 2)
                            continue
                        else:
                            logger.log(f"However, no imaginary frequency found for {molecule.name}. Resubmitting TS calculation")
                            molecule.converged = False
                            molecule.wrong_TS += 1
                            if molecule.wrong_TS >= 4:
                                continue
                            resubmit_job(molecule, logger, job_type)
                            time.sleep(interval / 2)
                            continue
                    else:
                        molecule.converged = True
                        logger.log(f"**{molecule.name} converged**")
                        xyz_coordinates = log2xyz(molecule)
                        if xyz_coordinates:
                            molecule.coordinates = xyz_coordinates
                            molecule.move_converged()
                        else:
                            logger.log(f"Error gathering XYZ coordinates from {log_file_name}. Check log file.")
                elif error_termination_detected:
                    molecule.error_termination_count += 1
                    molecule.converged = False
                    logger.log(f"Error termination found in {log_file_name}.")
                    xyz_coordinates = log2xyz(molecule)
                    if xyz_coordinates:
                        molecule.coordinates = xyz_coordinates
                        logger.log(f"Trying to resubmit job for {molecule.name} due to error termination")
                        resubmit_job(molecule, logger, job_type)
                        time.sleep(interval/2)
                        continue
                else:
                    # logger.log("No error or normal termination found")
                    molecule.converged = False
                    continue
            elif status == "Pending":
                if molecule.job_id not in pending:
                    logger.log(f"Job {molecule.job_id} is pending in the queue.")
                    pending.append(molecule.job_id)
                continue
            else:
                # logger.log(f"Status of job {molecule.job_id} is unknown. Checking log file.")
                molecule.converged = False # tmp solution
                
        if all(m.converged for m in molecules):
            all_converged = True

        attempts += 1
        if attempts % 10 == 0 or attempts == 1:
            logger.log(f"Log files of the {len(molecules)} conformers have been checked. Checking every {interval} seconds. Attempt: {attempts}/{max_attempts}")
        time.sleep(interval)

    if all_converged:
        if molecules:
            dir = molecules[0].directory
            input_files = molecules[0].input
            basename = os.path.basename(dir)
            pickle_path = os.path.join(dir, f'{basename}_{job_type}.pkl')
            Molecule.molecules_to_pickle(molecules, pickle_path)
            move_files(dir, input_files)
            for m in molecules: m.converged = False
            logger.log_with_stars(f"Yay! All conformer jobs have converged for job type: {job_type}.")
            if job_type == "DLPNO":
                for molecule in molecules: 
                    molecule.update_step()
                    global_molecules.append(molecule)
                return True
            elif args.auto:
                logger.log(f"Next step is: {molecules[0].next_step}")
                logger.log("----------------------------------------------------------------------")
                converged_job(molecules, logger, threads)
                return True
            else:
                logger.log(f"Automatic processing of calculations has been turned off and JKTS will not proceed with next step: {molecules[0].next_step}")
                logger.log(f"A pickle file {basename}_{job_type}.pkl has been created. Next step {molecules[0].next_step} for molecules can be started from this.")
                return True
        else:
            logger.log(f"No conformer has managed to converge for job type: {job_type}")
            return False


def check_convergence(molecule, logger, threads, time_seconds, max_attempts):
    log_file_name = f"{molecule.name}{molecule.output}"
    directory = molecule.directory
    log_file_path = os.path.join(directory, log_file_name)
    molecule.log_file_path = log_file_path
    job_type = molecule.current_step

    if args.initial_delay: initial_delay = args.initial_delay
    else: initial_delay = int(time_seconds*2)
    if args.interval: interval = args.interval
    else: interval = int(time_seconds)
    attempts = 0
    running = 1
    pending = 1
    crest_running = 1

    if args.freq_cutoff:
        freq_cutoff = args.freq_cutoff
    elif args.OH:
        freq_cutoff = -150
    elif args.CC:
        freq_cutoff = -50 # TEST
    else:
        freq_cutoff = -50 # TEST

    molecule.converged = False
    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)
    while attempts < max_attempts and molecule.converged is False and molecule.error_termination_count <= 3:
        status = check_job_status(molecule.job_id)

        if status in ["Running", "Completed or Not Found"] or molecule.converged is False:
            if running:
                logger.log(f"Job {job_type} for {molecule.name} with job id {molecule.job_id} is {status}.")
                running = 0

            if job_type == 'crest_sampling':
                if crest_running:
                    logger.log("CREST pickle file not generated yet")
                    crest_running = 0
                if os.path.exists(log_file_path) or os.path.exists(os.path.join(directory, f"collection{molecule.name}.pkl")):
                    logger.log("CREST log file and pickle file generated")
                    molecule.converged = True
                    file_names = ["constrain.inp", f"{molecule.name}.pkl", f"{molecule.name}.xyz"]
                    for file_name in file_names:
                        file_path = os.path.join(directory, file_name)
                        if os.path.exists(file_path):
                            shutil.move(file_path, os.path.join(directory, "input_files", file_name))
                    xyz_conformers = pkl_to_xyz(os.path.join(molecule.directory, f"collection{molecule.name}.pkl"))
                    if xyz_conformers:
                        conformer_molecules = []
                        for n, conformer_coords in enumerate(xyz_conformers, start=1):
                            conformer_molecule = copy.deepcopy(molecule)
                            conformer_molecule.name = f"{molecule.name}_conf{n}"
                            conformer_molecule.atoms = [atom[0] for atom in conformer_coords]
                            conformer_molecule.coordinates = [atom[1:] for atom in conformer_coords]
                            conformer_molecules.append(conformer_molecule)
                        logger.log_with_stars(f"Yay CREST job {molecule.name} converged! {len(xyz_conformers)} conformers generated.")
                        if args.auto:
                            logger.log(f"Proceeding with next step: {molecule.next_step}")
                            molecule.move_converged()
                            converged_job(conformer_molecules, logger, threads)
                        else:
                            logger.log(f"Automatic processing of calculations has been turned off and JKTS will not proceed with next step: {molecule.next_step}")
                            logger.log(f"Logging conformers to {molecule.name}_conformers.pkl")
                            Molecule.molecules_to_pickle(conformer_molecules, os.path.join(molecule.directory, f"{molecule.name}_conformers.pkl"))
                        return True
                    else:
                        logger.log(f"ERROR: could not extract conformers from CREST pickle file")
                        return False
                else:
                    attempts += 1
                    logger.log(f"attempts: {attempts}/{max_attempts}")
                    time.sleep(interval)
                    continue

            last_lines = read_last_lines(log_file_path, logger, 30) 
            if not last_lines:
                molecule.error_termination_count += 1
                time.sleep(interval)
                continue
            termination = termination_strings.get(molecule.program.lower(), "")
            error_termination = error_strings.get(molecule.program.lower(), "")
            termination_detected = any(termination in line for line in last_lines)
            error_termination_detected = any(error_termination in line for line in last_lines)
            link1_detected = any("Link1:  Proceeding to internal job step number  2" in line for line in last_lines)
            if link1_detected:
                time.sleep(interval)
                continue
            if termination_detected:
                if molecule.name in ["H2O_DLPNO", "OH_DLPNO"] and job_type == 'DLPNO':
                    molecule.converged = True
                    molecule.current_step = 'Done'
                    molecule.update_energy(logger, DLPNO=True)
                    global_molecules.append(molecule)
                    logger.log_with_stars(f"Yay! DLPNO single point for {molecule.name} molecule converged")
                    molecule.move_converged()
                    return True
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
                            molecule.move_converged()
                            if args.auto:
                                logger.log(f"Proccessing with next step: {molecule.next_step}")
                                converged_job(molecule, logger, threads)
                            else:
                                logger.log(f"Automatic processing of calculations has been turned off and JKTS will not proceed with next step: {molecule.next_step}")
                            return True
                        elif any(freq_cutoff <= freq < 0 for freq in molecule.vibrational_frequencies):
                            negative_freqs = ', '.join(str(freq) for freq in molecule.vibrational_frequencies if -freq_cutoff <= freq < 0)
                            logger.log(f"Small negative frequency found in the range {freq_cutoff} to 0: {negative_freqs}. This may be indicate wrong transition state. Resubmitting job")
                            molecule.wrong_TS += 1
                            if molecule.wrong_TS == 3:
                                molecule.set_active_site(perturb=True)
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
                        molecule.move_converged()
                        if args.auto:
                            logger.log(f"Proccesing with next step: {molecule.next_step}")
                            converged_job(molecule, logger, threads) 
                        else:
                            logger.log(f"Automatic processing of calculations has been turned off and JKTS will not proceed with next step: {molecule.next_step}")
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
            
        elif status == "Pending":
            if pending:
                logger.log(f"Job {molecule.job_id} is pending in the queue.")
                pending = 0
        attempts += 1
        logger.log(f"Attempt: {attempts}/{max_attempts}. Checking every {interval} seconds")
        time.sleep(interval)

    logger.log("Max attempts reached or job status is not running. Calculation may be stuck or completed.")
    return False

def check_crest_products(molecules, logger, threads, time_seconds, max_attempts):
    if args.initial_delay:
        initial_delay = args.initial_delay
    else:
        initial_delay = int(time_seconds*3)
    if args.interval:
        interval = args.interval
    else:
        interval = int(time_seconds) 

    attempts = 0
    all_conformers = []
    expected_files = {f"collection{molecule.name}.pkl" for molecule in molecules}
    
    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)

    while attempts < max_attempts:
        try:
            files_in_directory = set(os.listdir(molecules[0].directory))
        except FileNotFoundError:
            logger.log(f"Pickle files not generated yet. Retrying in {interval} seconds.")
            time.sleep(interval)
            attempts += 1
            continue

        if expected_files.issubset(files_in_directory):
            for molecule in molecules:
                try:
                    conformers = pkl_to_xyz(os.path.join(molecule.directory, f"collection{molecule.name}.pkl"))
                    for n, conformer_coords in enumerate(conformers, start=1):
                        conformer_molecule = copy.deepcopy(molecule)
                        conformer_molecule.name = f"{molecule.name}_conf{n}"
                        conformer_molecule.mult = 2
                        conformer_molecule.product = True
                        conformer_molecule.atoms = [atom[0] for atom in conformer_coords]
                        conformer_molecule.coordinates = [atom[1:] for atom in conformer_coords]
                        all_conformers.append(conformer_molecule)
                except Exception as e:
                    logger.log(f"Error processing molecule {molecule.name}: {e}")
                    return False
            converged_job(all_conformers, logger, threads)
            return True
        else:
            logger.log(f"Not all files found. Retrying in {interval} seconds.")
            time.sleep(interval)

        attempts += 1

    return False


def submit_and_monitor(molecule, logger, threads):
    if isinstance(molecule, list):
        if args.par == 'qtest':
            job_time = '1:00:00'
            time_seconds = 3600 / args.attempts
        else:
            time_seconds = get_time(molecule[0])
            total_seconds = time_seconds * (args.attempts + 5) # Taking into account initial_delay
            job_time = seconds_to_hours(total_seconds)
        job_id = submit_array_job(molecule, args.par, job_time, args.cpu, args.mem)
        logger.log(f"Submitted SLURM array job with job id {job_id} for conformers in {molecule[0].directory}")

        if job_id is not None:
            if molecule[0].current_step == 'crest_sampling':
                thread = threading.Thread(target=check_crest_products, args=(molecule, logger, threads, time_seconds, args.attempts))
            else:
                thread = threading.Thread(target=check_convergence_for_conformers, args=(molecule, logger, threads, time_seconds, args.attempts))
            threads.append(thread)
            thread.start()
        else: logger.log("Error getting job id")

    else:
        if args.par == 'qtest':
            job_time = '1:00:00'
            time_seconds = 3600 / args.attempts
        else:
            time_seconds = get_time(molecule)
            total_seconds = time_seconds * (args.attempts + 5)
            job_time = seconds_to_hours(total_seconds)
        job_id = submit_job(molecule, args.cpu, args.mem, args.par, job_time)
        logger.log(f"Submitting file {molecule.name}{molecule.input} for calculation in path {molecule.directory} with id {job_id}")

        if job_id is not None:
            thread = threading.Thread(target=check_convergence, args=(molecule, logger, threads, time_seconds, args.attempts))
            threads.append(thread)
            thread.start()
        else: logger.log("Error getting job id")



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
            molecule.file_path = output_file_path
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=True)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'DLPNO':
            molecule.name.replace("_TS", "")
            molecule.name += "_DLPNO"
            molecule.program = 'ORCA'
            QC_input(molecule, constrain=False, method='DLPNO', basis_set=high_basis, TS=False)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'optimization': #Usually for OH radical and H2O
            molecule.name.replace("_CREST", "")
            molecule.program = global_program
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.output}")
            molecule.file_path = output_file_path
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=False)
            submit_and_monitor(molecule, logger, threads)

        elif job_type == 'Done':
            logger.log(f"DLPNO single point calculation has been done for {molecule.name}")
        else:
            pass

    elif isinstance(molecule, list):
        conformer_molecules = molecule
        for conf in conformer_molecules:
            conf.update_step()
            job_type = conf.current_step
            conf.program = global_program

            if job_type == 'opt_constrain': 
                QC_input(conf, constrain=True, method=low_method, basis_set=low_basis, TS=False)

            elif job_type == 'opt':
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=False)
            
            elif job_type == 'TS_opt':
                conf.name += '_TS'
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=True)

            elif job_type == 'DLPNO':
                conf.program = 'ORCA'
                conf.name = conf.name.replace("_TS","")
                conf.name += '_DLPNO'
                QC_input(conf, constrain=False, method='DLPNO', basis_set=high_basis, TS=False)

            elif job_type  == 'Done':
                logger.log("DLPNO calculation has converged")
                break

            else:
                logger.log("Job type could not be determined for conformer list")
                break
        if conformer_molecules: 
            submit_and_monitor(conformer_molecules, logger, threads)

    else:
        logger.log(f"Error: next step could not be determined for molecule {molecule.name}")


def non_converged_job(molecule, logger, threads):
    if not isinstance(molecule, list):
        job_type = molecule.current_step

        if job_type == 'crest_sampling':
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}.xyz")
            molecule.write_xyz_file(output_file_path)

        elif job_type == 'TS_opt':
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.output}")
            molecule.file_path = output_file_path
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=True)

        elif job_type == 'DLPNO':
            molecule.program = 'ORCA'
            QC_input(molecule, constrain=False, method='DLPNO', basis_set=high_basis, TS=False)

        elif job_type == 'optimization': #Usually for OH radical and H2O
            molecule.program = global_program
            output_file_path = os.path.join(molecule.directory, f"{molecule.name}{molecule.output}")
            molecule.file_path = output_file_path
            QC_input(molecule, constrain=False, method=high_method, basis_set=high_basis, TS=False)

        if job_type:
            submit_and_monitor(molecule, logger, threads)

    elif isinstance(molecule, list):
        conformer_molecules = molecule
        job_type = conformer_molecules[0].current_step
        for conf in conformer_molecules:
            conf.program = global_program

            if job_type == 'opt_constrain': 
                QC_input(conf, constrain=True, method=low_method, basis_set=low_basis, TS=False)

            elif job_type == 'opt':
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=False)
            
            elif job_type == 'TS_opt':
                QC_input(conf, constrain=False, method=high_method, basis_set=high_basis, TS=True)

            elif job_type == 'DLPNO':
                conf.program = 'ORCA'
                QC_input(conf, constrain=False, method='DLPNO', basis_set=high_basis, TS=False)

            elif job_type  == "Done":
                logger.log("DLPNO calculation has converged")
                break

            else:
                logger.log("Job type could not be determined for conformer list")
                break
        if job_type: 
            submit_and_monitor(conformer_molecules, logger, threads)

    else:
        logger.log(f"Error: next step could not be determined for molecule {molecule.name}")


def determine_current_step(molecule):
    if molecule.program.lower() == 'crest' or 'collection' in molecule.name:
        return 'crest_sampling'

    if os.path.exists(molecule.log_file_path):
        with open(molecule.log_file_path, 'r') as f:
            for line in f:
                # Adjusted for ORCA log files which may have additional characters before the '!'
                if molecule.program.lower() == 'orca' and '|  1> !' in line:
                    line_content = line.split('! ', 1)[-1]  # Get content after '!'
                elif molecule.program.lower() == 'g16' and line.strip().startswith('#'):
                    line_content = line
                else:
                    continue  # Skip lines that don't match the criteria

                if 'ts' in line_content.lower() and molecule.program.lower() == 'g16':
                    return 'TS_opt'
                elif 'optts' in line_content.lower() and molecule.program.lower() == 'orca':
                    return 'TS_opt'
                elif 'opt' in line_content.lower():
                    return 'opt_constrain'
                elif 'dlpno' in line_content.lower() and molecule.program.lower() == 'orca':
                    return 'DLPNO'
                else:
                    return False

    else:
        if molecule.reactant or molecule.product:
            if 'OH' in molecule.name or 'H2O' in molecule.name:
                if 'DLPNO' in molecule.name:
                    return 'DLPNO'
                else:
                    return 'optimization'
            elif 'CREST' in molecule.name:
                return 'opt'
            elif 'DLPNO' in molecule.name:
                return 'DLPNO'
        # For TS molecules
        else:
            if 'TS' in molecule.name:
                return 'TS_opt'
            elif 'DLPNO' in molecule.name:
                return 'DLPNO'
            else:
                return 'opt_constrain'


def termination_status(molecule):
    termination = termination_strings.get(molecule.program, "")
    count = 0
    if molecule.current_step == 'TS_opt' and molecule.program.lower() == 'g16':
        required_count = 2
    else: required_count = 1
    with open(molecule.log_file_path, 'r') as f:
        for line in f:
            if termination in line:
                count += 1
                if count >= required_count:
                    return True
        return False


def generate_conformers(conformers, input_file):
    for n, conformer_coords in enumerate(conformers, start=1):
        name = input_file.split(".")[0].replace("collection", "")
        name += f"_conf{n}"
        conformer_molecule = Molecule(name=name,
        directory=start_dir,
        log_file_path=os.path.join(start_dir, f"{name}.log"),
        atoms=[atom[0] for atom in conformer_coords],
        coordinates=[atom[1:] for atom in conformer_coords])
        conformer_molecule.find_active_site()
        if 'reactant' in input_file:
            conformer_molecule.reactant = True
            conformer_molecule.mult = 1
        elif 'product' in input_file:
            conformer_molecule.product = True
        if 'collection' in input_file:
            conformer_molecule.current_step = 'crest_sampling'
        elif 'TS' in input_file:
            conformer_molecule.current_step = 'TS_opt'
        elif 'DLPNO' in input_file:
            conformer_molecule.current_step = 'DLPNO'
        else:
            conformer_molecule.current_step = 'opt_constrain'
    return conformers


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
    
def eckart(SP_TS, SP_reactant, SP_product, imag):
    tunneling_coefficient = 1.0
    return tunneling_coefficient


def rate_constant(TS_conformers, reactant_conformers, product_conformers, T=298.15):
    k_b = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck constant in J*s
    HtoJ = 4.3597447222e-18  # Conversion factor from Hartree to Joules
    Na = 6.022e23 # molecules/mol
    liters_to_cm3 = 1000 # 1 liter is 1000 cm^3
    
    # # Molecule lists
    # reactant_molecules = Molecule.load_molecules_from_pickle(reactant_molecules_path)
    # product_molecules = Molecule.load_molecules_from_pickle(product_molecules_path)
    # TS_molecules = Molecule.load_molecules_from_pickle(TS_molecules_path)

    # Calculation of rate constant
    if reactant_conformers and TS_conformers:
        # Seperate organic molecule and OH
        reactant_molecules = [mol for mol in reactant_conformers if 'OH' not in mol.name]
        OH = next((mol for mol in reactant_conformers if 'OH' in mol.name), None)

        if reactant_molecules:
            # Find the lowest single point energies for reactant (excluding OH) and TS
            lowest_reactant = min(reactant_molecules, key=lambda molecule: molecule.single_point)
            lowest_TS = min(TS_conformers, key=lambda molecule: molecule.single_point)

            # Lowest single point reactant and TS and convert from Hartree to Joules
            lowest_SP_TS_J = np.float64(lowest_TS.single_point * HtoJ)
            lowest_SP_reactant_J = np.float64(lowest_reactant.single_point * HtoJ)
            if OH:
                total_reactant_SP_J = lowest_SP_reactant_J + np.float64(OH.single_point*HtoJ)
                sum_reactant = np.sum([np.exp((lowest_SP_reactant_J - np.float64(mol.single_point * HtoJ)) / (k_b * T)) * np.float64(mol.Q) for mol in reactant_molecules])*OH.Q
            else:
                total_reactant_SP_J = lowest_SP_TS_J
                sum_reactant = np.sum([np.exp((lowest_SP_reactant_J - np.float64(mol.single_point * HtoJ)) / (k_b * T)) * np.float64(mol.Q) for mol in reactant_molecules])        

            sum_TS = np.sum([np.exp((lowest_SP_TS_J - np.float64(mol.single_point * HtoJ)) / (k_b * T)) * np.float64(mol.Q) for mol in TS_conformers])

            if product_conformers:
                # Seperate organic molecule and H2O
                product_molecules = [mol for mol in reactant_conformers if 'H2O' not in mol.name]
                H2O = next((mol for mol in product_conformers if 'H2O' in mol.name), None)

                if product_molecules:
                    lowest_product = min(product_molecules, key=lambda molecule: molecule.single_point)
                    lowest_SP_product_J = np.float64(lowest_product.single_point * HtoJ)

                    if H2O:
                        total_product_SP_J = lowest_SP_product_J + np.float64(H2O.single_point*HtoJ)
                    else:
                        total_product_SP_J = lowest_SP_product_J

                    imag = lowest_TS.vibrational_frequencies[0]
                    tunneling_coefficient = eckart(lowest_SP_TS_J, total_reactant_SP_J, total_product_SP_J, imag)
                    k = tunneling_coefficient * (k_b*T/h) * (sum_TS / (sum_reactant)) * np.exp(-total_reactant_SP_J / (k_b * T))
                else:
                    print("Error in product molecules")
                    return None
            else: # If products are not calculated assume tunneling coefficient is 1
                tunneling_coefficient = 1.0
                k = tunneling_coefficient * (k_b*T/h) * (sum_TS / (sum_reactant)) * np.exp(-(lowest_SP_TS_J-total_reactant_SP_J) / (k_b * T))
           
            k_mol_cm3_s = k/(Na/liters_to_cm3)
            return k_mol_cm3_s

    return None


def seconds_to_hours(seconds):
    seconds = round(seconds)
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    return f"{hours:02d}:{minutes:02d}:00"


def main():
    parser = argparse.ArgumentParser(description='''   -Dynamic Approach for Transition State-
    Automated tool for generating input files, primarily 
    for transition state geometry optimization. 
    Calculation of tunneling corrected multi-configurational 
    rate constants can also be calculated from log files.''',
                                     prog="JKTS",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
    Examples of use:
                JKTS pinonaldehyde.xyz -OH -auto false
                JKTS CH4.xyz -OH -auto --low_method "am1 3-21g" --high_method "B3LYP 6-31+g(d,p)"
                JKTS CH4_H1_opt_constrain.pkl -info
                JKTS benzene.xyz -OH -ORCA -par qtest -auto false
                JKTS *TS.log -time 5:00:00
                                     ''')


    parser.add_argument('input_files', metavar='reactant.xyz', nargs='+', help='Input XYZ files (e.g., pinonaldehyde.xyz or *.xyz)')

    reaction_options = parser.add_argument_group("Types of reactions")

    reaction_options.add_argument('-OH', action='store_true', help='Perform H abstraction with OH radical')
    reaction_options.add_argument('-CC', action='store_true', help='Perform addition to C=C bonds')
    reaction_options.add_argument('-OH_CC', action='store_true', help='Perform OH addition to C=C bonds')

    parser.add_argument('-G16', action='store_true', help='Use Gaussian16 for QC calculations (default)')
    parser.add_argument('-ORCA', action='store_true', help='Use ORCA for QC calculations')
    parser.add_argument('-constrain', action='store_true', default=True, help='Integrate constraints into relevant input files [def: True]')
    parser.add_argument('-reactants', type=str2bool, default=True, nargs='?',metavar='', const=True, help='Prepare folder for reactants [def: True]')
    parser.add_argument('-products', type=str2bool, default=True, nargs='?', const=True, metavar='', help='Prepare folder for products [def: True]')
    parser.add_argument('-NEB', type=str2bool, default=False, nargs='?', const=False, metavar='', help='Prepare input file for Nudged Elastic Band')
    parser.add_argument('-auto', type=str2bool, default=True, nargs='?', const=True, metavar='', help='Automated process with the following workflow:\n- CREST conformer sampling of xyz input file (GFN2-xTB -ewin=8 kcal/mol)\n- Preoptimization of geometry (Low Level of Theory)\n- Optimization towards transition state (High Level of Theory)\n- DLPNO-CCSD(T) SP energy calculations on top of TS conformers\n- Calculate rate constants and branching ratios for reaction type\n* Automatically resubmit failed calculations until convergence or wall-time limit is reached')

    additional_options = parser.add_argument_group("Additional arguments")
    additional_options.add_argument('-init', action='store_true', help=argparse.SUPPRESS)
    additional_options.add_argument('-k', action='store_true', default=True, help='Calculate Multiconformer Transition State rate constant [def: True]')
    additional_options.add_argument('-info', action='store_true', default=False, help='Print information of molecules in log files or .pkl file')
    additional_options.add_argument('--high_level', nargs='+', metavar='', help='Specify high-level theory for QC method TS optimization [def: uwB97X-D aug-cc-pVTZ]')
    additional_options.add_argument('--low_level', nargs='+', metavar='', help='Specify low-level theory for preoptimization [def: B3LYP 6-31+G(d,p)]')
    additional_options.add_argument('-cpu', metavar="int", nargs='?', const=1, type=int, default=4, help='Number of CPUs [def: 4]')
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=8000, help='Amount of memory allocated for the job [def: 8000MB]')
    additional_options.add_argument('-par', metavar="partition", nargs='?', const=1, default="q24,q28,q36,q40,q48,q64", help='Partition to use [def: qany]')
    additional_options.add_argument('-time', metavar="hh:mm:ss", type=str, default=None, help='Monitoring duration [def: 144 hours]')
    additional_options.add_argument('-interval', metavar="int", nargs='?', const=1, type=int, help='Time interval between log file checks [def: based on molecule size]')
    additional_options.add_argument('-initial_delay', metavar="int", nargs='?', const=1, type=int, help='Initial delay before checking log files [def: based on molecule size]')
    additional_options.add_argument('-attempts', metavar="int", nargs='?', const=1, type=int, default=100, help='Number of log file check attempts [def: 100]')
    additional_options.add_argument('-max_conformers', metavar="int", nargs='?', const=1, type=int, help='Maximum number of conformers from CREST [def: 50]')
    additional_options.add_argument('-freq_cutoff', metavar="int", nargs='?', const=1, type=int, help='TS imaginary frequency cutoff [def: -200 cm^-1]')
    additional_options.add_argument('-reaction_angle', metavar="float", nargs='?', default=175.0, const=1, type=float, help='Active site transition state angle [def: 175.0 degrees]')
    additional_options.add_argument('-ewin', metavar="int", nargs='?', const=1, default=8, type=int, help='Energy threshold for CREST conformer sampling [def: 8 kcal/mol]')
    additional_options.add_argument('-energy_cutoff', metavar="int", nargs='?', const=1, default=5, type=int, help='After preoptimization, remove conformers which are [int] kcal/mol higher in energy than the lowest conformer [def: 5 kcal/mol]')

    additional_options.add_argument('-test', action='store_true', default=False, help=argparse.SUPPRESS) # Run TEST section in main() and exit
    additional_options.add_argument('-num_molecules', type=int, default=None, help=argparse.SUPPRESS) # Set the number of directories created to [int]. Used when doing limited testing
    additional_options.add_argument('-XQC', action=SCFAction, nargs=0, dest='SCF')
    additional_options.add_argument('-YQC', action=SCFAction, nargs=0, dest='SCF')
    additional_options.add_argument('-QC', action=SCFAction, nargs=0, dest='SCF')
    parser.set_defaults(SCF="")

    global args, start_dir
    args = parser.parse_args()
    start_dir = os.getcwd()

    ################################ARGUMENT SPECIFICATIONS############################

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

    high_method, high_basis = extract_method_basis(args.high_level, ["uwb97xd", "aug-cc-pVTZ"])
    low_method, low_basis = extract_method_basis(args.low_level, ["uwb97xd", "6-31+g(d,p)"])

    global global_program
    if args.ORCA:
        global_program = "ORCA"
        if low_method in methods_no_basis:
            low_basis = ""
        if high_method in methods_no_basis:
            high_basis = ""
        if high_method.lower() == "uwb97xd":
            high_method = "WB97X-D3"
        if low_method.lower() == "uwb97xd":
            low_method = "WB97X-D3" 
    elif args.G16:
        global_program = "G16"
    else: 
        global_program = "G16" # Default case

    global termination_strings, error_strings
    termination_strings = {
        "g16": "Normal termination",
        "orca": "****ORCA TERMINATED NORMALLY****",
        "crest": " CREST terminated normally."
    }
    error_strings = {
        "g16": "Error termination",
        "orca": "error",
        "crest": "Find my error message"
    }


    #####################################################################################################
    # if args.test:
    #     print(args)
    #     exit()
    #     molecules = []
    #     for n, input_file in enumerate(args.input_files, start=1):
    #         input_file = args.input_files[0]
    #         file_name, file_type = os.path.splitext(input_file)
    #         input_file_path = os.path.join(start_dir, input_file)
    #         molecule = Molecule(input_file_path)
    #         molecule.print_items()
    #
    #     exit()
    #####################################################################################################

    threads = []
    input_molecules = []

    if args.init:
        input_file = args.input_files[0]
        file_name, file_type = os.path.splitext(input_file)
        input_file_path = os.path.join(start_dir, input_file)
        input_molecule = Molecule(input_file_path)

        if args.OH:
            reacted_molecules, product_molecules = input_molecule.H_abstraction(NEB=args.NEB, num_molecules=args.num_molecules)
        elif args.CC:
            other_molecule = args.input_files[1]
            reacted_molecules = input_molecule.addition(other_molecule)
        else:
            parser.error("Need to specify reaction type")

        for count, molecule in enumerate(reacted_molecules, start=1):
            molecule.name = f"{file_name}_H{count}" 
            molecule.directory = os.path.join(start_dir, molecule.name)
            molecule.current_step = 'crest_sampling'
            molecule.program = 'CREST'
            mkdir(molecule)
            logger = Logger(os.path.join(molecule.directory, "log"))
            # molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz")) Implement at some point just the writing of the xyz file instead of pkl file
            molecule.save_to_pickle(os.path.join(molecule.directory, f"{molecule.name}.pkl"))

        if args.reactants: 
            reactant_dir = os.path.join(start_dir, 'reactants')
            input_molecule.name = f"{file_name}_reactant"
            input_molecule.directory = reactant_dir
            input_molecule.reactant = True
            input_molecule.mult = 1
            input_molecule.program = 'CREST'
            input_molecule.current_step = 'crest_sampling'
            mkdir(input_molecule, CREST=False)
            input_molecule.save_to_pickle(os.path.join(input_molecule.directory, f"{input_molecule.name}.pkl"))

        if args.products and product_molecules:
            product_dir = os.path.join(start_dir, 'products')
            for count, molecule in enumerate(product_molecules, start=1):
                molecule.name = f"{file_name}_product_H{count}"
                molecule.directory = product_dir
                molecule.current_step = 'crest_sampling'
                molecule.program = 'CREST'
            mkdir(product_molecules[0], CREST=False)
            pickle_path = os.path.join(product_dir, "product_molecules.pkl")
            Molecule.molecules_to_pickle(product_molecules, pickle_path)
        else:
            parser.error("Error when handling product directory") # Tmp fix # Tmp fix

    elif args.info:
        for input_file in args.input_files:
            file_name, file_type = os.path.splitext(input_file)
            if file_type == '.pkl':
                molecules = Molecule.load_molecules_from_pickle(input_file)
                if isinstance(molecules, pd.DataFrame):
                    conformers = pkl_to_xyz(input_file)
                    conformers = generate_conformers(conformers, input_file)
                    for conf in conformers:
                        conf.print_items()
                else:
                    if isinstance(molecules, list):
                        for molecule in molecules:
                            molecule.print_items()
                    else: molecules.print_items()
            elif file_type in ['.log', '.out']:
                input_file_path = os.path.join(start_dir, input_file)
                molecule = Molecule(input_file_path)
                molecule.print_items()
            else:
                parser.error("Invalid input file format. The program expects input files with extensions '.pkl', '.log', or '.out'. Please ensure that your input file is in one of these formats and try again. If you provided a different type of file, convert it to a supported format or select an appropriate file for processing.")

    else: # We loop over all molecules given in the argument input and process them according to file type
        for n, input_file in enumerate(args.input_files, start=1):
            file_name, file_type = os.path.splitext(input_file)

            if file_type == '.xyz':
                if os.path.basename(start_dir) == 'reactants':
                    input_file_path = os.path.join(start_dir, f"{file_name}_reactant.pkl")
                    logger = Logger(os.path.join(start_dir, "log"))
                    OH = Molecule(name='OH')
                    if args.high_method is None and args.high_basis is None:
                        global_molecules.append(OH)
                    else:
                        QC_input(OH, constrain=False, method=high_method, basis_set=high_basis, TS=False)
                        submit_and_monitor(OH, logger, threads)

                    reactant = Molecule.load_from_pickle(input_file_path)
                    reactant.write_xyz_file(os.path.join(reactant.directory, f"{reactant.name}.xyz"))
                    submit_and_monitor(reactant, logger, threads)

                elif os.path.basename(start_dir) == 'products':
                    logger = Logger(os.path.join(start_dir, "log"))
                    H2O = Molecule(name='H2O')
                    if args.high_method is None and args.high_basis is None:
                        global_molecules.append(H2O)
                    else:
                        QC_input(H2O, constrain=False, method=high_method, basis_set=high_basis, TS=False)
                        submit_and_monitor(H2O, logger, threads)

                    pickle_path = os.path.join(start_dir, "product_molecules.pkl")
                    product_molecules = Molecule.load_molecules_from_pickle(pickle_path)
                    for product in product_molecules:
                        product.write_xyz_file(os.path.join(start_dir, f"{product.name}.xyz"))
                        product.product = True
                        product.current_step = 'crest_sampling'
                    submit_and_monitor(product_molecules, logger, threads)

                else:
                    input_file_path = os.path.join(start_dir, os.path.basename(start_dir)+'.pkl')
                    molecule = Molecule.load_from_pickle(input_file_path)
                    logger = Logger(os.path.join(start_dir, "log"))
                    molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz"))
                    submit_and_monitor(molecule, logger, threads)

    ### If input files are pickle or log files we collect them into input_molecules and process the list, by checking each input convergence status
            elif file_type == '.pkl':
                df = pd.read_pickle(input_file)
                if isinstance(df, pd.DataFrame): # If pkl file is from the CREST job it will come as a pandas dataframe, instead of a list
                    conformer_molecules = pkl_to_xyz(input_file)
                    conformer_molecules = generate_conformers(conformer_molecules, input_file)
                    input_molecules = conformer_molecules
                else:
                    input_molecules = Molecule.load_from_pickle(input_file)

            elif file_type in [".log", ".out"]: 
                # Initialize molecule from log file
                input_file_path = os.path.join(start_dir, input_file)
                molecule = Molecule(file_path=input_file_path)
                input_molecules.append(molecule)


        logger = Logger(os.path.join(start_dir, "log"))
        if input_molecules:
            current_step = input_molecules[0].current_step
            if all(m.current_step == current_step for m in input_molecules):
                logger.log(f"Detected molecules with current step: {current_step}")
                logger.log("Checking convergence status of molecules")
                for m in input_molecules:
                    m.directory = start_dir
                    if os.path.exists(m.log_file_path):
                        converge_status = termination_status(m)
                        if current_step == 'Done' and converge_status is True:
                            global_molecules.append(m)
                        elif converge_status is True:
                            m.converged = True
                        else:
                            m.converged = False
                    else: # in the case a single .pkl file has been given as input, and logfiles are not located in same directory:
                        # Check if the conformers in the pkl file have atoms and coordinates and send for calculation.
                        if m.atoms and m.coordinates and file_type == '.pkl':
                            m.converged = True
                if all(m.converged for m in input_molecules):
                    logger.log(f"All molecules converged for step {current_step}. Proceeding with next step: {input_molecules[0].next_step}")
                    converged_job(input_molecules, logger, threads)
                else:
                    logger.log(f"Not all molecules given have converged. Non-converged will be calculated and converged ones will be skipped.")
                    non_converged_job(input_molecules, logger, threads)
            else:
                logger.log(f"Not all molecules in the given files are of the same type. Resubmitting from the last common step.")
                # Implement way to find last common step and also log it so user can see
        else:
            logger.log("Error when generating input molecules")


            # if molecule.current_step == 'DLPNO':
            #     n = re.search(r'conf(\d+)', input_file)
            #     if n or 'OH' in molecule.name or 'H2O' in molecule.name:
            #         if molecule.reactant or molecule.product:
            #             if 'OH' in molecule.name: thermo_data_log = "OH.log"
            #             elif 'H2O' in molecule.name: thermo_data_log = "H2O.log"
            #             else: thermo_data_log = re.sub(r"_conf\d+_DLPNO", f"_conf{int(n.group(1))}", input_file) 
            #         else:
            #             thermo_data_log = re.sub(r"_conf\d+_DLPNO", f"_conf{int(n.group(1))}_TS", input_file)
            #         thermo_data_path = os.path.join(start_dir, thermo_data_log)            
            #         thermo_program = log2program(thermo_data_path)
            #         molecule.update_energy(log_logger, log_file_path=thermo_data_path, program=thermo_program) # first update thermochemistry
            #         molecule.update_energy(log_logger, DLPNO=True) # Then update single point energy for DLPNO
            #         global_molecules.append(molecule)
            #     continue
            

        
    # Monitor and handle convergence of submitted jobs
    while threads:
        for thread in list(threads):
            thread.join(timeout=0.1)
            if not thread.is_alive():
                threads.remove(thread)

    if global_molecules:
        logger = Logger(os.path.join(start_dir, "log"))
        molecules_logger = Logger(os.path.join(start_dir, "molecules.txt"))
        for molecule in global_molecules: 
            molecule.log_items(molecules_logger)

        if any(m.reactant for m in global_molecules):
            logger.log("Final DLPNO calculations for reactants is done. Logging properties to molecules_reactants.pkl")
            Molecule.molecules_to_pickle(global_molecules, os.path.join(start_dir, "molecules_reactants.pkl"))
        else:
            logger.log("Final DLPNO calculations for reactants is done. Logging properties to molecules_TS.pkl")
            TS_pkl_path = os.path.join(start_dir, "molecules_TS.pkl")
            Molecule.molecules_to_pickle(global_molecules, TS_pkl_path)

            if args.k:
                reactant_pkl_path = os.path.join(os.path.dirname(start_dir), 'reactants/molecules_reactants.pkl')
                product_pkl_path = os.path.join(os.path.dirname(start_dir), 'products/molecules_reactants.pkl')
                if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                    k = rate_constant(TS_pkl_path, reactant_pkl_path, product_pkl_path)
                    results_logger = Logger(os.path.join(os.path.dirname(start_dir), "results_log"))
                    results_logger.log_with_stars(f"{k} molecules cm^-3 s^-1")
                else:
                    reactant_pkl_path = os.path.join(start_dir, 'reactants/molecules_reactants.pkl')
                    product_pkl_path = os.path.join(start_dir, 'products/molecules_reactants.pkl')
                    if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                        k = rate_constant(TS_pkl_path, reactant_pkl_path, product_pkl_path)
                        results_logger = Logger(os.path.join(os.path.dirname(start_dir), "results_log"))
                        results_logger.log_with_stars(f"{k} molecules cm^-3 s^-1")
                    else:
                        logger.log("Done")



if __name__ == "__main__":
    main()
