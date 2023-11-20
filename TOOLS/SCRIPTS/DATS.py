#!/usr/bin/env python3
'''Dynamic Approach for Transition State'''
###############################LIBRARIES#####################################
import numpy as np
import pandas as pd
import argparse
import os
import shutil
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

maxtasks = 100

###########################VECTOR MANIPULATION################################
def calculate_vector(coord1, coord2):
    return np.array(coord2) - np.array(coord1)

def vector_length(vector):
    return np.linalg.norm(vector)

def atom_distance(atom1, atom2):
    return np.linalg.norm(np.array(atom2) - np.array(atom1))

def normalize_vector(vector):
    norm = np.linalg.norm(vector)
    if norm < 1e-8: 
        return np.zeros_like(vector)
    return vector / norm

def rotate_vector(vector, axis, angle):
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    cross_product = np.cross(axis, vector)
    rotated_vector = (vector * cos_theta +
                      cross_product * sin_theta +
                      axis * np.dot(axis, vector) * (1 - cos_theta))
    return rotated_vector

class Logger:
    def __init__(self, log_file):
        self.log_file = log_file

    def log(self, message):
        with open(self.log_file, 'a') as file:
            file.write(message + '\n')


#########################################SUBMIT JOBS############################
def submit_job(dir, input_file, job_program, ncpus, mem, partition, time, nnodes=1):
    job_name = os.path.splitext(os.path.basename(input_file))[0]
    file_extension = os.path.splitext(input_file)[1].lower()  # Extract the file extension and convert to lower case

    if job_program.lower() == "orca" or job_program.lower() == "g16":
        program_mem = mem + 2000
    else:
        program_mem = mem

    pwd = os.path.join(os.getcwd(), dir)
    path_submit_script = os.path.join(pwd, f"{job_name}_submit.sh")
    submit_file = os.path.join(pwd, "qsub.tmp")

    with open(path_submit_script, 'w') as file:
        file.write("#!/bin/bash\n")
        file.write(f"SUBMIT={submit_file}\n\n") 
        file.write("cat > $SUBMIT <<!EOF\n")
        file.write("#!/bin/sh\n")
        file.write(f"#SBATCH --job-name={job_name}\n")
        file.write(f"#SBATCH --nodes={nnodes}\n")
        file.write(f"#SBATCH --cpus-per-task={ncpus}\n")
        file.write(f"#SBATCH --ntasks={nnodes}\n")
        file.write(f"#SBATCH --error={pwd}/{job_name}_%j.err\n")
        file.write(f"#SBATCH --output={pwd}/{job_name}_%j.out\n")
        file.write(f"#SBATCH --time={time}\n")
        file.write(f"#SBATCH --partition={partition}\n")
        file.write(f"#SBATCH --no-requeue\n")
        file.write(f"#SBATCH --mem={program_mem}\n\n")

        if job_program.lower() == "g16" or file_extension == '.com':
            file.write("mkdir /scratch/\$SLURM_JOB_ID\n\n")
            file.write(f"cd {pwd}\n")
            file.write("export GAUSS_SCRDIR=/scratch/\$SLURM_JOB_ID\n\n")
            file.write(f"srun \$(which g16) {input_file} > {job_name}.log\n")
        elif job_program.lower() == "orca" or file_extension == '.inp':
            file.write("source /comm/groupstacks/chemistry/bin/modules.sh\n")
            file.write("ml orca/5.0.4\n")
            file.write("SCRATCH=/scratch/\$SLURM_JOB_ID\n")
            file.write("mkdir -p \$SCRATCH || exit $?\n")
            file.write("cd \$SCRATCH\n")
            file.write(f"cp {pwd}/{job_name}.inp .\n")
            file.write(f"\$(which orca) {job_name}.inp > {pwd}/{job_name}.log\n")
        elif job_program.lower() == "crest" or file_extension == '.xyz':
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

        file.write(f"NCPUS={ncpus}\n")
        file.write(f"NNODES={nnodes}\n\n")

        file.write(f"cat <<!EOF  | $submit\n")
        file.write("#!/bin/sh\n")
        file.write(f"#SBATCH --job-name=$JOB\n")
        file.write(f"#SBATCH --nodes=$NNODES\n")
        file.write(f"#SBATCH --cpus-per-task=$NCPUS\n")
        file.write(f"#SBATCH --ntasks=$NNODES\n")
        file.write(f"#SBATCH --error={dir}/{job_name}_%j.err\n")
        file.write(f"#SBATCH --output={dir}/{job_name}_%j.out\n")
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
            file.write("export GAUSS_SCRDIR=\$GSCR\n\n")

            file.write("source /comm/groupstacks/gaussian/bin/modules.sh\n")
            file.write("ml gaussian16/Rev.B.01\n")
            file.write("ml gcc/9.2.0\n")
            file.write("ml openmpi/4.0.1\n\n")

            file.write('GJ=\$(awk "NR == \$SLURM_ARRAY_TASK_ID" $IN)\n')
            file.write('LOG=\${GJ%.*}.log\n\n')
            file.write("srun $(which g16) \$GJ > \$LOG\n")
            file.write("#\n")
            file.write("!EOF")

        elif job_program.lower == "orca":
            file.write("  ulimit -c 0\n")
            file.write("source /comm/groupstacks/gaussian/bin/modules.sh\n")
            file.write("ml orca/5.0.4\n\n")


            file.write("SCRATCH=/scratch/\${SLURM_JOB_ID}/\${SLURM_ARRAY_TASK_ID}\n")
            file.write("mkdir -p \$SCRATCH\n\n")

            file.write("cd \$SCATCH\n")
            file.write("cp \$SLURM_SUBMIT_DIR/array.txt .\n")
            file.write("cp \$SLURM_SUBMIT_DIR/*.inp .\n\n")

            file.write('GJ=\$(awk "NR == \$SLURM_ARRAY_TASK_ID" $IN)\n')
            file.write('LOG=\${GJ%.*}.log\n\n')

            file.write("\$(which orca) \$GJ > \$SLURM_SUBMIT_DIR/\$LOG\n\n")
            file.write("!EOF")

    subprocess.run(['sh', path_submit_script, array_path])


def check_convergence(log_file_name, directory, logger, shared_data, job_type, job_program, initial_delay=300, interval=90, max_attempts=100):
    if job_program.lower() == "g16":
        termination = "Normal termination"
    elif job_program.lower() == "orca":
        termination = "****ORCA TERMINATED NORMALLY****"
    elif job_program.lower() == "crest":
        termination = "CREST terminated normally"
    else:
        logger.log("Invalid program specified. Unable to check convergence of log file")
        shared_data[log_file_name]['result'] = (None, None, False)
        return False

    log_file_path = os.path.join(directory, log_file_name)
    logger.log(f"Waiting for {initial_delay} seconds before first check.")
    time.sleep(initial_delay)
    attempts = 0
    while attempts < max_attempts:
        try:
            with open(log_file_path, 'r') as f:
                content = f.read()
                if termination in content:
                    logger.log(f"Calculation has converged. log file: {log_file_name}")
                    vibrations = re.findall(r"Frequencies --\s+(-?\d+\.\d+)", content)
                    [logger.log(vib) for vib in vibrations if float(vib) < 0]
                    if job_program.lower() == "crest":
                        xyz_coordinates = log2xyz(log_file_path.replace(".log", ".xyz"), job_program)
                    else:
                        xyz_coordinates = log2xyz(log_file_path, job_program)
                    if xyz_coordinates: 
                        logger.log(f"Extracting XYZ coordinates from {log_file_name}")
                        next_step = determine_next_step(log_file_path)

                        if next_step == 'transition_state_optimization':
                            new_input_file = log_file_name.replace(".log", "_TS")
                            shared_data[log_file_name]['result'] = (new_input_file, xyz_coordinates, next_step, True)

                        elif next_step == 'crest_sampling':
                            new_input_file = log_file_name.replace("_TS.log", "_CREST")
                            shared_data[log_file_name]['result'] = (new_input_file, xyz_coordinates, next_step, True)

                        elif next_step == 'ts_optimization_for_conformers':
                            conformers = []
                            new_input_file = log_file_name.replace("_CREST.log", "")
                            for conf in xyz_coordinates:
                                conformers.append(conf)
                            shared_data[log_file_name]['result'] = (new_input_file, conformers, next_step, True)
                        elif next_step == 'DLPNO_SP_for_conformers':
                            new_input_file = log_file_name.replace(".log", "_DLPNO")
                            shared_data[log_file_name]['result'] = (new_input_file, xyz_coordinates, next_step, True)
                        else:
                            logger.log("")
                    else:
                        logger.log(f"Normal termination of {job_program}. However, no XYZ coordinates found. Check log file")
                        shared_data[log_file_name]['result'] = (None, None, None, False)
                    return True
                elif "Error termination" in content:
                    logger.log(f"Error termination in {log_file_name}. Gathering last XYZ coordinates")
                    xyz_coordinates = log2xyz(log_file_path, job_program)
                    if xyz_coordinates:
                        logger.log(f"XYZ coordinates found in failed log file {log_file_name}. Trying to resubmit job")
                        new_input_file = log_file_name.replace(".log", "")
                        shared_data[log_file_name]['result'] = (new_input_file, xyz_coordinates, job_type, False) # None could be exchanged for job type being resubmitted
                    else:
                        logger.log(f"No XYZ coordinates found in {log_file_name}. Check log file for type of error.")
                    return False

                else:
                    attempts += 1
                    logger.log(f"No termination yet in {log_file_name}. Waiting for next check. Attempt: {attempts}/{max_attempts}")
                    time.sleep(interval)
        except FileNotFoundError:
            attempts += 1
            time.sleep(interval)
            logger.log(f"Log file {log_file_path} not found. Waiting for the next check. Attempt: {attempts}/{max_attempts}")
        
    logger.log("Max attempts reached. Calculation may be stock. Check for convergence.")

    shared_data[log_file_name]['result'] = (None, None, False)
    return False



def log2xyz(log_file_path, job_program, match_string='Standard orientation'):
    if job_program.lower() == "g16":
        atomic_number_to_symbol = {
            1: 'H',
            6: 'C',
            7: 'N',
            8: 'O',
            16: 'S'
        }

        with open(log_file_path, 'r') as file:
            coordinates = []
            start_reading = False
            for line in file:
                if start_reading:
                    parts = line.split()
                    if len(parts) >= 6 and parts[1].isdigit() and all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in parts[-3:]):
                        element_symbol = atomic_number_to_symbol.get(int(parts[1]), 'Unknown')
                        coords = [float(parts[3]), float(parts[4]), float(parts[5])]
                        coordinates.append([element_symbol] + coords)
                if match_string in line:
                    start_reading = True
                    coordinates = []
                if "Rotational" in line:
                    start_reading = False
            return coordinates

    elif job_program.lower() == "orca":
        pass # TODO

    elif job_program.lower() == "crest": 
        conformers_list = []
        with open(log_file_path, 'r') as file:
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
def read_xyz_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()[2:]
        coords = [[line.split()[0], float(line.split()[1]), float(line.split()[2]), float(line.split()[3])] for line in lines]
    return coords


def write_xyz_file(destination_path, output_file_name, updated_coords):
    file_path = os.path.join(destination_path, output_file_name)
    with open(file_path, 'w') as f:
        f.write(str(len(updated_coords)) + '\n' +'\n')
        for atom in updated_coords:
            f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')


def mkdir(file, index: list):
    file_name = file.split(".")[0]
    cwd = os.getcwd() # Current working directory
    dir_name = os.path.splitext(file)[0]
    new_dir = os.path.join(cwd, dir_name)
    if os.path.exists(new_dir):
        crest_constrain(new_dir, *index)
            
    else:
        os.mkdir(new_dir)
        crest_constrain(new_dir, *index)

    if args.NEB:
        manager_NEB(file, new_dir)
    else:
        manager_TS(file, new_dir)

    with open(new_dir + "/.constrain", "w") as c:
        c.write(f"{index[0]}, {index[1]}, {index[2]}") # C, H, O


def manager_NEB(file, dir):
    with open(dir + "/commands.txt", "w") as f:
        file_name = file.rsplit('.', 1)[0]
        f.write(f'sbatch -J NEB_TS_{file_name} -p {args.par} --mem={args.mem} -n {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_{program}  {dir}/NEB_TS.inp"\n')


def manager_TS(file, dir):
    with open(dir + "/commands.txt", "w") as f:
        file_name = file.rsplit('.', 1)[0]
        f.write(f'sbatch -J {file_name} -p {args.par} --mem={args.mem} -n {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_{program}  {dir}/{file_name}{dot_inputtype}"\n')
        f.write(f'sh /home/danayo/speciale/scripts/check_convergence.sh {file_name}{dot_outputtype}\n')
        f.write(f'if [ -e ".converged1" ]; then JKTS {file_name}.xyz -{program} -method {args.method} -basis "{args.basis}" -par {args.par} --no-xyz; else JKTS {file_name}.xyz -{program} -method {low_method} -basis "{low_basis}" -par {args.par} --no-xyz --no-TS -constrain; fi\n')
        f.write(f'sbatch -J {file_name} -p {args.par} --mem={args.mem} -n {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_{program}  {dir}/{file_name}{dot_inputtype}"\n')
        f.write(f'sh /home/danayo/speciale/scripts/check_convergence.sh {file_name}{dot_outputtype}\n')
        f.write(f'if [ -e ".TS_converged" ]; then sbatch -J {file_name}_CREST -p {args.par} --mem={args.mem} -n {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_CREST {dir}/{file_name}.xyz -gfn2 -ewin 2 -noreftopo -cinp {dir}/constrain.inp -uhf 1"; else JKTS {file_name}.xyz -{program} -method {args.method} -basis "{args.method}" -par {args.par} --no-xyz; fi; sbatch -J {file_name} -p {args.par} --mem={args.mem} -n {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_{program} {dir}/{file_name}{dot_inputtype}"\n')
        f.write(f'')
        f.write(f"rm {file_name}{dot_inputtype}\n")
        f.write(f'if [ -e ".TS_converged" ]; then JKTS collection{file_name}.pkl -{program} -method {args.method} -basis "{args.basis}" --no-xyz; else sbatch -J {file_name} -p {args.par} --mem={args.mem} -n {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_{program}  {dir}/{file_name}{dot_inputtype}"; fi\n')
        f.write(f'ls "$(pwd)"/*{dot_inputtype} > array.txt\n')
        f.write(f'JKCS3_run -p {program} -cpu {args.cpu} -mem {args.mem} -par {args.par} -mult 2 -maxtasks 100 -rf array.txt -nf LOWDFT_TS\n')
        f.write(f'rm *{dot_inputtype}\n')
        f.write(f'cp LOWDFT_TS/calc-LM/*.xyz .\n')
        f.write(f'JKTS *.xyz -{program} -method {args.method} -basis "{args.basis}" -par {args.par} --no-xyz\n')
        f.write(f'rm *.xyz\n')
        f.write(f'ls "$(pwd)"/*{dot_inputtype} > array.txt\n')
        f.write(f'JKCS3_run -p {program} -cpu {args.cpu} -mem {args.mem} -par {args.par} -mult 2 -maxtasks 100 -rf array.txt -nf DFT_TS\n')
        f.write(f'rm *{dot_inputtype}\n')
        f.write(f'rm *.xyz\n')
        f.write(f'cp DFT_TS/calc-LM/*.xyz .\n')
        f.write(f'sh /home/danayo/scripts/convert/xyz2orca_radical.sh\n')
        f.write(f'ls "$(pwd)"/{file_name}*.inp > array.txt\n')
        f.write(f'JKCS3_run -p ORCA -cpu {args.cpu} -mem {args.mem} -par {args.par} -mult 2 -maxtasks 100 -rf array.txt -nf DLPNO\n')
        f.write("\n")


def manager_reactant(file):
    cwd = os.getcwd()
    reactant_dir = os.path.join(cwd, "reactants")
    if os.path.exists(reactant_dir):
        shutil.copy2(cwd + "/" + file, reactant_dir + "/" + file)
    else:
        os.mkdir(reactant_dir)
        shutil.copy2(file, reactant_dir)
    with open(reactant_dir + "/commands.txt", "w") as f:
        f.write(f'sbatch -J {file} -p {args.par} --mem={args.mem} -c {args.cpu} JKsend "source ~/.JKCSusersetup.txt; program_CREST {reactant_dir}/{file} -gfn2 -ewin 2 -noreftopo"\n')
        f.write(f'rm {file}\n')
        f.write(f'JKTS collection{file.split(".")[0]}.pkl -{program} -method {args.method} -basis "{args.basis}" --no-xyz\n')
        f.write(f'ls "$(pwd)"/*.com > array.txt\n')
        f.write(f'JKCS3_run -p {program} -cpu {args.cpu} -mem {args.mem} -par {args.par} -maxtasks {maxtasks}  -rf array.txt -nf DFT_TS\n')

        
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
def crest_constrain(file_path, C_index, H_index, O_index, force_constant=0.95):
    '''Force constant tells how tight to constrain the atoms i.e. the magnitude of the oscillation between the constrained atoms'''
    with open (file_path + "/constrain.inp","w") as f:
        f.write("$constrain\n")
        f.write(f"  force constant={force_constant}\n") 
        f.write(f"  distance: {C_index}, {H_index}, auto\n")
        f.write(f"  distance: {H_index}, {O_index}, auto\n")
        f.write(f"  angle: {C_index}, {H_index}, {O_index}, auto\n")
        f.write("$end\n")


def QC_input(file_name, destination, coords,constrain, program, method, basis_set, TS, C_index=None, H_index=None, O_index=None):
    if constrain and C_index == None and H_index == None and O_index == None:
        with open(os.path.join(os.getcwd(), destination, "/.constrain"), "r") as f:
            content = f.read()
            C_index, H_index, O_index = [int(num) for num in content.split(",")]

    if program == "ORCA":
        file_name = file_name + ".inp"
        file_path = os.path.join(destination, file_name)
        with open(file_path, "w") as f:
            if method == "DLPNO" or basis_set == "DLPNO":
                f.write(f"! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF RI-JK aug-cc-pVTZ/JK\n")
                args.cpu = 1 # TEMP solution
            elif TS:
                f.write(f"! {method} {basis_set} OptTS freq\n")
            else:
                f.write(f"! {method} {basis_set} Opt\n")
            f.write(f"%pal nprocs {args.cpu} end\n")
            f.write(f"%maxcore {args.mem}\n")
            if constrain:
                f.write("%geom\n")
                f.write("Constraints\n")
                # ORCA start atom indexing from 0
                f.write(f"{{B {C_index-1} {H_index-1} C}}\n") # Bond being broken
                f.write(f"{{B {H_index-1} {O_index-1} C}}\n") # Bond being formed
                f.write("end\nend\n")
            f.write("\n")
            f.write("* xyz 0 2\n") # charge and spin multiplicity usually 0 and 2 for atmospheric radical reaction
            for atom in coords:
                f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')
            f.write("*")
    elif program == "G16":
        file_name = file_name + ".com"
        file_path = os.path.join(destination, file_name)
        with open(file_path, "w") as f:
            f.write(f"%nprocshared={args.cpu}\n")
            f.write(f"%mem={args.mem}mb\n")
            if TS and constrain:
                f.write(f"# {method} {basis_set} opt=(calcfc,ts,noeigen,modredundant) freq\n\n") # freq may be redundant
            elif TS and constrain is False:
                f.write(f"# {method} {basis_set} opt=(calcfc,ts,noeigen) freq\n\n")
            elif TS is False and constrain:
                f.write(f"# {method} {basis_set} opt=modredundant\n\n")
            else:
                f.write(f"# {method} {basis_set} opt\n\n")
            f.write("Title\n\n")
            f.write("0 2\n")
            for atom in coords:
                f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')
            f.write("\n")
            if constrain:
                f.write(f"B {C_index} {H_index} F\n")
                f.write(f"B {H_index} {O_index} F\n")
                f.write("\n")
    else:
        print("QC_input was called but no program was specified")


def NEP_input(file_path, file_name):
    if args.NEB:
        with open(file_path + "/NEB_TS.inp", "w") as f:
            f.write(f'! B3LYP 6-31+g(d,p)  NEB-TS FREQ\n')
            f.write(f'%NEB PREOPT_ENDS TRUE NEB_END_XYZFILE "{file_path + "/" + file_name}_product.xyz" END\n')
            f.write(f'* XYZfile 0 2 {file_path + "/" + file_name}_reactant.xyz\n')


##############################ADDITIONAL FUNCTIONS##################################
def get_terminal_O(coords, distance=1.5):
    O_index = []
    for i,atom in enumerate(coords):
        if atom[0] == "O":
            neighbors = sum(1 for j, other_atom in enumerate(coords) if i != j and atom_distance(atom[1:], other_atom[1:]) < distance)
            if neighbors == 1:
                O_index.append(i)
    return O_index[0]


#####################################MAIN FUNCTIONS####################################
def H_abstraction(file, distance=1.35, dist_OH=0.97, NEB=False):
    coords = read_xyz_file(file)
    num_atoms = len(coords)
    modified_coords = []
    sampled_hydrogens = set()

    for i in range(num_atoms):
        if coords[i][0] == "H":  # Atom label for hydrogen
            for j in range(num_atoms):
                if coords[j][0] == "C":  # Atom label for carbon
                    vector_CH = np.array(coords[i][1:]) - np.array(coords[j][1:])
                    dist_CH = vector_length(vector_CH)

                    if dist_CH < distance:
                        norm_vector_CH = normalize_vector(vector_CH)


                        H_perturb_axis = np.cross(norm_vector_CH, [0, 1, 0])
                        H_perturb_axis = normalize_vector(H_perturb_axis)
                        
                        # Update the initial hydrogen coordinates
                        if NEB:
                            new_coords = read_xyz_file(file) # Coordinates for TS guess
                            reactant_coords = read_xyz_file(file)
                            product_coords = read_xyz_file(file)
                            distance_reactant = 2
                            distance_product = 1.8
                        else:
                            new_coords = read_xyz_file(file)
                        
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + norm_vector_CH * (distance - dist_CH)

                        # Update oxygen in OH coordiantes
                        oxygen_coords = np.array(new_coords[i][1:]) + norm_vector_CH * distance

                        # Coordinates for hydrogen in OH
                        norm_vector_OH = normalize_vector(oxygen_coords - np.array(coords[i][1:]))
                        rotation_axis = np.cross(norm_vector_OH, H_perturb_axis)
                        rotation_axis = normalize_vector(rotation_axis)

                        rotation_angle_H = np.radians(75.5) # 180 - 75.5 = 104.5

                        rotated_vector_H = rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H)
                        hydrogen_coords = oxygen_coords + rotated_vector_H * dist_OH
                        
                        # Perturb inital hydrogen to avoid linear C--H--O
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + H_perturb_axis * 0.3

                        # Update XYZ file
                        new_coords.append(['O', *oxygen_coords])
                        new_coords.append(['H', *hydrogen_coords])

                        if NEB:
                            oxygen_reactant = np.array(reactant_coords[i][1:]) + norm_vector_CH *distance_reactant
                            hydrogen_reactant = oxygen_reactant + rotated_vector_H * dist_OH
                            reactant_coords.append(['O', *oxygen_reactant])
                            reactant_coords.append(['H', *hydrogen_reactant])

                            product_coords[i][1:] = np.array(product_coords[i][1:]) + norm_vector_CH * (distance_product - dist_CH)
                            oxygen_product = np.array(product_coords[i][1:]) + norm_vector_CH * dist_OH
                            hydrogen_product = oxygen_product + rotated_vector_H * dist_OH
                            product_coords.append(['O', *oxygen_product])
                            product_coords.append(['H', *hydrogen_product])

                        C_index = j+1
                        H_index = i+1
                        O_index = len(new_coords)-1
                        # OH_index = len(new_coords)
                        modified_coords.append((new_coords, (C_index, H_index, O_index)))

    return modified_coords

                        
def find_equivalent_hydrogens(coords, hydrogen_index, carbon_index, threshold=0.1):
    equivalent_hydrogens = []
    for index, atom in enumerate(coords):
        if atom[0] == "H" and index != hydrogen_index:
            vector_HC = np.array(coords[index][1:]) - np.array(coords[carbon_index][1:])
            dist_HC = vector_length(vector_HC)
            vector_HC_target = np.array(coords[hydrogen_index][1:]) - np.array(coords[carbon_index][1:])
            dist_HC_target = vector_length(vector_HC_target)

            if abs(dist_HC - dist_HC_target) < threshold:
                equivalent_hydrogens.append(index)

    return equivalent_hydrogens


def OH_addition(file, distance=1.45, double_bond_distance=1.36, dist_oh=0.97):
    coords = read_xyz_file(file)
    num_atoms = len(coords)
    count = 1
    modified_coords = []

    for i in range(num_atoms):
        if coords[i][0] == "c":
            for j in range(num_atoms):
                if coords[j][0] == "c" and i != j:
                    vector_cc = np.array(coords[i][1:]) - np.array(coords[j][1:])
                    dist_cc = vector_length(vector_cc)
                    if dist_cc <= double_bond_distance:
                        norm_vector_cc = normalize_vector(vector_cc)
                        new_coords = read_xyz_file(file)
                        # note: update the way perpendicular axis is calculated. copy from addition()
                        perpendicular_axis = np.cross(norm_vector_cc, [0,1,0])
                        perpendicular_axis = normalize_vector(perpendicular_axis)
                        # shift carbon
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + norm_vector_cc * 0.1
                        # update oxygen in oh coordiantes
                        oxygen_coords = np.array(new_coords[j][1:]) + perpendicular_axis * distance

                        rotation_axis = norm_vector_cc
                        rotation_angle_h = np.radians(45) 

                        rotated_vector_h = rotate_vector(perpendicular_axis, rotation_axis, rotation_angle_h)

                        hydrogen_coords = oxygen_coords + rotated_vector_h * dist_oh
                        
                        new_coords.append(['o', *oxygen_coords])
                        new_coords.append(['h', *hydrogen_coords])

                        modified_coords.append(new_coords)

                        base_file_name = os.path.splitext(os.path.basename(file))[0]
                        write_xyz_file(f"{base_file_name}_{count}.xyz", new_coords)

                        count += 1


def addition(file1, file2, method, basis_set, TS, no_xyz, program=None, constrain=False, distance = 1.55, double_bond_distance=1.36, elongation_factor=0.1):
    coords1 = read_xyz_file(file1)
    coords2 = read_xyz_file(file2)
    num_atoms1 = len(coords1)
    num_atoms2 = len(coords2)
    count = 1
    carbon_index = []
    vector_lst = []
    
    # gather carbon index
    for i in range(num_atoms1):
        if coords1[i][0] == "c":
            for j in range(num_atoms1):
                if atom_distance(coords1[i][1:], coords1[j][1:]) < 1.52:
                    perpendicular_vector = (coords1[i][1:], coords1[j][1:])
                    vector_lst.append(perpendicular_vector)
                if coords1[j][0] == "c" and i != j:
                    bond_cc = atom_distance(coords1[i][1:], coords1[j][1:])
                    if bond_cc <= double_bond_distance:
                        carbon_index.append(i)
                        carbon_index.append(j)

    for index in carbon_index[::2]:
        new_coords = read_xyz_file(file1)
        direction = calculate_vector(new_coords[index][1:], new_coords[index+1][1:]) # vector between c=c
        # perpendicular_vector = [calculate_vector(coords1[index][1:], coords1[k][1:]) for k in range(num_atoms1) if atom_distance(coords1[index][1:], coords1[k][1:]) < 1.52]

        for k in range(num_atoms1):
            if atom_distance(coords1[index][1:], coords1[k][1:]) < 1.52:
                perpendicular_vector = calculate_vector(coords1[index][1:], coords1[k][1:])
                break
            else: perpendicular_vector = np.array([0,0,1])

        new_coords[index][1:] -= direction * elongation_factor/2 # standard is 10% increase in bond length
        new_coords[index+1][1:] += direction * elongation_factor/2

        perpendicular_axis = np.cross(direction, perpendicular_vector)
        perpendicular_axis = normalize_vector(perpendicular_axis)

        o_index = get_terminal_O(coords2)

        first_pos = new_coords[index][1:] + perpendicular_axis * distance

        rotation_axis = normalize_vector(np.cross(direction, perpendicular_axis))

        # note: fix the rotation axis
        for i in range(num_atoms2):
            relative_vector = calculate_vector(coords2[o_index][1:], coords2[i][1:])
            dot = np.dot(perpendicular_axis, normalize_vector(relative_vector))
            if dot > 0:
                 rotation_angle = np.radians(90) # radians value should be function of dot product
                 rotated_vector = rotate_vector(relative_vector, rotation_axis, rotation_angle)
            else:
                rotation_angle = np.radians(-270)
                rotated_vector = rotate_vector(relative_vector, rotation_axis, rotation_angle)
            new_pos = first_pos + relative_vector + (rotated_vector * 0.1)
            new_coords.append([coords2[i][0], *new_pos])
        

        base_file_name1 = os.path.splitext(os.path.basename(file1))[0]
        base_file_name2 = os.path.splitext(os.path.basename(file2))[0]      
        write_xyz_file(f"{base_file_name1}_{base_file_name2}_{count}.xyz", new_coords)
        count += 1

                        

##################################thermochemistry##############################

def partition_function(vibrations: list, rot_constants, symmetry_num, mol_mass, multiplicity, T):
    # vibrational partition function
    qvib = 1 
    for vib in vibrations:
        vib = float(vib)
        if vib > 0:
            qvib *= 1 / (1 - np.exp(-(h * c * vib) / (k_b * T))) 

    # rotational partition function
    ghztowave = 29.978869
    rot_constants = [float(i) * ghztowave for i in rot_constants]
    if 0.0 in rot_constants:
        rot_constant = [e for e in rot_constants if e != 0.0]  
        qrot = ((k_b*T)/(int(symmetry_num)*h*c*rot_constant[0]))
    else:
        qrot = (1/int(symmetry_num)) * ((k_b*T)/(h*c))**1.5 * (np.pi / (rot_constants[0] * rot_constants[1] * rot_constants[2]))**0.5


    # translational partition function
    mol_mass = float(mol_mass) * 1.6605*(10**(-27))
    lmbda = h/(np.sqrt(2*np.pi*mol_mass*k_b*T))
    qtrans = 0.02479/(lmbda**3)

    qelec = float(multiplicity)

    q = qvib*qrot*qtrans*qelec

    return q


def rate_constant(files, T=293.15, program=None):
    reactants = []
    transition_states = []
    
    if program == "G16":
        for file in files:
            dic = {}
            with open (file, "r") as f:
                content = f.read()
                job_info = re.search(r"#(.*)", content)
                job_info = re.split('[, ()]', job_info.group(0)) 
                ee_zpe = re.search(r"Sum of electronic and zero-point Energies=\s+([-+]?\d*\.\d+|\d+)", content)
                if ee_zpe:
                    vibration = re.findall(r"Frequencies --\s+(-?\d+\.\d+)", content)
                    rot_constant = re.findall(r"Rotational constants \(GHZ\):\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", content)
                    rot_constant = rot_constant[-1]
                    symmetry_num = re.search(r"Rotational symmetry number\s*(\d+)", content)
                    if symmetry_num:
                        symmetry_num = symmetry_num.group(1)
                    else: 
                        print("No symmetry number found. assuming 1")
                        symmetry_num = 1
                    mol_mass = re.search(r"Molecular mass:\s+(-?\d+\.\d+)", content).group(1)
                    multiplicity = re.search(r"Multiplicity =\s*(\d+)", content).group(1)


                    Q = partition_function(vibration, rot_constant, symmetry_num, mol_mass, multiplicity, T)

                    dic[file] = ""
                    dic["ee_zpe"] = float(ee_zpe.group(1)) 
                    dic["Q"] = Q 

                    if 'ts' in job_info:
                        transition_states.append(dic)
                    else:
                        reactants.append(dic)

                else:
                    print(f"No energies in {file}. check if calculation has converged")


        lowest_zpec_ts = min(transition_states, key=lambda x: x['ee_zpe'])['ee_zpe']
        # lowest_zpec_r = min(reactants, key=lambda x: x['ee_zpe'])['ee_zpe']
        
        HtoJ = 4.3597447222*(10**(-18))
        sum_ts = sum([np.exp((lowest_zpec_ts*HtoJ - transition_states[i]['ee_zpe']*HtoJ) / (k_b*T)) *transition_states[i]['Q'] for i in range(len(transition_states))])
        # sum_r = [np.exp((lowest_zpec_r*htoj - reactants[i]['ee_zpe']*htoj) / (k_b*t)) *reactants[i]['q'] for i in range(len(transition_states))]
        print(sum_ts)

    elif program == "ORCA":
        for file in files:
            vibrations = []
            dic = {}
            with open(file, "r") as f:
                content = f.read()
                vibrations = []
                EE = re.search(r"Electronic energy\s*...\s*[-+]?\d*\.\d+", content)
                if EE:
                    EE = float(EE.group().split()[-1])
                
                    vib = re.search(r'[-+]?\d*\.\d+\s*cm\*\*-1', content)
                    if vib:
                        vibration = float(vib.group().split()[0])
                        vibrations.append(vibration)
                    else: print(f"No vibrations found in {file}")
        
                    rot = re.findall(r"Rotational constants in cm-1: \s*[-+]?\d*\.\d*  \s*[-+]?\d*\.\d* \s*[-+]?\d*\.\d*", content)
                    if rot:
                        rot_constants = rot[0].split()[-3:]
                    else: print(f"No rotational constants found in {file}")

                    symmetry_num = re.search(r'Symmetry Number:\s*(\d*)', content)
                    if symmetry_num:
                        symmetry_num = int(symmetry_num.group(1))
                    else: 
                        print(f"No symmetry number found in {file}. Assuming 1")
                        symmetry_num = 1

                    mol_mass = re.search(r'Total Mass\s*...\s*\d*\.\d+', content)
                    if mol_mass:
                        mol_mass = float(mol_mass.group().split()[-1])
                    else: print(f"No molecular mass found in {file}")

                    multiplicity = re.search(r'Mult\s* ....\s*(\d*)', content)
                    if multiplicity:
                        multiplicity = int(multiplicity.group().split()[-1])
                    else: print(f"No multiplicity found in {file}")

                else:
                    print(f"No energies in {file}. check if calculation has converged")

                Q = partition_function(vibrations, rot_constants, symmetry_num, mol_mass, multiplicity, T)
                print(Q)



def submit_and_monitor(dir, input_file_path, logger, convergence_info, threads, job_type, job_program):
    log_file_name = input_file_path.split("/")[-1].replace(f"{dot_inputtype}", ".log")
    logger.log(f"Submitting file {input_file_path.split('/')[-1]} for calculation in path {dir}")
    submit_job(dir, input_file_path, program, args.cpu, args.mem, args.par, args.time)
    convergence_info[log_file_name] = {'dir': dir, 'logger': logger, 'result': None, 'job_type': job_type, 'job_program': job_program}
    thread = threading.Thread(target=check_convergence, args=(log_file_name, dir, logger, convergence_info, job_type, job_program))
    threads.append(thread)
    thread.start()



def handle_convergence_result(convergence_info, threads):
    for log_file, job_info in list(convergence_info.items()):
        if job_info['result'] is not None:
            result = job_info['result']
            current_dir = job_info['dir']
            last_job = job_info['job_type']
            current_logger = job_info['logger']
            new_input_file, xyz_coordinates, job_type, converged = result
            log_file_path = os.path.join(current_dir, log_file)
            if xyz_coordinates:
                next_step = determine_next_step(log_file_path)
                print(next_step)

                if converged:
                    current_logger.log(f"Yay converged. Next step is: {job_type}")

                    if next_step == 'transition_state_optimization':
                        QC_input(file_name=new_input_file, destination=current_dir,coords=xyz_coordinates, constrain=False, method=high_method, basis_set=high_basis, program=program, TS=True)
                        submit_job(current_dir, new_input_file + dot_inputtype, program, args.cpu, args.mem, args.par, args.time)
                        current_logger.log(f'Submitted new {program} job with input file {new_input_file}{dot_inputtype}')
                        # Start monitoring the new job
                        new_log_file = new_input_file + ".log"
                        convergence_info[new_log_file] = {'dir': current_dir, 'logger': current_logger, 'job_type': next_step, 'result': None}
                        new_thread = threading.Thread(target=check_convergence, args=(new_log_file, current_dir, current_logger, convergence_info, 'transition_state_optimization', program))
                        threads.append(new_thread)
                        new_thread.start()
                        job_info['result'] = None


                    elif next_step == 'crest_sampling':
                        print("Doing CREST sampling")
                        write_xyz_file(current_dir, new_input_file + ".xyz", xyz_coordinates)
                        submit_job(current_dir, new_input_file + ".xyz", "CREST", args.cpu, args.mem, args.par, args.time)
                        current_logger.log(f'Submitted CREST job with input file {new_input_file}.xyz')
                        new_log_file = new_input_file + ".log"
                        convergence_info[new_log_file] = {'dir': current_dir, 'logger': current_logger, 'job_type': next_step, 'result': None}
                        new_thread = threading.Thread(target=check_convergence, args=(new_log_file, current_dir, current_logger, convergence_info, 'crest_sampling', "CREST"))
                        threads.append(new_thread)
                        new_thread.start()
                        job_info['result'] = None

                    elif next_step == 'ts_optimization_for_conformers':
                        print("Doing TS optimization on conformers")
                        conformers_dir = os.path.join(current_dir, "conformers_TS")
                        os.makedirs(conformers_dir, exist_ok=True)
                        job_files = []
                        for n, conf in enumerate(xyz_coordinates, start=1):
                            conf_file_name = f"{new_input_file}_conf{n}"
                            conf_file_path = os.path.join(conformers_dir, conf_file_name)
                            QC_input(file_name=conf_file_name, destination=conformers_dir, coords=conf, constrain=False, method=high_method, basis_set=high_basis, program=program, TS=True)
                            job_files.append(conf_file_name + dot_inputtype)
                            new_log_file = conf_file_name + ".log"
                            convergence_info[new_log_file] = {'dir': conformers_dir, 'logger': current_logger,'job_type': next_step, 'result': None}
                            new_thread = threading.Thread(target=check_convergence, args=(new_log_file, conformers_dir, current_logger, convergence_info, 'transition_state_optimization_for_conformers', program))
                            threads.append(new_thread)
                            new_thread.start()
                                                    
                        input_array_list_name = "array.txt"
                        submit_array_job(conformers_dir, job_files, input_array_list_name, f"{new_input_file}_array", program, args.par, args.time, args.cpu, args.mem)
                        current_logger.log(f'Submitting TS optimization on CREST conformers with array file: {input_array_list_name}')


                    elif next_step == 'DLPNO_SP_for_conformers':
                        print("Doing Couple Cluster")
                        DLPNO_dir = os.path.join(current_dir, "DLPNO")
                        os.makedirs(DLPNO_dir, exist_ok=True)
                        QC_input(file_name=new_input_file, destination=DLPNO_dir, coords=xyz_coordinates, method="DLPNO", basis_set="DLPNO", program="ORCA", TS=False, constrain=False)
                        submit_job(DLPNO_dir, new_input_file+".inp", "ORCA", 1, args.mem, args.par, args.time)
                        current_logger.log(f"DLPNO single point calculation with input file: {new_input_file}.inp in folder confomers_TS/DLPNO")
                        new_log_file = new_input_file + ".log"
                        convergence_info[new_log_file] = {'dir': DLPNO_dir, 'logger': current_logger, 'job_type': next_step, 'result': None}
                        new_thread = threading.Thread(target=check_convergence, args=(new_log_file, DLPNO_dir, current_logger, convergence_info, 'DLPNO_SP_for_conformers', "ORCA"))
                        threads.append(new_thread)
                        new_thread.start()
                        job_info['result'] = None

                    else:
                        current_logger.log(f"Error: next step could not be determined for log file: {log_file_path}")


                else:
                    print(f"failed converged log file: {log_file_path}. Resubmiting calculation {job_type} {next_step}")
                    if next_step == "transition_state_optimization":
                        current_logger.log(f"Failed preoptimization for log file: {log_file_path}. Redoing calculation {job_type} {next_step}")
                        QC_input(file_name=new_input_file, destination=current_dir, coords=xyz_coordinates, constrain=True, method=low_method, basis_set=low_basis, program=program, TS=False)
                        submit_job(current_dir, new_input_file + dot_inputtype, program, args.cpu, args.mem, args.par, args.time)
                        current_logger.log(f"Resubmitted {program} job with input file {new_input_file}{dot_inputtype}")
                        new_log_file = new_input_file + ".log"
                        convergence_info[new_log_file] = {'dir': current_dir, 'logger': current_logger, 'job_type': next_step, 'result': None}
                        new_thread = threading.Thread(target=check_convergence, args=(new_log_file, current_dir, current_logger, convergence_info, 'transition_state_optimization', program))
                        threads.append(new_thread)
                        new_thread.start()
                        job_info['result'] = None

                    elif next_step == "crest_sampling":
                        print("Redoing TS optimization")
                        current_logger.log(f"Failed optimization towards first order saddle point. Log file: {log_file_path}. Redoing calculation")
                        QC_input(file_name=new_input_file, destination=current_dir, coords=xyz_coordinates, constrain=False, method=high_method, basis_set=high_basis, program=program, TS=True)
                        submit_job(current_dir, new_input_file + dot_inputtype, program, args.cpu, args.mem, args.par, args.time)
                        current_logger.log(f"Resubmitted {program} job with input file {new_input_file}{dot_inputtype}")
                        new_log_file = new_input_file + ".log"
                        convergence_info[new_log_file] = {'dir': current_dir, 'logger': current_logger, 'job_type': next_step, 'result': None}
                        new_thread = threading.Thread(target=check_convergence, args=(new_log_file, current_dir, current_logger, convergence_info, 'transition_state_optimization', program))
                        threads.append(new_thread)
                        new_thread.start()
                        job_info['result'] = None

                    elif next_step == "ts_optimization_for_conformers":
                        print("Redoing CREST sampling")
                        write_xyz_file(current_dir, new_input_file + ".xyz", xyz_coordinates)
                        submit_job(current_dir, new_input_file + ".xyz", "CREST", args.cpu, args.mem, args.par, args.time)
                        current_logger.log(f'Submitted CREST job with input file {new_input_file}.xyz')
                        new_log_file = new_input_file + ".log"
                        convergence_info[new_log_file] = {'dir': current_dir, 'logger': current_logger, 'job_type': next_step, 'result': None}
                        new_thread = threading.Thread(target=check_convergence, args=(new_log_file, current_dir, current_logger, convergence_info, 'crest_sampling', "CREST"))
                        threads.append(new_thread)
                        new_thread.start()
                        job_info['result'] = None


                job_info['result'] = None

def determine_next_step(log_file_path):
    log_file_name = log_file_path.split("/")[-1]
    try:
        with open(log_file_path, 'r') as file:
            for i in range(30, 200): # Check for method section between line 50 and 200. Adjust accordingly
                try:
                    line = next(file).strip()
                    if line.startswith('#'):
                        components = line.split(",")
                        if 'ts' not in components:
                            return 'transition_state_optimization'
                        elif 'ts' in components: 
                            if 'conf' in log_file_name:
                                vibrations = log2vib(log_file_path, program)
                                print(vibrations)
                                if vibrations:
                                    return 'DLPNO_SP_for_conformers'
                            else: return 'crest_sampling' 
                    elif line.startswith('xTB'):
                        return 'ts_optimization_for_conformers'
                except StopIteration:
                    break  
        return 'Method section not found'
    except FileNotFoundError:
        print(f"Error: Log file {log_file_name} not found.")
        return 'error'


def log2vib(log_file_path, program):
    with open(log_file_path, 'r') as file:
        content = file.read()
        if program.lower == "g16":
            vibrations = re.findall(r"Frequencies --\s+(-?\d+\.\d+)", content)
        elif program.lower == "orca":
            vibrations = []
            vib = re.search(r'[-+]?\d*\.\d+\s*cm\*\*-1', content)
            if vib:
                vibration = float(vib.group().split()[0])
                vibrations.append(vibration)
        else:
            return 'No vibrations found'
    return vibrations
        
    

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
    parser.add_argument('-auto', action='store_true', help='Automated process with the workflow: \n-> Preoptimization as low level of theory \n-> TS optimization as high level of theory \n-> CREST TS conformer sampling \n-> DLPNO-CCSD(T) SP energy calculations on top of TS conformers \n-> calculate rate constants and branching rations for reaction type \n- Resubmission of failed calculations is automatically done until convergence')

    additional_options = parser.add_argument_group("Additional arguments")

    additional_options.add_argument('--no-TS', action='store_true', help='Input files for geometry relaxation are generated')
    additional_options.add_argument('--no-xyz', action='store_true', help='No XYZ files generated')
    additional_options.add_argument('-k', action='store_true', help='Calculate Multiconformer Transition State rate constant')
    additional_options.add_argument('--high_level', nargs=2, metavar='', default=['wb97xd', 'aug-cc-pVTZ'],  help='Specify the high level of theory for QC method TS optimization [def = wB97X-D aug-cc-pVTZ]')
    additional_options.add_argument('--low_level', nargs=2, metavar='', default=['B3LYP', '6-31+G(d,p)'],  help='Specify the low level of theory for preoptimization [def = B3LYP 6-31+G(d,p)]')
    additional_options.add_argument('-method', nargs="?", default='wb97xd',  help='Specify the QC method [def = wB97X-D]') # redundant
    additional_options.add_argument('-basis',  nargs='?', default="6-31+g(d,p)", help='Specify the basis set used [def = 6-31+G(d,p)]') # redundant
    additional_options.add_argument('-cpu', metavar="int", nargs='?', const=1, type=int, default=4, help='Number of CPUs [def = 4]')
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=4000, help='Amount of memory allocated for job [def = 400mb]')
    additional_options.add_argument('-par', metavar="partition", nargs='?', const=1, default="qany", help='Partition to use [def = qany]')
    additional_options.add_argument('-time', metavar="hours:minutes:seconds", nargs='?', const=1, default="72:00:00", help='Specify total time for calculation [def = 72 Hours]')

    global args
    args = parser.parse_args()
    start_dir = os.getcwd()

    ################################ARGUMENT SPECIFICATIONS############################
    if args.auto:
        args.constrain=True; args.no_TS=True; args.G16=True; args.reactants=True; args.no_xyz=True; args.time="18:00:00";

    # Program and method 
    global low_basis, low_method, high_basis, high_method
    low_basis, low_method = args.low_level
    high_basis, high_method = args.high_level

    methods = ["B97-3c", "r2scan-3c", "pm3", "am1", "pm6", "pm7"] 

    global program
    global dot_inputtype
    global dot_outputtype
    if args.G16:
        program = "G16"
        dot_inputtype = ".com"
        dot_outputtype = ".log"

    elif args.ORCA:
        program = "ORCA"
        dot_inputtype = ".inp"
        dot_outputtype = ".out"
        if low_method in methods:
            low_basis = ""
        if high_method in methods:
            high_basis = ""
        if high_method == "wb97xd":
            high_method = "WB97X-D3"
        if low_method == "wb97xd":
            low_method = "WB97X-D3" 
    else: 
        program = None # Just produce XYZ-file
        dot_inputtype = ".xyz"
        dot_outputtype = ".log"

    if args.no_TS:
        TS = False
    else: TS = True

    
    ####################################################################################################
    threads = []
    convergence_info = {}
    global conformer_info
    conformer_info = {}

    for n, input_file in enumerate(args.input_files):
        file_name, file_type = os.path.splitext(input_file)

        if file_type == ".xyz":
            # logger = Logger(os.path.join(start_dir, file_name, "log"))

            if args.H:
                # Hydrogen abstraction process
                modified_coords = H_abstraction(input_file, NEB=args.NEB)
                for count, (coords, (C_index, H_index, O_index)) in enumerate(modified_coords, start=1):
                    input_file_count = f"{file_name}_H{count}"
                    mkdir(input_file_count, [C_index, H_index, O_index])
                    logger = Logger(os.path.join(start_dir, input_file_count, "log")) # Create file for keeping log
                    dir_for_each_H = os.path.join(start_dir, input_file_count)
                    if args.no_xyz is False:
                        write_xyz_file(dir_for_each_H, f"{input_file_count}.xyz", coords) # Modify to move xyz file to new dir

                    if program is not None:
                        QC_input(file_name=input_file_count, destination=dir_for_each_H, coords=coords, C_index=C_index, H_index=H_index, O_index=O_index, constrain=args.constrain, method=low_method, basis_set=low_basis, program=program, TS=False)
                        input_file_path = os.path.join(dir_for_each_H, f"{input_file_count}{dot_inputtype}")
                        submit_and_monitor(dir_for_each_H, input_file_path, logger, convergence_info, threads, job_type='preoptimization', job_program=program)

            elif args.crest:
                executing_path = os.getcwd()
                logger = Logger(os.path.join(executing_path, "log"))
                input_file_path = os.path.join(executing_path, input_file)
                submit_and_monitor(executing_path, input_file_path, logger, convergence_info, threads, job_type='crest', job_program=program)

            elif args.CC and len(args.input_files) == 2:
                # C=C addition
                addition(args.input_files[n], args.input_files[n+1], constrain=args.constrain, program=program, method=low_method, TS=False, no_xyz=args.no_xyz, basis_set=low_basis)

            elif not args.H and not args.CC:
                # Default case for XYZ files without specific reactions
                coords = read_xyz_file(input_file)
                QC_input(file_name=input_file.split(".")[0], destination=os.getcwd(), coords=coords, constrain=args.constrain, program=program, method=args.method, basis_set=args.basis, TS=TS)

            else:
                print("Unsupported reaction type or combination of arguments.")

        elif file_type in [".log", ".out"]:
            executing_path = os.getcwd()
            logger = Logger(os.path.join(executing_path, "log"))
            next_step = determine_next_step(input_file)

            if next_step == 'transition_state_optimization':
                logger.log("Detected preoptimized structure. Submitting calculation for first order saddle point optimization")
                xyz_coordinates = log2xyz(input_file, program)
                if xyz_coordinates:
                    new_input_file = f"{file_name}_TS{dot_inputtype}"
                    input_file_path = os.path.join(executing_path, new_input_file)
                    submit_and_monitor(executing_path, input_file_path, logger, convergence_info, threads, job_type='transition_state_optimization', job_program=program)

            elif next_step == 'crest_sampling':
                logger.log("Detected log file with transition state structure. Submitting for CREST conformer sampling of transition state")
                xyz_coordinates = log2xyz(input_file + "_CREST.xyz", program)
                if xyz_coordinates:
                    new_input_file = f"{file_name}_CREST.xyz" # This is prob wrong. TODO: _CREST.xyz is not generated
                    input_file_path = os.path.join(executing_path, new_input_file)
                    submit_and_monitor(executing_path, input_file_path, logger, convergence_info, threads, job_type='crest_sampling', job_program=program)

            elif next_step == 'ts_optimization_for_conformers':
                logger.log("Detected CREST log file. Assuming conformers are of transition state. Submitting calculation for first order saddle point optimization for all conformers")
                conformers = log2xyz(input_file, program)
                if conformers:
                    for n, conf in enumerate(conformers):
                        QC_input(file_name=f"{file_name}_conf{n}", destination=executing_path, coords=conf, constrain=False, program=program, method=high_method, basis_set=high_basis, TS=True)
                        submit_and_monitor(executing_path, f"{file_name}_conf{n}{dot_inputtype}", logger, convergence_info, threads, job_type='ts_optimization_for_conformers', job_program=program)


            elif next_step == 'cc_single_point_energy':
                logger.log("Detected the need for performing DLPNO-CCSD(T) calculations on top of TS conformers")
                print("DLPNO to be done")

            elif next_step == 'default_action':
                print("Doing nothing")


        else:
            print(f"Unsupported file type for input: {input_file}")

    # Monitor and handle convergence of submitted jobs
    while threads:
        for thread in list(threads):
            thread.join(timeout=0.1)
            if not thread.is_alive():
                threads.remove(thread)
                handle_convergence_result(convergence_info, threads)

    print("You are now finished")


if __name__ == "__main__":
    main()
