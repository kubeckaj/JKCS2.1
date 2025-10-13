import os
import time
import re
import subprocess
from threading import Thread

from .config import config
from .slurm_submit import submit_job, submit_array_job
from .ArbAlign import energy_cutoff, filter_molecules
from .utils import is_aldehyde, read_last_lines
from .QC import QC_input
from .transition_states import check_transition_state, good_active_site
from .input import initiate_conformers
from .molecule import Molecule
from .molecule_manager import global_molecules

termination_strings = {
"g16": ["normal termination"],
"orca": ["****orca terminated normally****", "ORCA TERMINATED NORMALLY", "TOTAL RUN TIME", "total run time"],
"crest": ["crest done", "crest terminated normally."]  # Multiple possible termination messages
}
error_strings = {
"g16": ["error termination", "another termination example"],
"orca": ["aborting the run", "this wavefunction is not fully converged", "not enough memory", "orca finished by error termination", "error", "this wavefunction is not converged"],
"crest": ["some crest error message"]  # Multiple possible error messages
}

def crest_constrain(molecule, args, force_constant=1):
    if molecule.reactant or molecule.product:
        pass
    else:
        if not molecule.constrained_indexes:
            molecule.find_active_site(indexes=args.CHO)
        C_index = molecule.constrained_indexes['C']
        H_index = molecule.constrained_indexes['H']
        O_index = molecule.constrained_indexes['O']

        aldehyde, aldehyde_O = is_aldehyde(molecule, C_index, H_index)
        CHO_angle = molecule.calculate_angle(molecule.coordinates[C_index-1], molecule.coordinates[H_index-1], molecule.coordinates[O_index-1])

        with open (molecule.directory + "/constrain.inp","w") as f:
            f.write("$constrain\n")
            f.write(f"  force constant={force_constant}\n") 
            f.write(f"  distance: {C_index}, {H_index}, auto\n")
            f.write(f"  distance: {H_index}, {O_index}, auto\n")
            f.write(f"  angle: {C_index}, {H_index}, {O_index}, {CHO_angle}\n")
            if args.OH and aldehyde: 
                f.write(f"  dihedral: {aldehyde_O}, {C_index}, {O_index}, {O_index+1}, 0\n")
                f.write(f"  dihedral: {C_index}, {H_index}, {O_index}, {O_index+1}, 0\n")
            f.write("$end\n")


def check_crest(molecules, logger, threads, interval, max_attempts, args):
    initial_delay = int(interval * 2)
    attempts = 0
    sleeping = False
    pending = []
    all_conformers = []
    expected_files = {f"collection{molecule.name}.pkl" for molecule in molecules}
    constrained_indexes = molecules[0].constrained_indexes
    mult = molecules[0].mult
    charge = molecules[0].charge
    current_step = molecules[0].current_step
    dir = molecules[0].directory
    method = molecules[0].method
    
    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)

    while attempts < max_attempts:
        update_molecules_status(molecules)
        for molecule in molecules:
            if molecule.status == 'pending':
                if molecule.job_id not in pending:
                    logger.log(f"Job {molecule.job_id} is pending in the queue.")
                    pending.append(molecule.job_id)
                continue
            else:
                if molecule.job_id in pending:
                    pending.remove(molecule.job_id)
        if len(pending) >= max(1, int(len(molecules) / 1.5)):
            if len(pending) == len(molecules):
                msg = "All the submitted jobs are pending. Sleeping for now."
                sleep_time = 2*interval
            else: 
                msg = "Majority of submitted jobs are pending. Sleeping for now."
                sleep_time = interval
            if not sleeping:
                logger.log(msg)
                time.sleep(sleep_time)
                sleeping = True
            time.sleep(interval)
            continue

        try:
            files_in_directory = set(os.listdir(molecules[0].directory))
        except FileNotFoundError:
            logger.log(f"Pickle file(s) not generated yet. Retrying in {interval} seconds.")
            time.sleep(interval)
            attempts += 1
            continue

        if expected_files.issubset(files_in_directory):
            for molecule in molecules:
                try:
                    pickle_file_path = os.path.join(molecule.directory, f"collection{molecule.name}.pkl")
                    conformers = initiate_conformers(pickle_file_path)
                    logger.log_with_stars(f"{len(conformers)} conformers generated for {molecule.name.replace('_CREST', '')}")
                    for conf in conformers:
                        conf.constrained_indexes = constrained_indexes
                        conf.mult = mult
                        conf.charge = charge
                        conf.current_step = current_step
                        conf.directory = dir
                        conf.method = method
                        all_conformers.append(conf)
                    molecule.move_inputfile() #NOTE: test
                    molecule.move_converged()
                except Exception as e:
                    logger.log(f"Error processing molecule {molecule.name}: {e}")
                    return False
            logger.log(f"CREST conformers generated. Proceeding with next step: {all_conformers[0].next_step}")
            handle_termination(all_conformers, logger, threads, converged=True, args=args)
            return True
        else:
            if attempts == 1:
                logger.log(f"Not all files found. Retrying every {interval} seconds.")
            time.sleep(interval)

        attempts += 1

    return False


def update_molecules_status(molecules):
    # Gather all unique job IDs
    unique_job_ids = set(molecule.job_id.split("_")[0] for molecule in molecules)
    job_statuses = {}

    # Query the status of each unique job ID
    for job_id in unique_job_ids:
        try:
            result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True, check=True)
            job_statuses.update(parse_job_statuses(result.stdout, job_id))
        except subprocess.CalledProcessError:
            # If the command fails, assume 'unknown' status for all jobs with this prefix
            for molecule in molecules:
                if molecule.job_id.startswith(job_id):
                    molecule.status = 'unknown'
            continue  # Skip to the next job_id if there's an error

    # Update molecule status based on job_statuses dictionary
    for molecule in molecules:
        molecule.status = job_statuses.get(molecule.job_id, 'completed or not found')


def parse_job_statuses(output, main_job_id):
    job_statuses = {}
    lines = output.strip().split("\n")

    for line in lines[1:]:  # Skip the header line
        parts = line.split()
        if len(parts) > 4:
            job_id_field = parts[0]
            status = 'running' if 'R' in parts else 'pending' if 'PD' in parts else 'completed or not found'

            # Handle range of job IDs
            if '[' in job_id_field:
                range_match = re.search(r'\[(\d+)-(\d+)', job_id_field)
                if range_match:
                    start, end = range_match.groups(1)
                    for i in range(int(start), int(end) + 1):
                        job_statuses[f"{main_job_id}_{i}"] = status
            else:
                job_statuses[job_id_field] = status

    return job_statuses


def submit_and_monitor(molecules, logger, threads, args):
    if not isinstance(molecules, list):
        molecules = [molecules]

    if len(molecules) == 1:
        job_id, interval = submit_job(molecules[0], args)
        logger.log(f"Submitting file {molecules[0].name}{molecules[0].input} for calculation in path {molecules[0].directory} with job id {job_id}")
    else:
        job_id, interval = submit_array_job(molecules, args)
        logger.log(f"Submitted SLURM array job with job id {job_id} for conformers in {molecules[0].directory}")

    if job_id:
        if molecules[0].current_step == 'crest_sampling':
            # pass args into check_crest so it can call back into handle_termination with args
            thread = Thread(target=check_crest, args=(molecules, logger, threads, interval, args.attempts, args))
        else:
            # Pass args explicitly into the check_convergence thread
            thread = Thread(target=check_convergence, args=(molecules, logger, threads, interval, args.attempts, args))
        threads.append(thread)
        thread.start()
    else: 
        logger.log("Error getting job id")

        
def handle_termination(molecules, logger, threads, converged, args):
    if not isinstance(molecules, list):
        try:
            molecules = list(molecules)
        except TypeError:
            raise ValueError("molecules must be a list or convertible to a list")
    if converged:
        for m in molecules: 
            m.converged = False
            m.update_step()
    current_step = molecules[0].current_step
    QC_program = molecules[0].program
    if args.skip_preopt and current_step in ["opt_constrain", "opt_constrain_conf"]: # Works as long as only TS molecules are associated with the opt_constrain(_conf)
        for m in molecules:
            m.update_step() # Update step to skip preoptimization
        current_step = molecules[0].current_step
    logger.log(f"Job to be performed: {current_step}")
    if current_step in ['opt_constrain_conf','DLPNO'] and converged and 'H2O' not in molecules[0].name:
        if molecules[0].product: # For products ArbAlign is needed to be done on every individual H
            conformer_molecules = []
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in molecules if "_H" in m.name))
            grouped_lists = [[m for m in molecules if f"_H{h_num}_" in m.name] for h_num in h_numbers]
            for h, group in zip(h_numbers, grouped_lists):
                logger.log(f"Filtering product molecules for H{h} using ArbAlign alogrithm")
                filtered_group = filter_molecules(group, logger)
                for molecule in filtered_group:
                    molecule.converged = False
                    conformer_molecules.append(molecule)
        else:
            logger.log("Filtering molecules using ArbAlign algorithm")
            conformer_molecules = filter_molecules(molecules, logger)
            if len(conformer_molecules) >= 50:
                conformer_molecules = energy_cutoff(conformer_molecules, energy_cutoff=args.energy_cutoff)
    else:
        conformer_molecules = molecules

    for conf in conformer_molecules:
        if conf.converged is False:
            conf.name = conf.name.replace("_TS", "").replace("_CREST", "").replace("_DLPNO", "")
            job_type = conformer_molecules[0].current_step

            if job_type == 'crest_sampling':
                conf.name += '_CREST'
                conf.program = 'CREST'
                output_file_path = os.path.join(conf.directory, f"{conf.name}.xyz")
                conf.write_xyz_file(output_file_path)
                crest_constrain(conf, args)
            elif job_type in ['opt_constrain', 'opt_constrain_conf']:
                conf.program = QC_program
                QC_input(conf, constrain=True, TS=False, args=args)

            elif job_type in ['optimization', 'optimization_conf']:
                conf.program = QC_program
                QC_input(conf, constrain=False, TS=False, args=args)
            
            elif job_type in ['TS_opt', 'TS_opt_conf']:
                conf.program = QC_program
                conf.name += '_TS'
                QC_input(conf, constrain=False, TS=True, args=args)

            elif job_type == 'DLPNO':
                conf.program = 'ORCA'
                conf.method = 'DLPNO-CCSD(T)'
                conf.name += '_DLPNO'
                QC_input(conf, constrain=False, TS=False, args=args)

            elif job_type  == 'Done':
                logger.log("DLPNO calculation has converged")
                break

            else:
                logger.log("Job type could not be determined for conformer list")
                break
    if conformer_molecules: 
        submit_and_monitor(conformer_molecules, logger, threads, args)

def resubmit_job(molecule, logger, error=None, args=None):
    molecule.move_failed()
    job_type = molecule.current_step
    if job_type in ['opt_constrain', 'opt_constrain_conf']:
        QC_input(molecule, constrain=True, TS=False, args=args)

    elif job_type in ['optimization', 'optimization_conf']:
        QC_input(molecule, constrain=False, TS=False, args=args)

    elif job_type in ['TS_opt', 'TS_opt_conf']:
        QC_input(molecule, constrain=False, TS=True, args=args)

    elif job_type == 'DLPNO':
        QC_input(molecule, constrain=False, TS=False, args=args)

    else:
        logger.log(f"Error determining job type for resubmission of {molecule.name}")
    job_id, _ = submit_job(molecule, args)
    molecule.job_id = f"{job_id}"
    logger.log(f"submitted file {molecule.name}{molecule.input} with job type: {job_type} and job id {molecule.job_id} in directory {molecule.directory}")


def check_convergence(molecules, logger, threads, interval, max_attempts, args, all_converged=False):
    initial_delay = int(interval * 2)
    attempts = 0
    sleeping = False
    max_terminations_allowed = 2
    pending, running = [], []
    job_type = molecules[0].current_step

    for m in molecules:  # Initialize with all molecules not being converged and no terminations counted
        m.converged = False
        m.error_termination_count = 0
        m.log_file_path = os.path.join(m.directory, f"{m.name}{m.output}")

    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)

    while attempts < max_attempts and not all_converged:
        i = 0
        while i < len(molecules):
            molecule = molecules[i]
            if molecule.error_termination_count >= max_terminations_allowed:
                logger.log(f"!!! Dropping molecule conformer {molecule.name} due to repeated error terminations!!!")
                molecules.pop(i)
                molecule.move_failed()
                # Check if all remaining molecules are converged
                if all(m.converged for m in molecules):
                    all_converged = True
                    break
            else:
                i += 1
        if all_converged:
            break

        update_molecules_status(molecules)

        for molecule in molecules:
            if molecule.converged:
                continue
            if molecule.status == 'pending':
                if molecule.job_id not in pending:
                    logger.log(f"Job {molecule.job_id} is pending in the queue.")
                    pending.append(molecule.job_id)
                continue

            elif molecule.status in ['running', 'completed or not found'] or not molecule.converged:
                if molecule.job_id not in running:
                    logger.log(f"Job {job_type} for {molecule.name} with job id {molecule.job_id} is running.") if molecule.status == 'running' else None
                    molecule.move_inputfile()
                    running.append(molecule.job_id)
                if molecule.job_id in pending:
                    pending.remove(molecule.job_id)

                normal_termination_detected, termination_string = termination_status(molecule, logger, args)

                if normal_termination_detected:
                    logger.log_with_stars(termination_string) if 'Yay' in termination_string else logger.log(termination_string)
                    molecule.converged = True
                    molecule.error_termination_count = 0
                elif normal_termination_detected is False:
                    if termination_string in ['Error string not found', 'Could not read molecule log file']:
                        continue
                    else:
                        molecule.error_termination_count += 1
                        if molecule.error_termination_count >= max_terminations_allowed:
                            continue
                        handle_error_termination(molecule, logger, termination_string, args)
            else:
                logger.log(f"Status {molecule} could not be determined. Ensure it is running. Job id: {molecule.job_id}")
                continue
            time.sleep(10)
                
        if all(m.converged for m in molecules):
            all_converged = True
            break

        if len(pending) >= max(1, int(len(molecules) / 1.5)):
            if len(pending) == len(molecules):
                msg = "All the submitted jobs are pending. Sleeping for now."
                sleep_time = 2*interval
            else: 
                msg = "Majority of submitted jobs are pending. Sleeping for now."
                sleep_time = interval
            if not sleeping:
                logger.log(msg)
                time.sleep(sleep_time)
                sleeping = True
            time.sleep(interval)
        else:
            attempts += 1
            if attempts % 10 == 0 or attempts == 1:
                logger.log(f"Log files of the {len(molecules)} conformers have been checked. Checking every {interval} seconds. Attempt: {attempts}/{max_attempts}")
            time.sleep(interval)

    if all_converged:
        if molecules: # Check if not all molecules from list has been dropped
            dir = molecules[0].directory
            basename = os.path.basename(dir)
            pickle_path = os.path.join(dir, f'{basename}_{job_type}.pkl')
            molecules[0].move_files()
            if len(molecules) > 1:
                Molecule.molecules_to_pickle(molecules, pickle_path)
            logger.log_with_stars(f"Yay! All conformer jobs have converged for job type: {job_type}.")
            if job_type == "DLPNO":
                for molecule in molecules: 
                    global_molecules.append(molecule)
                gm_list = global_molecules.as_list()
                # Determine QC program from args if provided, otherwise fallback to config
                qc_program = getattr(args, 'QC_program', None) or config.QC_program
                if molecules[0].product and not any('H2O' in mol.name for mol in gm_list):
                    H2O = Molecule.create_H2O()
                    H2O.program = qc_program
                    QC_input(H2O, constrain=False, TS=False, args=args)
                    submit_and_monitor(H2O, logger, threads, args)
                elif molecules[0].reactant and not any('OH' in mol.name for mol in gm_list):
                    OH = Molecule.create_OH()
                    OH.program = qc_program
                    QC_input(OH, constrain=False, TS=False, args=args)
                    submit_and_monitor(OH, logger, threads, args)
                return True
            elif getattr(args, 'auto', False):
                handle_termination(molecules, logger, threads, converged=True, args=args)
                return True
            else:
                logger.log(f"Automatic processing of calculations has been turned off and JKTS will not proceed with next step: {molecules[0].next_step}")
                logger.log(f"A pickle file {basename}_{job_type}.pkl has been created. Next step {molecules[0].next_step} for molecules can be started from this.")
                return True
        else:
            logger.log(f"No conformer has managed to converge for job type: {job_type}")
            print("Program terminated")
            exit()


def handle_error_termination(molecule, logger, error_termination_string, args=None):
    # Ensure args is available; fall back to global config
    args = args or config.args
    if molecule.program.lower() == 'orca':
        logger.log(f"Error termination was found in {molecule.name}. Resubmitting")
        xyz_coordinates = molecule.log2xyz()
        molecule.coordinates = xyz_coordinates
        resubmit_job(molecule, logger, args=args)
    else:
        convergence_errors = ["l9999", "l508"]
        intervention_errors = ["l301"]
        G16_common_errors = "https://wongzit.github.io/gaussian-common-errors-and-solutions/"
        last_lines = read_last_lines(molecule.log_file_path, 30)

        if last_lines:
            detected_convergence_errors = [error for error in convergence_errors if any(error in line for line in last_lines)]
            detected_intervention_errors = [error for error in intervention_errors if any(error in line for line in last_lines)]

            if detected_convergence_errors:
                error = detected_convergence_errors[0]
                logger.log(f"Convergence error '{error}' termination found in {molecule.name}")
                xyz_coordinates = molecule.log2xyz()
                molecule.coordinates = xyz_coordinates
                logger.log(f"Trying to resubmit job for {molecule.name} due to error termination")
                resubmit_job(molecule, logger, error, args=args)
            elif detected_intervention_errors:
                error = detected_intervention_errors[0]
                logger.log(f"Error '{error}' detected in {molecule.name} which needs taken care of manually")
                logger.log(f"Removing the conformer {molecule.name} for now so user can inspect error")
                molecule.error_termination_count = 3
                if molecule.program.lower() == 'g16':
                    logger.log(f"Common G16 errors can be found in: {G16_common_errors}")
            elif 'TS' in molecule.current_step and good_active_site(molecule) or "l801" in last_lines:
                logger.log("SCF didnt seem to converge due to offset of OH radical. Trying to correct and resubmit.")
                molecule.set_active_site(indexes=args.CHO)
                molecule.perturb_active_site(indexes=args.CHO)
                resubmit_job(molecule, logger, args=args)
            else:
                logger.log(f"Error termination found in {molecule.name}. Trying to resubmit")
                xyz_coordinates = molecule.log2xyz()
                if xyz_coordinates:
                    molecule.coordinates = xyz_coordinates
                    resubmit_job(molecule, logger, args=args)
                else:
                    logger.log(f"Geometry coordiantes could not be found in {molecule.log_file_path}. Dropping molecule.")
                    molecule.error_termination_count = 3
        else:
            logger.log(f"Could not find log file in path: {molecule.log_file_path}")
            logger.log(f"Dropping molecule {molecule.name}")
        

def termination_status(molecule, logger, args=None):
    last_lines = read_last_lines(molecule.log_file_path, 30)
    if not last_lines:
        molecule.error_termination_count += 1
        return False, 'Could not read molecule log file'

    termination_detected = any(
        termination in line.lower() for termination in termination_strings[molecule.program] for line in last_lines)
    
    if termination_detected:
        # Update geometry and energy
        xyz_coordinates = molecule.log2xyz()
        molecule.coordinates = xyz_coordinates
        molecule.update_energy(logger)
        if 'TS' in molecule.current_step:
            ts_check_passed, msg = check_transition_state(molecule, args=args)
            return ts_check_passed, msg
        else:
            return True, f"***{molecule.name} converged***"

    for error_termination in error_strings[molecule.program]:
        if any(error_termination in line.lower() for line in last_lines):
            return False, error_termination

    return False, 'Error string not found'