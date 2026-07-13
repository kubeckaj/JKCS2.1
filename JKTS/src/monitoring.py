import os
import re
import sys
import time
from threading import Thread

from classes import Molecule, Step
from slurm_submit import submit_job, submit_array_job, update_molecules_status
from qc_input import QC_input, crest_constrain
from conformer_tools import filter_molecules, energy_cutoff, initiate_conformers
from ts_validation import check_transition_state, good_active_site
import runtime

# Workflow steps after which duplicate conformers are filtered out
FILTER_STEPS = ('opt_constrain_conf', 'DLPNO')


def read_last_lines(filename, num_lines, interval=200):
    attempts = 0
    max_attempts = 5

    while attempts < max_attempts:
        try:
            with open(filename, 'rb') as f:
                f.seek(0, os.SEEK_END)
                buffer_size = 8192
                offset = f.tell()
                content = ''
                while len(content.splitlines()) < num_lines + 1 and offset > 0:
                    offset = max(0, offset - buffer_size)
                    f.seek(offset)
                    content = f.read().decode('utf-8', errors='replace')
                return content.splitlines()[-num_lines:]
        except FileNotFoundError:
            attempts += 1
            time.sleep(interval)

    return []


def resubmit_job(molecule, logger, error=None):
    molecule.move_failed()
    job_type = molecule.current_step
    if job_type in ['opt_constrain', 'opt_constrain_conf']:
        QC_input(molecule, constrain=True, TS=False)
    elif job_type in ['optimization', 'optimization_conf']:
        QC_input(molecule, constrain=False, TS=False)
    elif job_type in ['TS_opt', 'TS_opt_conf']:
        QC_input(molecule, constrain=False, TS=True)
    elif job_type == 'DLPNO':
        QC_input(molecule, constrain=False, TS=False)
    else:
        logger.warning(f"Could not determine job type for resubmission of {molecule.name}")

    job_id, _ = submit_job(molecule, runtime.args)
    molecule.job_id = f"{job_id}"
    logger.event(f"Resubmitted {molecule.name}{molecule.input} ({job_type}, job id {molecule.job_id})")


def termination_status(molecule, logger):
    last_lines = read_last_lines(molecule.log_file_path, 30)
    if not last_lines:
        molecule.error_termination_count += 1
        return False, 'Could not read molecule log file'

    termination_detected = any(
        termination in line.lower()
        for termination in runtime.termination_strings[molecule.program]
        for line in last_lines
    )

    if termination_detected:
        xyz_coordinates = molecule.log2xyz()
        if not xyz_coordinates:
            return False, 'Termination detected but no geometry found in log file'
        molecule.coordinates = xyz_coordinates
        molecule.update_energy(logger)
        if 'TS' in molecule.current_step:
            ts_check_passed, msg = check_transition_state(molecule)
            if ts_check_passed:
                return True, msg
            return None, msg  # converged normally but failed TS validation
        else:
            return True, f"{molecule.name} converged"

    for error_termination in runtime.error_strings[molecule.program]:
        if any(error_termination in line.lower() for line in last_lines):
            return False, error_termination

    return False, 'Error string not found'


def handle_error_termination(molecule, logger, error_termination_string):
    if molecule.program.lower() == 'orca':
        logger.warning(f"Error termination found in {molecule.name}. Resubmitting")
        xyz_coordinates = molecule.log2xyz()
        if xyz_coordinates:
            molecule.coordinates = xyz_coordinates
            resubmit_job(molecule, logger)
        else:
            logger.error(f"Geometry coordinates could not be found in {molecule.log_file_path}. Dropping molecule.")
            molecule.error_termination_count = 3
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
                logger.warning(f"Convergence error '{error}' found in {molecule.name}. Resubmitting")
                xyz_coordinates = molecule.log2xyz()
                if xyz_coordinates:
                    molecule.coordinates = xyz_coordinates
                    resubmit_job(molecule, logger, error)
                else:
                    logger.error(f"Geometry coordinates could not be found in {molecule.log_file_path}. Dropping molecule.")
                    molecule.error_termination_count = 3
            elif detected_intervention_errors:
                error = detected_intervention_errors[0]
                logger.error(f"Error '{error}' detected in {molecule.name} needs manual attention. Removing the conformer so the error can be inspected")
                molecule.error_termination_count = 3
                if molecule.program.lower() == 'g16':
                    logger.info(f"Common G16 errors can be found in: {G16_common_errors}")
            elif 'TS' in molecule.current_step and good_active_site(molecule) or "l801" in last_lines:
                logger.warning(f"SCF did not converge for {molecule.name}, likely due to offset of the abstracting radical. Correcting geometry and resubmitting")
                molecule.set_active_site(indexes=runtime.args.CHO)
                molecule.perturb_active_site(indexes=runtime.args.CHO)
                resubmit_job(molecule, logger)
            else:
                logger.warning(f"Error termination found in {molecule.name}. Resubmitting")
                xyz_coordinates = molecule.log2xyz()
                if xyz_coordinates:
                    molecule.coordinates = xyz_coordinates
                    resubmit_job(molecule, logger)
                else:
                    logger.error(f"Geometry coordinates could not be found in {molecule.log_file_path}. Dropping molecule.")
                    molecule.error_termination_count = 3
        else:
            logger.error(f"Could not find log file {molecule.log_file_path}. Dropping molecule {molecule.name}")


def check_convergence(molecules, logger, threads, interval, max_attempts, all_converged=False):
    initial_delay = runtime.args.initial_delay if runtime.args.initial_delay else int(interval * 2)
    interval = runtime.args.interval if runtime.args.interval else int(interval)
    attempts = 0
    sleeping = False
    max_terminations_allowed = 2
    last_status = {}
    job_type = molecules[0].current_step

    for m in molecules:
        m.converged = False
        m.error_termination_count = 0
        m.log_file_path = os.path.join(m.directory, f"{m.name}{m.output}")

    logger.info(f"Monitoring {len(molecules)} {job_type} job(s); first check in {initial_delay} seconds, then every {interval} seconds")
    time.sleep(initial_delay)

    while attempts < max_attempts and not all_converged:
        i = 0
        while i < len(molecules):
            molecule = molecules[i]
            if molecule.error_termination_count >= max_terminations_allowed:
                logger.error(f"Dropping conformer {molecule.name} due to repeated error terminations")
                molecules.pop(i)
                molecule.move_failed()
                if all(m.converged for m in molecules):
                    all_converged = True
                    break
            else:
                i += 1
        if all_converged:
            break

        update_molecules_status(molecules)

        # Report queue-status transitions once, aggregated over the batch
        active = [m for m in molecules if not m.converged]
        newly_running = [m for m in active if m.status == 'running' and last_status.get(m.name) != 'running']
        if not last_status and active:
            n_pending = sum(1 for m in active if m.status == 'pending')
            logger.info(f"Queue status: {n_pending} pending, {len(active) - n_pending} running/finished of {len(active)} job(s) (job id {active[0].job_id})")
        if newly_running:
            if len(newly_running) == 1:
                logger.event(f"{job_type} running for {newly_running[0].name} (job id {newly_running[0].job_id})")
            else:
                logger.event(f"{job_type} running for {len(newly_running)} conformers (job id {newly_running[0].job_id})")
        for molecule in active:
            if last_status.get(molecule.name, 'pending') == 'pending' and molecule.status != 'pending':
                molecule.move_inputfile()
        last_status = {m.name: m.status for m in molecules}

        pending_count = sum(1 for m in active if m.status == 'pending')

        for molecule in active:
            if molecule.status == 'pending':
                continue

            normal_termination_detected, termination_string = termination_status(molecule, logger)

            if normal_termination_detected is True:
                logger.success(termination_string)
                molecule.converged = True
                molecule.error_termination_count = 0
            elif normal_termination_detected is None:
                # Job converged normally but failed TS geometry/frequency validation
                if "Trying to correct" in termination_string and molecule.error_termination_count == 0:
                    logger.warning(termination_string)
                    molecule.error_termination_count += 1
                    resubmit_job(molecule, logger)
                else:
                    logger.error(f"TS validation failed for {molecule.name}: {termination_string} Dropping molecule.")
                    molecule.error_termination_count = max_terminations_allowed
            elif normal_termination_detected is False:
                if termination_string in ['Error string not found', 'Could not read molecule log file']:
                    continue
                else:
                    molecule.error_termination_count += 1
                    if molecule.error_termination_count >= max_terminations_allowed:
                        continue
                    handle_error_termination(molecule, logger, termination_string)
            time.sleep(10)

        if all(m.converged for m in molecules):
            all_converged = True
            break

        if pending_count >= max(1, int(len(molecules) / 1.5)):
            if pending_count == len(molecules):
                msg = "All submitted jobs are pending; waiting for the queue."
                sleep_time = 2*interval
            else:
                msg = "Most submitted jobs are pending; waiting for the queue."
                sleep_time = interval
            if not sleeping:
                logger.info(msg)
                time.sleep(sleep_time)
                sleeping = True
            time.sleep(interval)
        else:
            sleeping = False
            attempts += 1
            if attempts % 10 == 0 or attempts == 1:
                logger.info(f"Checked log files of {len(molecules)} conformer(s); polling every {interval} seconds (attempt {attempts}/{max_attempts})")
            time.sleep(interval)

    if all_converged:
        if molecules:
            dir = molecules[0].directory
            basename = os.path.basename(dir)
            pickle_path = os.path.join(dir, f'{basename}_{job_type}.pkl')
            molecules[0].move_files()
            if len(molecules) > 1:
                Molecule.molecules_to_pickle(molecules, pickle_path)
            logger.success(f"All {len(molecules)} job(s) converged for step {job_type}")
            if job_type == "DLPNO":
                with runtime.global_molecules_lock:
                    for molecule in molecules:
                        runtime.global_molecules.append(molecule)
                    if runtime.args.Cl:
                        needs_small_molecule = molecules[0].product and not any(mol.name in ('HCl', 'HCl_DLPNO') for mol in runtime.global_molecules)
                        needs_radical = molecules[0].reactant and not any(mol.name in ('Cl', 'Cl_DLPNO') for mol in runtime.global_molecules)
                    elif runtime.args.NO3:
                        needs_small_molecule = molecules[0].product and not any(mol.name in ('HNO3', 'HNO3_DLPNO') for mol in runtime.global_molecules)
                        needs_radical = molecules[0].reactant and not any(mol.name in ('NO3', 'NO3_DLPNO') for mol in runtime.global_molecules)
                    else:
                        needs_small_molecule = molecules[0].product and not any('H2O' in mol.name for mol in runtime.global_molecules)
                        needs_radical = molecules[0].reactant and not any('OH' in mol.name for mol in runtime.global_molecules)
                if needs_small_molecule:
                    if runtime.args.Cl:
                        small_mol = Molecule.create_HCl()
                    elif runtime.args.NO3:
                        small_mol = Molecule.create_HNO3()
                    else:
                        small_mol = Molecule.create_H2O()
                    small_mol.program = runtime.QC_program
                    QC_input(small_mol, constrain=False, TS=False)
                    submit_and_monitor(small_mol, logger, threads)
                elif needs_radical:
                    if runtime.args.Cl:
                        radical = Molecule.create_Cl()
                    elif runtime.args.NO3:
                        radical = Molecule.create_NO3()
                    else:
                        radical = Molecule.create_OH()
                    radical.program = runtime.QC_program
                    QC_input(radical, constrain=False, TS=False)
                    submit_and_monitor(radical, logger, threads)
                return True
            elif not runtime.args.stop:
                handle_termination(molecules, logger, threads, converged=True)
                return True
            else:
                logger.event(f"Stopping as requested (-stop); next step would be {molecules[0].next_step}. Resume from {basename}_{job_type}.pkl with -restart.")
                return True
        else:
            logger.error(f"No conformer managed to converge for step {job_type}. Terminating.")
            sys.exit(1)


def check_crest(molecules, logger, threads, interval, max_attempts):
    initial_delay = runtime.args.initial_delay if runtime.args.initial_delay else int(interval * 2)
    interval = runtime.args.interval if runtime.args.interval else int(interval)
    attempts = 0
    sleeping = False
    last_pending = None
    all_conformers = []
    expected_files = {f"collection{molecule.name}.pkl" for molecule in molecules}
    constrained_indexes = molecules[0].constrained_indexes
    mult = molecules[0].mult
    charge = molecules[0].charge
    current_step = molecules[0].current_step
    dir = molecules[0].directory
    method = molecules[0].method

    logger.info(f"Monitoring CREST sampling; first check in {initial_delay} seconds, then every {interval} seconds")
    time.sleep(initial_delay)

    while attempts < max_attempts:
        update_molecules_status(molecules)
        pending_count = sum(1 for m in molecules if m.status == 'pending')
        if last_pending is None and pending_count:
            logger.info(f"Queue status: {pending_count} of {len(molecules)} CREST job(s) pending (job id {molecules[0].job_id})")
        elif last_pending and pending_count < last_pending:
            logger.event(f"CREST sampling running (job id {molecules[0].job_id})")
        last_pending = pending_count
        if pending_count >= max(1, int(len(molecules) / 1.5)):
            if pending_count == len(molecules):
                msg = "All submitted jobs are pending; waiting for the queue."
                sleep_time = 2*interval
            else:
                msg = "Most submitted jobs are pending; waiting for the queue."
                sleep_time = interval
            if not sleeping:
                logger.info(msg)
                time.sleep(sleep_time)
                sleeping = True
            time.sleep(interval)
            continue
        sleeping = False

        try:
            files_in_directory = set(os.listdir(molecules[0].directory))
        except FileNotFoundError:
            logger.info(f"Pickle file(s) not generated yet. Retrying in {interval} seconds.")
            time.sleep(interval)
            attempts += 1
            continue

        if expected_files.issubset(files_in_directory):
            for molecule in molecules:
                try:
                    pickle_file_path = os.path.join(molecule.directory, f"collection{molecule.name}.pkl")
                    conformers = initiate_conformers(pickle_file_path)
                    logger.success(f"CREST sampling done: {len(conformers)} conformers generated for {molecule.name.replace('_CREST', '')}")
                    for conf in conformers:
                        conf.constrained_indexes = constrained_indexes
                        conf.mult = mult
                        conf.charge = charge
                        conf.current_step = current_step
                        conf.directory = dir
                        conf.method = method
                        all_conformers.append(conf)
                    molecule.move_inputfile()
                    molecule.move_converged()
                except Exception as e:
                    logger.error(f"Error processing molecule {molecule.name}: {e}")
                    return False
            handle_termination(all_conformers, logger, threads, converged=True)
            return True
        else:
            if attempts == 1:
                logger.info(f"CREST results not complete yet. Retrying every {interval} seconds.")
            time.sleep(interval)

        attempts += 1

    return False


def TS_search(molecule):
    import numpy as np
    from copy import deepcopy
    increment = 0.1
    min_CH, max_CH = 1.2, 1.6
    min_HO, max_HO = 1.1, 1.6
    min_angle, max_angle = 170, 179
    i = 1
    for distance_CH in np.arange(min_CH, max_CH, increment):
        for distance_HO in np.arange(min_HO, max_HO, increment):
            for reaction_angle in range(max_angle, min_angle-1, -5):
                m_copy = deepcopy(molecule)
                m_copy.name += f"_{i}"
                i += 1
                m_copy.set_active_site(indexes=runtime.args.CHO, distance_CH=distance_CH, distance_HO=distance_HO, reaction_angle=reaction_angle)
                QC_input(m_copy, constrain=False, TS=True)


# ---------------------------------------------------------------------------
# Step handlers — one function per workflow step, called by handle_termination
# ---------------------------------------------------------------------------
def _step_crest(conf):
    conf.name += '_CREST'
    conf.program = 'CREST'
    conf.write_xyz_file(os.path.join(conf.directory, f"{conf.name}.xyz"))
    crest_constrain(conf)


def _step_opt_constrain(conf):
    conf.program = runtime.QC_program
    QC_input(conf, constrain=True, TS=False)


def _step_optimization(conf):
    conf.program = runtime.QC_program
    QC_input(conf, constrain=False, TS=False)


def _step_ts_opt(conf):
    conf.program = runtime.QC_program
    conf.name += '_TS'
    QC_input(conf, constrain=False, TS=True)


def _step_dlpno(conf):
    conf.program = 'ORCA'
    conf.method = 'DLPNO-CCSD(T)'
    conf.name += '_DLPNO'
    QC_input(conf, constrain=False, TS=False)


STEP_HANDLERS = {
    Step.CREST_SAMPLING:     _step_crest,
    Step.OPT_CONSTRAIN:      _step_opt_constrain,
    Step.OPT_CONSTRAIN_CONF: _step_opt_constrain,
    Step.OPTIMIZATION:       _step_optimization,
    Step.TS_OPT:             _step_ts_opt,
    Step.TS_OPT_CONF:        _step_ts_opt,
    Step.DLPNO:              _step_dlpno,
}
# ---------------------------------------------------------------------------


def submit_and_monitor(molecules, logger, threads):
    if isinstance(molecules, Molecule):
        molecules = [molecules]

    if len(molecules) == 1:
        job_id, interval = submit_job(molecules[0], runtime.args)
        logger.event(f"Submitted {molecules[0].name}{molecules[0].input} ({molecules[0].current_step}, job id {job_id})")
    else:
        job_id, interval = submit_array_job(molecules, runtime.args)
        logger.event(f"Submitted SLURM array job for {len(molecules)} conformers ({molecules[0].current_step}, job id {job_id})")

    if job_id:
        if molecules[0].current_step == 'crest_sampling':
            thread = Thread(target=check_crest, args=(molecules, logger, threads, interval, runtime.args.attempts))
        else:
            thread = Thread(target=check_convergence, args=(molecules, logger, threads, interval, runtime.args.attempts))
        threads.append(thread)
        thread.start()
    else:
        logger.error("Could not get a job id from SLURM; the job was not submitted correctly")


def handle_termination(molecules, logger, threads, converged):
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
    logger.event(f"Starting step: {current_step}")
    if current_step in FILTER_STEPS and converged and not molecules[0].small_molecule:
        if molecules[0].product:
            conformer_molecules = []
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in molecules if "_H" in m.name))
            grouped_lists = [[m for m in molecules if f"_H{h_num}_" in m.name] for h_num in h_numbers]
            for h, group in zip(h_numbers, grouped_lists):
                logger.info(f"Filtering product molecules for H{h} using the ArbAlign algorithm")
                filtered_group = filter_molecules(group, logger)
                for molecule in filtered_group:
                    molecule.converged = False
                    conformer_molecules.append(molecule)
        else:
            logger.info("Filtering molecules using the ArbAlign algorithm")
            conformer_molecules = filter_molecules(molecules, logger)
            if len(conformer_molecules) >= 50:
                conformer_molecules = energy_cutoff(conformer_molecules)
    else:
        conformer_molecules = molecules

    job_type = conformer_molecules[0].current_step
    handler = STEP_HANDLERS.get(job_type)
    if handler is None:
        if job_type == 'Done':
            logger.success("DLPNO calculation has converged")
        else:
            logger.warning(f"Job type '{job_type}' could not be determined for conformer list")
        return

    for conf in conformer_molecules:
        if conf.converged is False:
            conf.name = conf.name.replace("_TS", "").replace("_CREST", "").replace("_DLPNO", "")
            handler(conf)
    if conformer_molecules:
        submit_and_monitor(conformer_molecules, logger, threads)


def handle_input_molecules(molecules, logger, threads):
    current_step = molecules[0].current_step
    if all(m.current_step == current_step for m in molecules):
        logger.event(f"Detected {len(molecules)} molecule(s) at step {current_step}; checking convergence status")
        if current_step == 'crest_sampling':
            all_converged = True
        else:
            for m in molecules:
                if os.path.exists(m.log_file_path):
                    converge_status, _ = termination_status(m, logger)
                    if converge_status is True:
                        m.converged = True
                    else:
                        m.converged = False
                else:
                    if len(molecules) == 1:
                        m.converged = False
                    elif m.atoms and m.coordinates:
                        m.converged = True
            all_converged = all(m.converged for m in molecules)

        if all_converged:
            if all(mol.current_step == 'DLPNO' for mol in molecules):
                if all(mol.product for mol in molecules):
                    if runtime.args.Cl:
                        product_missing = not any(mol.name in ('HCl', 'HCl_DLPNO') for mol in molecules)
                    elif runtime.args.NO3:
                        product_missing = not any(mol.name in ('HNO3', 'HNO3_DLPNO') for mol in molecules)
                    else:
                        product_missing = not any('H2O' in mol.name for mol in molecules)
                    if product_missing:
                        if runtime.args.Cl:
                            small_mol_name = "HCl"
                        elif runtime.args.NO3:
                            small_mol_name = "HNO3"
                        else:
                            small_mol_name = "H2O"
                        logger.event(f"Converged DLPNOs. However, {small_mol_name} is needed for product energies")
                        check_convergence(molecules, logger, threads, 30, runtime.args.attempts, all_converged=True)
                elif all(mol.reactant for mol in molecules):
                    if runtime.args.Cl:
                        reactant_missing = not any(mol.name in ('Cl', 'Cl_DLPNO') for mol in molecules)
                    elif runtime.args.NO3:
                        reactant_missing = not any(mol.name in ('NO3', 'NO3_DLPNO') for mol in molecules)
                    else:
                        reactant_missing = not any('OH' in mol.name for mol in molecules)
                    if reactant_missing:
                        if runtime.args.Cl:
                            radical_name = "Cl"
                        elif runtime.args.NO3:
                            radical_name = "NO3"
                        else:
                            radical_name = "OH"
                        logger.event(f"Converged DLPNOs. However, {radical_name} is needed for reactant energies")
                        check_convergence(molecules, logger, threads, 30, runtime.args.attempts, all_converged=True)
                else:
                    logger.success("All given molecules are converged DLPNOs")
                    for m in molecules:
                        runtime.global_molecules.append(m)
            else:
                logger.success(f"All molecules converged for step {current_step}. Proceeding with next step: {molecules[0].next_step}")
                molecules[0].move_files()
                handle_termination(molecules, logger, threads, converged=True)
        else:
            logger.event("Not all molecules have converged. Non-converged ones will be calculated; converged ones will be skipped.")
            handle_termination(molecules, logger, threads, converged=False)
    else:
        logger.error("Not all molecules in the given files are at the same workflow step. Check that the correct input is provided.")
        sys.exit(1)
