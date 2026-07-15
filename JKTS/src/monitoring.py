import os
import re
import sys
import time
from threading import Thread

from classes import Molecule, Step
from slurm_submit import submit_job, submit_array_job, update_molecules_status, get_interval_seconds
from qc_input import QC_input, crest_constrain
from conformer_tools import filter_molecules, energy_cutoff, initiate_conformers
from ts_validation import check_transition_state, good_active_site
import checkpoint
import runtime

# Workflow steps after which duplicate conformers are filtered out
FILTER_STEPS = ('opt_constrain_conf', 'DLPNO')

VANISHED_GRACE_POLLS = 3
MAX_NODE_FAILURES = 3


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


def termination_status(molecule, logger, quick=False):
    # quick: don't wait minutes for a log that may never appear (vanished job,
    # restart reconciliation) — read_last_lines retries 5x at `interval` seconds.
    last_lines = read_last_lines(molecule.log_file_path, 30, interval=5 if quick else 200)
    if not last_lines:
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
        molecule.update_energy()
        if 'TS' in molecule.current_step:
            ts_check_passed, msg = check_transition_state(molecule)
            if ts_check_passed:
                return True, msg
            return None, msg  # converged normally but failed TS validation
        else:
            return True, f"*{molecule.name} converged*"

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
    missing_polls = {}
    job_type = molecules[0].current_step

    # Only reset state for molecules that still have work to do: on a restart
    # reattach, already-converged molecules ride along untouched (their logs
    # may have been moved to log_files/ already).
    for m in molecules:
        if m.converged:
            continue
        m.error_termination_count = 0
        m.node_failure_count = getattr(m, 'node_failure_count', 0)
        m.log_file_path = os.path.join(m.directory, f"{m.name}{m.output}")

    logger.info(f"Monitoring {len(molecules)} {job_type} job(s); first check in {initial_delay} seconds, then every {interval} seconds")
    time.sleep(initial_delay)

    while attempts < max_attempts and not all_converged:
        dirty = False  # any state change this sweep is checkpointed at the end
        i = 0
        while i < len(molecules):
            molecule = molecules[i]
            if molecule.error_termination_count >= max_terminations_allowed:
                logger.error(f"Dropping conformer {molecule.name} due to repeated error terminations")
                molecules.pop(i)
                molecule.move_failed()
                checkpoint.save_checkpoint(molecules)
                if all(m.converged for m in molecules):
                    all_converged = True
                    break
            else:
                i += 1
        if all_converged:
            break

        active = [m for m in molecules if not m.converged]
        update_molecules_status(active)
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

            job_gone = molecule.status == 'completed or not found'
            normal_termination_detected, termination_string = termination_status(molecule, logger, quick=job_gone)

            if normal_termination_detected is True:
                logger.success(termination_string)
                molecule.converged = True
                molecule.error_termination_count = 0
                missing_polls.pop(molecule.name, None)
                dirty = True
            elif normal_termination_detected is None:
                # Job converged normally but failed TS geometry/frequency validation
                if "Trying to correct" in termination_string and molecule.error_termination_count == 0:
                    logger.warning(termination_string)
                    molecule.error_termination_count += 1
                    resubmit_job(molecule, logger)
                else:
                    logger.error(f"TS validation failed for {molecule.name}: {termination_string} Dropping molecule.")
                    molecule.error_termination_count = max_terminations_allowed
                dirty = True
            elif normal_termination_detected is False:
                if termination_string in ['Error string not found', 'Could not read molecule log file']:
                    # No termination string yet: either still running, or the job
                    # died without a trace (node failure). Only the latter — job
                    # gone from squeue for several consecutive polls — triggers
                    # a resubmission; 'unknown' (squeue outage) never does.
                    if job_gone:
                        missing_polls[molecule.name] = missing_polls.get(molecule.name, 0) + 1
                        if missing_polls[molecule.name] >= VANISHED_GRACE_POLLS:
                            missing_polls[molecule.name] = 0
                            molecule.node_failure_count = getattr(molecule, 'node_failure_count', 0) + 1
                            if molecule.node_failure_count > MAX_NODE_FAILURES:
                                logger.error(f"{molecule.name}: job vanished {MAX_NODE_FAILURES} times without finishing; dropping conformer")
                                molecule.error_termination_count = max_terminations_allowed
                            else:
                                logger.warning(f"{molecule.name}: job {molecule.job_id} disappeared from the queue without finishing (node failure?). Resubmitting.")
                                resubmit_job(molecule, logger)
                            dirty = True
                    else:
                        missing_polls.pop(molecule.name, None)
                    continue
                else:
                    molecule.error_termination_count += 1
                    dirty = True
                    if molecule.error_termination_count >= max_terminations_allowed:
                        continue
                    handle_error_termination(molecule, logger, termination_string)
            time.sleep(10)

        if all(m.converged for m in molecules):
            all_converged = True
            break

        if dirty:
            checkpoint.save_checkpoint(molecules)

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
            Molecule.molecules_to_pickle(molecules, pickle_path)
            logger.success(f"All {len(molecules)} job(s) converged for step {job_type}")
            checkpoint.save_checkpoint(molecules)
            if job_type == "DLPNO":
                with runtime.global_molecules_lock:
                    for molecule in molecules:
                        runtime.global_molecules.append(molecule)
                checkpoint.save_checkpoint(molecules)
                submit_missing_partner(molecules, logger, threads)
                return True
            elif not runtime.args.stop:
                handle_termination(molecules, logger, threads, converged=True)
                return True
            else:
                logger.event(f"Stopping as requested (-stop); next step would be {molecules[0].next_step}. Resume from {basename}_checkpoint.pkl (or {basename}_{job_type}.pkl) with -restart.")
                return True
        else:
            logger.error(f"No conformer managed to converge for step {job_type}. Terminating.")
            sys.exit(1)


def submit_missing_partner(batch, logger, threads, known=None):
    if known is None:
        with runtime.global_molecules_lock:
            known = list(runtime.global_molecules)
    if runtime.args.Cl:
        needs_small_molecule = batch[0].product and not any(mol.name in ('HCl', 'HCl_DLPNO') for mol in known)
        needs_radical = batch[0].reactant and not any(mol.name in ('Cl', 'Cl_DLPNO') for mol in known)
    elif runtime.args.NO3:
        needs_small_molecule = batch[0].product and not any(mol.name in ('HNO3', 'HNO3_DLPNO') for mol in known)
        needs_radical = batch[0].reactant and not any(mol.name in ('NO3', 'NO3_DLPNO') for mol in known)
    else:
        needs_small_molecule = batch[0].product and not any('H2O' in mol.name for mol in known)
        needs_radical = batch[0].reactant and not any('OH' in mol.name for mol in known)

    if needs_small_molecule:
        if runtime.args.Cl:
            partner = Molecule.create_HCl()
        elif runtime.args.NO3:
            partner = Molecule.create_HNO3()
        else:
            partner = Molecule.create_H2O()
    elif needs_radical:
        if runtime.args.Cl:
            partner = Molecule.create_Cl()
        elif runtime.args.NO3:
            partner = Molecule.create_NO3()
        else:
            partner = Molecule.create_OH()
    else:
        return False

    partner.program = runtime.QC_program
    QC_input(partner, constrain=False, TS=False)
    submit_and_monitor(partner, logger, threads)
    return True


def crest_collection_path(molecule):
    name = molecule.name if molecule.name.endswith('_CREST') else f"{molecule.name}_CREST"
    return os.path.join(molecule.directory, f"collection{name}.pkl")


def _process_crest_output(molecules, logger, threads):
    all_conformers = []
    constrained_indexes = molecules[0].constrained_indexes
    mult = molecules[0].mult
    charge = molecules[0].charge
    current_step = molecules[0].current_step
    dir = molecules[0].directory
    method = molecules[0].method

    for molecule in molecules:
        try:
            conformers = initiate_conformers(crest_collection_path(molecule))
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


def check_crest(molecules, logger, threads, interval, max_attempts):
    initial_delay = runtime.args.initial_delay if runtime.args.initial_delay else int(interval * 2)
    interval = runtime.args.interval if runtime.args.interval else int(interval)
    attempts = 0
    sleeping = False
    last_pending = None
    crest_missing_polls = 0
    expected_files = {os.path.basename(crest_collection_path(molecule)) for molecule in molecules}

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
            return _process_crest_output(molecules, logger, threads)

        # All jobs gone from the queue but the CREST output never appeared:
        # the node likely died. Resubmit the affected jobs after a grace period.
        if all(m.status == 'completed or not found' for m in molecules):
            crest_missing_polls += 1
            if crest_missing_polls >= VANISHED_GRACE_POLLS:
                crest_missing_polls = 0
                for molecule in molecules:
                    if os.path.basename(crest_collection_path(molecule)) in files_in_directory:
                        continue
                    molecule.node_failure_count = getattr(molecule, 'node_failure_count', 0) + 1
                    if molecule.node_failure_count > MAX_NODE_FAILURES:
                        logger.error(f"{molecule.name}: CREST job vanished {MAX_NODE_FAILURES} times without producing results. Giving up.")
                        return False
                    xyz_path = os.path.join(molecule.directory, f"{molecule.name}.xyz")
                    if not os.path.exists(xyz_path):
                        molecule.write_xyz_file(xyz_path)
                    job_id, _ = submit_job(molecule, runtime.args)
                    molecule.job_id = f"{job_id}"
                    logger.warning(f"{molecule.name}: CREST job disappeared from the queue without producing results (node failure?). Resubmitted as job {job_id}.")
                checkpoint.save_checkpoint(molecules)
        else:
            crest_missing_polls = 0
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
        # Persist job id + current step before monitoring starts, so a crash
        # from here on can reattach to the submitted job instead of losing it.
        checkpoint.save_checkpoint(molecules)
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
            conf.job_id = ""  # no job submitted yet for this step
            handler(conf)
    if conformer_molecules:
        # Record the step advance (and any conformer expansion/filtering) so a
        # crash before submission resumes at this step with a clean resubmit.
        checkpoint.save_checkpoint(conformer_molecules)
        submit_and_monitor(conformer_molecules, logger, threads)


def handle_input_molecules(molecules, logger, threads):
    groups = {}
    for m in molecules:
        # current_step may be a Step enum member or a plain string depending on
        # how the molecule was created; normalize to the step name.
        step = m.current_step.value if isinstance(m.current_step, Step) else str(m.current_step)
        groups.setdefault(step, []).append(m)
    if len(groups) > 1:
        logger.event(f"Molecules span {len(groups)} workflow steps ({', '.join(groups)}); reconciling each group separately")
    for step, group in groups.items():
        reconcile_group(group, step, logger, threads, all_molecules=molecules)


def reconcile_group(group, step, logger, threads, all_molecules=None):
    known = list(all_molecules) if all_molecules else list(group)
    logger.event(f"Reconciling {len(group)} molecule(s) at step {step}")

    if step == 'Done':
        with runtime.global_molecules_lock:
            for m in group:
                if m not in runtime.global_molecules:
                    runtime.global_molecules.append(m)
        logger.success(f"{len(group)} molecule(s) have already finished the workflow")
        return

    if step == 'crest_sampling':
        _reconcile_crest_group(group, logger, threads)
        return

    # ----- QC steps: classify each molecule as converged / running / lost -----
    update_molecules_status(group)

    resubmit_list = []
    for m in group:
        if m.converged:
            continue
        # Logs of finished jobs may already have been moved to log_files/
        candidates = [os.path.join(m.directory, f"{m.name}{m.output}"),
                      os.path.join(m.directory, "log_files", f"{m.name}{m.output}")]
        m.log_file_path = next((p for p in candidates if os.path.exists(p)), candidates[0])
        if os.path.exists(m.log_file_path):
            converge_status, _ = termination_status(m, logger, quick=True)
            if converge_status is True:
                m.converged = True
                continue
        if m.status in ('running', 'pending', 'unknown'):
            continue  # still in the queue (or squeue unreachable): reattach below
        resubmit_list.append(m)

    converged = [m for m in group if m.converged]
    running = [m for m in group if not m.converged and m not in resubmit_list]

    if len(converged) == len(group):
        if step == 'DLPNO':
            with runtime.global_molecules_lock:
                for m in group:
                    if m not in runtime.global_molecules:
                        runtime.global_molecules.append(m)
                known_all = list(runtime.global_molecules)
            known_all += [m for m in known if m not in known_all]
            checkpoint.save_checkpoint(group)
            logger.success(f"All {len(group)} DLPNO calculation(s) are converged")
            if submit_missing_partner(group, logger, threads, known=known_all):
                logger.event("Partner molecule needed for reaction energies was missing; submitted it")
        else:
            logger.success(f"All molecules converged for step {step}. Proceeding with next step: {group[0].next_step}")
            group[0].move_files()
            handle_termination(group, logger, threads, converged=True)
        return

    if resubmit_list and not running and not converged:
        # Nothing queued and nothing finished: cleanly re-run the whole step
        # (regenerates the input files and submits one job/array).
        logger.event(f"No jobs for step {step} are queued and none finished; re-running the step for {len(group)} molecule(s)")
        handle_termination(group, logger, threads, converged=False)
        return

    for m in resubmit_list:
        logger.warning(f"{m.name}: job is no longer queued and its log shows no completion; resubmitting")
        resubmit_job(m, logger)
    if running:
        logger.event(f"Reattaching to {len(running)} queued/running job(s) for step {step} (job id {running[0].job_id})")

    checkpoint.save_checkpoint(group)
    thread = Thread(target=check_convergence, args=(group, logger, threads, get_interval_seconds(group[0]), runtime.args.attempts))
    threads.append(thread)
    thread.start()


def _reconcile_crest_group(group, logger, threads):
    # 1) CREST already produced its conformer collections → process them
    if all(os.path.exists(crest_collection_path(m)) for m in group):
        logger.event("CREST output found on disk; processing conformers")
        _process_crest_output(group, logger, threads)
        return
    # 2) CREST job(s) still queued/running → reattach the CREST monitor
    if any(m.job_id for m in group):
        update_molecules_status(group)
        if any(m.status in ('running', 'pending', 'unknown') for m in group):
            logger.event(f"Reattaching to CREST job (job id {group[0].job_id})")
            thread = Thread(target=check_crest, args=(group, logger, threads, get_interval_seconds(group[0]), runtime.args.attempts))
            threads.append(thread)
            thread.start()
            return
    # 3) Nothing queued and no output → (re)submit CREST sampling
    logger.event("No CREST job queued and no CREST output found; submitting CREST sampling")
    for m in group:
        m.converged = False
    handle_termination(group, logger, threads, converged=False)
