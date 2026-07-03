import os
import subprocess
import re
import shlex
import math

_JKSEND_PATH = None

def _find_jksend():
    global _JKSEND_PATH
    if _JKSEND_PATH:
        return _JKSEND_PATH
    start = os.path.dirname(os.path.realpath(__file__))
    d = start
    while True:
        candidate = os.path.join(d, 'TOOLS', 'SCRIPTS', 'JKsend')
        if os.path.isfile(candidate):
            _JKSEND_PATH = candidate
            return candidate
        parent = os.path.dirname(d)
        if parent == d:
            raise FileNotFoundError(f"Cannot locate TOOLS/SCRIPTS/JKsend above {start}; is JKTS installed inside the JKCS2.1 tree?")
        d = parent


def _get_sbatch_prefix(job_name, partition, slurm_time, nnodes, mem_per_cpu_mb, ncpus=None):
    command = (f"source ~/.JKCSusersetup.txt; "
               f"program_SBATCH '{job_name}' {partition} {slurm_time} {nnodes} {mem_per_cpu_mb}mb")
    if ncpus is not None:
        command += f" {ncpus}"
    result = subprocess.run(['bash', '-c', command], capture_output=True, text=True)
    lines = [line for line in result.stdout.strip().splitlines() if line.strip()]
    prefix = lines[-1] if lines else ''
    if not prefix.startswith('sbatch'):
        raise RuntimeError(f"program_SBATCH from ~/.JKCSusersetup.txt did not produce an sbatch command: stdout={result.stdout!r} stderr={result.stderr!r}")
    return shlex.split(prefix)[1:]


def _compute_time_and_interval(molecule, args):
    if args.par == 'qtest':
        interval = min(get_interval_seconds(molecule), 3600/args.attempts)
        slurm_time = '1:00:00'
    else:
        if 'aug' in args.basis_set:
            interval = 8640
        else:
            interval = get_interval_seconds(molecule)
            if args.F12:
                interval *= 5
        total_seconds = interval * (args.attempts + 3)
        slurm_time = seconds_to_hours(total_seconds)
    return slurm_time, interval


def _compute_program_mem(program, current_step, mem, error_count):
    if program in ('orca', 'g16'):
        if current_step == 'DLPNO':
            if error_count == 1:
                return round((mem+16000) * 1.25)
            elif error_count >= 2:
                return round((mem+20000) * 1.25)
            else:
                return round((mem+12000) * 1.25)
        return round(mem * 1.25) # Headspace for program
    return mem


def _build_submit_command(job_name, program, args, nnodes, program_mem, slurm_time, payload, output_pattern, array_range=None):
    ncpus = args.cpu
    mem_per_cpu = max(1, math.ceil(program_mem / ncpus))
    if program == 'g16':
        prefix = _get_sbatch_prefix(job_name, args.par, slurm_time, nnodes, mem_per_cpu)
        cpu_flags = [f'--cpus-per-task={ncpus}', '--ntasks-per-node=1', '-n', str(nnodes)]
    else:
        prefix = _get_sbatch_prefix(job_name, args.par, slurm_time, nnodes, mem_per_cpu, ncpus)
        cpu_flags = []
    extras = [f'--output=./slurm_output/{job_name}_{output_pattern}',
              f'--error=./slurm_output/{job_name}_{output_pattern}',
              '--no-requeue']
    if array_range:
        extras.append(f'--array={array_range}')
    if args.account:
        extras.append(f'--account={args.account}')
    return ['sbatch'] + prefix + cpu_flags + extras + [_find_jksend(), payload]


def _submit(cmd_argv, cwd):
    if os.environ.get('JKTS_DRYRUN'):
        print(f"JKTS_DRYRUN (cwd={cwd}): {shlex.join(cmd_argv)}")
        return 0
    result = subprocess.run(cmd_argv, capture_output=True, text=True, check=True, cwd=cwd)
    return extract_job_id(result.stdout)


def submit_array_job(molecules, args, nnodes=1):
    dir = molecules[0].directory
    job_files = [f"{conf.name}{conf.input}" for conf in molecules if conf.converged is False]
    job_name = f"{molecules[0].current_step.split('_')[0]}_{os.path.basename(dir)}"

    program = molecules[0].program.lower()
    slurm_time, interval = _compute_time_and_interval(molecules[0], args)
    max_error_count = max((m.error_termination_count for m in molecules), default=0)
    program_mem = _compute_program_mem(program, molecules[0].current_step, args.mem, max_error_count)

    with open(os.path.join(dir, 'array.txt'), 'w') as f:
        for file in job_files:
            f.write(f"{file}\n")

    d = shlex.quote(dir)
    pick_input = 'IN=$(awk "NR == $SLURM_ARRAY_TASK_ID" array.txt)'
    if program == 'g16':
        payload = f'source ~/.JKCSusersetup.txt; cd {d}; {pick_input}; program_G16 $IN'
    elif program == 'orca':
        payload = f'source ~/.JKCSusersetup.txt; cd {d}; {pick_input}; program_ORCA $IN'
    else:
        payload = f'source ~/.JKCSusersetup.txt; cd {d}; {pick_input}; program_CREST $IN -gfn{args.gfn} -ewin {args.ewin} --tstep 1 -noreftopo'

    n_jobs = len(job_files)
    try:
        submit_command = _build_submit_command(job_name, program, args, nnodes, program_mem, slurm_time, payload, "%A_%a.out", array_range=f"1-{n_jobs}%{n_jobs}")
        job_id = _submit(submit_command, dir)
        for n, molecule in enumerate(molecules, start=1):
            molecule.job_id = f"{job_id}_{n}"
        return job_id, interval
    except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError) as e:
        stderr = getattr(e, 'stderr', str(e))
        print(f"Error submitting array job (exit {getattr(e, 'returncode', 1)}): {stderr}\nCommand: {getattr(e, 'cmd', '')}")
        exit()


def submit_job(molecule, args, nnodes=1):
    dir = molecule.directory
    program = molecule.program.lower()
    slurm_time, interval = _compute_time_and_interval(molecule, args)
    program_mem = _compute_program_mem(program, molecule.current_step, args.mem, molecule.error_termination_count)

    d = shlex.quote(dir)
    if program == 'crest':
        job_name = f"{molecule.name}_crest"
        if molecule.reactant:
            crest_input = f"{molecule.name}.xyz --gfn{args.gfn} --ewin {args.ewin} --noreftopo --tstep 1"
        elif molecule.product:
            crest_input = f"{molecule.name}.xyz --gfn{args.gfn} --ewin {args.ewin} --noreftopo --uhf 1 --tstep 1"
        else:
            crest_input = f"{molecule.name}.xyz --gfn{args.gfn} --ewin {args.ewin} --noreftopo --uhf 1 --tstep 1 --cinp {dir}/constrain.inp"
        payload = f"source ~/.JKCSusersetup.txt; cd {d}; program_CREST {crest_input}"
    elif program == 'g16':
        job_name = molecule.name
        payload = f"source ~/.JKCSusersetup.txt; cd {d}; program_G16 {molecule.name}{molecule.input}"
    elif program == 'orca':
        job_name = molecule.name
        payload = f"source ~/.JKCSusersetup.txt; cd {d}; program_ORCA {molecule.name}{molecule.input}"
    else:
        print(f"Error: program '{molecule.program}' is not supported for job submission")
        exit()

    try:
        submit_command = _build_submit_command(job_name, program, args, nnodes, program_mem, slurm_time, payload, "%j.output")
        job_id = _submit(submit_command, dir)
        molecule.job_id = f"{job_id}"
        return job_id, interval
    except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError) as e:
        stderr = getattr(e, 'stderr', str(e))
        print(f"Error submitting job (exit {getattr(e, 'returncode', 1)}): {stderr}\nCommand: {getattr(e, 'cmd', '')}")
        exit()


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
                range_match = re.search(r'\[(\d+)-(\d+)\]', job_id_field)
                if range_match:
                    start, end = range_match.groups()
                    for i in range(int(start), int(end) + 1):
                        job_statuses[f"{main_job_id}_{i}"] = status
            else:
                job_statuses[job_id_field] = status

    return job_statuses


def extract_job_id(output):
    match = re.search(r'Submitted batch job (\d+)', output)
    if match:
        return int(match.group(1))
    else:
        print("Could not extract job ID.")
        return None


def get_interval_seconds(molecule):
    job_type = molecule.current_step
    atoms = molecule.atoms
    if molecule.reactant or molecule.product:
        heavy_count = 0
    else:
        heavy_count = -1
    H_count = 0
    for atom in atoms:
        if atom == 'H':
            H_count += 1
        else:
            heavy_count += 1

    a = 5.90; b = -19.79; c = 190.33; d = -90.7

    if 'H2O' in molecule.name or 'OH' in molecule.name:
        return 20
    elif 'conf' in job_type:
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    elif job_type == 'crest_sampling':
        a = 27; b = 15
        interval = (a*heavy_count) + b
    elif 'DLPNO' in job_type:
        factor = 0.66
        a*=factor; b*=factor; c*=factor; d*=factor
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    elif 'TS' in job_type:
        factor = 0.80
        a*=factor; b*=factor; c*=factor; d*=factor
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    else: # Geometry optimization
        factor = 0.6
        a*=factor; b*=factor; c*=factor; d*=factor
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d

    return max(interval, 30)


def seconds_to_hours(seconds):
    seconds = round(seconds)
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    return f"{hours:02d}:{minutes:02d}:00"
