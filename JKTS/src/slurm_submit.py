import os
import subprocess
import re
import time
from .config import config


def submit_array_job(molecules, args=None, nnodes=1):
    args = args or config.args
    dir = molecules[0].directory
    job_files = [f"{conf.name}{conf.input}" for conf in molecules if conf.converged is False]
    job_name = f"{molecules[0].current_step.split('_')[0]}_{os.path.basename(dir)}"

    job_program = molecules[0].program
    submit_name = f"{molecules[0].current_step}_submit.sh"
    path_submit_script = os.path.join(dir, submit_name)
    array_txt = "array.txt"
    array_path = os.path.join(dir, array_txt)

    ncpus = args.cpu
    mem = args.mem
    partition = args.par
    if partition == 'qtest':
        interval = min(get_interval_seconds(molecules[0]), 3600 / args.attempts)
        slurm_time = '1:00:00'
    else:
        if 'aug' in args.basis_set:
            interval = 8640
        else:
            interval = get_interval_seconds(molecules[0])
            if args.F12:
                interval *= 5
        total_seconds = interval * (args.attempts + 3)
        slurm_time = seconds_to_hours(total_seconds)

    if job_program.lower() in ("orca", "g16"):
        if molecules[0].current_step == 'DLPNO':
            program_mem = round((mem + 12000) * 1.25)
        else:
            program_mem = round(mem * 1.25)
    else:
        program_mem = mem

    with open(array_path, 'w') as f:
        for jf in job_files:
            f.write(f"{jf}\n")

    with open(path_submit_script, 'w') as file:
        file.write("#!/bin/bash\n\n")
        file.write("submit=sbatch\n\n")
        file.write("IN=$1\n")
        file.write("[ `cut -c1 <<< $IN` == '-' ] && { submit=cat; IN=`cut -c2- <<< $IN`; }\n\n")

        file.write('[ -f "${IN:-##}" ] || { echo "File not found. Good Bye"; exit 1; }\n')
        file.write("PAR=${2:-" + str(len(job_files)) + "}\n")
        file.write("egrep -q '^X[0-9]*$' <<< \"X${PAR}\" || { echo 'Illegal number, $PAR. Good Bye'; exit 1; }\n")

        file.write("MAX=$(wc -l < $IN)\n")
        file.write("[ $PAR -gt $MAX ] && PAR=$MAX\n")
        file.write('ARRAY="1-${MAX}%${PAR}"\n\n')

        file.write('JOB=${IN%.*}\n\n')
        file.write('SUBMIT=qsub.tmp\n')
        file.write("PWD=`pwd`\n\n")

        file.write("cat > $SUBMIT <<!EOF\n")
        file.write("#!/bin/sh\n")
        file.write(f"#SBATCH --job-name={job_name}\n")
        if job_program.lower() == 'g16':
            file.write(f"#SBATCH --nodes={nnodes}\n")
            file.write(f"#SBATCH --cpus-per-task={ncpus}\n")
            file.write(f"#SBATCH --ntasks={nnodes}\n")
        else:
            file.write(f"#SBATCH --nodes={nnodes}\n")
            file.write("#SBATCH --cpus-per-task=1\n")
            file.write(f"#SBATCH --ntasks={ncpus}\n")
        file.write(f"#SBATCH --output=./slurm_output/{job_name}_%A_%a.out\n")
        file.write(f"#SBATCH --time={slurm_time}\n")
        file.write(f"#SBATCH --partition={partition}\n")
        file.write("#SBATCH --no-requeue\n")
        file.write(f"#SBATCH --mem={program_mem}\n")
        file.write("#SBATCH --array=$ARRAY\n\n")
        if getattr(args, 'account', None):
            file.write(f"#SBATCH --account={args.account}\n")

        file.write("source ~/.JKCSusersetup.txt\n")

        if job_program.lower() == 'g16':
            file.write("eval \\${MODULE_G16}\n")
            file.write("export GAUSS_EXEDIR=\\${PATH_G16}g16/\n")
            file.write("# create and set scratch dir\n")
            file.write("mkdir -p \\${WRKDIR} || exit $?\n")
            file.write("cd $PWD\n")
            file.write("export GAUSS_SCRDIR=\\${WRKDIR}\n\n")
            file.write('GJ=\\$(awk "NR == \\${SLURM_ARRAY_TASK_ID}" $IN)\n')
            file.write('LOG=\\${GJ%.*}.log\n\n')
            file.write('\\${PATH_G16}g16/g16 \\${GJ} > \\${LOG}\n')
            file.write("!EOF\n\n")
        elif job_program.lower() == 'orca':
            file.write("ulimit -c 0\n")
            file.write("export LD_LIBRARY_PATH=\\${PATH_ORCA}:\\${LD_LIBRARY_PATH}\n")
            file.write("eval \\${MODULE_ORCA}\n")
            file.write("# create and set scratch dir\n")
            file.write("mkdir -p \\${WRKDIR} || exit $?\n")
            file.write("cd \\${WRKDIR}\n")
            file.write(f"cp \\${{SLURM_SUBMIT_DIR}}/{array_txt} .\n")
            file.write("cp \\${SLURM_SUBMIT_DIR}/*.inp .\n\n")
            file.write('GJ=\\$(awk "NR == \\${SLURM_ARRAY_TASK_ID}" $IN)\n')
            file.write('LOG=\\${GJ%.*}.out\n\n')
            file.write('\\${PATH_ORCA}/orca \\${GJ} > \\${SLURM_SUBMIT_DIR}/\\${LOG}\n\n')
            if getattr(args, 'IRC', False):
                file.write('cp *xyz \\${SLURM_SUBMIT_DIR}/.\n')
            file.write("!EOF\n\n")
        else:
            file.write('MOLECULE_NAME=\\$(awk "NR == \\${SLURM_ARRAY_TASK_ID}" "$IN")\n')
            file.write('cp \\${SLURM_SUBMIT_DIR}/\\${MOLECULE_NAME} \\${WRKDIR}/\n')
            file.write('cd \\${WRKDIR}\n')
            file.write(f"program_CREST \\${{MOLECULE_NAME}} -gfn{args.gfn} -ewin {args.ewin} --tstep 1 -noreftopo\n")
            file.write('cp *pkl \\${SLURM_SUBMIT_DIR}/.\n')
            file.write('cp *output \\${SLURM_SUBMIT_DIR}/.\n')
            file.write("!EOF\n\n")

        file.write("mkdir -p slurm_output\n")
        file.write("sbatch $SUBMIT")

    submit_command = ['sh', path_submit_script, array_txt]
    try:
        result = _run_submit_with_retries(submit_command, cwd=dir)
        job_id = extract_job_id(result.stdout)
        for n, molecule in enumerate(molecules, start=1):
            molecule.job_id = f"{job_id}_{n}"
        return job_id, interval
    except subprocess.CalledProcessError as e:
        print(f"Error in job submission: {e}!! Check g16/orca_submit.sh")
        if e.stdout:
            print("STDOUT:\n" + e.stdout)
        if e.stderr:
            print("STDERR:\n" + e.stderr)
        exit()


def submit_job(molecule, args=None, nnodes=1):
    args = args or config.args
    input_file_name = f"{molecule.name}{molecule.input}"
    dir = molecule.directory    
    submit_name = f"{molecule.program}_submit.sh"
    path_submit_script = os.path.join(dir, submit_name)
    ncpus = args.cpu
    mem = args.mem
    partition = args.par
    if partition == 'qtest':
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

    if molecule.reactant:
        crest_input = f"{molecule.name}.xyz --gfn{args.gfn} --ewin {args.ewin} --noreftopo --tstep 1"
    elif molecule.product:
        crest_input = f"{molecule.name}.xyz --gfn{args.gfn} --ewin {args.ewin} --noreftopo --uhf 1 --tstep 1"
    else:
        crest_input = f"{molecule.name}.xyz --gfn{args.gfn} --ewin {args.ewin} --noreftopo --uhf 1 --tstep 1 --cinp {dir}/constrain.inp"

    if molecule.program.lower() == "orca" or molecule.program.lower() == "g16":
        if molecule.current_step == 'DLPNO':
            if molecule.error_termination_count == 1:
                program_mem = round((mem+16000) * 1.25)
            elif molecule.error_termination_count == 2:
                program_mem = round((mem+20000) * 1.25)
            else:
                program_mem = round((mem+12000) * 1.25)
        else:
            program_mem = round(mem * 1.25) # Headspace for program
    else:
        program_mem = mem

    if args.account:
        account_string = f'#SBATCH --account={args.account}'
    else:
        account_string = ''


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
#SBATCH --output=./slurm_output/${{JOB}}_%j.output
#SBATCH --time={slurm_time}
#SBATCH --partition={partition}
#SBATCH --no-requeue
#SBATCH --mem={program_mem}
{account_string}

source ~/.JKCSusersetup.txt

mkdir -p \\$WRKDIR || exit $?
cd \\$WRKDIR
cp \\$SLURM_SUBMIT_DIR/{molecule.name}.xyz .
program_CREST {crest_input}
cp *pkl \\$SLURM_SUBMIT_DIR/.
cp *output \\$SLURM_SUBMIT_DIR/.
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
#SBATCH --output=./slurm_output/${{JOB}}_%j.output
#SBATCH --time={slurm_time}
#SBATCH --partition={partition}
#SBATCH --no-requeue
#SBATCH --mem={program_mem}
{account_string}

source ~/.JKCSusersetup.txt
eval \\$MODULE_G16
export GAUSS_EXEDIR=\\${{PATH_G16}}/g16/

# Create scratch folder
mkdir -p \\$WRKDIR || exit $?
cd \\$PWD
export GAUSS_SCRDIR=\\$WRKDIR

\\${{PATH_G16}}/g16/g16 $JOB.com > \\$SLURM_SUBMIT_DIR/$JOB.log

!EOF

mkdir -p slurm_output
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
#SBATCH --output=./slurm_output/${{JOB}}_%j.output
#SBATCH --time={slurm_time}
#SBATCH --partition={partition}
#SBATCH --no-requeue
#SBATCH --mem={program_mem}
{account_string}

source ~/.JKCSusersetup.txt
export LD_LIBRARY_PATH=\\$PATH_ORCA:\\$LD_LIBRARY_PATH
eval \\$MODULE_ORCA

# create and set scratch dir
mkdir -p \\$WRKDIR || exit $?

cd \\$WRKDIR
cp \\$SLURM_SUBMIT_DIR/$JOB.inp .
\\${{PATH_ORCA}}/orca $JOB.inp > \\$SLURM_SUBMIT_DIR/$JOB.out
cp *Full_trj.xyz \\$SLURM_SUBMIT_DIR/.

!EOF

mkdir -p slurm_output
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
        result = _run_submit_with_retries(submit_command, cwd=dir)
        job_id = extract_job_id(result.stdout)
        molecule.job_id = f"{job_id}"
        return job_id, interval
    except subprocess.CalledProcessError as e:
        print(f"Error in job submission: {e}!! Check g16/orca_submit.sh")
        if e.stdout:
            print("STDOUT:\n" + e.stdout)
        if e.stderr:
            print("STDERR:\n" + e.stderr)
        exit()


def _run_submit_with_retries(cmd, cwd=None, retries=10, backoff_seconds=30):
    """Run sbatch wrapper script with retries when QOS submit limit is hit."""
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=cwd)
        except subprocess.CalledProcessError as e:
            stdout = e.stdout or ""
            stderr = e.stderr or ""
            qos_limit = (
                "QOSMaxSubmitJobPerUserLimit" in stderr
                or "Job violates accounting/QOS policy" in stderr
                or "QOSMaxSubmitJobPerUserLimit" in stdout
                or "Job violates accounting/QOS policy" in stdout
            )
            if qos_limit and attempt < retries:
                wait = backoff_seconds * attempt
                print(f"sbatch submit limit reached (attempt {attempt}/{retries}). Waiting {wait}s before retry…")
                time.sleep(wait)
                last_err = e
                continue
            # Not a QOS limit error or retries exhausted: re-raise with context
            if stdout:
                print("STDOUT:\n" + stdout)
            if stderr:
                print("STDERR:\n" + stderr)
            raise
    # If we somehow exit loop without returning
    if last_err:
        raise last_err
    raise RuntimeError("Unknown submission failure without exception")


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
    
    a = 5.90
    b = -19.79
    c = 190.33
    d = -90.7

    if 'H2O' in molecule.name or 'OH' in molecule.name:
        return 20
    elif 'conf' in job_type:
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    elif job_type == 'crest_sampling':
        a = 27
        b = 15
        interval = (a*heavy_count) + b
    elif 'DLPNO' in job_type:
        factor = 0.66
        a *= factor
        b *= factor
        c *= factor
        d *= factor
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    elif 'TS' in job_type:
        factor = 0.80
        a *= factor
        b *= factor
        c *= factor
        d *= factor
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d
    else: # Geometry optimization
        factor = 0.6
        a *= factor
        b *= factor
        c *= factor
        d *= factor
        interval = (a*heavy_count**3) + (b*heavy_count**2) + (c*heavy_count) + d

    return max(interval, 30)


def seconds_to_hours(seconds):
    seconds = round(seconds)
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    return f"{hours:02d}:{minutes:02d}:00"
