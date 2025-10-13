import os 
import re
import time


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


def is_aldehyde(molecule, C_index, H_index):
    _, aldehyde_groups, _ = molecule.identify_functional_groups()
    for aldehyde in aldehyde_groups:
        if aldehyde['C'] == C_index-1 and aldehyde['H'] == H_index-1:
            return True, aldehyde['O']+1
        elif aldehyde['C'] == C_index and aldehyde['H'] == H_index:
            return True, aldehyde['O']
    return False, None

    

def mkdir(molecule, method, crest_constrain=True):
    if not os.path.exists(molecule.directory):
        os.makedirs(molecule.directory, exist_ok=True)
        if not os.path.exists(os.path.join(molecule.directory, "input_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "input_files")))
        if not os.path.exists(os.path.join(molecule.directory, "log_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "log_files")))
        if not os.path.exists(os.path.join(molecule.directory, "slurm_output")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "slurm_output")))
    with open(molecule.directory + "/.method", "w") as f:
        f.write(f"{method}")

    if crest_constrain:
        force_constant = 1
        C_index = molecule.constrained_indexes['C']
        H_index = molecule.constrained_indexes['H']
        O_index = molecule.constrained_indexes['O']
        # with open (molecule.directory + "/constrain.inp","w") as f:
        #     f.write("$constrain\n")
        #     f.write(f"  force constant={force_constant}\n") 
        #     f.write(f"  distance: {C_index}, {H_index}, auto\n")
        #     f.write(f"  distance: {H_index}, {O_index}, auto\n")
        #     f.write("$end\n")
        #
        with open (molecule.directory + "/.constrain","w") as f:
            f.write(f"C: {C_index}\n")
            f.write(f"H: {H_index}\n")
            f.write(f"O: {O_index}\n")

        with open (molecule.directory + "/log_files/.constrain","w") as f:
            f.write(f"C: {C_index}\n")
            f.write(f"H: {H_index}\n")
            f.write(f"O: {O_index}\n")


def read_last_lines(filename, num_lines, interval=200):
    attempts = 0
    max_attempts = 5

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
            time.sleep(interval)

    # logger.log(f"Failed to find the file {filename}")
    return []

    
def collect_DFT_and_DLPNO(molecules):
    collected_molecules = []
    
    for m in molecules:
        process_conditions = (m.current_step == 'optimization' and (m.reactant or m.product)) or ('TS_opt_conf' in m.current_step)
        
        if process_conditions:
            if 'OH' in m.name or 'H2O' in m.name:
                identifier = 'OH' if 'OH' in m.name else 'H2O'
            else:
                match = re.search(r'conf(\d+)', m.name)
                identifier = match.group() if match else None

            if identifier:
                pattern = re.compile(re.escape(identifier) + r'(?!\d)')
                for m_DLPNO in molecules:
                    if pattern.search(m_DLPNO.name) and m_DLPNO.current_step == 'DLPNO':
                        m.electronic_energy = m_DLPNO.electronic_energy
                        m.update_step()
                        m.name = m.name.replace('TS', 'DLPNO')
                        collected_molecules.append(m)
                        break

    return collected_molecules