from threading import Lock

args = None
QC_program = None
termination_strings = {}
error_strings = {}
global_molecules = []
global_molecules_lock = Lock()
start_dir = None
DFT_method = None
