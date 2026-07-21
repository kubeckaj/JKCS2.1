import sys
import argparse
import re
import os
import runtime
from output import logger
from classes import Molecule
from conformer_tools import initiate_conformers


def read_input():
    molecules = []
    max_conformers = runtime.args.max_conformers

    def f(molecule):
        if runtime.args.info:
            molecule.print_items()
        else:
            molecule.directory = runtime.start_dir
            molecules.append(molecule)

    def extract_conf_number(molecule):
        match = re.search(r'conf(\d+)', molecule.name)
        return int(match.group(1)) if match else -1

    for input_file in runtime.args.input_files:
        file_name, file_type = os.path.splitext(input_file)
        if file_type == '.pkl':
            from pandas import DataFrame
            pandas_molecules = Molecule.load_molecules_from_pickle(input_file)
            if isinstance(pandas_molecules, DataFrame):
                conformers = initiate_conformers(input_file)
                for conf in conformers:
                    f(conf)
            else:
                if isinstance(pandas_molecules, list):
                    for molecule in pandas_molecules:
                        f(molecule)
                else:
                    f(pandas_molecules)
        elif file_type in ['.log', '.out', '.com', '.inp', '.xyz']:
            QC_program = 'ORCA' if file_type in ['.out', '.inp'] else 'G16'
            input_file_path = os.path.join(runtime.start_dir, input_file)
            molecule = Molecule(input_file_path, indexes=runtime.args.CHO, program=QC_program, method=runtime.args.method)
            molecule.name = file_name
            f(molecule)
        else:
            logger.error(f"Unsupported input file '{input_file}'. Without a specified reaction type, JKTS expects '.pkl', '.log', '.out', '.com', '.inp', or '.xyz' files.")
            sys.exit(1)

    molecules.sort(key=extract_conf_number)
    return molecules[:max_conformers]


def slurm_time(value):
    if not re.fullmatch(r'\d+|\d+:\d{1,2}|\d+:\d{1,2}:\d{1,2}|\d+-\d{1,2}(:\d{1,2}(:\d{1,2})?)?', value):
        raise argparse.ArgumentTypeError(f"'{value}' is not a valid SLURM time (e.g. 72:00:00, 3-00:00:00)")
    return value


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    raise argparse.ArgumentTypeError('Boolean value expected.')