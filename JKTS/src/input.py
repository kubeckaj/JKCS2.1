import os
import re
from pandas import read_pickle

from .molecule import Molecule
# Avoid importing monitor at module import time (circular import with monitor).
# Import required monitor functions inside the functions that need them.
from .molecule_manager import global_molecules

def initiate_conformers(input_file=None):
    conformer_molecules = []
    pickle_file_path = input_file

    if not pickle_file_path:
        print("No input file specified.")
        return []

    with open(pickle_file_path, 'rb') as f:
        df = read_pickle(f)

    for index, row in df.iterrows():
        el_energy = row[('log', 'electronic_energy')]
        xyz = row[('xyz', 'structure')]
        file_basename = row[('info', 'file_basename')]  # Accessing the file basename
        
        coords = xyz.get_positions()
        atoms = xyz.get_chemical_symbols()
        coordinates_list = [[atom] + list(coord) for atom, coord in zip(atoms, coords)]

        name = file_basename.replace("-str", "_conf")
        conformer_molecule = Molecule(name=name, directory=os.path.dirname(input_file), electronic_energy=el_energy, atoms=[], coordinates=[], program='CREST')

        if 'reactant' in input_file:
            conformer_molecule.reactant = True
            conformer_molecule.mult = 1
        elif 'product' in input_file:
            conformer_molecule.product = True
        else:
            conformer_molecule.reactant = False
            conformer_molecule.product = False

        conformer_molecule.electronic_energy = el_energy
        conformer_molecule.atoms = [atom[0] for atom in coordinates_list]
        conformer_molecule.coordinates = [atom[1:] for atom in coordinates_list]
        conformer_molecule.program = 'CREST'
        
        # Workflow management
        conformer_molecule.workflow = conformer_molecule.determine_workflow()
        conformer_molecule.set_current_step()
        conformer_molecule.method = conformer_molecule.log2method()
        
        conformer_molecules.append(conformer_molecule)

    return sorted(conformer_molecules, key=lambda m: m.electronic_energy)


def handle_input_molecules(molecules, logger, threads, attempts, args):
    # Local import to avoid circular import (monitor imports this module).
    from .monitor import termination_status, check_convergence, handle_termination
    current_step = molecules[0].current_step
    if all(m.current_step == current_step for m in molecules):
        logger.log(f"Detected molecules with current step: {current_step}")
        logger.log("Checking convergence status of molecules")
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
                else: # in the case a single .pkl file has been given as input, and logfiles are not located in same directory:
                    if len(molecules) == 1:
                        m.converged = False
                    # Check if the conformers in the pkl file have atoms and coordinates and send for calculation.
                    elif m.atoms and m.coordinates:
                        m.converged = True
            all_converged = all(m.converged for m in molecules)

        if all_converged:
            if all(mol.current_step == 'DLPNO' for mol in molecules):
                if all(mol.product for mol in molecules) and not any('H2O' in mol.name for mol in molecules):
                    logger.log("Converged DLPNOs. However, H2O needed for product energies")
                    check_convergence(molecules, logger, threads, 30, attempts, all_converged=True)
                elif all(mol.reactant for mol in molecules) and not any('OH' in mol.name for mol in molecules):
                    logger.log("Converged DLPNOs. However, OH needed for reactant energies")
                    check_convergence(molecules, logger, threads, 30, attempts, all_converged=True)
                else:
                    logger.log("All given molecules are converged DLPNOs")
                    for m in molecules:
                        global_molecules.append(m)
            else:
                logger.log(f"All molecules converged for step {current_step}. Proceeding with next step: {molecules[0].next_step}")
                molecules[0].move_files()
                handle_termination(molecules, logger, threads, converged=True, args=args)
        else:
            logger.log("Not all molecules given have converged. Non-converged will be calculated and converged ones will be skipped.")
            handle_termination(molecules, logger, threads, converged=False, args=args)
    else:
        logger.log("Not all molecules in the given files are of the same type. Check if correct input is provided. Exiting for now.")
        exit()


def read_input(args, start_dir):
    molecules = []
    max_conformers = args.max_conformers if args.max_conformers else 1000

    def f(molecule):
        if args.info:
            molecule.print_items()
        else:
            molecule.directory = start_dir
            molecules.append(molecule)
    
    def extract_conf_number(molecule):
        match = re.search(r'conf(\d+)', molecule.name)
        return int(match.group(1)) if match else -1

    for input_file in args.input_files:
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
            input_file_path = os.path.join(start_dir, input_file)
            molecule = Molecule(input_file_path, indexes=args.CHO, program=QC_program, method=args.method)
            molecule.name = file_name
            f(molecule)
        else:
            print("Invalid input file format. Without a specified reaction type the program expects input files with extensions '.pkl', '.log', '.out', '.com', or '.inp'\nPlease ensure that your input file is in one of these formats and try again. If you provided a different type of file, convert it to a supported format or select an appropriate file for processing.")
            exit()

    molecules.sort(key=extract_conf_number)
   
    methods = {molecule.method for molecule in molecules}

    if len(methods) == 1:
        return molecules[:max_conformers]
    elif len(methods) == 2:
        molecules_level1 = []
        molecules_level2 = []
        for m in molecules:
            if m.method == list(methods)[0]:
                molecules_level1.append(m)
            else:
                molecules_level2.append(m)
        return (molecules_level1[:max_conformers], molecules_level2[:max_conformers])
    elif len(methods) == 3:
        molecules_level1 = []
        molecules_level2 = []
        molecules_level3 = []
        for m in molecules:
            if m.method == list(methods)[0]:
                molecules_level1.append(m)
            elif m.method == list(methods)[1]:
                molecules_level2.append(m)
            else:
                molecules_level3.append(m)
        return (molecules_level1[:max_conformers], molecules_level2[:max_conformers], molecules_level3[:max_conformers])
        