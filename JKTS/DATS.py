#!/usr/bin/env python3
import os
import shlex
import re

from src.molecule import Molecule, Logger
from src.molecule_manager import global_molecules
from src.cli import get_parser
from src.utils import mkdir, collect_DFT_and_DLPNO
from src.monitor import handle_termination
from src.input import read_input, initiate_conformers, handle_input_molecules
from src.kinetics import rate_constant
from src import config


def main():
    parser = get_parser()

    global start_dir
    # Parse args and store them in a shared config object so other modules
    # can import `from src.config import config` and access `config.args`.
    config.args = parser.parse_args()
    # Local convenience names for backward compatibility inside this module
    args = config.args
    # Normalize incorrectly packed arguments passed as a single positional string (e.g., "-smiles C -OH -par qtest")
    if args.input_files and len(args.input_files) == 1:
        single = args.input_files[0]
        if isinstance(single, str) and single.startswith('-') and ' ' in single:
            reparsed = parser.parse_args(shlex.split(single))
            config.args = reparsed
            args = reparsed
            # Ensure input_files is a list (argparse will set [] if none provided)
            if not hasattr(args, 'input_files') or args.input_files is None:
                args.input_files = []
    start_dir = os.getcwd()

    # Determine QC program and store it
    if args.ORCA:
        args.G16 = False
        config.QC_program = "ORCA"
    else:
        args.G16 = True
        config.QC_program = "G16"

    # The chosen QC program is stored in config.QC_program

    ###########################################TEST######################################################
    if args.test:
        molecules = []
        print("TEST")
        
        exit()

    threads = []
    input_molecules = []
    final_TS = []
    final_reactants = []
    final_products = []
    molecules = []

    if args.init: # only executes if input file is an xyz file
        if args.smiles:
            input_molecule = Molecule(smiles=args.smiles, reactant=True, method=args.method)
            file_name, file_type = input_molecule.name, 'xyz'
        else:
            input_file = args.input_files[0]
            file_name, file_type = os.path.splitext(input_file)
            input_file_path = os.path.join(start_dir, input_file)
            input_molecule = Molecule(input_file_path, reactant=True, method=args.method) # Reuse input molecule as template for reactant

        if args.OH or args.Cl:
            reacted_molecules, product_molecules = input_molecule.H_abstraction(Cl=args.Cl, products=args.products, num_molecules=args.num_molecules)
        elif args.CC:
            parser.error("Reaction type not supported yet.")
            other_molecule = args.input_files[1]
            reacted_molecules = input_molecule.addition(other_molecule)
        elif args.OH_CC:
            # parser.error("Reaction type not supported yet.")
            reacted_molecules = input_molecule.OH_addition()
        else:
            parser.error("Need to specify reaction type")

        if args.TS:
            for count, molecule in enumerate(reacted_molecules, start=1):
                molecule.name = f"{file_name}_H{count}" 
                molecule.directory = os.path.join(start_dir, molecule.name)
                mkdir(molecule, args.method)
                logger = Logger(os.path.join(molecule.directory, "log"))
                molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz")) 
                molecule.save_to_pickle(os.path.join(molecule.directory, f"{molecule.name}.pkl"))

        if args.reactants: 
            reactant_dir = os.path.join(start_dir, 'reactants')
            input_molecule.name = f"{file_name}_reactant"
            input_molecule.directory = reactant_dir
            input_molecule.set_current_step('crest_sampling')
            mkdir(input_molecule, args.method, crest_constrain=False)
            input_molecule.save_to_pickle(os.path.join(input_molecule.directory, f"{input_molecule.name}.pkl"))

        if args.products and product_molecules:
            product_dir = os.path.join(start_dir, 'products')
            for count, molecule in enumerate(product_molecules, start=1):
                molecule.name = f"{file_name}_product_H{count}"
                molecule.directory = product_dir
            mkdir(product_molecules[0], args.method, crest_constrain=False)
            pickle_path = os.path.join(product_dir, "product_molecules.pkl")
            Molecule.molecules_to_pickle(product_molecules, pickle_path)

    elif args.info:
        if args.smiles:
            smiles_molecule = Molecule(smiles=args.smiles, reactant=True, method=args.method)
            smiles_molecule.print_items()
        else:
            read_input(args, start_dir)


    else: # We loop over all molecules given in the argument input and process them according to file type
        if args.smiles:
            input_molecule = Molecule(smiles=args.smiles, reactant=True, method=args.method)
            args.input_files = [f"{input_molecule.name}.xyz"]
        for n, input_file in enumerate(args.input_files, start=1):
            file_name, file_type = os.path.splitext(input_file)

            if file_type == '.xyz':
                if os.path.basename(start_dir) == 'reactants':
                    if args.reactants is False:
                        exit()
                    input_file_path = os.path.join(start_dir, f"{file_name}_reactant.pkl")
                    logger = Logger(os.path.join(start_dir, "log"))
                    reactant = Molecule.load_from_pickle(input_file_path)
                    molecules.append(reactant)

                elif os.path.basename(start_dir) == 'products':
                    if args.products is False:
                        exit()
                    logger = Logger(os.path.join(start_dir, "log"))
                    pickle_path = os.path.join(start_dir, "product_molecules.pkl")
                    product_molecules = Molecule.load_molecules_from_pickle(pickle_path)
                    for product in product_molecules:
                        product.product = True
                        product.current_step = 'crest_sampling'
                        molecules.append(product)
                else:
                    input_file_path = os.path.join(start_dir, os.path.basename(start_dir)+'.pkl')
                    molecule = Molecule.load_from_pickle(input_file_path)
                    logger = Logger(os.path.join(start_dir, "log"))
                    molecules.append(molecule)

                handle_termination(molecules, logger, threads, converged=False, args=args)

    ### If input files are pickle or log files we collect them into input_molecules and process the list, by checking each input convergence status
            elif file_type == '.pkl':
                from pandas import read_pickle, DataFrame
                df = read_pickle(input_file)
                if isinstance(df, DataFrame): # If pkl file is from the CREST job it will come as a pandas dataframe, instead of a list
                    conformer_molecules = initiate_conformers(input_file)
                    for m in conformer_molecules:
                        input_molecules.append(m)
                else:
                    molecules = Molecule.load_molecules_from_pickle(input_file)
                    if not isinstance(molecules, list):
                        molecules = [molecules]
                    for m in molecules: 
                        if m.current_step == 'Done' or m.current_step == 'DLPNO':
                            if m.reactant:
                                final_reactants.append(m)
                            elif m.product:
                                final_products.append(m)
                            else:
                                final_TS.append(m)
                        else:
                            input_molecules.append(m)

            elif file_type in [".log", ".out", ".com", ".inp"]: 
                # Initialize molecule from log file
                input_file_path = os.path.join(start_dir, input_file)
                molecule = Molecule(file_path=input_file_path, indexes=args.CHO, method=args.method)
                input_molecules.append(molecule)

            else:
                print(f"File type {file_type} is not supported")
                exit()

        if final_TS and final_reactants: # In the case final_products is empty symmetrical eckart tunneling used
            k, kappa = rate_constant(final_TS, final_reactants, final_products)
            print(f"{k} cm^3 molecules^-1 s^-1")
            print(f"Tunneling coefficient: {kappa}")
            exit()

        # FIX file_type variable
        if input_molecules and file_type != '.xyz':
            logger = Logger(os.path.join(start_dir, "log"))
            if args.collect:
                collected_molecules = collect_DFT_and_DLPNO(input_molecules)
                Molecule.molecules_to_pickle(collected_molecules, os.path.join(start_dir, "collected_molecules.pkl"))
            elif args.restart:
                logger.log("\nJKTS restarted")
                with open(os.path.join(start_dir, '.method'), 'w') as f:
                    f.write(f"{args.method}")
                # handle_input_molecules will add converged molecules to the shared manager
                handle_input_molecules(input_molecules, logger, threads, attempts=args.attempts)
            elif args.rerun:
                logger.log(f"\nJKTS - rerunning calculations of type: {input_molecules[0].current_step}")
                handle_termination(input_molecules, logger, threads, converged=False, args=args)
            elif args.pickle:
                filename = re.sub("_conf\d{1,2}", "",input_molecules[0].name)
                Molecule.molecules_to_pickle(input_molecules, os.path.join(start_dir, f"collection{filename}.pkl"))
            else:
                parser.error("Error")

        if final_TS:
            for m in final_TS:
                global_molecules.append(m)
                    
        elif not input_molecules and file_type != '.xyz':
            logger = Logger(os.path.join(start_dir, "log"))
            logger.log("Error when generating input molecules. Could not create list from given input")
            print("Error when generating input molecules. Could not create list from given input")
            if file_type in [".log", ".out"]:
                logger.log(".log or .out extension detected. Make sure input files are from ORCA or G16")
            elif file_type == '.pkl':
                logger.log("Detected .pkl file. Make sure the structure of the pickle file is either a python list, set, tuple or pandas.DataFrame")


    # Monitor and handle convergence of submitted jobs
    while threads:
        for thread in list(threads):
            thread.join(timeout=0.1)
            if not thread.is_alive():
                threads.remove(thread)

    
    if global_molecules:
        logger = Logger(os.path.join(start_dir, "log"))
        molecules_logger = Logger(os.path.join(start_dir, "molecules.txt"))
        for molecule in global_molecules: 
            molecule.print_items(molecules_logger)
    
        # Use snapshot lists for decision-making and pickling
        gm_list = global_molecules.as_list()
        if all(m.reactant for m in gm_list):
            molecule_name = gm_list[0].name.split("_")[0]
            logger.log(f"Final DLPNO calculations for reactants is done. Logging molecules to Final_reactants_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(gm_list, os.path.join(start_dir, f"Final_reactants_{molecule_name}.pkl"))
        elif all(m.product for m in gm_list):
            # Identify unique H numbers
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in gm_list if "_H" in m.name))
            # Group molecules by H numbers and include H2O in each group
            grouped_lists = [[m for m in gm_list if f"_H{h_num}_" in m.name] + 
                             [m for m in gm_list if "H2O" in m.name] 
                             for h_num in h_numbers]

            # Create a pickle file for each group
            for h, molecules_group in zip(h_numbers, grouped_lists):
                molecule_name = f"{molecules_group[0].name.split('_')[0]}_H{h}"
                pickle_path = os.path.join(start_dir, f"Final_products_{molecule_name}.pkl")
                logger.log(f"Final DLPNO calculations for products are done. Logging properties to Final_products_{molecule_name}.pkl")
                Molecule.molecules_to_pickle(molecules_group, pickle_path)
        else:
            molecule_name = os.path.basename(start_dir)
            logger.log(f"Final DLPNO calculations for transition state molecules is done. Logging properties to Final_TS_{molecule_name}.pkl")
            TS_pkl_path = os.path.join(start_dir, f"Final_TS_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(gm_list, TS_pkl_path)

            if args.k:
                reactant_pkl_name = os.path.basename(start_dir).split("_")[0]
                product_pkl_name = os.path.basename(start_dir)
                reactant_pkl_path = os.path.join(os.path.dirname(start_dir), f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                product_pkl_path = os.path.join(os.path.dirname(start_dir), f'products/Final_products_{product_pkl_name}.pkl')
                logger.log(f"{product_pkl_path}, {reactant_pkl_path}")
                if os.path.exists(reactant_pkl_path):
                    final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                    if os.path.exists(product_pkl_path):
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                    k, kappa = rate_constant(gm_list, final_reactants, final_products)
                    results_logger = Logger(os.path.join(os.path.dirname(start_dir), "Rate_constants.txt"))
                    results_logger.log_with_stars(f"{molecule_name}: {k} cm^3 molecules^-1 s^-1 with tunneling coefficient {kappa}")
                else:
                    reactant_pkl_path = os.path.join(start_dir, f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                    product_pkl_path = os.path.join(start_dir, f'products/Final_products_{product_pkl_name}.pkl')
                    logger.log(f"{product_pkl_path}, {reactant_pkl_path}")
                    if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                        final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                        k, kappa = rate_constant(gm_list, final_reactants, final_products)
                        results_logger = Logger(os.path.join(start_dir, "Rate_constants.txt"))
                        results_logger.log_with_stars(f"1 {molecule_name}: {k} cm^3 molecules^-1 s^-1 with tunneling coefficient {kappa}")
                    else:
                        logger.log(f"Could not find pickle files in path: {product_pkl_name} and {reactant_pkl_path}")


if __name__ == "__main__":
    main()
