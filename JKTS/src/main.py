#!/usr/bin/env python3
import argparse
import os
import re

import runtime
from classes import Molecule, Logger
from qc_input import mkdir
from monitoring import handle_termination, handle_input_molecules
from conformer_tools import initiate_conformers, collect_DFT_and_DLPNO
from rate_constant import rate_constant


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def read_input():
    molecules = []
    max_conformers = runtime.args.max_conformers if runtime.args.max_conformers else 1000

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


def main():
    parser = argparse.ArgumentParser(description='''    Dynamic Approach for Transition State-
    Automated tool for generating input files, primarily
    for transition state geometry optimization.
    Calculation of tunneling corrected multi-configurational
    rate constants can also be calculated from log and pickle files.''',
                                     prog="JKTS",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
    Examples of use:
                JKTS pinonaldehyde.xyz -OH
                JKTS -smiles CO -Cl  # Methanol hydrogen abstraction with Cl radical
                JKTS CH4_H1_opt_constrain.pkl -info
                JKTS benzene.xyz -OH -ORCA -par qtest -auto false
                JKTS *TS.log -time 5:00:00
                                     ''')

    parser.add_argument('input_files', metavar='reactant.xyz', nargs='*', help='Input XYZ files (e.g., pinonaldehyde.xyz)')
    parser.add_argument('-G16', action='store_true', default=True, help='Use Gaussian16 for QC calculations')
    parser.add_argument('-ORCA', action='store_true', default=False, help='Use ORCA for QC calculations (default)')
    parser.add_argument('-reactants', type=str2bool, default=True, metavar='<boolean>', help='Prepare folder for reactants [def: True]')
    parser.add_argument('-products', type=str2bool, default=True, metavar='<boolean>', help='Prepare folder for products [def: True]')
    parser.add_argument('-TS', type=str2bool, default=True, metavar='<boolean>', help=argparse.SUPPRESS)
    parser.add_argument('-auto', type=str2bool, default=True, metavar='<boolean>', help='Automated process with the following workflow:\n- CREST conformer sampling of xyz input file (GFN2-xTB -ewin=5 kcal/mol)\n- Preoptimization of geometry (Low Level of Theory)\n- Optimization towards transition state (High Level of Theory)\n- DLPNO-CCSD(T) SP energy calculations on top of TS conformers\n- Calculate rate constants and branching ratios for reaction type\n* Automatically resubmit failed calculations until convergence or wall-time limit is reached')

    # Argparse reactions
    reaction_options = parser.add_argument_group("Types of reactions")
    reaction_options.add_argument('-OH', action='store_true', help='Perform H abstraction with OH radical')
    reaction_options.add_argument('-Cl', action='store_true', help='Perform H abstraction with Cl radical')
    reaction_options.add_argument('-CC', action='store_true', help='(TBA) Perform addition to C=C bonds')
    reaction_options.add_argument('-OH_CC', action='store_true', help='(TBA) Perform OH addition to C=C bonds')

    # Argparse additional
    additional_options = parser.add_argument_group("Additional arguments")
    additional_options.add_argument('-k', type=str2bool, metavar='<boolean>', default=True, help='Calculate Multiconformer Transition State rate constant [def: True]')
    additional_options.add_argument('-restart', action='store_true', default=False, help='Need: .log, .out, .pkl - Restart from molecule current step')
    additional_options.add_argument('-smiles', metavar='string', type=str, help='Input molecule as a SMILES string')
    additional_options.add_argument('-movie', action='store_true', default=False, help='Produce movie.xyz for Molden viewing')
    additional_options.add_argument('-info', action='store_true', default=False, help='Print information of molecules in log files or .pkl file')
    additional_options.add_argument('-CHO', dest='CHO', nargs='*', type=int, help="Set indexes of atoms for active site. Indexing starting from 1")
    additional_options.add_argument('-collect', action='store_true', default=False, help='Collect thermochemical data from TS structures and single point correction from DLPNO')
    additional_options.add_argument('-method', type=str, default='wB97XD', help='Specify QC method to use for optimization and TS search [def: WB97X-D3BJ]')
    additional_options.add_argument('-basis_set', type=str, default='6-31++g(d,p)', help='Specify basis set to use with QC method [def: 6-31++g(d,p)]')
    additional_options.add_argument('-skip_preopt', action='store_true', default=False, help='Skip the preoptimization of the structures before the TS search')
    additional_options.add_argument('-filter_step', dest='filter_step', default=['opt_constrain_conf', 'DLPNO'], nargs='*', help="Steps at which to perform filtering of conformers")
    additional_options.add_argument('-cpu', metavar="int", nargs='?', const=1, type=int, default=4, help='Number of CPUs [def: 4]')
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=8000, help='Amount of memory allocated for the job [def: 8000MB]')
    additional_options.add_argument('-par', metavar="partition", type=str, default="q64", help='Partition to use [def: q64,q48,q28,q24]')
    additional_options.add_argument('-time', metavar="hh:mm:ss", type=str, default=None, help='Monitoring duration [def: 144 hours]')
    additional_options.add_argument('-interval', metavar="int", nargs='?', const=1, type=int, help='Time interval between log file checks [def: based on molecule size]')
    additional_options.add_argument('-initial_delay', metavar="int", nargs='?', const=1, type=int, help='Initial delay before checking log files [def: based on molecule size]')
    additional_options.add_argument('-attempts', metavar="int", nargs='?', const=1, type=int, default=100, help='Number of log file check attempts [def: 100]')
    additional_options.add_argument('-max_conformers', metavar="int", nargs='?', const=1, type=int, default=None, help='Maximum number of conformers from CREST [def: 50]')
    additional_options.add_argument('-freq_cutoff', metavar="int", nargs='?', const=1, type=int, default=-100, help='TS imaginary frequency cutoff [def: -100 cm^-1]')
    additional_options.add_argument('-F12', action='store_true', help='Use CCSD(T)-F12-pVTZ instead of DLPNO')
    additional_options.add_argument('-energy_cutoff', metavar="digit", nargs='?', const=1, default=5, type=float, help='After preoptimization, remove conformers which are [int] kcal/mol higher in energy than the lowest conformer [def: 5 kcal/mol]')
    additional_options.add_argument('-pickle', action='store_true', default=False, help='Store given log files into pickle file')
    additional_options.add_argument('--gfn', default='2', choices=['1', '2'], help='Specify the GFN version for CREST (1 or 2, default: 2)')
    additional_options.add_argument('-ewin', metavar="int", nargs='?', const=1, default=5, type=int, help='Energy threshold for CREST conformer sampling [def: 5 kcal/mol]')
    additional_options.add_argument('--account', type=str, help=argparse.SUPPRESS)
    additional_options.add_argument('-test', action='store_true', default=False, help=argparse.SUPPRESS)
    additional_options.add_argument('-num_molecules', type=int, default=None, help=argparse.SUPPRESS)
    additional_options.add_argument('-rerun', action='store_true', help=argparse.SUPPRESS)
    additional_options.add_argument('-plot', type=str, default=None, help='Generate plot')

    hidden_options = parser.add_argument_group()
    hidden_options.add_argument('-init', action='store_true', help=argparse.SUPPRESS)
    hidden_options.add_argument('-hybrid', action='store_true', default=False, help=argparse.SUPPRESS)
    hidden_options.add_argument('-dispersion', action='store_true', default=False, help=argparse.SUPPRESS)

    runtime.args = parser.parse_args()
    runtime.start_dir = os.getcwd()

    ################################ARGUMENT SPECIFICATIONS############################

    if runtime.args.ORCA:
        runtime.args.G16 = False
        runtime.QC_program = "ORCA"
    else:
        runtime.args.G16 = True
        runtime.QC_program = "G16"

    runtime.termination_strings = {
        "g16": ["normal termination"],
        "orca": ["****orca terminated normally****", "ORCA TERMINATED NORMALLY", "TOTAL RUN TIME", "total run time"],
        "crest": ["crest done", "crest terminated normally."]
    }
    runtime.error_strings = {
        "g16": ["error termination", "another termination example"],
        "orca": ["aborting the run", "this wavefunction is not fully converged", "not enough memory", "orca finished by error termination", "error", "this wavefunction is not converged"],
        "crest": ["some crest error message"]
    }

    #####################################################################################################
    if runtime.args.test:
        from ts_validation import check_transition_state
        threads = []
        logger = Logger(os.path.join(runtime.start_dir, "test.txt"))
        molecules = read_input()
        if molecules:
            for molecule in molecules:
                print(molecule.name, check_transition_state(molecule))
        exit()
    ####################################################################################################

    threads = []
    input_molecules = []
    final_TS = []
    final_reactants = []
    final_products = []
    molecules = []

    if runtime.args.init:
        if runtime.args.smiles:
            input_molecule = Molecule(smiles=runtime.args.smiles, reactant=True, method=runtime.args.method)
            file_name, file_type = input_molecule.name, 'xyz'
        else:
            input_file = runtime.args.input_files[0]
            file_name, file_type = os.path.splitext(input_file)
            input_file_path = os.path.join(runtime.start_dir, input_file)
            input_molecule = Molecule(input_file_path, reactant=True, method=runtime.args.method)

        if runtime.args.OH or runtime.args.Cl:
            reacted_molecules, product_molecules = input_molecule.H_abstraction(Cl=runtime.args.Cl, products=runtime.args.products, num_molecules=runtime.args.num_molecules)
        elif runtime.args.CC:
            parser.error("Reaction type not supported yet.")
            other_molecule = runtime.args.input_files[1]
            reacted_molecules = input_molecule.addition(other_molecule)
        elif runtime.args.OH_CC:
            reacted_molecules = input_molecule.OH_addition()
        else:
            parser.error("Need to specify reaction type")

        if runtime.args.TS:
            for count, molecule in enumerate(reacted_molecules, start=1):
                molecule.name = f"{file_name}_H{count}"
                molecule.directory = os.path.join(runtime.start_dir, molecule.name)
                mkdir(molecule)
                logger = Logger(os.path.join(molecule.directory, "log"))
                molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz"))
                molecule.save_to_pickle(os.path.join(molecule.directory, f"{molecule.name}.pkl"))

        if runtime.args.reactants:
            reactant_dir = os.path.join(runtime.start_dir, 'reactants')
            input_molecule.name = f"{file_name}_reactant"
            input_molecule.directory = reactant_dir
            input_molecule.set_current_step('crest_sampling')
            mkdir(input_molecule, crest_constrain_flag=False)
            input_molecule.save_to_pickle(os.path.join(input_molecule.directory, f"{input_molecule.name}.pkl"))

        if runtime.args.products and product_molecules:
            product_dir = os.path.join(runtime.start_dir, 'products')
            for count, molecule in enumerate(product_molecules, start=1):
                molecule.name = f"{file_name}_product_H{count}"
                molecule.directory = product_dir
            mkdir(product_molecules[0], crest_constrain_flag=False)
            pickle_path = os.path.join(product_dir, "product_molecules.pkl")
            Molecule.molecules_to_pickle(product_molecules, pickle_path)

    elif runtime.args.info:
        if runtime.args.smiles:
            smiles_molecule = Molecule(smiles=runtime.args.smiles, reactant=True, method=runtime.args.method)
            smiles_molecule.print_items()
        else:
            read_input()

    elif runtime.args.movie:
        molecules = read_input()
        if molecules is not None:
            with open(os.path.join(runtime.start_dir, 'movie.xyz'), 'w') as f:
                for m in molecules:
                    f.write(f"{len(m.atoms)}\n")
                    f.write(f"{m.name}\n")
                    for atom, coord in zip(m.atoms, m.coordinates):
                        f.write(f"{atom} {coord[0]} {coord[1]} {coord[2]}\n")
        print("movie.xyz generated")

    elif runtime.args.plot:
        print("Plotting is not available (plotting module has been removed). Use -info to inspect molecule data.")

    else:
        if runtime.args.smiles:
            input_molecule = Molecule(smiles=runtime.args.smiles, reactant=True, method=runtime.args.method)
            runtime.args.input_files = [f"{input_molecule.name}.xyz"]
        for n, input_file in enumerate(runtime.args.input_files, start=1):
            file_name, file_type = os.path.splitext(input_file)

            if file_type == '.xyz':
                if os.path.basename(runtime.start_dir) == 'reactants':
                    if runtime.args.reactants is False:
                        exit()
                    input_file_path = os.path.join(runtime.start_dir, f"{file_name}_reactant.pkl")
                    logger = Logger(os.path.join(runtime.start_dir, "log"))
                    reactant = Molecule.load_from_pickle(input_file_path)
                    molecules.append(reactant)

                elif os.path.basename(runtime.start_dir) == 'products':
                    if runtime.args.products is False:
                        exit()
                    logger = Logger(os.path.join(runtime.start_dir, "log"))
                    pickle_path = os.path.join(runtime.start_dir, "product_molecules.pkl")
                    product_molecules = Molecule.load_molecules_from_pickle(pickle_path)
                    for product in product_molecules:
                        product.product = True
                        product.current_step = 'crest_sampling'
                        molecules.append(product)
                else:
                    input_file_path = os.path.join(runtime.start_dir, os.path.basename(runtime.start_dir)+'.pkl')
                    molecule = Molecule.load_from_pickle(input_file_path)
                    logger = Logger(os.path.join(runtime.start_dir, "log"))
                    molecules.append(molecule)

                handle_termination(molecules, logger, threads, converged=False)

            elif file_type == '.pkl':
                from pandas import read_pickle, DataFrame
                df = read_pickle(input_file)
                if isinstance(df, DataFrame):
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
                input_file_path = os.path.join(runtime.start_dir, input_file)
                molecule = Molecule(file_path=input_file_path, indexes=runtime.args.CHO, method=runtime.args.method)
                input_molecules.append(molecule)

            else:
                print(f"File type {file_type} is not supported")
                exit()

        if final_TS and final_reactants:
            k, kappa = rate_constant(final_TS, final_reactants, final_products)
            print(f"{k} cm^3 molecules^-1 s^-1")
            print(f"Tunneling coefficient: {kappa}")
            exit()

        if input_molecules and file_type != '.xyz':
            logger = Logger(os.path.join(runtime.start_dir, "log"))
            if runtime.args.collect:
                collected_molecules = collect_DFT_and_DLPNO(input_molecules)
                Molecule.molecules_to_pickle(collected_molecules, os.path.join(runtime.start_dir, "collected_molecules.pkl"))
            elif runtime.args.restart:
                logger.log("\nJKTS restarted")
                with open(os.path.join(runtime.start_dir, '.method'), 'w') as f:
                    f.write(f"{runtime.args.method}")
                handle_input_molecules(input_molecules, logger, threads)
            elif runtime.args.rerun:
                logger.log(f"\nJKTS - rerunning calculations of type: {input_molecules[0].current_step}")
                handle_termination(input_molecules, logger, threads, converged=False)
            elif runtime.args.pickle:
                filename = re.sub("_conf\d{1,2}", "", input_molecules[0].name)
                Molecule.molecules_to_pickle(input_molecules, os.path.join(runtime.start_dir, f"collection{filename}.pkl"))
            else:
                parser.error("Error")

        if final_TS:
            for m in final_TS:
                runtime.global_molecules.append(m)

        elif not input_molecules and file_type != '.xyz':
            logger = Logger(os.path.join(runtime.start_dir, "log"))
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

    if runtime.global_molecules:
        logger = Logger(os.path.join(runtime.start_dir, "log"))
        molecules_logger = Logger(os.path.join(runtime.start_dir, "molecules.txt"))
        for molecule in runtime.global_molecules:
            molecule.print_items(molecules_logger)

        if all(m.reactant for m in runtime.global_molecules):
            molecule_name = runtime.global_molecules[0].name.split("_")[0]
            logger.log(f"Final DLPNO calculations for reactants is done. Logging molecules to Final_reactants_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(runtime.global_molecules, os.path.join(runtime.start_dir, f"Final_reactants_{molecule_name}.pkl"))
        elif all(m.product for m in runtime.global_molecules):
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in molecules if "_H" in m.name))
            grouped_lists = [[m for m in runtime.global_molecules if f"_H{h_num}_" in m.name] +
                             [m for m in runtime.global_molecules if "H2O" in m.name]
                             for h_num in h_numbers]

            for h, molecules_group in zip(h_numbers, grouped_lists):
                molecule_name = f"{molecules_group[0].name.split('_')[0]}_H{h}"
                pickle_path = os.path.join(runtime.start_dir, f"Final_products_{molecule_name}.pkl")
                logger.log(f"Final DLPNO calculations for products are done. Logging properties to Final_products_{molecule_name}.pkl")
                Molecule.molecules_to_pickle(molecules_group, pickle_path)
        else:
            molecule_name = os.path.basename(runtime.start_dir)
            logger.log(f"Final DLPNO calculations for transition state molecules is done. Logging properties to Final_TS_{molecule_name}.pkl")
            TS_pkl_path = os.path.join(runtime.start_dir, f"Final_TS_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(runtime.global_molecules, TS_pkl_path)

            if runtime.args.k:
                reactant_pkl_name = os.path.basename(runtime.start_dir).split("_")[0]
                product_pkl_name = os.path.basename(runtime.start_dir)
                reactant_pkl_path = os.path.join(os.path.dirname(runtime.start_dir), f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                product_pkl_path = os.path.join(os.path.dirname(runtime.start_dir), f'products/Final_products_{product_pkl_name}.pkl')
                logger.log(f"{product_pkl_path}, {reactant_pkl_path}")
                if os.path.exists(reactant_pkl_path):
                    final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                    if os.path.exists(product_pkl_path):
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                    k, kappa = rate_constant(runtime.global_molecules, final_reactants, final_products)
                    results_logger = Logger(os.path.join(os.path.dirname(runtime.start_dir), "Rate_constants.txt"))
                    results_logger.log_with_stars(f"{molecule_name}: {k} cm^3 molecules^-1 s^-1 with tunneling coefficient {kappa}")
                else:
                    reactant_pkl_path = os.path.join(runtime.start_dir, f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                    product_pkl_path = os.path.join(runtime.start_dir, f'products/Final_products_{product_pkl_name}.pkl')
                    logger.log(f"{product_pkl_path}, {reactant_pkl_path}")
                    if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                        final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                        k, kappa = rate_constant(runtime.global_molecules, final_reactants, final_products)
                        results_logger = Logger(os.path.join(runtime.start_dir, "Rate_constants.txt"))
                        results_logger.log_with_stars(f"1 {molecule_name}: {k} cm^3 molecules^-1 s^-1 with tunneling coefficient {kappa}")
                    else:
                        logger.log(f"Could not find pickle files in path: {product_pkl_name} and {reactant_pkl_path}")


if __name__ == "__main__":
    main()
