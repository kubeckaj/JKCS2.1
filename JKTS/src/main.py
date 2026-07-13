#!/usr/bin/env python3
import argparse
import os
import re
import sys

import runtime
from classes import Molecule
from output import Logger, console, banner
from qc_input import mkdir
from monitoring import handle_termination, handle_input_molecules
from conformer_tools import initiate_conformers, collect_DFT_and_DLPNO
from rate_constant import rate_constant
from results import record_rate, write_molecule_summary, format_rate


def read_reaction_path_degeneracy(directory):
    # σᵢ written per TS directory by mkdir(); defaults to 1
    for path in (os.path.join(directory, ".symmetry"), os.path.join(directory, "log_files", ".symmetry")):
        try:
            with open(path) as f:
                return int(f.read().strip())
        except (OSError, ValueError):
            continue
    return 1


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
            console.error(f"Unsupported input file '{input_file}'. Without a specified reaction type, JKTS expects '.pkl', '.log', '.out', '.com', '.inp', or '.xyz' files.")
            sys.exit(1)

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


def build_parser():
    parser = argparse.ArgumentParser(description='''    JKTS - automated transition state search and MC-TST rate constants.
    Workflow for an .xyz input: CREST conformer sampling, constrained
    preoptimization, TS optimization, DLPNO-CCSD(T) single points, and
    tunneling-corrected multiconformer rate constants. Failed jobs are
    resubmitted automatically until convergence or the wall-time limit.''',
                                     prog="JKTS",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
    Examples of use:
                JKTS pinonaldehyde.xyz -OH
                JKTS -smiles CO -Cl  # Methanol hydrogen abstraction with Cl radical
                JKTS molecule.xyz -NO3  # Nighttime H abstraction with NO3 radical
                JKTS CH4_H1_opt_constrain.pkl -info
                JKTS benzene.xyz -OH -ORCA -par qtest -stop
                JKTS *TS.log -time 5:00:00
                                     ''')

    parser.add_argument('input_files', metavar='reactant.xyz', nargs='*', help='Input files (.xyz, .pkl, .log, .out, .com, .inp)')

    reaction = parser.add_argument_group("Reaction type")
    reaction.add_argument('-OH', action='store_true', help='H abstraction by OH radical')
    reaction.add_argument('-Cl', action='store_true', help='H abstraction by Cl radical')
    reaction.add_argument('-NO3', action='store_true', help='H abstraction by NO3 radical (nighttime chemistry)')

    workflow = parser.add_argument_group("Workflow control")
    workflow.add_argument('-reactants', type=str2bool, nargs='?', const=False, default=True, metavar='<bool>', help='Skip the reactants workflow (bare flag or explicit false) [def: run it]')
    workflow.add_argument('-products', type=str2bool, nargs='?', const=False, default=True, metavar='<bool>', help='Skip the products workflow (bare flag or explicit false) [def: run it]')
    workflow.add_argument('-TS', type=str2bool, nargs='?', const=False, default=True, metavar='<bool>', help=argparse.SUPPRESS)
    workflow.add_argument('-stop', action='store_true', help='Stop after the current workflow step completes instead of continuing automatically (resume with -restart)')
    workflow.add_argument('-restart', action='store_true', help='Resume the workflow from the current step of the given .pkl/.log/.out files')
    workflow.add_argument('-k', type=str2bool, metavar='<bool>', default=True, help='Calculate the MC-TST rate constant at the end [def: true]')
    workflow.add_argument('-time', metavar="hh:mm:ss", type=str, default=None, help='SLURM wall time for the monitoring job (used by the JKTS wrapper) [def: 240:00:00]')

    qc = parser.add_argument_group("QC settings")
    qc.add_argument('-G16', action='store_true', help='Use Gaussian16 (default)')
    qc.add_argument('-ORCA', action='store_true', help='Use ORCA instead of Gaussian16')
    qc.add_argument('-method', type=str, default='wB97XD', help='QC method for optimization and TS search [def: wB97XD]')
    qc.add_argument('-basis_set', type=str, default='6-31++g(d,p)', help='Basis set for the QC method [def: 6-31++g(d,p)]')
    qc.add_argument('-F12', action='store_true', help='Use CCSD(T)-F12/cc-pVTZ-F12 instead of DLPNO for single points')
    qc.add_argument('--gfn', default='2', choices=['1', '2'], help='GFN version for CREST [def: 2]')
    qc.add_argument('-ewin', metavar="int", type=int, default=5, help='CREST conformer sampling energy window [def: 5 kcal/mol]')
    qc.add_argument('-energy_cutoff', metavar="float", type=float, default=5, help='Drop conformers this much above the lowest after preoptimization [def: 5 kcal/mol]')
    qc.add_argument('-max_conformers', metavar="int", type=int, default=1000, help='Maximum number of conformers taken from CREST [def: 1000]')
    qc.add_argument('-freq_cutoff', metavar="int", type=int, default=-100, help='TS imaginary frequency cutoff [def: -100 cm^-1]')

    slurm = parser.add_argument_group("SLURM and monitoring")
    slurm.add_argument('-par', metavar="partition", type=str, default="q64", help='SLURM partition [def: q64]')
    slurm.add_argument('-cpu', metavar="int", type=int, default=4, help='Number of CPUs per job [def: 4]')
    slurm.add_argument('-mem', metavar="int", type=int, default=8000, help='Memory per job in MB [def: 8000]')
    slurm.add_argument('-interval', metavar="int", type=int, help='Seconds between log file checks [def: based on molecule size]')
    slurm.add_argument('-initial_delay', metavar="int", type=int, help='Seconds before the first log file check [def: based on molecule size]')
    slurm.add_argument('-attempts', metavar="int", type=int, default=100, help='Maximum number of log file checks [def: 100]')

    tools = parser.add_argument_group("Utilities")
    tools.add_argument('-smiles', metavar='string', type=str, help='Input molecule as a SMILES string instead of an .xyz file')
    tools.add_argument('-info', action='store_true', help='Print molecule information from log or .pkl files')
    tools.add_argument('-movie', action='store_true', help='Write movie.xyz of the given structures for viewing')
    tools.add_argument('-collect', action='store_true', help='Collect DFT thermochemistry and DLPNO single points into a .pkl file')
    tools.add_argument('-pickle', action='store_true', help='Store the given log files into a .pkl file')
    tools.add_argument('-CHO', nargs='*', type=int, help='Atom indexes of the active site (C H O), 1-indexed')

    hidden = parser.add_argument_group()
    hidden.add_argument('-init', action='store_true', help=argparse.SUPPRESS)
    hidden.add_argument('-test', action='store_true', help=argparse.SUPPRESS)
    hidden.add_argument('-rerun', action='store_true', help=argparse.SUPPRESS)
    hidden.add_argument('-num_molecules', type=int, default=None, help=argparse.SUPPRESS)
    hidden.add_argument('-hybrid', action='store_true', help=argparse.SUPPRESS)
    hidden.add_argument('-dispersion', action='store_true', help=argparse.SUPPRESS)
    hidden.add_argument('--account', type=str, help=argparse.SUPPRESS)
    # Accepted for compatibility with the JKTS bash wrapper, which passes these through
    hidden.add_argument('-loc', dest='collect', action='store_true', help=argparse.SUPPRESS)
    return parser


def main():
    parser = build_parser()
    runtime.args = parser.parse_args()
    runtime.start_dir = os.getcwd()

    runtime.QC_program = "ORCA" if runtime.args.ORCA else "G16"

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
    logger = Logger(os.path.join(runtime.start_dir, "log"))

    if runtime.args.init:
        if runtime.args.smiles:
            input_molecule = Molecule(smiles=runtime.args.smiles, reactant=True, method=runtime.args.method)
            file_name, file_type = input_molecule.name, 'xyz'
        else:
            input_file = runtime.args.input_files[0]
            file_name, file_type = os.path.splitext(input_file)
            input_file_path = os.path.join(runtime.start_dir, input_file)
            input_molecule = Molecule(input_file_path, reactant=True, method=runtime.args.method)

        if runtime.args.OH or runtime.args.Cl or runtime.args.NO3:
            reacted_molecules, product_molecules = input_molecule.H_abstraction(Cl=runtime.args.Cl, NO3=runtime.args.NO3, products=runtime.args.products, num_molecules=runtime.args.num_molecules)
        else:
            parser.error("Need to specify a reaction type: -OH, -Cl or -NO3")

        if runtime.args.TS:
            for count, molecule in enumerate(reacted_molecules, start=1):
                molecule.name = f"{file_name}_H{count}"
                molecule.directory = os.path.join(runtime.start_dir, molecule.name)
                mkdir(molecule)
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

        radical = 'Cl' if runtime.args.Cl else 'NO3' if runtime.args.NO3 else 'OH'
        workflow_dirs = [d for d, wanted in (('reactants', runtime.args.reactants),
                                             ('products', runtime.args.products),
                                             ('TS', runtime.args.TS)) if wanted]
        banner([("Molecule", file_name),
                ("Reaction", f"H abstraction by {radical}"),
                ("Sites", f"{len(reacted_molecules)} unique abstraction site(s)"),
                ("Method", f"{runtime.args.method} {runtime.args.basis_set}   backend: {runtime.QC_program}"),
                ("SLURM", f"partition {runtime.args.par}, {runtime.args.cpu} CPUs, {runtime.args.mem} MB"),
                ("Workflows", ', '.join(workflow_dirs))])

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
        console.success("movie.xyz generated")

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
                    reactant = Molecule.load_from_pickle(input_file_path)
                    molecules.append(reactant)

                elif os.path.basename(runtime.start_dir) == 'products':
                    if runtime.args.products is False:
                        exit()
                    pickle_path = os.path.join(runtime.start_dir, "product_molecules.pkl")
                    product_molecules = Molecule.load_molecules_from_pickle(pickle_path)
                    for product in product_molecules:
                        product.product = True
                        product.current_step = 'crest_sampling'
                        molecules.append(product)
                else:
                    input_file_path = os.path.join(runtime.start_dir, os.path.basename(runtime.start_dir)+'.pkl')
                    molecule = Molecule.load_from_pickle(input_file_path)
                    molecules.append(molecule)

                logger.info(f"JKTS run started (pid {os.getpid()}, method {runtime.args.method} {runtime.args.basis_set}, program {runtime.QC_program})")
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
                console.error(f"File type {file_type} is not supported")
                sys.exit(1)

        if final_TS and final_reactants:
            symmetry = read_reaction_path_degeneracy(runtime.start_dir)
            result = rate_constant(final_TS, final_reactants, final_products, symmetry=symmetry)
            console.results(format_rate(result))
            exit()

        if input_molecules and file_type != '.xyz':
            if runtime.args.collect:
                collected_molecules = collect_DFT_and_DLPNO(input_molecules)
                Molecule.molecules_to_pickle(collected_molecules, os.path.join(runtime.start_dir, "collected_molecules.pkl"))
            elif runtime.args.restart:
                logger.event(f"JKTS restarted (method {runtime.args.method}, program {runtime.QC_program})")
                with open(os.path.join(runtime.start_dir, '.method'), 'w') as f:
                    f.write(f"{runtime.args.method}")
                handle_input_molecules(input_molecules, logger, threads)
            elif runtime.args.rerun:
                logger.event(f"JKTS rerunning calculations of type: {input_molecules[0].current_step}")
                handle_termination(input_molecules, logger, threads, converged=False)
            elif runtime.args.pickle:
                filename = re.sub("_conf\d{1,2}", "", input_molecules[0].name)
                Molecule.molecules_to_pickle(input_molecules, os.path.join(runtime.start_dir, f"collection{filename}.pkl"))
            else:
                parser.error("Input files were loaded but no action was specified. Use -restart, -rerun, -collect or -pickle.")

        if final_TS:
            for m in final_TS:
                runtime.global_molecules.append(m)

        elif not input_molecules and file_type != '.xyz':
            logger.error("Could not create molecule list from the given input files")
            if file_type in [".log", ".out"]:
                logger.info(".log or .out extension detected. Make sure input files are from ORCA or G16")
            elif file_type == '.pkl':
                logger.info("Detected .pkl file. Make sure the structure of the pickle file is either a python list, set, tuple or pandas.DataFrame")

    # Monitor and handle convergence of submitted jobs
    while threads:
        for thread in list(threads):
            thread.join(timeout=0.1)
            if not thread.is_alive():
                threads.remove(thread)

    if runtime.global_molecules:
        write_molecule_summary(os.path.join(runtime.start_dir, "molecules.txt"),
                               runtime.global_molecules, title=os.path.basename(runtime.start_dir))
        logger.info("Molecule summary written to molecules.txt")

        if all(m.reactant for m in runtime.global_molecules):
            molecule_name = runtime.global_molecules[0].name.split("_")[0]
            logger.success(f"Final DLPNO calculations for reactants are done. Results saved to Final_reactants_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(runtime.global_molecules, os.path.join(runtime.start_dir, f"Final_reactants_{molecule_name}.pkl"))
        elif all(m.product for m in runtime.global_molecules):
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in molecules if "_H" in m.name))
            small_product_names = ('H2O', 'H2O_DLPNO', 'HCl', 'HCl_DLPNO', 'HNO3', 'HNO3_DLPNO')
            grouped_lists = [[m for m in runtime.global_molecules if f"_H{h_num}_" in m.name] +
                             [m for m in runtime.global_molecules if m.name in small_product_names]
                             for h_num in h_numbers]

            for h, molecules_group in zip(h_numbers, grouped_lists):
                molecule_name = f"{molecules_group[0].name.split('_')[0]}_H{h}"
                pickle_path = os.path.join(runtime.start_dir, f"Final_products_{molecule_name}.pkl")
                logger.success(f"Final DLPNO calculations for products are done. Results saved to Final_products_{molecule_name}.pkl")
                Molecule.molecules_to_pickle(molecules_group, pickle_path)
        else:
            molecule_name = os.path.basename(runtime.start_dir)
            logger.success(f"Final DLPNO calculations for transition state molecules are done. Results saved to Final_TS_{molecule_name}.pkl")
            TS_pkl_path = os.path.join(runtime.start_dir, f"Final_TS_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(runtime.global_molecules, TS_pkl_path)

            if runtime.args.k:
                symmetry = read_reaction_path_degeneracy(runtime.start_dir)
                reactant_pkl_name = os.path.basename(runtime.start_dir).split("_")[0]
                product_pkl_name = os.path.basename(runtime.start_dir)
                reactant_pkl_path = os.path.join(os.path.dirname(runtime.start_dir), f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                product_pkl_path = os.path.join(os.path.dirname(runtime.start_dir), f'products/Final_products_{product_pkl_name}.pkl')
                logger.info(f"Looking for reactant and product pickles: {reactant_pkl_path}, {product_pkl_path}")
                if os.path.exists(reactant_pkl_path):
                    final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                    if os.path.exists(product_pkl_path):
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                    result = rate_constant(runtime.global_molecules, final_reactants, final_products, symmetry=symmetry)
                    record_rate(os.path.dirname(runtime.start_dir), molecule_name, result, method=runtime.args.method)
                    logger.results(format_rate(result))
                else:
                    reactant_pkl_path = os.path.join(runtime.start_dir, f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                    product_pkl_path = os.path.join(runtime.start_dir, f'products/Final_products_{product_pkl_name}.pkl')
                    logger.info(f"Looking for reactant and product pickles: {reactant_pkl_path}, {product_pkl_path}")
                    if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                        final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                        result = rate_constant(runtime.global_molecules, final_reactants, final_products, symmetry=symmetry)
                        record_rate(runtime.start_dir, molecule_name, result, method=runtime.args.method)
                        logger.results(format_rate(result))
                    else:
                        logger.error(f"Could not find pickle files: {reactant_pkl_path} and {product_pkl_path}")


if __name__ == "__main__":
    main()
