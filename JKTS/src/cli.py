#!/usr/bin/env python3
"""Shared CLI helpers for JKTS.

This module exports ParseList (an argparse.Action) and str2bool
so they can be reused by other modules (for example DATS.py).
"""
import argparse

__all__ = ["ParseList", "str2bool"]


class ParseList(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        processed_list = []
        if isinstance(values, list):
            flat_list = [item for sublist in values for item in (sublist.split(',') if ',' in sublist else [sublist])]
        else:
            flat_list = [item for item in values.split(',')]

        for item in flat_list:
            try:
                processed_item = int(item)
            except ValueError:
                processed_item = item  # Keep as string if it can't be converted to int

            processed_list.append(processed_item)

        setattr(namespace, self.dest, processed_list)


def str2bool(v):
    """Convert common truthy/falsy strings to boolean for argparse.

    Accepts: yes/no, true/false, t/f, y/n, 1/0 (case-insensitive) and
    will return the corresponding boolean. If a boolean is passed in,
    it is returned unchanged.
    """
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_parser():
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
    parser.add_argument('-constrain', type=str2bool, metavar='<boolean>', default=True, help='Integrate constraints into relevant input files [def: True]')
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
    additional_options.add_argument('-IRC', action='store_true', default=False, help='Perform IRC calcuation on ORCA or G16 output file')
    additional_options.add_argument('-CHO', dest='CHO', action=ParseList, nargs='*', help="Set indexes of atoms for active site. Indexing starting from 1")
    additional_options.add_argument('-collect', action='store_true', default=False, help='Collect thermochemical data from TS structures and single point correction from DLPNO')
    additional_options.add_argument('-method', type=str, default='wB97XD', help='Specify QC method to use for optimization and TS search [def: WB97X-D3BJ]')
    additional_options.add_argument('-basis_set', type=str, default='6-31++g(d,p)', help='Specify basis set to use with QC method [def: 6-31++g(d,p)]')
    # additional_options.add_argument('--low_level', nargs='+', metavar='', help='Specify low-level theory for preoptimization [def: B3LYP 6-31+G(d,p)]')
    additional_options.add_argument('--gfn', default='2', choices=['1','2'], help='Specify the GFN version (1 or 2, default: 2)')
    additional_options.add_argument('-skip_preopt', action='store_true', default=False, help='Skip the preoptimization of the structures before the TS search')
    additional_options.add_argument('-cpu', metavar="int", nargs='?', const=1, type=int, default=4, help='Number of CPUs [def: 4]')
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=8000, help='Amount of memory allocated for the job [def: 8000MB]')
    additional_options.add_argument('-par', metavar="partition", type=str, default="q64", help='Partition to use [def: q64,q48,q28,q24]')
    additional_options.add_argument('-time', metavar="hh:mm:ss", type=str, default=None, help='Monitoring duration [def: 144 hours]')
    additional_options.add_argument('-interval', metavar="int", nargs='?', const=1, type=int, help='Time interval between log file checks [def: based on molecule size]')
    additional_options.add_argument('-initial_delay', metavar="int", nargs='?', const=1, type=int, help='Initial delay before checking log files [def: based on molecule size]')
    additional_options.add_argument('-attempts', metavar="int", nargs='?', const=1, type=int, default=100, help='Number of log file check attempts [def: 100]')
    additional_options.add_argument('-max_conformers', metavar="int", nargs='?', const=1, type=int, default=None, help='Maximum number of conformers from CREST [def: 50]')
    additional_options.add_argument('-freq_cutoff', metavar="int", nargs='?', const=1, type=int, default=-100, help='TS imaginary frequency cutoff [def: -100 cm^-1]')
    additional_options.add_argument('-ewin', metavar="int", nargs='?', const=1, default=5, type=int, help='Energy threshold for CREST conformer sampling [def: 5 kcal/mol]')
    additional_options.add_argument('-F12', action='store_true', help='Use CCSD(T)-F12-pVTZ instead of DLPNO')
    additional_options.add_argument('-energy_cutoff', metavar="digit", nargs='?', const=1, default=5, type=float, help='After preoptimization, remove conformers which are [int] kcal/mol higher in energy than the lowest conformer [def: 5 kcal/mol]')
    additional_options.add_argument('-pickle', action='store_true', default=False, help='Store given log files into pickle file')
    additional_options.add_argument('-filter', type=str2bool, metavar='<boolean>', default=True, help='Filter identical conformers after transition state optimization [def: True]')
    additional_options.add_argument('--account', type=str, help=argparse.SUPPRESS) # For Finland Puhti
    additional_options.add_argument('-test', action='store_true', default=False, help=argparse.SUPPRESS) # Run TEST section in main() and exit
    additional_options.add_argument('-num_molecules', type=int, default=None, help=argparse.SUPPRESS) # Set the number of directories created to [int]. Used when doing limited testing
    additional_options.add_argument('-rerun', action='store_true', help=argparse.SUPPRESS) 


    hidden_options = parser.add_argument_group()
    hidden_options.add_argument('-init', action='store_true', help=argparse.SUPPRESS)
    hidden_options.add_argument('-loc', action='store_true', help=argparse.SUPPRESS)
    hidden_options.add_argument('-reaction_angle', metavar="float", nargs='?', default=169.5, const=1, type=float, help=argparse.SUPPRESS) 
    hidden_options.add_argument('-hybrid', action='store_true', default=False, help=argparse.SUPPRESS) #'Hybrid Hessian. Full Hessian for active site and approximate for rest of atoms'
    hidden_options.add_argument('-dispersion', action='store_true', default=False, help=argparse.SUPPRESS)
    additional_options.add_argument('-skip', type=str2bool, default=False, help=argparse.SUPPRESS)

    return parser