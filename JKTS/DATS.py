#!/usr/bin/env python3
'''Dynamic Approach for Transition States'''
###############################LIBRARIES#####################################
import argparse
import os
import re
import time
from threading import Thread
import json
from classes import Molecule, Logger
from slurm_submit import submit_array_job, submit_job, update_molecules_status

global_molecules = []
###############################################################################
class SCFAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        algorithm = option_string.strip('-').upper()
        scf_string = f"SCF=({algorithm})"
        if values is not None:
            additional_args = ', '.join(values)
            if additional_args:
                scf_string = f"SCF=({algorithm}, maxcycle={additional_args})"
            else:
                scf_string = f"SCF={algorithm}"
        setattr(namespace, self.dest, scf_string)

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


def load_last_run_args():
    try:
        with open('.last_run_args', 'r') as file:
            return json.load(file)
    except FileNotFoundError:
        return None

def save_last_run_args(args):
    with open('.last_run_args', 'w') as file:
        # Convert args to a dictionary to easily serialize it
        args_dict = vars(args)
        json.dump(args_dict, file)


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def mkdir(molecule, crest_constrain=True):
    if not os.path.exists(molecule.directory):
        os.makedirs(molecule.directory, exist_ok=True)
        if not os.path.exists(os.path.join(molecule.directory, "input_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "input_files")))
        if not os.path.exists(os.path.join(molecule.directory, "log_files")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "log_files")))
        if not os.path.exists(os.path.join(molecule.directory, "slurm_output")):
            os.makedirs(os.path.join(os.path.join(molecule.directory, "slurm_output")))

    if crest_constrain:
        force_constant = 1
        C_index = molecule.constrained_indexes['C']
        H_index = molecule.constrained_indexes['H']
        O_index = molecule.constrained_indexes['O']
        with open (molecule.directory + "/constrain.inp","w") as f:
            f.write("$constrain\n")
            f.write(f"  force constant={force_constant}\n") 
            f.write(f"  distance: {C_index}, {H_index}, auto\n")
            f.write(f"  distance: {H_index}, {O_index}, auto\n")
            f.write("$end\n")

        with open (molecule.directory + "/.constrain","w") as f:
            f.write(f"C: {C_index}\n")
            f.write(f"H: {H_index}\n")
            f.write(f"O: {O_index}\n")

        with open (molecule.directory + "/log_files/.constrain","w") as f:
            f.write(f"C: {C_index}\n")
            f.write(f"H: {H_index}\n")
            f.write(f"O: {O_index}\n")


def pkl_to_xyz(pkl_file_path, max_conformers=50):
    from pandas import read_pickle
    if args.max_conformers:
        max_conformers = args.max_conformers
    all_conformers = []
    with open(pkl_file_path, 'rb') as f:
        df = read_pickle(f) 
        for xyz, el in zip(df[('xyz', 'structure')], df[('log', 'electronic_energy')]):
            coordinates_list = []
            coords = xyz.get_positions()
            atom = xyz.get_chemical_symbols()
            for atom, coord in zip(atom, coords):
                coordinates = [atom, el, coord[0], coord[1], coord[2]]
                coordinates_list.append(coordinates)
            all_conformers.append(coordinates_list)
    return all_conformers[:max_conformers]


def extract_normal_coordinates(molecule):
    '''Return format: [[0.02, -0.02, 0.06], [0.07, -0.11, -0.02], ...]''' 
    if molecule.program.lower() == 'g16':
        with open(molecule.log_file_path, 'r') as file:
            lines = file.readlines()
        
        negative_freq_found = False
        read_xyz = False  # Flag to start reading XYZ coordinates
        xyz_coordinates = []
        for line in lines:
            if "Frequencies --" in line and '-' in line:
                negative_freq_found = True
                continue
            if negative_freq_found and not read_xyz:
                if "Atom  AN      X      Y      Z" in line:  # This indicates the start of the XYZ section
                    read_xyz = True  # Now start reading XYZ coordinates on subsequent lines
                continue  # Skip all lines until XYZ section is reached
            if read_xyz:
                # Assuming the coordinates section ends before a line not matching the XYZ format
                match = re.match(r'\s*\d+\s+\d+\s+(-?\d+\.\d+\s+)*', line)
                if match:
                    xyz = re.findall(r'-?\d+\.\d+', line)[:3]  # Extract first three floats
                    xyz = [float(coord) for coord in xyz]
                    xyz_coordinates.append(xyz)
                else:
                    break 
        return xyz_coordinates[:-1]
    else:
        print("extract_normal_coordinates() was called, but only works for G16 log files")


def check_normal_mode_displacement(molecule):
    if molecule.program.lower() == 'orca':
        imag_index = 0
        imag = -100
        with open(molecule.log_file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if "imaginary mode" in line:
                    line_split = line.split()
                    if abs(float(line_split[1])) > abs(imag):
                        imag = float(line_split[1])
                        imag_index = line_split[0][0]


def normal_mode_displacement_significant(file_content):
    lines = file_content.split('\n')
    
    in_vibrational_frequencies_section = False
    in_normal_modes_section = False
    lowest_imaginary_index = None
    lowest_imaginary_value = None
    
    for line in lines:
        if 'VIBRATIONAL FREQUENCIES' in line:
            in_vibrational_frequencies_section = True
        elif 'NORMAL MODES' in line:
            in_vibrational_frequencies_section = False
            in_normal_modes_section = True
            normal_modes_start = False
            mode_values = []
        elif 'IR SPECTRUM' in line:
            in_normal_modes_section = False
        
        if in_vibrational_frequencies_section:
            if line.strip().startswith(tuple(str(i) for i in range(10))) and 'cm**-1' in line and 'imaginary mode' in line:
                # This line contains a frequency
                parts = line.split()
                index = int(parts[0].rstrip(':'))
                value = float(parts[1])
                if lowest_imaginary_value is None or value < lowest_imaginary_value:
                    lowest_imaginary_value = value
                    lowest_imaginary_index = index
        
        # If we are in the normal modes section, parse it
        if in_normal_modes_section:
            if str(lowest_imaginary_index) in line:
                normal_modes_start = True
            elif normal_modes_start:
                if line.strip() == '':
                    break  # End of the relevant normal mode values
                parts = line.split()
                # Check if the first part is an integer, indicating we're still reading mode values
                try:
                    int(parts[0])  # Just to check if it's an integer
                    mode_values.append(float(parts[1]))  # Assume the value of interest is the second column
                except (ValueError, IndexError):
                    break  # Not a mode value line, end of the section

    # Check if any of the mode values exceed the absolute threshold
    value_exceeds_threshold = any(abs(val) > 0.9 for val in mode_values)

    return value_exceeds_threshold


def bad_active_site(molecule, CH_threshold=1.1, HO_threshold=1.95, CO_threshold=3.2, angle_threshold=160):
    C_index = molecule.constrained_indexes['C']-1
    H_index = molecule.constrained_indexes['H']-1
    O_index = molecule.constrained_indexes['O']-1

    distance_CH = molecule.atom_distance(molecule.coordinates[C_index], molecule.coordinates[H_index])
    distance_HO = molecule.atom_distance(molecule.coordinates[H_index], molecule.coordinates[O_index])
    distance_CO = molecule.atom_distance(molecule.coordinates[C_index], molecule.coordinates[O_index])
    angle_CHO = molecule.calculate_angle(molecule.coordinates[C_index], molecule.coordinates[H_index], molecule.coordinates[O_index])
    if distance_CH < CH_threshold or distance_HO > HO_threshold or distance_CO > CO_threshold:
        if angle_CHO < angle_CHO:
            return True
    return False


def check_transition_state(molecule, logger, threshold=0.6):
    from numpy import array
    from numpy.linalg import norm
    freq_cutoff = -abs(args.freq_cutoff) if args.freq_cutoff else -100
    CH_threshold = 1.1
    HO_threshold = 1.95
    CO_threshold = 3.2
    angle_threshold = 160
    try:
        H_index = molecule.constrained_indexes['H']-1 # python 0-based indexing
        C_index = molecule.constrained_indexes['C']-1
        O_index = molecule.constrained_indexes['O']-1
    except:
        molecule.find_active_site()
        H_index = molecule.constrained_indexes['H']-1 # python 0-based indexing
        C_index = molecule.constrained_indexes['C']-1
        O_index = molecule.constrained_indexes['O']-1

    methyl_C_indexes, aldehyde_groups, ketone_methyl_groups = molecule.identify_functional_groups()
    for aldehyde in aldehyde_groups:
        if aldehyde['C'] == C_index and aldehyde['H'] == H_index:
            freq_cutoff = -80
            threshold = 0.4
            CH_threshold = 1.1
            HO_threshold = 1.98
            CO_threshold = 3.2
            angle_threshold = 110

    bad_TS = bad_active_site(molecule, CH_threshold, HO_threshold, CO_threshold, angle_threshold)

    logger.log(f"Checking transition state of {molecule.name}")
    sorted_negative_freqs = sorted((freq for freq in molecule.vibrational_frequencies if freq < 0))
    if sorted_negative_freqs:
        imag = sorted_negative_freqs[0]
        if imag < freq_cutoff:
            logger.log("Imaginary frequency found under cutoff. Proceeding with check of normal coordinates of transition state")
        elif freq_cutoff < imag < 0:
            logger.log(f"Small negative frequency between cutoff {freq_cutoff} and 0: {imag} for molecule {molecule.name}")
            if molecule.error_termination_count == 0:
                logger.log(f"This may indicate wrong transition state. Checking geometry")
                if bad_TS:
                    logger.log("Unusual bond lengths and angles detected. Trying to correct")
                    molecule.set_active_site(indexes=args.CHO)
                    return False
                else:
                    logger.log(f"No unusual bond angles or lengths detected. However, dropping molecule {molecule.name}, due to small imaginary frequency.")
                    return False
            else:
                logger.log(f"Already tried to correct geoemtry for {molecule.name}. Dropping molecule from conformers. However, check by IRC should be done")
                return False

    else:
        logger.log("No negative frequency found for transition state. Trying to correct geometry and resubmit.")
        molecule.set_active_site(indexes=args.CHO)  # testing this in pinic_acid/bad_TS_test
        return False

    if molecule.program.lower() == 'g16':
        normal_coords = extract_normal_coordinates(molecule)
        
        if len(normal_coords) != len(molecule.coordinates):
            print("Error: The number of normal mode displacements does not match the number of atoms in the molecule.")
            return 

        # Copy original coordinates and apply displacement
        displaced_coordinates_plus = [array(original) + array(displacement) for original, displacement in zip(molecule.coordinates, normal_coords)]
        displaced_coordinates_minus = [array(original) - array(displacement) for original, displacement in zip(molecule.coordinates, normal_coords)]

        # Calculate original and new distances between H and O
        original_distance_HO = norm(array(molecule.coordinates[H_index]) - array(molecule.coordinates[O_index]))
        new_distance_HO_plus = norm(displaced_coordinates_plus[H_index] - displaced_coordinates_plus[O_index])
        new_distance_HO_minus = norm(displaced_coordinates_minus[H_index] - displaced_coordinates_minus[O_index])
        relative_change_HO_minus = abs((new_distance_HO_minus - original_distance_HO) / original_distance_HO)
        relative_change_HO_plus = abs((new_distance_HO_plus - original_distance_HO) / original_distance_HO)

        original_distance_CH = norm(array(molecule.coordinates[C_index]) - array(molecule.coordinates[H_index]))
        new_distance_CH_plus = norm(displaced_coordinates_plus[C_index] - displaced_coordinates_plus[H_index])
        new_distance_CH_minus = norm(displaced_coordinates_minus[C_index] - displaced_coordinates_minus[H_index])
        relative_change_CH_minus = abs((new_distance_CH_minus - original_distance_CH) / original_distance_CH)
        relative_change_CH_plus = abs((new_distance_CH_plus - original_distance_CH) / original_distance_CH)

        combined_values_plus_minus = relative_change_CH_plus + relative_change_HO_minus
        combined_values_minus_plus = relative_change_CH_minus + relative_change_HO_plus
        combined_values_plus = relative_change_CH_plus + relative_change_HO_plus
        combined_values_minus = relative_change_CH_minus + relative_change_HO_minus
        # print(f"{molecule.name}, Bad geometry: {bad_TS}, {combined_values_plus:.3f}, {combined_values_minus:.3f}, {combined_values_plus_minus:.3f},  {combined_values_minus_plus:.3f}")


        if any(value > threshold for value in [combined_values_plus, combined_values_minus, combined_values_plus_minus, combined_values_minus_plus]) and not bad_TS:
            logger.log_with_stars(f"Yay! Normal mode analysis indicate correct TS for {molecule.name} with imaginary frequency: {imag}")
            return True
        else:
            if imag < freq_cutoff:
                logger.log(f"Change in bond length does not meet threshold. However, magnitude of imaginiary frequency fulfills criteria. Dropping {molecule.name} for now, but check geometry")
                return False
            else:
                logger.log("The change in bond length does not meet the threshold for significant shortening or elongation and neither did the magnitude of imaginary frequency.")
                return False

    elif molecule.program.lower() == 'orca':
        with open(molecule.log_file_path, 'r') as f:
            content = f.read()
        if normal_mode_displacement_significant(content) and not bad_active_site(molecule):
            logger.log_with_stars(f"Yay! Normal mode analysis indicate correct TS for {molecule.name} with imaginary frequency: {imag}")
            return True
        else:
            if imag < freq_cutoff:
                logger.log(f"Change in bond length does not meet threshold. However, magnitude of imaginiary frequency fulfills criteria. Dropping {molecule.name} for now, but check geometry")
                return False
            else:
                logger.log("The change in bond length does not meet the threshold for significant shortening or elongation and neither did the magnitude of imaginary frequency.")
                return False


def ArbAlign_compare_molecules(molecules, logger=None, pickle=True, RMSD_threshold=0.38):
    from ArbAlign import compare
    initial_len = len(molecules)
    filtered_molecules = []

    while molecules:
        reference = molecules.pop(0) # Take first molecule and assume its unique for now
        filtered_molecules.append(reference)

        non_similar_molecules = []

        for molecule in molecules:
            rmsd_value = compare(reference, molecule)
            if not logger:
                print(f"Molecule: {molecule.name}, Reference: {reference.name}, RMSD: {rmsd_value}")
            if rmsd_value <= RMSD_threshold:
                if molecule.electronic_energy and reference.electronic_energy:
                    # Compare their energies and keep the one with lower energy
                    if molecule.electronic_energy <= reference.electronic_energy:
                        filtered_molecules.pop()  # Remove the previous reference
                        filtered_molecules.append(molecule)  # Add the new one with lower energy
                        reference = molecule  # Update reference to the new molecule
            else:
                non_similar_molecules.append(molecule)

        molecules = non_similar_molecules

    if logger:
        logger.log(f"Filtered {initial_len} conformers to {len(filtered_molecules)} conformers using RMSD threshold: {RMSD_threshold}")
    else:
        print(f"Filtered {initial_len} conformers to {len(filtered_molecules)} conformers using RMSD threshold: {RMSD_threshold}")

    if pickle:
        Molecule.molecules_to_pickle(filtered_molecules, os.path.join(start_dir, "filtered_molecules.pkl"))
    return filtered_molecules


def ArbAlign_pair(TS_conformers, reactants=None, products=None):
    from copy import deepcopy
    from ArbAlign import compare

    TS_reactant_pairs = []
    TS_product_pairs = []

    try:
        H_index = TS_conformers[0].constrained_indexes['H'] - 1
    except KeyError:
        # If extraction fails, attempt to find the active site and retry
        TS_conformers[0].find_active_site()
        H_index = TS_conformers[0].constrained_indexes['H'] - 1

    excluded_indexes_TS = {H_index, len(TS_conformers[0].atoms) - 1, len(TS_conformers[0].atoms) - 2}
    excluded_indexes_reactant = {H_index}

    for TS in TS_conformers:
        TS_copy = deepcopy(TS)
        TS_copy.atoms = [atom for idx, atom in enumerate(TS_copy.atoms) if idx not in excluded_indexes_TS]
        TS_copy.coordinates = [coord for idx, coord in enumerate(TS_copy.coordinates) if idx not in excluded_indexes_TS]

        if reactants:
            reactants = [mol for mol in reactants if 'OH' not in mol.name]
            min_rmsd_reactant = float('inf')
            best_match_reactant = None  # Now storing the reactant object

            for reactant in reactants:
                if not 'OH' in reactant.name:
                    reactant_copy = deepcopy(reactant)
                    reactant_copy.atoms = [atom for idx, atom in enumerate(reactant_copy.atoms) if idx not in excluded_indexes_reactant]
                    reactant_copy.coordinates = [coord for idx, coord in enumerate(reactant_copy.coordinates) if idx not in excluded_indexes_reactant]

                    rmsd = compare(TS_copy, reactant_copy)
                    if rmsd < min_rmsd_reactant:
                        min_rmsd_reactant = rmsd
                        best_match_reactant = reactant  # Storing the reactant object

            if best_match_reactant:
                print(f"{TS.name}, {best_match_reactant.name}, {min_rmsd_reactant}")
                TS_reactant_pairs.append((TS, best_match_reactant, min_rmsd_reactant))

        if products:
            products = [mol for mol in products if 'H2O' not in mol.name]
            min_rmsd_product = float('inf')
            best_match_product = None  # Preparing to store the product object

            for product in products:
                if not 'H2O' in product.name:
                    product_copy = deepcopy(product)

                    rmsd = compare(TS_copy, product_copy)
                    if rmsd < min_rmsd_product:
                        min_rmsd_product = rmsd
                        best_match_product = product  # Storing the product object

            if best_match_product:
                print(f"Best: {TS.name}, {best_match_product.name}, {min_rmsd_product}")
                TS_product_pairs.append((TS, best_match_product, min_rmsd_product))

    return TS_reactant_pairs, TS_product_pairs


def QC_input(molecule, constrain, method, basis_set, TS):
    file_name = f"{molecule.name}{molecule.input}"
    file_path = os.path.join(molecule.directory, file_name)
    atoms = molecule.atoms
    coords = molecule.coordinates
    maxiter = 150 # Maximum iterations for geoemtry optimization
    if molecule.product or method == 'DLPNO' or 'opt_constrain' in molecule.current_step:
        freq = ''
    else:
        freq = 'freq'
    SCF = 'NoTrah'
    disp = 'EmpiricalDispersion=GD3BJ' if args.dispersion else ''


    if constrain or TS:
        if not molecule.constrained_indexes:
            molecule.find_active_site(indexes=args.CHO)
        if molecule.program.lower() == 'orca': # ORCA indexes from 0
            C_index = molecule.constrained_indexes['C']-1
            H_index = molecule.constrained_indexes['H']-1
            O_index = molecule.constrained_indexes['O']-1
        else:
            C_index = molecule.constrained_indexes['C']
            H_index = molecule.constrained_indexes['H']
            O_index = molecule.constrained_indexes['O']

    if molecule.program.lower() == "orca" and molecule.converged is False:
        with open(file_path, "w") as f:
            if method == "DLPNO":
                # Consider scaling of the memory with the molecule size
                if molecule.error_termination_count == 1:
                    f.write(f"! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) VeryTightSCF RI-JK aug-cc-pVTZ/JK {SCF}\n")
                    # f.write(f"! cc-pVTZ-F12 cc-pVTZ/C cc-pVTZ-F12-CABS DLPNO-CCSD(T)-F12 TightPNO TightSCF\n")
                    f.write(f"%pal nprocs {args.cpu} end\n")
                    f.write(f"%maxcore {round((args.mem+16000)/args.cpu)}\n") 
                elif molecule.error_termination_count == 2:
                    f.write(f"! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) VeryTightSCF RI-JK aug-cc-pVTZ/JK {SCF}\n")
                    f.write(f"%pal nprocs {args.cpu} end\n")
                    f.write(f"%maxcore {round((args.mem+20000)/args.cpu)}\n") 
                else:
                    f.write(f"! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF RI-JK aug-cc-pVTZ/JK {SCF}\n")
                    f.write(f"%pal nprocs {args.cpu} end\n")
                    f.write(f"%maxcore {round((args.mem+12000)/args.cpu)}\n") 
            elif TS:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OptTS defgrid2 {freq}\n")
                f.write(f"%pal nprocs {args.cpu} end\n")
                f.write(f"%maxcore {round(args.mem/args.cpu)}\n")
            else:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OPT defgrid2 {freq}\n")
                f.write(f"%pal nprocs {args.cpu} end\n")
                f.write(f"%maxcore {round(args.mem/args.cpu)}\n")
            if constrain:
                f.write("%geom\n")
                f.write("Constraints\n")
                f.write(f"{{B {C_index} {H_index} C}}\n") # Bond being broken
                f.write(f"{{B {H_index} {O_index} C}}\n") # Bond being formed
                # f.write(f"{{A {C_index-1} {H_index-1} {O_index-1} C}}\n") # CHO angle
                f.write("end\nend\n")
            elif TS:
                f.write("%geom\n")
                f.write(f"maxiter {maxiter}\n")
                f.write("Calc_Hess true\n")
                if molecule.error_termination_count == 0:
                    f.write("Recalc_Hess 5\n")
                else:
                    f.write("Recalc_Hess 1\n")
                    f.write("Trust -0.2\n") # Maybe -0.1?
                if args.hybrid:
                    f.write(f"Hybrid_Hess {{ {C_index} {H_index} {O_index} }} end\n")
                f.write(f"TS_Mode {{ B {H_index} {O_index} }} end\n")
                f.write(f"TS_Active_Atoms {{ {C_index} {H_index} {O_index} }} end\n")
                f.write("TS_Active_Atoms_Factor 3\n")
                f.write("end\n")
            f.write("\n")
            f.write(f"* xyz {molecule.charge} {molecule.mult}\n") 
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("*")
########################################################G16################################################################
    elif molecule.program.lower() == "g16" and molecule.converged is False:
        with open(file_path, "w") as f:
            f.write(f"%nprocshared={args.cpu}\n")
            f.write(f"%mem={args.mem}mb\n")
            if TS and constrain: # should only be in the case of hard convergence problems where some flexible parts should be constrained.
                f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,modredundant,MaxCycles={maxiter}) {freq} {disp} {args.SCF}\n\n')
            elif TS and constrain is False:
                if molecule.error_termination_count == 0:
                    f.write(f'# {method} {basis_set} opt=(calcfc,ts,noeigen,MaxCycles={maxiter},ReCalcFC=5) {freq} {disp} {args.SCF}\n\n') # RecalcFC=N also option, recalc Hessian every N iteration
                else:
                    f.write(f'# {method} {basis_set} opt=(calcall,ts,noeigen,MaxCycles={maxiter}) {freq} {disp} {args.SCF}\n\n')
            elif TS is False and constrain:
                f.write(f'# {method} {basis_set} opt=modredundant {freq} {disp} {args.SCF}\n\n')
            else:
                f.write(f"# {method} {basis_set} opt {freq} {disp} {args.SCF}\n\n")
            f.write("Title\n\n")
            f.write(f"{molecule.charge} {molecule.mult}\n")
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("\n")
            if constrain:
                f.write(f"B {C_index} {H_index} F\n")
                f.write(f"B {H_index} {O_index} F\n")
                # f.write(f"A {C_index} {H_index} {O_index} F\n") # CHO angle
                f.write("\n")
    else:
        print(f"QC_input was called but no program was specified for {molecule.name}")


def IRC_submit(molecules):
    maxiter = 50
    step_size = 10
    path = "both"
    for molecule in molecules:
        molecule.program = QC_program
        molecule.name += '_IRC'
        file_name = f"{molecule.name}{molecule.input}"
        file_path = os.path.join(molecule.directory, file_name)
        atoms = molecule.atoms
        coords = molecule.coordinates
        if molecule.program.lower() == "orca":
            with open(file_path, "w") as f:
                f.write(f"! {args.method} {args.basis_set} TightSCF defgrid3 freq IRC\n")
                f.write(f"%pal nprocs {args.cpu} end\n")
                f.write(f"%maxcore {round(args.mem/args.cpu)}\n")
                f.write("%irc\n")
                f.write(f"maxiter {maxiter}\n")
                f.write(f"Direction {path}\n")
                f.write(f"PrintLevel 1\n")
                f.write(f"InitHess calc_anfreq\n")
                f.write("end\n\n")
                f.write(f"* xyz {molecule.charge} {molecule.mult}\n") 
                for atom, coord in zip(atoms, coords):
                    f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
                f.write("*")
        elif molecule.program.lower() == 'g16':
            with open(file_path, "w") as f:
                f.write(f"%nprocshared={args.cpu}\n")
                f.write(f"%mem={args.mem}mb\n")
                f.write(f'# {args.method} {args.basis_set} irc=(calcall, maxpoints=300, stepsize=10, maxcycle={maxiter}) Guess(Mix,Always) {args.SCF}\n\n')
                f.write("Title\n\n")
                f.write(f"{molecule.charge} {molecule.mult}\n")
                for atom, coord in zip(atoms, coords):
                    f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
                f.write("\n")

    
    if len(molecules) == 1:
        submit_job(molecules[0], args)
    else:
        submit_array_job(molecules, args)

    
# def NEP_input(file_path, file_name):
#     if args.NEB:
#         with open(file_path + "/NEB_TS.inp", "w") as f:
#             f.write(f'! B3LYP 6-31+g(d,p)  NEB-TS FREQ\n')
#             f.write(f'%NEB PREOPT_ENDS TRUE NEB_END_XYZFILE "{file_path + "/" + file_name}_product.xyz" END\n')
#             f.write(f'* XYZfile 0 2 {file_path + "/" + file_name}_reactant.xyz\n')


def energy_cutoff(molecules, logger=None):
    cutoff = args.energy_cutoff if args.energy_cutoff else 5
    Htokcalmol = 627.509

    lowest_energy = sorted(molecules, key=lambda molecule: molecule.electronic_energy)[0].electronic_energy * Htokcalmol
    filtered_molecules = [m for m in molecules if m.electronic_energy*Htokcalmol - lowest_energy <= cutoff]

    return filtered_molecules


def read_last_lines(filename, logger, num_lines, interval=200):
    attempts = 0
    max_attempts = 10

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
            # logger.log(f"Attempt {attempts}: File {filename} not found. Sleeping for {interval} seconds before retrying.")
            time.sleep(interval)

    logger.log(f"Failed to find the file {filename} after {max_attempts} attempts.")
    return []


def resubmit_job(molecule, logger, error=None):
    molecule.move_failed()
    job_type = molecule.current_step
    if job_type in ['opt_constrain', 'opt_constrain_conf']:
        QC_input(molecule, constrain=True, method=args.method, basis_set=args.basis_set, TS=False)

    elif job_type in ['optimization', 'optimization_conf']:
        QC_input(molecule, constrain=False, method=args.method, basis_set=args.basis_set, TS=False)

    elif job_type in ['TS_opt', 'TS_opt_conf']:
        QC_input(molecule, constrain=False,  method=args.method, basis_set=args.basis_set, TS=True)

    elif job_type == 'DLPNO':
        QC_input(molecule, constrain=False, method="DLPNO", basis_set=args.basis_set, TS=False)

    elif job_type == 'optimization':
        QC_input(molecule, constrain=False, method=args.method, basis_set=args.basis_set, TS=False)

    else:
        logger.log(f"Error determining job type for resubmission of {molecule.name}")

    job_id, _ = submit_job(molecule, args)
    molecule.job_id = f"{job_id}"
    logger.log(f"submitted file {molecule.name}{molecule.input} with job type: {job_type} and job id {molecule.job_id} in directory {molecule.directory}")


def check_convergence(molecules, logger, threads, interval, max_attempts):
    initial_delay = args.initial_delay if args.initial_delay else int(interval * 2)
    interval = args.interval if args.interval else int(interval)
    attempts = 0
    sleeping = False
    max_terminations_allowed = 2
    pending, running = [], []
    job_type = molecules[0].current_step

    all_converged = False
    for m in molecules:  # Initialize with all molecules not being converged and no terminations counted
        m.converged = False
        m.error_termination_count = 0
        m.log_file_path = os.path.join(m.directory, f"{m.name}{m.output}")

    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)

    while attempts < max_attempts and not all_converged:
        i = 0
        while i < len(molecules):
            molecule = molecules[i]
            if molecule.error_termination_count >= max_terminations_allowed:
                logger.log(f"!!! Dropping molecule conformer {molecule.name} due to repeated error terminations!!!")
                molecules.pop(i)
                molecule.move_failed()
                # Check if all remaining molecules are converged
                if all(m.converged for m in molecules):
                    all_converged = True
                    break
            else:
                i += 1
        if all_converged:
            break

        update_molecules_status(molecules)

        for molecule in molecules:
            if molecule.converged:
                continue
            if molecule.status == 'pending':
                if molecule.job_id not in pending:
                    logger.log(f"Job {molecule.job_id} is pending in the queue.")
                    pending.append(molecule.job_id)
                continue

            elif molecule.status in ['running', 'completed or not found'] or not molecule.converged:
                if molecule.job_id not in running:
                    logger.log(f"Job {job_type} for {molecule.name} with job id {molecule.job_id} is running.") if molecule.status == 'running' else None
                    molecule.move_inputfile()
                    running.append(molecule.job_id)
                if molecule.job_id in pending:
                    pending.remove(molecule.job_id)

                normal_termination_detected, termination_string = termination_status(molecule, logger)

                if normal_termination_detected:
                    logger.log(termination_string) if termination_string else None
                    molecule.converged = True
                    molecule.error_termination_count = 0
                elif normal_termination_detected is False:
                    if termination_string in ['Error string not found', 'Could not read molecule log file']:
                        continue
                    else:
                        molecule.error_termination_count += 1
                        if molecule.error_termination_count >= max_terminations_allowed:
                            continue
                        handle_error_termination(molecule, logger, termination_string)
            else:
                logger.log(f"Status {molecule} could not be determined. Ensure it is running. Job id: {molecule.job_id}")
                continue
            time.sleep(10)
                
        if all(m.converged for m in molecules):
            all_converged = True
            break

        if len(pending) >= max(1, int(len(molecules) / 1.5)):
            if len(pending) == len(molecules):
                msg = "All the submitted jobs are pending. Sleeping for now."
                sleep_time = 2*interval
            else: 
                msg = "Majority of submitted jobs are pending. Sleeping for now."
                sleep_time = interval
            if not sleeping:
                logger.log(msg)
                time.sleep(sleep_time)
                sleeping = True
            time.sleep(interval)
        else:
            attempts += 1
            if attempts % 10 == 0 or attempts == 1:
                logger.log(f"Log files of the {len(molecules)} conformers have been checked. Checking every {interval} seconds. Attempt: {attempts}/{max_attempts}")
            time.sleep(interval)

    if all_converged:
        if molecules: # Check if not all molecules from list has been dropped
            dir = molecules[0].directory
            basename = os.path.basename(dir)
            pickle_path = os.path.join(dir, f'{basename}_{job_type}.pkl')
            molecules[0].move_files()
            if job_type == 'crest_sampling':
                logger.log_with_stars(f"Yay! CREST job converged. {len(molecules)} conformers generated from CREST")
            elif len(molecules) > 1:
                Molecule.molecules_to_pickle(molecules, pickle_path)
                logger.log_with_stars(f"Yay! All conformer jobs have converged for job type: {job_type}.")
            if job_type == "DLPNO":
                for molecule in molecules: 
                    global_molecules.append(molecule)
                if molecules[0].product and any(mol for mol in global_molecules if 'H2O' in mol.name) is False:
                    H2O = Molecule.create_H2O()
                    H2O.program = QC_program
                    QC_input(H2O, constrain=False, method=args.method, basis_set=args.basis_set, TS=False)
                    submit_and_monitor(H2O, logger, threads)
                elif molecules[0].reactant and any(mol for mol in global_molecules if 'OH' in mol.name) is False:
                    OH = Molecule.create_OH()
                    OH.program = QC_program
                    QC_input(OH, constrain=False, method=args.method, basis_set=args.basis_set, TS=False)
                    submit_and_monitor(OH, logger, threads)
                return True
            elif args.auto:
                handle_termination(molecules, logger, threads, converged=True)
                return True
            else:
                logger.log(f"Automatic processing of calculations has been turned off and JKTS will not proceed with next step: {molecules[0].next_step}")
                logger.log(f"A pickle file {basename}_{job_type}.pkl has been created. Next step {molecules[0].next_step} for molecules can be started from this.")
                return True
        else:
            logger.log(f"No conformer has managed to converge for job type: {job_type}")
            print("Program terminated")
            exit()


def handle_error_termination(molecule, logger, error_termination_string):
    if molecule.program.lower() == 'orca':
        logger.log(f"Error termination was found in {molecule.name}. Resubmitting")
        xyz_coordinates = molecule.log2xyz()
        molecule.coordinates = xyz_coordinates
        resubmit_job(molecule, logger)
    else:
        convergence_errors = ["l9999", "l508"]
        intervention_errors = ["l301"]
        G16_common_errors = "https://wongzit.github.io/gaussian-common-errors-and-solutions/"
        last_lines = read_last_lines(molecule.log_file_path, logger, 30)

        detected_convergence_errors = [error for error in convergence_errors if any(error in line for line in last_lines)]
        detected_intervention_errors = [error for error in intervention_errors if any(error in line for line in last_lines)]

        if detected_convergence_errors:
            error = detected_convergence_errors[0]
            logger.log(f"Convergence error '{error}' termination found in {molecule.name}")
            xyz_coordinates = molecule.log2xyz()
            molecule.coordinates = xyz_coordinates
            logger.log(f"Trying to resubmit job for {molecule.name} due to error termination")
            resubmit_job(molecule, logger, error)
        elif detected_intervention_errors:
            error = detected_intervention_errors[0]
            logger.log(f"Error '{error}' detected in {molecule.name} which needs taken care of manually")
            logger.log(f"Removing the conformer {molecule.name} for now so user can inspect error")
            molecule.error_termination_count = 3
            if molecule.program.lower() == 'g16':
                logger.log(f"Common G16 errors can be found in: {G16_common_errors}")
        else:
            logger.log(f"Error termination found in {molecule.name}. Trying to resubmit")
            xyz_coordinates = molecule.log2xyz()
            molecule.coordinates = xyz_coordinates
            resubmit_job(molecule, logger)


def check_crest(molecules, logger, threads, interval, max_attempts):
    initial_delay = args.initial_delay if args.initial_delay else int(interval * 2)
    interval = args.interval if args.interval else int(interval)
    attempts = 0
    sleeping = False
    pending = []
    all_conformers = []
    expected_files = {f"collection{molecule.name}.pkl" for molecule in molecules}
    constrained_indexes = molecules[0].constrained_indexes
    mult = molecules[0].mult
    charge = molecules[0].charge
    current_step = molecules[0].current_step
    dir = molecules[0].directory
    
    logger.log(f"Waiting {initial_delay} seconds before first check")
    time.sleep(initial_delay)

    while attempts < max_attempts:
        update_molecules_status(molecules)
        for molecule in molecules:
            if molecule.status == 'pending':
                if molecule.job_id not in pending:
                    logger.log(f"Job {molecule.job_id} is pending in the queue.")
                    pending.append(molecule.job_id)
                continue
            else:
                if molecule.job_id in pending:
                    pending.remove(molecule.job_id)
        if len(pending) >= max(1, int(len(molecules) / 1.5)):
            if len(pending) == len(molecules):
                msg = "All the submitted jobs are pending. Sleeping for now."
                sleep_time = 2*interval
            else: 
                msg = "Majority of submitted jobs are pending. Sleeping for now."
                sleep_time = interval
            if not sleeping:
                logger.log(msg)
                time.sleep(sleep_time)
                sleeping = True
            time.sleep(interval)
            continue

        try:
            files_in_directory = set(os.listdir(molecules[0].directory))
        except FileNotFoundError:
            logger.log(f"Pickle file(s) not generated yet. Retrying in {interval} seconds.")
            time.sleep(interval)
            attempts += 1
            continue

        if expected_files.issubset(files_in_directory):
            for molecule in molecules:
                try:
                    pickle_file_path = os.path.join(molecule.directory, f"collection{molecule.name}.pkl")
                    conformers = initiate_conformers(pickle_file_path)
                    logger.log_with_stars(f"{len(conformers)} conformers generated for {molecule.name.replace('_CREST', '')}")
                    for conf in conformers:
                        conf.constrained_indexes = constrained_indexes
                        conf.mult = mult
                        conf.charge = charge
                        conf.current_step = current_step
                        conf.directory = dir
                        all_conformers.append(conf)
                    molecule.move_inputfile() #NOTE: test
                    molecule.move_converged()
                except Exception as e:
                    logger.log(f"Error processing molecule {molecule.name}: {e}")
                    return False
            logger.log(f"CREST conformers generated. Proceeding with next step: {all_conformers[0].next_step}")
            handle_termination(all_conformers, logger, threads, converged=True)
            return True
        else:
            if attempts == 1:
                logger.log(f"Not all files found. Retrying every {interval} seconds.")
            time.sleep(interval)

        attempts += 1

    return False


def submit_and_monitor(molecules, logger, threads):
    if isinstance(molecules, Molecule):
        molecules = [molecules]

    if len(molecules) == 1:
        job_id, interval = submit_job(molecules[0], args)
        logger.log(f"Submitting file {molecules[0].name}{molecules[0].input} for calculation in path {molecules[0].directory} with job id {job_id}")
    else:
        job_id, interval = submit_array_job(molecules, args)
        logger.log(f"Submitted SLURM array job with job id {job_id} for conformers in {molecules[0].directory}")

    if job_id:
        if molecules[0].current_step == 'crest_sampling':
            thread = Thread(target=check_crest, args=(molecules, logger, threads, interval, args.attempts))
        else:
            thread = Thread(target=check_convergence, args=(molecules, logger, threads, interval, args.attempts))
        threads.append(thread)
        thread.start()
    else: logger.log("Error getting job id")


def handle_termination(molecules, logger, threads, converged):
    if not isinstance(molecules, list):
        try:
            molecules = list(molecules)
        except TypeError:
            raise ValueError("molecules must be a list or convertible to a list")
    if converged:
        for m in molecules: 
            m.converged = False
            m.update_step()
    current_step = molecules[0].current_step
    if args.skip_preopt and current_step in ["opt_constrain", "opt_constrain_conf"]: # Works as long as only TS molecules are associated with the opt_constrain(_conf)
        for m in molecules:
            m.update_step() # Update step to skip preoptimization
        current_step = molecules[0].current_step
    logger.log(f"Job to be performed: {current_step}")
    if current_step in args.filter_step and converged and 'H2O' not in molecules[0].name:
        if molecules[0].product: # For products ArbAlign is needed to be done on every individual H
            conformer_molecules = []
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in molecules if "_H" in m.name))
            grouped_lists = [[m for m in molecules if f"_H{h_num}_" in m.name] for h_num in h_numbers]
            for h, group in zip(h_numbers, grouped_lists):
                if len(group) > 5:
                    logger.log(f"Filtering product molecules for H{h} using ArbAlign alogrithm")
                    filtered_group = ArbAlign_compare_molecules(group, logger)
                    for molecule in filtered_group:
                        conformer_molecules.append(molecule)
                else:
                    logger.log(f"Skipping filtering for H{h} as the number of conformers is less than 5")
                    for molecule in group:
                        conformer_molecules.append(molecule)
        else:
            if len(molecules) > 5:
                logger.log("Filtering molecules using ArbAlign algorithm")
                conformer_molecules = ArbAlign_compare_molecules(molecules, logger)
            else:
                conformer_molecules = molecules
                logger.log("Skipping filtering as list of molecules is too short")
    else:
        conformer_molecules = molecules

    for conf in conformer_molecules:
        if conf.converged is False:
            conf.name = conf.name.replace("_TS", "").replace("_CREST", "").replace("_DLPNO", "")
            job_type = conformer_molecules[0].current_step

            if job_type == 'crest_sampling':
                conf.name += '_CREST'
                conf.program = 'CREST'
                output_file_path = os.path.join(conf.directory, f"{conf.name}.xyz")
                conf.write_xyz_file(output_file_path)
            elif job_type in ['opt_constrain', 'opt_constrain_conf']:
                conf.program = QC_program
                QC_input(conf, constrain=True, method=args.method, basis_set=args.basis_set, TS=False)

            elif job_type in ['optimization', 'optimization_conf']:
                conf.program = QC_program
                QC_input(conf, constrain=False, method=args.method, basis_set=args.basis_set, TS=False)
            
            elif job_type in ['TS_opt', 'TS_opt_conf']:
                conf.program = QC_program
                conf.name += '_TS'
                QC_input(conf, constrain=False, method=args.method, basis_set=args.basis_set, TS=True)

            elif job_type == 'DLPNO':
                conf.program = 'ORCA'
                conf.name += '_DLPNO'
                QC_input(conf, constrain=False, method='DLPNO', basis_set=args.basis_set, TS=False)

            elif job_type  == 'Done':
                logger.log("DLPNO calculation has converged")
                break

            else:
                logger.log("Job type could not be determined for conformer list")
                break
    if conformer_molecules: 
        submit_and_monitor(conformer_molecules, logger, threads)


def termination_status(molecule, logger):
    last_lines = read_last_lines(molecule.log_file_path, logger, 30)
    if not last_lines:
        molecule.error_termination_count += 1
        return False, 'Could not read molecule log file'

    termination_detected = any(
        termination in line.lower() for termination in termination_strings[molecule.program] for line in last_lines)
    
    if termination_detected:
        # Update geometry and energy
        xyz_coordinates = molecule.log2xyz()
        molecule.coordinates = xyz_coordinates
        molecule.update_energy(logger)
        if 'TS' in molecule.current_step:
            ts_check_passed = check_transition_state(molecule, logger)
            msg = "" if ts_check_passed else 'Normal termination, but transition state check didn’t pass'
            return ts_check_passed, msg
        else:
            return True, f"***{molecule.name} converged***"

    for error_termination in error_strings[molecule.program]:
        if any(error_termination in line.lower() for line in last_lines):
            return False, error_termination

    return False, 'Error string not found'


def initiate_conformers(input_file=None):
    from pandas import read_pickle
    import os

    conformer_molecules = []
    pickle_file_path = input_file
    max_conformers = args.max_conformers # default 50

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

        conformer_molecule.electronic_energy = el_energy
        conformer_molecule.atoms = [atom[0] for atom in coordinates_list]
        conformer_molecule.coordinates = [atom[1:] for atom in coordinates_list]
        conformer_molecule.program = 'CREST'
        
        # Workflow management
        conformer_molecule.workflow = conformer_molecule.determine_workflow()
        conformer_molecule.set_current_step()
        
        conformer_molecules.append(conformer_molecule)

    filtered_molecules = energy_cutoff(conformer_molecules)

    return filtered_molecules[:max_conformers]

def collect_DFT_and_DLPNO(molecules):
    collected_molecules = []
    
    for m in molecules:
        process_conditions = (m.current_step == 'optimization' and (m.reactant or m.product)) or ('TS_opt' in m.current_step)
        
        if process_conditions:
            if 'OH' in m.name or 'H2O' in m.name:
                identifier = 'OH' if 'OH' in m.name else 'H2O'
            else:
                match = re.search(r'conf(\d+)', m.name)
                identifier = match.group() if match else None

            if identifier:
                # Adjusted pattern to ensure exact match
                pattern = re.compile(re.escape(identifier) + r'(?!\d)')
                for m_DLPNO in molecules:
                    if pattern.search(m_DLPNO.name) and m_DLPNO.current_step == 'DLPNO':
                        m.electronic_energy = m_DLPNO.electronic_energy
                        m.update_step()
                        m.name = m.name.replace('TS', 'DLPNO')
                        collected_molecules.append(m)
                        break

    return collected_molecules


def handle_input_molecules(molecules, logger, threads):
    current_step = molecules[0].current_step
    if all(m.current_step == current_step for m in molecules):
        logger.log(f"\nDetected molecules with current step: {current_step}")
        logger.log("Checking convergence status of molecules")
        for m in molecules:
            m.directory = start_dir
            if os.path.exists(m.log_file_path):
                converge_status, _ = termination_status(m, logger)
                if current_step in ['Done', 'DLPNO'] and converge_status is True:
                    global_molecules.append(m)
                elif converge_status is True:
                    m.converged = True
                else:
                    m.converged = False
            else: # in the case a single .pkl file has been given as input, and logfiles are not located in same directory:
                if len(molecules) == 1:
                    m.converged = False
                # Check if the conformers in the pkl file have atoms and coordinates and send for calculation.
                elif m.atoms and m.coordinates:
                    m.converged = True
        if all(m.converged for m in molecules):
            logger.log(f"All molecules converged for step {current_step}. Proceeding with next step: {molecules[0].next_step}")
            molecules[0].move_files()
            handle_termination(molecules, logger, threads, converged=True)
        else:
            logger.log(f"Not all molecules given have converged. Non-converged will be calculated and converged ones will be skipped.")
            handle_termination(molecules, logger, threads, converged=False)
    else:
        logger.log(f"Not all molecules in the given files are of the same type. Check if correct input is provided. Exiting for now.")
        exit()


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
        

def eckart_test(TS_E, reactant_E, product_E, imag, reduced_mass=1, temperature=298.15):
    try:
        import numpy as np

        # Define physical constants
        c = 2.99792458e10
        AMU_TO_KG = 1.66053886E-27
        PI = 3.14159265359
        reduced_mass_amu = 1

        def barrier(transition_state_energy, reactant_energy, product_energy):
            return (transition_state_energy ** 0.5 + (transition_state_energy - (product_energy - reactant_energy)) ** 0.5) ** 2
        
        def alpha(barrier_parameter, force_constant, transition_state_energy, reactant_energy, product_energy):
            return (barrier_parameter * force_constant / (2 * transition_state_energy * (transition_state_energy - (product_energy - reactant_energy)))) ** 0.5
        
        def A(energy, reduced_mass_kg, alpha):
            return 2 * PI * (2 * reduced_mass_kg * energy)**0.5 / alpha
        
        def B(energy, reduced_mass_kg, product_energy, reactant_energy, alpha):
            return 2 * PI * (2 * reduced_mass_kg * (energy - (product_energy - reactant_energy)))**0.5 / alpha
        
        def delta(barrier_parameter, reduced_mass_kg, alpha):
            return 2 * PI * abs((2 * reduced_mass_kg * barrier_parameter - (alpha/2)**2))**0.5 / alpha
        
        def P(a_argument, b_argument, d_argument):
            return (np.cosh(a_argument + b_argument) - np.cosh(a_argument - b_argument)) / (np.cosh(a_argument + b_argument) + np.cosh(d_argument))

        def force_constant(imag):
            return (imag*2*np.pi*c)**2
        
        # Convert mass to kg and calculate parameters
        reduced_mass_kg = reduced_mass_amu * AMU_TO_KG
        force_constant = force_constant(imag)

        barrier_parameter = barrier(TS_E, reactant_E, product_E)
        alpha = alpha(barrier_parameter, force_constant, TS_E, reactant_E, product_E)

        # Energy range for integration
        energies = np.linspace(0, TS_E, 1000)
        transmission_probabilities = np.zeros_like(energies)

        # Calculate transmission probability for each energy
        for i, energy in enumerate(energies):
            a_argument = A(energy, reduced_mass_kg, alpha)
            b_argument = B(energy, reduced_mass_kg, product_E, reactant_E, alpha)
            d_argument = delta(barrier_parameter, reduced_mass_kg, alpha)
            transmission_probabilities[i] = P(a_argument, b_argument, d_argument)

        # Integrate transmission probability over the energy range using the trapezoidal rule
        tunneling_correction_factor = np.trapz(transmission_probabilities, energies)

        return tunneling_correction_factor
    except Exception as e:
        print("Error in calculating the Eckart tunneling. Returning tunneling coefficient 1")
        return 1


def gcalc_optimized(emin, gsize, v, a, b, l, h, kb, m, t, v1):
    from numpy import pi, arange, sqrt, cosh, exp, zeros, trapz
    pi = pi
    e = arange(emin, 2 * max(v), gsize)
    k = zeros(len(e))
    ek = e.copy()
    eka = e * 627.509

    c = (h ** 2) / (8 * m * (l ** 2)) 
    for i, ei in enumerate(e):
        a = 0.5 * sqrt(ei / c)
        b = 0.5 * sqrt((ei - a) / c)
        d = 0.5 * sqrt((b - c) / c)
        k[i] = 1 - ((cosh(2 * pi * (a - b)) + cosh(2 * pi * d)) / (cosh(2 * pi * (a + b)) + cosh(2 * pi * d)))

    g = zeros(len(t))
    for q, ti in enumerate(t):
        gk = k * exp(-ek / (kb * ti))
        gi = trapz(gk, ek)
        gi = gi * (exp(v1 / (kb * ti)) / (kb * ti))

        # For energies above the interval where transmission is less than 1, use
        # analytical integration to obtain the final result
        g[q] = gi + exp(v1 / (kb * ti)) * exp(-ek[-1] / (kb * ti))

    return g, eka, k, gk


def eckart_optimized(sp_ts, sp_reactant, sp_product, imag, t=[298.15]):
    from numpy import pi, sqrt, arange, exp
    try:
        # Constants
        c = 2.99792458e+8  # speed of light (m s-1)
        pi = 3.141592653589
        kb = 3.1668152e-6
        h = 2 * pi
        na = 6.0221409e+23  # Avogadro's number (mol-1)
        mu = 1  # Assuming this is a constant as it's not defined how it changes

        # Energy conversions
        e1 = (sp_ts - sp_reactant) * 4184 / na / 4.3597447222071e-18
        e2 = (sp_ts - sp_product) * 4184 / na / 4.3597447222071e-18
        wau = imag * 100 * c * 2.418884326509e-17
        m = mu * 1822.888479

        # Force constants and reaction coordinate calculations
        f = -4 * (pi ** 2) * (wau ** 2) * m
        a = e1 - e2
        b = (sqrt(e2) + sqrt(e1)) ** 2
        l = -pi * (a - b) * (b + a) / (sqrt(-2 * f * b) * b)

        # Defining reaction coordinate
        x = arange(-3, 3, 0.01) / sqrt(mu)  # Adjusted for reduced mass
        
        # Calculating potential
        y = -exp((2 * pi * x) / l)
        v = (-(y * a) / (1 - y)) - ((y * b) / ((1 - y) ** 2))
        xa = 0.529177 * x * sqrt(mu)  # Reduced mass re-inserted for plotting
        va = v * 627.509  # Potential converted to kcal/mol for plotting

        # Calculation for tunneling correction factors
        vb = [v[0], v[-1]]
        emin = max(vb)
        gsize = max(v) / 50
        gold, eka, k, gk = gcalc_optimized(emin, gsize, v, a, b, l, h, kb, m, t, e1)
        gdiff = 1
        while gdiff >= 0.001:
            gsize /= 10
            gnew, eka, k, gk = gcalc_optimized(emin, gsize, v, a, b, l, h, kb, m, t, e1)
            gdiff = max(abs(gnew - gold) / gold)
            gold = gnew

        kappa = gold[0]
        return kappa
    except Exception as e:
        print("Error in calculating the Eckart tunneling. Returning tunneling coefficient 1")
        return 1

    
def gcalc(emin,gsize,v,a,b,l,h,kb,m,t,v1):
    from numpy import pi, arange, sqrt, cosh, exp
    pi = pi
    e=arange(emin,2*max(v),gsize)
    k=[0 for x in range(len(e))]
    ek=[0 for x in range(len(e))]
    eka=[0 for x in range(len(e))]
    for i in range(len(e)):
        c=(h**2)/(8*m*(l**2))     # calculation of c factor
        a=0.5*sqrt(e[i]/c)
        b=0.5*sqrt((e[i]-a)/c)
        d=0.5*sqrt((b-c)/c)
        k[i]=1-((cosh(2*pi*(a-b))+cosh(2*pi*d))/(cosh(2*pi*(a+b))+cosh(2*pi*d)))
        ek[i]=e[i]
        eka[i]=e[i]*627.509
    
    g=[0 for x in range(len(t))]
    for q in range(len(t)):                                # temperature dependent operations begin
        gk=[0 for x in range(len(k))]                                              # clean variable
        for j in range(len(k)):                            
            gk[j]=k[j]*exp(-ek[j]/(kb*t[q]))           # calculate integrand
        gi=[0 for x in range(len(k))]
        for l in range(len(ek)-1):
            gi[l]=(0.5*(gk[l]+gk[l+1])*abs(ek[l]-ek[l+1])) # numerical integration by using the area of squares
        
        gi=sum(gi)
        gi=gi*(exp(v1/(kb*t[q]))/(kb*t[q]))                     # multiplying integral with prefactor
        
        # for energies above the interval where transmission is less than 1, we use
        # analytical integration to obtain the final result
        
        g[q]=gi+exp(v1/(kb*t[q]))*exp(-ek[len(ek)-1]/(kb*t[q]))
        
    return g, eka, k, gk


def eckart(sp_ts, sp_reactant, sp_product, imag, t=[298.15]):
    from numpy import pi, sqrt, arange, exp
    try:
        c=2.99792458e+8          # speed of light (m s-1)
        pi=pi               # π
        kb=3.1668152e-6
        h=2*pi
        na=6.0221409e+23         # avogadro's number (mol-1)

        e1 = sp_ts - sp_reactant
        e2 = sp_ts - sp_product
        mu = 1
        v1=((e1*4184)/na)/4.3597447222071e-18
        v2=((e2*4184)/na)/4.3597447222071e-18
        wau=(imag*100)*c*2.418884326509e-17
        m=mu*1822.888479

        # calculate force constant, a, b and l
        f=-4*(pi**2)*(wau**2)*m;
        f2=-4*(pi**2)*(wau**2)*1;
        a=v1-v2;
        b=(sqrt(v2)+sqrt(v1))**2;
        l=-pi*(a-b)*(b+a)/(sqrt(-2*f*b)*b);

        # defining reaction coordinate
        x = arange(-3, 3, 0.01)
        x = x/(sqrt(mu))      # removing reduced mass from reaction coordinate in order to get a well defined potential

        # calculating potential
        y=[0 for i in range(len(x))]
        v=[0 for i in range(len(x))]
        xa=[0 for i in range(len(x))]
        va=[0 for i in range(len(x))]

        for i in range(len(x)):
            y[i]=-exp( (2*pi*x[i])/l )
            v[i]=( (-(y[i]*a)/(1-y[i]) ) - ( (y[i]*b)/((1-y[i])**2)) )
            xa[i]=0.529177*x[i]*sqrt(mu)         # reduced mass re-inserted for plotting
            va[i]=v[i]*627.509                        # potential converted to kcal/mol for plotting

        # calculating the correction factors for all t's
        vb=[0,0]
        vb[0]=v[0]                                                # value of potential at reactants
        vb[1]=v[len(x)-1]                                           # value of potential at products
        emin=max(vb)                                           # minimum energy at which tunnelling can occur
        gdiff=1                                                   # initial convergence control set to 1
        gsize=max(v)/50                                        # initial integration stepsize 
        [gold,eka,k,gk]=gcalc(emin,gsize,v,a,b,l,h,kb,m,t,v1)                      # calculate g
        gsize=gsize/10                                            # reduce integration stepsize
        runs=0                                                    # initial number of runs
        while gdiff >= 0.001:                                        # convergence criteria
            [gnew,eka,k,gk]=gcalc(emin,gsize,v,a,b,l,h,kb,m,t,v1)  # new g
        #    print("tunneling factor", gnew)                                         # display new correction factor
            gsize=gsize/10                                        # reduce integration stepsize
            gdiffcalc=[0 for x in range(len(t))]
            for j in range(len(t)):
                gdiffcalc[j]=abs(gnew[j]-gold[j])/gold[j]         # calculate convergence
            gdiff=max(gdiffcalc)                                  # max convergence control value
        #    print("convergence control value", gdiff)                                        # display convergence control value
            gold=gnew                                             # replace old correction factor with new
            runs=runs+1                                           # a run completed
        #    print("runs done", runs)                                         # display run number

        [g,eka,k,gk]=gcalc(emin,gsize,v,a,b,l,h,kb,m,t,v1)        #final g

        kappa = g[0]
        return kappa
    except Exception as e:
        print("Error in calculating the eckart tunneling. returning tunneling coefficient 1")
        return 1

def rate_constant_test(TS_conformers, reactant_conformers, product_conformers, T=298.15):
    from numpy import  exp, sum, float64, NaN
    k_b = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck constant in J*s
    HtoJ = 43.597447222e-19  # Conversion factor from Hartree to Joules
    Htokcalmol = 627.509
    Na = 6.022e23 # molecules/mol
    liters_to_cm3 = 1000 # 1 liter is 1000 cm^3
    mol_per_liter = 1/(0.0821*T)
    p_ref = (mol_per_liter*Na)/liters_to_cm3
    kappa = 1
    
    # Calculation of rate constant
    if reactant_conformers and TS_conformers:
        # Seperate organic molecule and OH
        reactant_molecules = [mol for mol in reactant_conformers if 'OH' not in mol.name]
        OH = next((mol for mol in reactant_conformers if 'OH' in mol.name))
        
        if product_conformers:
            product_molecules = [mol for mol in product_conformers if 'H2O' not in mol.name]
            H2O = next((mol for mol in product_conformers if 'H2O' in mol.name), None)
        else:
            product_molecules = []

        # imag = abs(lowest_TS.vibrational_frequencies[0])
        # kappa = eckart(lowest_EE_TS_kcalmol, lowest_EE_reactant_kcalmol, lowest_EE_product_kcalmol, imag)

        # Pair up reactants, TS, and products
        TS_reactant_pairs, TS_products_pairs = ArbAlign_pair(TS_conformers, reactant_molecules, product_molecules)

        sum_TS = sum([exp(-TS.free_energy*HtoJ/(k_b*T)) for TS, _ in TS_reactant_pairs])
        Q_TS = sum([exp(-(lowest_ZP_TS_J - float64(mol.zero_point_corrected * HtoJ)) / (k_b * T)) * float64(mol.Q) for mol in TS_conformers])
        # sum_reactant = sum([exp(-(reactant.free_energy+OH.free_energy)*HtoJ/(k_b*T)) for _, reactant in TS_reactant_pairs])
        # k = (k_b*T)/(h*p_ref) * sum_TS/sum_reactant
        k = kappa * (k_b*T)/(h*p_ref) *  sum([exp(-(mol_TS.free_energy*HtoJ - (mol_reac.free_energy+OH.free_energy)*HtoJ)/(k_b*T)) for mol_TS, mol_reac in TS_reactant_pairs])

        return k, kappa


def rate_constant(TS_conformers, reactant_conformers, product_conformers, T=298.15):
    from numpy import  exp, sum, float64, NaN
    k_b = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck constant in J*s
    HtoJ = 43.597447222e-19  # Conversion factor from Hartree to Joules
    Htokcalmol = 627.509
    Na = 6.022e23 # molecules/mol
    liters_to_cm3 = 1000 # 1 liter is 1000 cm^3
    mol_per_liter = 1/(0.0821*T)
    p_ref = (mol_per_liter*Na)/liters_to_cm3
    kappa = 1
    
    # Calculation of rate constant
    if reactant_conformers and TS_conformers:
        # Seperate organic molecule and OH
        reactant_molecules = [mol for mol in reactant_conformers if 'OH' not in mol.name]
        OH = next((mol for mol in reactant_conformers if 'OH' in mol.name))

        # Find the lowest single point energies for reactant (excluding OH) and TS
        lowest_reactant = min(reactant_molecules, key=lambda molecule: molecule.zero_point_corrected)
        lowest_TS = min(TS_conformers, key=lambda molecule: molecule.zero_point_corrected)

        # Lowest single point reactant and TS and convert from Hartree to joules and kcal/mol
        lowest_ZP_TS_J = float64(lowest_TS.zero_point_corrected * HtoJ)
        lowest_ZP_reactant_J = float64(lowest_reactant.zero_point_corrected * HtoJ)
        lowest_ZP_TS_kcalmol = float64(lowest_TS.zero_point_corrected * Htokcalmol)
        lowest_ZP_reactant_kcalmol = float64((lowest_reactant.zero_point_corrected + OH.zero_point_corrected) * Htokcalmol)

        if OH:
            sum_reactant_ZP_J = lowest_ZP_reactant_J + (float64(OH.zero_point_corrected * HtoJ))
            Q_reactant = sum([exp(-(lowest_ZP_reactant_J - float64(mol.zero_point_corrected * HtoJ)) / (k_b * T)) * float64(mol.Q) for mol in reactant_molecules]) * float64(OH.Q)
        else:
            sum_reactant_ZP_J = lowest_ZP_reactant_J
            Q_reactant = sum([exp(-(lowest_ZP_reactant_J - float64(mol.zero_point_corrected * HtoJ)) / (k_b * T)) * float64(mol.Q) for mol in reactant_molecules])        

        Q_TS = sum([exp(-(lowest_ZP_TS_J - float64(mol.zero_point_corrected * HtoJ)) / (k_b * T)) * float64(mol.Q) for mol in TS_conformers])

        if product_conformers:
            # Seperate organic molecule and H2O
            product_molecules = [mol for mol in product_conformers if 'H2O' not in mol.name]
            H2O = next((mol for mol in product_conformers if 'H2O' in mol.name), None)

            if product_molecules and H2O:
                lowest_product = min(product_molecules, key=lambda molecule: molecule.electronic_energy)
                lowest_EE_product_kcalmol = float64((lowest_product.electronic_energy + H2O.electronic_energy) * Htokcalmol)
                lowest_EE_TS_kcalmol = float64(lowest_TS.electronic_energy * Htokcalmol)
                lowest_EE_reactant_kcalmol = float64((lowest_reactant.electronic_energy + OH.electronic_energy) * Htokcalmol)

                imag = abs(lowest_TS.vibrational_frequencies[0])
                kappa = eckart(lowest_EE_TS_kcalmol, lowest_EE_reactant_kcalmol, lowest_EE_product_kcalmol, imag)
                k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))
            else:
                print("Error in product molecules")
                k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))
        else: # If products are not calculated assume tunneling coefficient is 1
            k = kappa * (k_b*T)/(h*p_ref) * (Q_TS/Q_reactant) * exp(-(lowest_ZP_TS_J - sum_reactant_ZP_J) / (k_b * T))

        print(f"Ea: {lowest_ZP_TS_kcalmol - lowest_ZP_reactant_kcalmol:.4f} kcal/mol  Q_TS: {Q_TS}  k: {k}")
        return k, kappa # cm^3 molecules^-1 s^-1

    return NaN, NaN


def main():
    last_run_args = load_last_run_args()
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
    parser.add_argument('-transition_states', type=str2bool, default=True, metavar='<boolean>', help=argparse.SUPPRESS)
    # parser.add_argument('-NEB', type=str2bool, default=False, metavar='<boolean>', help='Prepare input file for Nudged Elastic Band')
    parser.add_argument('-auto', type=str2bool, default=True, metavar='<boolean>', help='Automated process with the following workflow:\n- CREST conformer sampling of xyz input file (GFN2-xTB -ewin=10 kcal/mol)\n- Preoptimization of geometry (Low Level of Theory)\n- Optimization towards transition state (High Level of Theory)\n- DLPNO-CCSD(T) SP energy calculations on top of TS conformers\n- Calculate rate constants and branching ratios for reaction type\n* Automatically resubmit failed calculations until convergence or wall-time limit is reached')

    # Argparse reactions
    reaction_options = parser.add_argument_group("Types of reactions")
    reaction_options.add_argument('-OH', action='store_true', help='Perform H abstraction with OH radical')
    reaction_options.add_argument('-Cl', action='store_true', help='Perform H abstraction with Cl radical')
    reaction_options.add_argument('-CC', action='store_true', help='(TBA) Perform addition to C=C bonds')
    reaction_options.add_argument('-OH_CC', action='store_true', help='(TBA) Perform OH addition to C=C bonds')

    # Argparse additional
    additional_options = parser.add_argument_group("Additional arguments")
    # additional_options.add_argument('-restart', action='store_true', default=False, help='Restart and resume log files')
    additional_options.add_argument('-k', type=str2bool, metavar='<boolean>', default=True, help='Calculate Multiconformer Transition State rate constant [def: True]')
    additional_options.add_argument('-restart', action='store_true', default=False, help='Need: .log, .out, .pkl - Restart from molecule current step')
    additional_options.add_argument('-smiles', metavar='string', type=str, help='Input molecule as a SMILES string')
    additional_options.add_argument('-info', action='store_true', default=False, help='Print information of molecules in log files or .pkl file')
    additional_options.add_argument('-IRC', action='store_true', default=False, help='Perform IRC calcuation on ORCA or G16 output file')
    additional_options.add_argument('-CHO', dest='CHO', action=ParseList, nargs='*', help="Set indexes of atoms for active site. Indexing starting from 1")
    additional_options.add_argument('-collect', action='store_true', default=False, help='Collect thermochemical data from TS structures and single point correction from DLPNO')
    additional_options.add_argument('-method', type=str, default='wB97XD', help='Specify QC method to use for optimization and TS search [def: WB97X-D3BJ]')
    additional_options.add_argument('-basis_set', type=str, default='6-31++g(d,p)', help='Specify basis set to use with QC method [def: 6-31++g(d,p)]')
    # additional_options.add_argument('--low_level', nargs='+', metavar='', help='Specify low-level theory for preoptimization [def: B3LYP 6-31+G(d,p)]')
    additional_options.add_argument('--gfn', default='2', choices=['1','2'], help='Specify the GFN version (1 or 2, default: 2)')
    additional_options.add_argument('-skip_preopt', action='store_true', default=False, help='Skip the preoptimization of the structures before the TS search')
    additional_options.add_argument('-filter_step', dest='filter_step', default=['DLPNO'], action=ParseList, nargs='*', help="Steps at which to perform filtering of conformers")
    additional_options.add_argument('-cpu', metavar="int", nargs='?', const=1, type=int, default=4, help='Number of CPUs [def: 4]')
    additional_options.add_argument('-mem', metavar="int", nargs='?', const=1, type=int, default=8000, help='Amount of memory allocated for the job [def: 8000MB]')
    additional_options.add_argument('-par', metavar="partition", type=str, default="q64,q48,q28,q24", help='Partition to use [def: q64,q48,q28,q24]')
    additional_options.add_argument('-time', metavar="hh:mm:ss", type=str, default=None, help='Monitoring duration [def: 144 hours]')
    additional_options.add_argument('-interval', metavar="int", nargs='?', const=1, type=int, help='Time interval between log file checks [def: based on molecule size]')
    additional_options.add_argument('-initial_delay', metavar="int", nargs='?', const=1, type=int, help='Initial delay before checking log files [def: based on molecule size]')
    additional_options.add_argument('-attempts', metavar="int", nargs='?', const=1, type=int, default=100, help='Number of log file check attempts [def: 100]')
    additional_options.add_argument('-max_conformers', metavar="int", nargs='?', const=1, type=int, default=50, help='Maximum number of conformers from CREST [def: 50]')
    additional_options.add_argument('-freq_cutoff', metavar="int", nargs='?', const=1, type=int, default=-100, help='TS imaginary frequency cutoff [def: -100 cm^-1]')
    additional_options.add_argument('-ewin', metavar="int", nargs='?', const=1, default=10, type=int, help='Energy threshold for CREST conformer sampling [def: 10 kcal/mol]')
    additional_options.add_argument('-energy_cutoff', metavar="digit", nargs='?', const=1, default=10, type=float, help='After preoptimization, remove conformers which are [int] kcal/mol higher in energy than the lowest conformer [def: 10 kcal/mol]')
    additional_options.add_argument('-pickle', action='store_true', default=False, help='Store given log files into pickle file')
    additional_options.add_argument('-filter', type=str2bool, metavar='<boolean>', default=True, help='Filter identical conformers after transition state optimization [def: True]')
    additional_options.add_argument('--account', type=str, help=argparse.SUPPRESS) # For Finland Puhti
    additional_options.add_argument('-test', action='store_true', default=False, help=argparse.SUPPRESS) # Run TEST section in main() and exit
    additional_options.add_argument('-num_molecules', type=int, default=None, help=argparse.SUPPRESS) # Set the number of directories created to [int]. Used when doing limited testing
    additional_options.add_argument('-XQC', '-YQC', '-QC', action=SCFAction, nargs='*', const=[], default='SCF=(XQC)', dest='SCF', help='Use G16 SCF=(X,Y)QC algorithm for SCF procedure - https://gaussian.com/scf/')
    parser.set_defaults(SCF="")


    hidden_options = parser.add_argument_group()
    hidden_options.add_argument('-init', action='store_true', help=argparse.SUPPRESS)
    hidden_options.add_argument('-loc', action='store_true', help=argparse.SUPPRESS)
    hidden_options.add_argument('-reaction_angle', metavar="float", nargs='?', default=169.5, const=1, type=float, help=argparse.SUPPRESS) 
    hidden_options.add_argument('-hybrid', action='store_true', default=False, help=argparse.SUPPRESS) #'Hybrid Hessian. Full Hessian for active site and approximate for rest of atoms'
    hidden_options.add_argument('-dispersion', action='store_true', default=False, help=argparse.SUPPRESS)
    additional_options.add_argument('-skip', type=str2bool, default=False, help=argparse.SUPPRESS)

    global args, start_dir
    args = parser.parse_args()
    start_dir = os.getcwd()

    ################################ARGUMENT SPECIFICATIONS############################

    # Program and method 
    global method, basis_set

    methods_no_basis = {"b97-3c", "r2scan-3c", "pm3", "am1", "pm6", "pm7", 'g3mp2', 'g3'}
    if args.method.lower() in methods_no_basis:
        args.basis_set = ""

    global QC_program
    if args.ORCA:
        args.G16 = False
        QC_program = "ORCA"
        if args.method.lower() == "wb97xd":
            args.method = "WB97X-D3BJ"
    elif args.G16:
        args.ORCA = False
        QC_program = "G16"
        if args.method.lower() == "wb97xd":
            args.method = "wb97xd"

    global termination_strings, error_strings  # Keep in lower case!!
    termination_strings = {
    "g16": ["normal termination"],
    "orca": ["****orca terminated normally****"],
    "crest": ["crest done", "crest terminated normally."]  # Multiple possible termination messages
    }
    error_strings = {
    "g16": ["error termination", "another termination example"],
    "orca": ["aborting the run", "this wavefunction is not fully converged", "not enough memory", "orca finished by error termination", "error", "this wavefunction is not converged"],
    "crest": ["some crest error message"]  # Multiple possible error messages
    }

    #####################################################################################################
    if args.test:
        from ArbAlign import compare
        molecules = []
        TS = []
        reactants = []
        products = []
        threads = []
        logger = Logger(os.path.join(start_dir, "test.txt"))

        for n, input_file in enumerate(args.input_files, start=1):
            file_name, file_type = os.path.splitext(input_file)
            if file_type == '.pkl':
                from pandas import DataFrame
                pandas_molecules = Molecule.load_molecules_from_pickle(input_file)
                if isinstance(pandas_molecules, DataFrame):
                    conformers = initiate_conformers(input_file) 
                    for conf in conformers:
                        molecules.append(conf)
                else:
                    if isinstance(pandas_molecules, list):
                        for molecule in pandas_molecules:
                            if not 'H2O' in molecule.name:
                                if not 'OH' in molecule.name:
                                        molecules.append(molecule)
            elif file_type in ['.log', '.out', '.com', '.inp', '.xyz']:
                QC_program = 'ORCA' if file_type in ['.out', '.inp'] else 'G16'
                input_file_path = os.path.join(start_dir, input_file)
                molecule = Molecule(input_file_path, indexes=args.CHO, program=QC_program)
                molecule.name = file_name
                molecules.append(molecule)
            else:
                molecule = Molecule(input_file, indexes=args.CHO)
                molecules.append(molecule)
        # 

        # for n, molecule in enumerate(molecules,start=1):
        #     check_transition_state(molecule, logger)
            # molecule.perturb_active_site()
            # molecule.write_xyz_file(os.path.join(os.getcwd(), f"perturb{n}.xyz"))
            # print(molecule.name, check_transition_state(molecule,logger))
            # if 'CREST' in molecule.name:
            #     TS.append(molecule)
            # else:
            #     products.append(molecule)
        # #
        # ArbAlign_pair(TS, reactants=None, products=products)
        # for i in TS:
        #     for j in products:
        #         rmsd = compare(i,j)
        #         print(i.name, j.name, rmsd)
        new_list = ArbAlign_compare_molecules(molecules, pickle=True, RMSD_threshold=0.34)
        # energy_cutoff(new_list) 




        exit()
    ####################################################################################################

    threads = []
    input_molecules = []
    final_TS = []
    final_reactants = []
    final_products = []
    molecules = []

    if args.init: # only executes if input file is an xyz file
        if args.smiles:
            input_molecule = Molecule(smiles=args.smiles, reactant=True)
            file_name, file_type = input_molecule.name, 'xyz'
        else:
            input_file = args.input_files[0]
            file_name, file_type = os.path.splitext(input_file)
            input_file_path = os.path.join(start_dir, input_file)
            input_molecule = Molecule(input_file_path, reactant=True) # Reuse input molecule as template for reactant

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

        if args.transition_states:
            for count, molecule in enumerate(reacted_molecules, start=1):
                molecule.name = f"{file_name}_H{count}" 
                molecule.directory = os.path.join(start_dir, molecule.name)
                mkdir(molecule)
                logger = Logger(os.path.join(molecule.directory, "log"))
                molecule.write_xyz_file(os.path.join(molecule.directory, f"{molecule.name}.xyz")) 
                molecule.save_to_pickle(os.path.join(molecule.directory, f"{molecule.name}.pkl"))

        if args.reactants: 
            reactant_dir = os.path.join(start_dir, 'reactants')
            input_molecule.name = f"{file_name}_reactant"
            input_molecule.directory = reactant_dir
            input_molecule.set_current_step('crest_sampling')
            mkdir(input_molecule, crest_constrain=False)
            input_molecule.save_to_pickle(os.path.join(input_molecule.directory, f"{input_molecule.name}.pkl"))

        if args.products and product_molecules:
            product_dir = os.path.join(start_dir, 'products')
            for count, molecule in enumerate(product_molecules, start=1):
                molecule.name = f"{file_name}_product_H{count}"
                molecule.directory = product_dir
            mkdir(product_molecules[0], crest_constrain=False)
            pickle_path = os.path.join(product_dir, "product_molecules.pkl")
            Molecule.molecules_to_pickle(product_molecules, pickle_path)

    elif args.info:
        if args.smiles:
            smiles_molecule = Molecule(smiles=args.smiles, reactant=True)
            smiles_molecule.print_items()
        for input_file in args.input_files[:args.max_conformers]:
            file_name, file_type = os.path.splitext(input_file)
            if file_type == '.pkl':
                from pandas import DataFrame
                molecules = Molecule.load_molecules_from_pickle(input_file)
                if isinstance(molecules, DataFrame):
                    conformers = initiate_conformers(input_file) 
                    for conf in conformers:
                        conf.print_items()
                else:
                    if isinstance(molecules, list):
                        molecules = energy_cutoff(molecules)
                        for molecule in molecules[:args.max_conformers]:
                            molecule.print_items()
                    else: molecules.print_items()
            elif file_type in ['.log', '.out', '.com', '.inp', '.xyz']:
                QC_program = 'ORCA' if file_type in ['.out', '.inp'] else 'G16'
                input_file_path = os.path.join(start_dir, input_file)
                molecule = Molecule(input_file_path, indexes=args.CHO, program=QC_program)
                molecule.print_items()
            else:
                parser.error("Invalid input file format. Without a specified reaction type the program expects input files with extensions '.pkl', '.log', '.out', '.com', or '.inp'\nPlease ensure that your input file is in one of these formats and try again. If you provided a different type of file, convert it to a supported format or select an appropriate file for processing.")

    else: # We loop over all molecules given in the argument input and process them according to file type
        if args.smiles:
            input_molecule = Molecule(smiles=args.smiles, reactant=True)
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

                handle_termination(molecules, logger, threads, converged=False)

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
                molecule = Molecule(file_path=input_file_path, indexes=args.CHO)
                input_molecules.append(molecule)

            else:
                print(f"File type {file_type} is not supported")
                exit()

        if final_TS and final_reactants: # In the case final_products is empty symmetrical eckart tunneling used
            k, kappa = rate_constant(final_TS, final_reactants, final_products)
            # print(f"{k} cm^3 molecules^-1 s^-1")
            # print(f"Tunneling coefficient: {kappa}")
            exit()

        # FIX file_type variable
        if input_molecules and file_type != '.xyz':
            logger = Logger(os.path.join(start_dir, "log"))
            if args.collect:
                collected_molecules = collect_DFT_and_DLPNO(input_molecules)
                Molecule.molecules_to_pickle(collected_molecules, os.path.join(start_dir, "collected_molecules.pkl"))
            elif args.IRC:
                IRC_submit(input_molecules)
            elif args.restart:
                handle_input_molecules(input_molecules, logger, threads)
            elif args.pickle:
                filename = re.sub("_conf\d{1,2}", "",input_molecules[0].name)
                Molecule.molecules_to_pickle(input_molecules, os.path.join(start_dir, f"collection{filename}.pkl"))
            else:
                parser.error("Error")
                    
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

        if all(m.reactant for m in global_molecules):
            molecule_name = global_molecules[0].name.split("_")[0]
            logger.log(f"Final DLPNO calculations for reactants is done. Logging molecules to Final_reactants_{molecule_name}.pkl")
            Molecule.molecules_to_pickle(global_molecules, os.path.join(start_dir, f"Final_reactants_{molecule_name}.pkl"))
        elif all(m.product for m in global_molecules):
            # Identify unique H numbers
            h_numbers = sorted(set(re.search(r'_H(\d+)[._]', m.name + '_').group(1) for m in molecules if "_H" in m.name))
            # Group molecules by H numbers and include H2O in each group
            grouped_lists = [[m for m in global_molecules if f"_H{h_num}_" in m.name] + 
                             [m for m in global_molecules if "H2O" in m.name] 
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
            Molecule.molecules_to_pickle(global_molecules, TS_pkl_path)

            if args.k:
                reactant_pkl_name = os.path.basename(start_dir).split("_")[0]
                product_pkl_name = os.path.basename(start_dir)
                reactant_pkl_path = os.path.join(os.path.dirname(start_dir), f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                product_pkl_path = os.path.join(os.path.dirname(start_dir), f'products/Final_products_{product_pkl_name}.pkl')
                logger.log(f"{product_pkl_path}, {reactant_pkl_path}")
                if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                    final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                    final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                    k, kappa = rate_constant(global_molecules, final_reactants, final_products)
                    results_logger = Logger(os.path.join(os.path.dirname(start_dir), "Rate_constants.txt"))
                    results_logger.log_with_stars(f"{molecule_name}: {k} molecules cm^-3 s^-1 with tunneling coefficient {kappa}")
                else:
                    reactant_pkl_path = os.path.join(start_dir, f'reactants/Final_reactants_{reactant_pkl_name}.pkl')
                    product_pkl_path = os.path.join(start_dir, f'products/Final_products_{product_pkl_name}.pkl')
                    logger.log(f"{product_pkl_path}, {reactant_pkl_path}")
                    if os.path.exists(reactant_pkl_path) and os.path.exists(product_pkl_path):
                        final_reactants = Molecule.load_molecules_from_pickle(reactant_pkl_path)
                        final_products = Molecule.load_molecules_from_pickle(product_pkl_path)
                        k, kappa = rate_constant(global_molecules, final_reactants, final_products)
                        results_logger = Logger(os.path.join(start_dir, "results_log"))
                        results_logger.log_with_stars(f"1 {molecule_name}: {k} molecules cm^-3 s^-1 with tunneling coefficient {kappa}")
                    else:
                        logger.log(f"Could not find pickle files in path: {product_pkl_name} and {reactant_pkl_path}")


if __name__ == "__main__":
    main()
