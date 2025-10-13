import numpy as np
from copy import deepcopy

from .QC import extract_normal_coordinates
from .QC import QC_input
from .config import config

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

    
    
def good_active_site(molecule, aldehyde=False):
    if aldehyde:
        CH_threshold_low, CH_threshold_high = 1.1, 1.4
        HO_threshold_low, HO_threshold_high = 1.5, 2.0
        CO_threshold_low, CO_threshold_high = 2.4, 3.2
        angle_threshold_low, angle_threshold_high = 135, 170
    else:
        CH_threshold_low, CH_threshold_high = 1.1, 1.4
        HO_threshold_low, HO_threshold_high = 1.2, 1.6
        CO_threshold_low, CO_threshold_high = 2.2, 2.8
        angle_threshold_low, angle_threshold_high = 135, 179


    C_index = molecule.constrained_indexes['C']-1
    H_index = molecule.constrained_indexes['H']-1
    O_index = molecule.constrained_indexes['O']-1

    distance_CH = molecule.atom_distance(molecule.coordinates[C_index], molecule.coordinates[H_index])
    distance_HO = molecule.atom_distance(molecule.coordinates[H_index], molecule.coordinates[O_index])
    distance_CO = molecule.atom_distance(molecule.coordinates[C_index], molecule.coordinates[O_index])
    angle_CHO = molecule.calculate_angle(molecule.coordinates[C_index], molecule.coordinates[H_index], molecule.coordinates[O_index])
    print(distance_CH, distance_HO, distance_CO, angle_CHO)

    if (CH_threshold_low <= distance_CH <= CH_threshold_high and
        HO_threshold_low <= distance_HO <= HO_threshold_high and
        CO_threshold_low <= distance_CO <= CO_threshold_high and
        angle_threshold_low <= angle_CHO <= angle_threshold_high):
        return True
    return False

    
def check_transition_state(molecule, threshold=0.5, freq_cutoff=80, args=None):
    # Allow callers to pass args explicitly or rely on the shared config
    args = args or config.args
    from numpy import array
    from numpy.linalg import norm
    aldehyde_TS = False
    msg = 'error'
    try:
        H_index = molecule.constrained_indexes['H']-1 # python 0-based indexing
        C_index = molecule.constrained_indexes['C']-1
        O_index = molecule.constrained_indexes['O']-1
    except Exception:
        molecule.find_active_site()
        H_index = molecule.constrained_indexes['H']-1 # python 0-based indexing
        C_index = molecule.constrained_indexes['C']-1
        O_index = molecule.constrained_indexes['O']-1

    methyl_C_indexes, aldehyde_groups, ketone_methyl_groups = molecule.identify_functional_groups()
    for aldehyde in aldehyde_groups:
        if aldehyde['C'] == C_index and aldehyde['H'] == H_index:
            aldehyde_TS = True
            threshold = 0.3

    good_TS = good_active_site(molecule, aldehyde=aldehyde_TS)

    sorted_freqs = sorted((freq for freq in molecule.vibrational_frequencies))
    if sorted_freqs:
        imag = sorted_freqs[0]
        if imag > freq_cutoff:
            if molecule.error_termination_count == 0 and not good_TS:
                msg = f"No frequency under threshold found for {molecule.name}. Trying to correct geometry and resubmit."
                molecule.set_active_site(indexes=args.CHO)
            else:
                msg = f"No frequency under threshold found for {molecule.name}"
            return False, msg

    if molecule.program.lower() == 'g16':
        normal_coords = extract_normal_coordinates(molecule)
        
        if len(normal_coords) != len(molecule.coordinates):
            print("Error: The number of normal mode displacements does not match the number of atoms in the molecule.")
            return False, ""

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
        # print(f"{molecule.name}  {'!GOOD Geometry!' if good_TS else '':<2}   Over threshold:  {any(value > threshold for value in [combined_values_plus, combined_values_minus, combined_values_plus_minus, combined_values_minus_plus])}     imag: {imag:.2f}")


        if any(value > threshold for value in [combined_values_plus, combined_values_minus, combined_values_plus_minus, combined_values_minus_plus]) and good_TS:
            msg = f"Yay! Normal mode analysis indicate correct TS for {molecule.name} with imaginary frequency: {imag}"
            return True, msg

        else:
            if imag < freq_cutoff:
                msg = f"Change in bond length for {molecule.name} does not meet threshold."
            else:
                msg = "The change in bond length does not meet the threshold for significant shortening or elongation and neither did the magnitude of imaginary frequency."
            return False, msg

    elif molecule.program.lower() == 'orca':
        with open(molecule.log_file_path, 'r') as f:
            content = f.read()
        if normal_mode_displacement_significant(content) and not good_active_site(molecule):
            msg = f"Yay! Normal mode analysis indicate correct TS for {molecule.name} with imaginary frequency: {imag}"
            return True, msg
        else:
            if imag < freq_cutoff:
                msg = f"Change in bond length for {molecule.name} does not meet threshold."
            else:
                msg = "The change in bond length does not meet the threshold for significant shortening or elongation and neither did the magnitude of imaginary frequency."
            return False, msg

            
def TS_search(molecule, indexes=None):
    increment = 0.1
    min_CH, max_CH = 1.2, 1.6
    min_HO, max_HO = 1.1, 1.6
    min_angle, max_angle = 170, 179
    i = 1
    for distance_CH in np.arange(min_CH, max_CH, increment):  
        for distance_HO in np.arange(min_HO, max_HO, increment): 
            for reaction_angle in range(max_angle, min_angle-1, -5):
                m_copy = deepcopy(molecule)
                m_copy.name += f"_{i}"
                i += 1
                m_copy.set_active_site(indexes=indexes, distance_CH=distance_CH, distance_HO=distance_HO, reaction_angle=reaction_angle)
                # Use global config args if available for consistency
                QC_input(m_copy, constrain=False, TS=True, args=config.args)