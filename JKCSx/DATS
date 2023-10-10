#!/usr/bin/env python3


###############################LIBRARIES#####################################
import numpy as np
import sys
import argparse

###########################VECTOR MANIPULATION################################
def calculate_vector(coord1, coord2):
    return np.array(coord2) - np.array(coord1)

def vector_length(vector):
    return np.linalg.norm(vector)

def atom_distance(atom1, atom2):
    return np.linalg.norm(np.array(atom2) - np.array(atom1))

def normalize_vector(vector):
    norm = np.linalg.norm(vector)
    if norm < 1e-8:  # A small threshold to handle very small or zero norms
        return np.zeros_like(vector)  # Return a zero vector if the norm is very small
    return vector / norm

def rotate_vector(vector, axis, angle):
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    cross_product = np.cross(axis, vector)
    rotated_vector = (vector * cos_theta +
                      cross_product * sin_theta +
                      axis * np.dot(axis, vector) * (1 - cos_theta))
    return rotated_vector

#########################################FILES MANIPULATION#############################
# if len(sys.argv) > 6:
#     print("Usage: python3 TS_generation.py -reactant1 structure.xyz C=C(bond length) -reactant2 peroxy.xyz ROO()")
#     sys.exit(1)

# file_in = sys.argv[2]
# file_in_strip_xyz = file_in.split(".")[0]

def read_xyz_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()[2:]
        coords = [[line.split()[0], float(line.split()[1]), float(line.split()[2]), float(line.split()[3])] for line in lines]
    return coords

def write_xyz_file(file_path, updated_coords):
    with open(file_path, 'w') as f:
        f.write(str(len(updated_coords)) + '\n' +'\n')
        for atom in updated_coords:
            f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')


def write_constrain(file_path, C_index, H_index, O_index, OH_index, force_constant=0.25):
    '''Force constant tells how tight to constrain the atoms i.e. the magnitude of the oscillation between the constrained atoms'''
    with open (file_path,"w") as f:
        f.write("$constrain\n")
        f.write(f"  force constant={force_constant}\n") 
        f.write(f"  distance: {C_index}, {H_index}, auto\n")
        f.write(f"  distance: {H_index}, {O_index}, auto\n")
        f.write(f"  distance: {O_index}, {OH_index}, auto\n")
        f.write(f"  angle: {C_index}, {H_index}, {O_index}, auto\n")
        f.write("$end\n")


def find_atom_to_constrain(file, desired_distance=1.25):
    coords = read_xyz_file(file)

    c_index = None
    h_index = None
    o_index = None

    for i, carbon in enumerate(coords):
        if carbon[0] == 'O':
            c_index = i
            for j, hydrogen in enumerate(coords):
                if hydrogen[0] == "H":
                    distance_CH = atom_distance(carbon[1:], hydrogen[1:])

    return c_index, h_index, o_index


######################################################################################
def OH_addition(file, distance=1.45, double_bond_distance=1.36):
    coords = read_xyz_file(file)
    num_atoms = len(coords)
    dist_OH = 0.97
    count = 1

    for i in range(num_atoms):
        if coords[i][0] == "C":
            for j in range(num_atoms):
                if coords[j][0] == "C" and i != j:
                    vector_CC = np.array(coords[i][1:]) - np.array(coords[j][1:])
                    dist_CC = vector_length(vector_CC)
                    if dist_CC <= double_bond_distance:
                        norm_vector_CC = normalize_vector(vector_CC)
                        new_coords = read_xyz_file(file)
                        perpendicular_axis = np.cross(norm_vector_CC, [1,0,0])
                        perpendicular_axis = normalize_vector(perpendicular_axis)
                        # Shift carbon
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + norm_vector_CC * 0.1
                        # Update oxygen in OH coordiantes
                        oxygen_coords = np.array(new_coords[j][1:]) + perpendicular_axis * distance

                        rotation_axis = norm_vector_CC
                        rotation_angle_H = np.radians(45) # 180 - 75.5 = 104.5

                        rotated_vector_H = rotate_vector(perpendicular_axis, rotation_axis, rotation_angle_H)

                        hydrogen_coords = oxygen_coords + rotated_vector_H * dist_OH
                        
                        new_coords.append(['O', *oxygen_coords])
                        new_coords.append(['H', *hydrogen_coords])

                        write_xyz_file(f"{file_in_strip_xyz}_H{count}.xyz", new_coords)
                        count += 1


def OH_abstraction(file, distance=1.35):
    coords = read_xyz_file(file)
    num_atoms = len(coords)
    dist_OH = 0.97
    count = 1

    for i in range(num_atoms):
        if coords[i][0] == "H":  # Atom label for hydrogen
            for j in range(num_atoms):
                if coords[j][0] == "C":  # Atom label for carbon
                    vector_CH = np.array(coords[i][1:]) - np.array(coords[j][1:])
                    dist_CH = vector_length(vector_CH)

                    if dist_CH < distance:
                        norm_vector_CH = normalize_vector(vector_CH)

                        H_perturb_axis = np.cross(norm_vector_CH, [0, 1, 0])
                        H_perturb_axis = normalize_vector(H_perturb_axis)
                        
                        # Update the initial hydrogen coordinates
                        new_coords = read_xyz_file(file)
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + norm_vector_CH * (distance - dist_CH)

                        # Update oxygen in OH coordiantes
                        oxygen_coords = np.array(new_coords[i][1:]) + norm_vector_CH * distance

                        # Coordinates for hydrogen in OH
                        norm_vector_OH = normalize_vector(oxygen_coords - np.array(coords[i][1:]))
                        rotation_axis = np.cross(norm_vector_OH, H_perturb_axis)
                        rotation_axis = normalize_vector(rotation_axis)

                        rotation_angle_H = np.radians(75.5) # 180 - 75.5 = 104.5

                        rotated_vector_H = rotate_vector(norm_vector_CH, rotation_axis, rotation_angle_H)
                        hydrogen_coords = oxygen_coords + rotated_vector_H * dist_OH
                        
                        # Perturb inital hydrogen to avoid linear C--H--O
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + H_perturb_axis * 0.3

                        # Update XYZ file
                        new_coords.append(['O', *oxygen_coords])
                        new_coords.append(['H', *hydrogen_coords])
                        write_xyz_file(f"{file_in_strip_xyz}_H{count}.xyz", new_coords)
                        C_index = j+1
                        H_index = i+1
                        O_index = len(new_coords)-1
                        OH_index = len(new_coords)
                        write_constrain(f"constraints_{count}.inp",C_index, H_index, O_index, OH_index)
                        count += 1
                       



def main():
    parser = argparse.ArgumentParser(description='TS Generation Tool')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-reactant1', help='First reactant XYZ file')
    group.add_argument('-reactant2', help='Second reactant XYZ file')
    group.add_argument('-constrain', action='store_true', help='Find atoms for constraining')

    args = parser.parse_args()

    if args.reactant1:
        file_in = args.reactant1
    elif args.reactant2:
        file_in = args.reactant2

    file_in_strip_xyz = file_in.split(".")[0]

    if args.reactant1 or args.reactant2:
        if 'C=C' in sys.argv:
            OH_addition(file_in)
        elif 'OH' in sys.argv:
            OH_abstraction(file_in)
    elif args.constrain:
        find_atom_to_constrain(file_in)

if __name__ == "__main__":
    main()
