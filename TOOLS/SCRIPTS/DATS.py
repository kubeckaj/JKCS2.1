#!/usr/bin/env python3
'''Dynamic Approach for Transition State'''
###############################LIBRARIES#####################################
import numpy as np
import argparse
import os
import re

###########################VECTOR MANIPULATION################################
def calculate_vector(coord1, coord2):
    return np.array(coord2) - np.array(coord1)

def vector_length(vector):
    return np.linalg.norm(vector)

def atom_distance(atom1, atom2):
    return np.linalg.norm(np.array(atom2) - np.array(atom1))

def normalize_vector(vector):
    norm = np.linalg.norm(vector)
    if norm < 1e-8: 
        return np.zeros_like(vector)
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


def crest_constrain(file_path, C_index, H_index, O_index, OH_index, force_constant=0.25):
    '''Force constant tells how tight to constrain the atoms i.e. the magnitude of the oscillation between the constrained atoms'''
    with open (file_path,"w") as f:
        f.write("$constrain\n")
        f.write(f"  force constant={force_constant}\n") 
        f.write(f"  distance: {C_index}, {H_index}, auto\n")
        f.write(f"  distance: {H_index}, {O_index}, auto\n")
        f.write(f"  distance: {O_index}, {OH_index}, auto\n")
        f.write(f"  angle: {C_index}, {H_index}, {O_index}, auto\n")
        f.write("$end\n")


def orca_constrain(file_path, coords, C_index, H_index, O_index, OH_index, method="wB97X-D3", basis_set="6-31+G(d,p)"):
    with open(file_path, "w") as f:
        f.write(f"! {method} {basis_set} Opt\n")
        f.write("%pal nprocs 8 end\n")
        f.write("%maxcore 8000\n")
        f.write("%geom\n")
        f.write("Constraints\n")
        # ORCA start atom indexing from 0
        f.write(f"{{B {C_index-1} {H_index-1} C}}\n") # Bond being broken
        f.write(f"{{B {H_index-1} {O_index-1} C}}\n") # Bond being formed
        f.write(f"{{B {O_index-1} {OH_index-1} C}}\n") # Bond being formed
        f.write("end\nend\n \n * xyz 0 2\n") # charge and spin multiplicity usually 0 and 2 for atmospheric radical reaction
        for atom in coords:
            f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')
        f.write("*")


def gauss_constrain(file_path, coords):
    pass


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


def get_terminal_O(coords, distance=1.5):
    O_index = []
    for i,atom in enumerate(coords):
        if atom[0] == "O":
            neighbors = sum(1 for j, other_atom in enumerate(coords) if i != j and atom_distance(atom[1:], other_atom[1:]) < distance)
            if neighbors == 1:
                O_index.append(i)
    return O_index[0]

def molecule_controid(coords):
    coords = np.array(coords)
    
    # Extract x, y, z coordinates
    x_coords = coords[:, 1].astype(float)
    y_coords = coords[:, 2].astype(float)
    z_coords = coords[:, 3].astype(float)

    # Calculate centroid (average) for each coordinate
    centroid_x = np.mean(x_coords)
    centroid_y = np.mean(y_coords)
    centroid_z = np.mean(z_coords)

    return (centroid_x, centroid_y, centroid_z)
######################################################################################

def H_abstraction(file, distance=1.35, dist_OH=0.97, constraint=False, constrained_opt=False):
    coords = read_xyz_file(file)
    num_atoms = len(coords)
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

                        base_file_name = os.path.splitext(os.path.basename(file))[0]
                        write_xyz_file(f"{base_file_name}_H{count}.xyz", new_coords)


                        if constrained_opt:
                            C_index = j+1
                            H_index = i+1
                            O_index = len(new_coords)-1
                            OH_index = len(new_coords)
                            orca_constrain(file_path=f"{base_file_name}_H{count}.inp", coords=new_coords,C_index=C_index, H_index=H_index, O_index=O_index, OH_index=OH_index)

                        count += 1


def OH_addition(file, distance=1.45, double_bond_distance=1.36, dist_OH=0.97):
    coords = read_xyz_file(file)
    num_atoms = len(coords)
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
                        # NOTE: update the way perpendicular axis is calculated. Copy from addition()
                        perpendicular_axis = np.cross(norm_vector_CC, [0,1,0])
                        perpendicular_axis = normalize_vector(perpendicular_axis)
                        # Shift carbon
                        new_coords[i][1:] = np.array(new_coords[i][1:]) + norm_vector_CC * 0.1
                        # Update oxygen in OH coordiantes
                        oxygen_coords = np.array(new_coords[j][1:]) + perpendicular_axis * distance

                        rotation_axis = norm_vector_CC
                        rotation_angle_H = np.radians(45) 

                        rotated_vector_H = rotate_vector(perpendicular_axis, rotation_axis, rotation_angle_H)

                        hydrogen_coords = oxygen_coords + rotated_vector_H * dist_OH
                        
                        new_coords.append(['O', *oxygen_coords])
                        new_coords.append(['H', *hydrogen_coords])

                        base_file_name = os.path.splitext(os.path.basename(file))[0]
                        write_xyz_file(f"{base_file_name}_{count}.xyz", new_coords)

                        count += 1


def addition(file1, file2, distance = 1.55, double_bond_distance=1.36, elongation_factor=0.1):
    coords1 = read_xyz_file(file1)
    coords2 = read_xyz_file(file2)
    num_atoms1 = len(coords1)
    num_atoms2 = len(coords2)
    count = 1
    carbon_index = []
    vector_lst = []
    
    # Gather carbon index
    for i in range(num_atoms1):
        if coords1[i][0] == "C":
            for j in range(num_atoms1):
                if atom_distance(coords1[i][1:], coords1[j][1:]) < 1.52:
                    perpendicular_vector = (coords1[i][1:], coords1[j][1:])
                    vector_lst.append(perpendicular_vector)
                if coords1[j][0] == "C" and i != j:
                    bond_CC = atom_distance(coords1[i][1:], coords1[j][1:])
                    if bond_CC <= double_bond_distance:
                        carbon_index.append(i)
                        carbon_index.append(j)

    for index in carbon_index[::2]:
        new_coords = read_xyz_file(file1)
        direction = calculate_vector(new_coords[index][1:], new_coords[index+1][1:]) # Vector between C=C
        for k in range(num_atoms1):
            if atom_distance(coords1[index][1:], coords1[k][1:]) < 1.52:
                perpendicular_vector = calculate_vector(coords1[index][1:], coords1[k][1:])
                break

        new_coords[index][1:] -= direction * elongation_factor/2 # Standard is 10% increase in bond length
        new_coords[index+1][1:] += direction * elongation_factor/2

        perpendicular_axis = np.cross(direction, perpendicular_vector)
        perpendicular_axis = normalize_vector(perpendicular_axis)

        O_index = get_terminal_O(coords2)

        first_pos = new_coords[index][1:] + perpendicular_axis * distance

        rotation_axis = normalize_vector(np.cross(direction, perpendicular_axis))

        # NOTE: fix the rotation axis
        for i in range(num_atoms2):
            relative_vector = calculate_vector(coords2[O_index][1:], coords2[i][1:])
            dot = np.dot(perpendicular_axis, normalize_vector(relative_vector))
            if dot > 0:
                 rotation_angle = np.radians(90) # radians value should be function of dot product
                 rotated_vector = rotate_vector(relative_vector, rotation_axis, rotation_angle)
            else:
                rotation_angle = np.radians(-270)
                rotated_vector = rotate_vector(relative_vector, rotation_axis, rotation_angle)
            new_pos = first_pos + relative_vector + (rotated_vector * 0.1)
            new_coords.append([coords2[i][0], *new_pos])
        

        base_file_name1 = os.path.splitext(os.path.basename(file1))[0]
        base_file_name2 = os.path.splitext(os.path.basename(file2))[0]      
        write_xyz_file(f"{base_file_name1}_{base_file_name2}_{count}.xyz", new_coords)
        count += 1

                        

def rate_constant_gauss(file1):
    with open (file1, "r") as file:
        content = file.read()
        match = re.search(r"Sum of electronic and zero-point Energies=\s+([-+]?\d*\.\d+|\d+)", content) 
    if match:
        float_value = float(match.group(1))
        print(float_value)




def main():
    parser = argparse.ArgumentParser(description='TS Generation Tool',
                                     prog="JKTS",
                                     epilog="")

    parser.add_argument('reactant1', metavar='file1.xyz', help='First reactant XYZ file (e.g., pinonaldehyde.xyz)')
    parser.add_argument('reactant2', nargs='?', metavar='file2.xyz', help='Second reactant XYZ file (e.g., peroxy.xyz)')
    parser.add_argument('-CC', action='store_true', help='Perform addition to C=C')
    parser.add_argument('-H', action='store_true', help='Perform H abstraction')
    parser.add_argument('Program', action='store_true', help='Gaussian or ORCA')
    parser.add_argument('-constrained_opt', action='store_true', help='Create constraint.inp file')
    parser.add_argument('-k', action='store_true', help='Write me later')

    args = parser.parse_args()

    if args.reactant1 and args.reactant2 == None:
        file_in = args.reactant1
    elif args.reactant1 and args.reactant2:
        file1 = args.reactant1
        file2 = args.reactant2
    else:
        parser.error('Please provide xyz file(s) of reactant(s)')


    if args.reactant1:
        if args.H:
            if args.constrained_opt:
                H_abstraction(file_in, constrained_opt=True)
            else:
                H_abstraction(file_in)
        elif args.CC:
            OH_addition(file_in)
        elif args.reactant1 and args.reactant2:
            addition(file1, file2)
        elif args.k:
            rate_constant_gauss(file_in)
        else:
            parser.error('Please specify either -H or -CC or another xyz-file.')

if __name__ == "__main__":
    main()

