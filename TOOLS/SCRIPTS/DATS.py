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


#########################################GENERATE INPUT FILES#############################
def crest_constrain(file_path, C_index, H_index, O_index, force_constant=0.25):
    '''Force constant tells how tight to constrain the atoms i.e. the magnitude of the oscillation between the constrained atoms'''
    with open (file_path,"w") as f:
        f.write("$constrain\n")
        f.write(f"  force constant={force_constant}\n") 
        f.write(f"  distance: {C_index}, {H_index}, auto\n")
        f.write(f"  distance: {H_index}, {O_index}, auto\n")
        f.write(f"  angle: {C_index}, {H_index}, {O_index}, auto\n")
        f.write("$end\n")


def orca_input(file_path, coords, TS,constrain, C_index=None, H_index=None, O_index=None, method="wB97X-D3", basis_set="6-31+G(d,p)"):
    with open(file_path, "w") as f:
        if TS:
            f.write(f"! {method} {basis_set} OptTS freq\n")
        elif TS is False:
            f.write(f"! {method} {basis_set} Opt\n")
        f.write("%pal nprocs 8 end\n")
        f.write("%maxcore 8000\n")
        if constrain:
            f.write("%geom\n")
            f.write("Constraints\n")
            # ORCA start atom indexing from 0
            f.write(f"{{B {C_index-1} {H_index-1} C}}\n") # Bond being broken
            f.write(f"{{B {H_index-1} {O_index-1} C}}\n") # Bond being formed
            f.write("end\nend\n")
        f.write("\n")
        f.write("* xyz 0 2\n") # charge and spin multiplicity usually 0 and 2 for atmospheric radical reaction
        for atom in coords:
            f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')
        f.write("*")


def gauss_input(file_path, coords, TS, constrain, C_index=None, H_index=None, O_index=None, method="wb97xd", basis_set="6-31+G(d,p)"):
    with open(file_path, "w") as f:
        f.write("%nprocshared=4\n")
        f.write("%mem=4GB\n")
        if TS and constrain:
            f.write(f"# {method} {basis_set} opt=(calcfc,ts,noeigen,modredundant) freq\n\n") # freq may be redundant
        elif TS and constrain is False:
            f.write(f"# {method} {basis_set} opt=(calcfc,ts,noeigen) freq\n\n")
        elif TS is False and constrain:
            f.write(f"# {method} {basis_set} opt=modredundant\n\n")
        else:
            f.write(f"# {method} {basis_set} opt\n\n")
        f.write("Title\n\n")
        f.write("0 2\n")
        for atom in coords:
            f.write(f'{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n')
        f.write("\n")
        if constrain:
            f.write(f"B {C_index} {H_index} F\n")
            f.write(f"B {H_index} {O_index} F\n")
            f.write("\n")


##############################ADDITIONAL FUNCTIONS##################################
def get_terminal_O(coords, distance=1.5):
    O_index = []
    for i,atom in enumerate(coords):
        if atom[0] == "O":
            neighbors = sum(1 for j, other_atom in enumerate(coords) if i != j and atom_distance(atom[1:], other_atom[1:]) < distance)
            if neighbors == 1:
                O_index.append(i)
    return O_index[0]

def molecule_centroid(coords):
    coords = np.array(coords)
    
    x_coords = coords[:, 1].astype(float)
    y_coords = coords[:, 2].astype(float)
    z_coords = coords[:, 3].astype(float)

    centroid_x = np.mean(x_coords)
    centroid_y = np.mean(y_coords)
    centroid_z = np.mean(z_coords)

    return (centroid_x, centroid_y, centroid_z)


#####################################MAIN FUNCTIONS####################################
def prepare_file(file, method, basis_set, TS, XYZ=False, program=None, constrain=False):
    coords = read_xyz_file(file)
    file_name = file.split(".")[0] 

    if program == "gauss" or program == None: 
        gauss_input(file_path=f"{file_name}.com", coords=coords, constrain=constrain, TS=TS, method=method, basis_set=basis_set)
    elif program == "orca":
        orca_input(file_path=f"{file_name}.inp", coords=coords, constrain=constrain, TS=TS, method=method, basis_set=basis_set)
    # elif program == "crest":
    #     crest_constrain(file_path=f"{file}.inp")




def H_abstraction(file, method, basis_set, TS, XYZ, program=None, distance=1.35, dist_OH=0.97, constrain=False):
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
                        if XYZ is True:
                            write_xyz_file(f"{base_file_name}_H{count}.xyz", new_coords)


                        if constrain or method is not None or basis_set is not None:
                            C_index = j+1
                            H_index = i+1
                            O_index = len(new_coords)-1
                            OH_index = len(new_coords)
                            # Default program if nothing specified is Gaussian
                            if program == "gauss" or program == None: 
                                gauss_input(file_path=f"{base_file_name}_H{count}.com", coords=new_coords,C_index=C_index, H_index=H_index, O_index=O_index, constrain=constrain, TS=TS, method=method, basis_set=basis_set)
                            elif program == "orca":
                                orca_input(file_path=f"{base_file_name}_H{count}.inp", coords=new_coords, C_index=C_index, H_index=H_index, O_index=O_index, constrain=constrain, TS=TS, method=method, basis_set=basis_set)
                            elif program == "crest":
                                crest_constrain(file_path=f"{base_file_name}_H{count}.inp", C_index=C_index, H_index=H_index, O_index=O_index)

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


def addition(file1, file2, method, basis_set, TS, XYZ, program=None, constrain=False, distance = 1.55, double_bond_distance=1.36, elongation_factor=0.1):
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
        # perpendicular_vector = [calculate_vector(coords1[index][1:], coords1[k][1:]) for k in range(num_atoms1) if atom_distance(coords1[index][1:], coords1[k][1:]) < 1.52]

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

                        

##################################THERMOCHEMISTRY##############################
def rate_constant(files, T=293.15, program=None):
    k_b = 1.380649*10**-23 # J/K
    h = 6.62607*10**-34 # J*s
    R = 8.314 # J/mol*K 

    files_energies = []
    count = 0
# Loop over inp file and determine whether reactant or TS. Append to respective Bolztmann sumation. Exit loop and start list comprehension to calculate the Bolztmann summation factor. Lastly a list comprehension for the rate constant looping over the E_TS - E_R in the exponential.
    for file in files:
        dic = {}
        with open (file, "r") as f:
            content = f.read()
            ee_zpe = re.search(r"Sum of electronic and zero-point Energies=\s+([-+]?\d*\.\d+|\d+)", content) 

        if ee_zpe:
            dic[file] = count
            dic["Sum of eletronic and ZPE"] = float(ee_zpe.group(1))
            files_energies.append(dic)
            count += 1
        else: print(f"No energies in {file}. Check if calculation has converged")
    print(files_energies)



def main():
    ascii_art = '''
     ______   _______ _________ _______ 
    (  __  \\ (  ___  )\\__   __/(  ____ \\
    | (  \\  )| (   ) |   ) (   | (    \\/
    | |   ) || (___) |   | |   | (_____ 
    | |   | ||  ___  |   | |   (_____  )
    | |   ) || (   ) |   | |         ) |
    | (__/  )| )   ( |   | |   /\\____) |
    (______/ |/     \\|   )_(   \\_______)
'''
    
    parser = argparse.ArgumentParser(description='''    'Dynamic Approach for Transition State'
Automated tool for generating input files, primarily 
for transition state geometry optimization. 
Calculation of tunneling corrected multi-configurational 
rate constants can also be calculated from log files.''',
                                     prog="JKTS",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='''
Examples of use:
                JKTS file.xyz -H -orca -method R2SCAN-3C
                JKTS *.xyz -basis def2-TZVPPD --no-xyz
                JKTS reactant.log product.log -k
                                     ''')


    parser.add_argument('input_files', metavar='file.xyz', nargs='+', help='XYZ files (e.g., pinonaldehyde.xyz or even *.xyz)')

    parser.add_argument('-CC', action='store_true', help='Perform addition to C=C')
    parser.add_argument('-H', action='store_true', help='Perform H abstraction')
    parser.add_argument('-gauss', action='store_true', help='Create Gaussian16 input file [defualt]')
    parser.add_argument('-orca', action='store_true', help='Create ORCA input file')
    parser.add_argument('-crest', action='store_true', help='Create crest input file')
    parser.add_argument('-constrain', action='store_true', help='Constrain is integrated into input file')

    additional_options = parser.add_argument_group("Additional arguments")

    additional_options.add_argument('--no-TS',action='store_true', help='Input files for normal geometry relaxation are generated')
    additional_options.add_argument('--no-xyz',action='store_true', help='No XYZ files generated')
    additional_options.add_argument('-k', action='store_true', help='Write me later')
    additional_options.add_argument('-method', dest='method_name', help='Specify the QC method [def = wB97X-D]')
    additional_options.add_argument('-basis', dest='basis_name', help='Specify the basis set used [def = 6-311+G(d,p)]')

    args = parser.parse_args()
    

    ################################ARGUMENT SPECIFICATIONS############################
    # Program 
    if args.gauss:
        program = "gauss"
    elif args.orca:
        program = "orca"
    elif args.crest:
        program = "crest"
    else: program = None # Just produce XYZ-file

    # Method
    if args.method_name:
        method = args.method_name
    elif args.method_name not in args and program == "orca":
        method = "WB97X-D3" # ORCA 
    else: method = "wb97xd" # Gaussian

    # Basis set
    if args.basis_name:
        basis_set = args.basis_name
    else: basis_set = "6-31+G(d,p)"

    # Constrain 
    if args.constrain:
        constrain = True
    else: constrain = False

    if args.no_TS is True:
        TS = False
    else: TS = True

    if args.no_xyz is True:
        XYZ = False
    else: XYZ = True


    for n, input_file in enumerate(args.input_files):
        if len(args.input_files) == 2 and args.CC:
            addition(args.input_files[n], args.input_files[n+1], constrain=constrain, program=program, method=method, TS=TS, XYZ=XYZ, basis_set=basis_set)
            break
        elif len(args.input_files) > 0 and args.H:
            H_abstraction(input_file, constrain=constrain, TS=TS, program=program, method=method, basis_set=basis_set, XYZ=XYZ)
        elif len(args.input_files) > 0 and args.H==False and args.CC==False and args.k==False:
            prepare_file(input_file, constrain=constrain, TS=TS, program=program, method=method, basis_set=basis_set, XYZ=XYZ)
        elif len(args.input_files) > 0 and args.k:
            rate_constant(args.input_files)
            break
        else:
            parser.error("No input file given")
   

if __name__ == "__main__":
    main()

