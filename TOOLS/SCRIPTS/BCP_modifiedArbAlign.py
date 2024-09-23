#!/usr/bin/python

import numpy as np
import lapjv
from collections import Counter 
import operator


# ArbAlign Lite Edit
# New stuff: 
#  - Takes ASE atoms as input i.e. import ArbAlign and call ArbAlign.compare(ase.atom,ase.atom)
#  - Converted from python2 to python3
#  - Uses the labjv Jonker–Volgenant algorithm instead of the python2 hungarian library
#  - Disabled options: Will always do swap/reflections and will always include H.¨
#  - Simplified output: Will always just return the lowest RMSD as a float


def kabsch(A, B):
   """
   Kabsch Algorithm as implemented by Jimmy Charnley Kromann

   Calculate RMSD between two XYZ files

   by: Jimmy Charnley Kromann <jimmy@charnley.dk> and 
   Lars Andersen Bratholm <larsbratholm@gmail.com>
   project: https://github.com/charnley/rmsd
   license: https://github.com/charnley/rmsd/blob/master/LICENSE

   A - set of coordinates
   B - set of coordinates

   Performs the kabsch algorithm to calculate the RMSD between A and B

   Returns an RMSD
   """
   A_new = np.array(A)
   A_new = A_new - sum(A_new) / len(A_new)
   A = A_new
   B_new = np.array(B)
   B_new = B_new - sum(B_new) / len(B_new)
   B = B_new

   # Compute covariance matrix
   C = np.dot(np.transpose(A), B)

   # Compute singular value decomposition (SVD)
   V, S, W = np.linalg.svd(C)
   d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

   if d:
      S[-1] = -S[-1]
      V[:, -1] = -V[:, -1]

   # Compute rotation matrix
   U = np.dot(V, W)

   # Rotate A
   A = np.dot(A, U)

   return rmsd(A, B)

def rmsd(V, W):
    """
    V - set of coordinates
    W - set of coordinates

    Returns root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def sorted_xyz(filename, noHydrogens):
   qq = [list(x) for x in filename.get_positions()]
   test = sorted(zip(filename.get_chemical_symbols(),qq,[x for x in range(len(filename.get_positions()))]),key=lambda x: (x[0],x[1][0],x[1][1],x[1][2]))
   sortedlabels = [x[0] for x in test]
   sortedcoords = [x[1] for x in test]
   sortedorder = [x[2] for x in test]
   NA = len(sortedlabels)
   return (sortedlabels, sortedcoords, NA, sortedorder)
   """
   Reads an xyz file into lists
   filename - name of xyz file
   noHydrogens - if true, hydrogens are ignored; if not, hydrogens are included

   Returns a tuple (a, b) sorted by atom labels and coordinates
    where a is a list of coordinate labels and b is a set of coordinates
   (i.e) a = ["O", "H", "H"], b = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]
   Sorts the file by atom labels first and coordinates second 
   such that atoms of the same label/type are grouped together
   """
   xyz = open(filename, "r")
   num_atoms = int(xyz.readline().strip())
   xyz.readline()

   sortedlabels = []
   sortedcoords = []
   sortedorder = []
   sortedlines = []
   atomcount = 0
   for i in range(num_atoms):
      line = xyz.readline().strip().split()
      if noHydrogens and line[0].upper().startswith("H") : 
         continue
      else:
         sortedlines.append([line[0].upper(),float(line[1]), float(line[2]), float(line[3]), atomcount])
         atomcount += 1
   xyz.close()

   # sort by element followed by first coordinate (x) then second coordinate (y)
   sortedlines.sort(key=lambda x: (x[0],x[1],x[2]))
   for i in range(atomcount):
      sortedlabels.append(sortedlines[i][0])
      sortedcoords.append([float(sortedlines[i][1]), float(sortedlines[i][2]), float(sortedlines[i][3])])
      sortedorder.append(sortedlines[i][4])
   xyz.close()
   #print sortedlabels
   NA = len(sortedlabels)
   return sortedlabels, sortedcoords, NA, sortedorder

def parse_for_atom(labels, coords, atom):
   """
   labels - a list of coordinate labels
   coords - a set of coordinates
   atom - the atom that is to be parsed

   Returns a set of coordinates corresponding to parsed atom
   """
   atom_coords = []
   for i in range(len(labels)):
      if labels[i] == atom:
           atom_coords.append(coords[i])
   return atom_coords

def transform_coords(coords, swap, reflect):
   """
   coords - a set of coordinates
   swap - the swap transformation (i.e. (0, 2, 1) --> (x, z, y))
   reflect - the reflection transformation (i.e. (1, -1, 1) --> (x, -y, z))

   Returns the transformed coordinates
   """
   new_coords = []
   for i in range(len(coords)):
      new_coords.append([coords[i][swap[0]]*reflect[0],
                          coords[i][swap[1]]*reflect[1],
                          coords[i][swap[2]]*reflect[2]])
   return new_coords

def transform_atoms(coords, swap, reflect, atom_indices):
   """
   coords - a set of coordinates
   swap - the swap transformation (i.e. (0, 2, 1) --> (x, z, y))
   reflect - the reflection transformation (i.e. (1, -1, 1) --> (x, -y, z))
   atom_indices - indices of all desired atoms in [coords]

   Returns coordinates after transforming specific atoms
   """
   new_coords = [x[:] for x in coords]
   for i in atom_indices:
      new_coords[i][0] = coords[i][swap[0]]*reflect[0]
      new_coords[i][1] = coords[i][swap[1]]*reflect[1]
      new_coords[i][2] = coords[i][swap[2]]*reflect[2]
   return new_coords

def permute_coords(coords, permutation):
   """
   UNUSED at the moment 

   coords - a set of coordinates
   permutation - permutation of atoms (i.e. [0, 2, 3, 1])

   Returns the permuted coordinates
   """
   new_coords = []
   for i in permutation:
      new_coords.append(coords[i])
   return new_coords

def permute_atoms(coords, permutation, atom_indices):
   """
   coords - a set of coordinates
   permuation - a permutation of atoms
   atom_indices - indices of all desired atoms in [coords]

   Returns the coordinates after permuting just the specified atom
   """
   new_coords = coords[:]
   for i in range(len(permutation)):
      j = atom_indices[permutation[i]]
      k = atom_indices[i]
      new_coords[k] = coords[j] 
   return new_coords

def permute_all_atoms(labels, coords, permutation):
   """
   labels - atom labels 
   coords - a set of coordinates
   permuation - a permutation of atoms

   Returns the permuted labels and coordinates
   """
   new_coords = coords[:]
   new_labels = labels[:]
   for i in range(len(permutation)):
      new_coords[permutation[i]] = coords[i]
      new_labels[permutation[i]] = labels[i]
   return new_labels, new_coords

def get_atom_indices(labels, atom):
   """
   labels - a list of coordinate labels ("Elements")
   atom - the atom whose indices in labels are sought
   Returns a list of all location of [atom] in [labels]
   """
   indices = []
   for i in range(len(labels)):
      if labels[i] == atom:
         indices.append(i)
   return indices

def coords_to_xyz(labels, coords):
   """
   Displays coordinates
   """
   s = ""
   for i in range(len(coords)):
      s += labels[i] + " " + str(coords[i][0]) + " " + str(coords[i][1]) + " " + str(coords[i][2]) + "\n"
   return s

def write_to_xyz(num_atoms, name, labels, coords):
   """
   num_atoms - number of atoms 
   name - name of file to write coordinates to
   labels - a list of coordinate labels ("Elements")
   coords - a list of XYZ coordinates

   Writes the Cartesian coordinates to a file called 'name'
   """
   xyz = open(name, "w")
   xyz.write(str(num_atoms) + "\n")
   xyz.write(name + "\n")
   xyz.write(coords_to_xyz(labels, coords))
   xyz.close()

def read_xyz_from_molecule(molecule, noHydrogens=False):
    """
    Reads molecule data into lists

    molecule - Molecule object
    noHydrogens - if true, hydrogens are ignored; if not, hydrogens are included

    Returns a tuple (a, b)
    where a is a list of coordinate labels and b is a set of coordinates
    (i.e) a = ["O", "H", "H"], b = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]
    """
    if not molecule.atoms or not molecule.coordinates:
        raise ValueError("Molecule object does not have atoms or coordinates data")

    unsorted_labels = []
    unsorted_coords = []
    
    for atom, coord in zip(molecule.atoms, molecule.coordinates):
        if noHydrogens and atom.upper().startswith("H"):
            continue
        else:
            unsorted_labels.append(atom.upper())
            unsorted_coords.append(coord)

    NA = len(unsorted_labels)
    return unsorted_labels, unsorted_coords, NA

def sorted_xyz_from_molecule(molecule, noHydrogens=False):
    """
    Reads molecule data into lists and sorts them

    molecule - Molecule object
    noHydrogens - if true, hydrogens are ignored; if not, hydrogens are included

    Returns a tuple (a, b, c, d) sorted by atom labels and coordinates
    where a is a list of coordinate labels, b is a set of coordinates, c is the number of atoms, and d is the sorted order
    (i.e) a = ["O", "H", "H"], b = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2]]
    Sorts by atom labels first and coordinates second 
    such that atoms of the same label/type are grouped together
    """
    if not molecule.atoms or not molecule.coordinates:
        raise ValueError("Molecule object does not have atoms or coordinates data")

    sortedlines = []
    atomcount = 0

    for atom, coord in zip(molecule.atoms, molecule.coordinates):
        if noHydrogens and atom.upper().startswith("H"):
            continue
        else:
            sortedlines.append([atom.upper()] + coord + [atomcount])
            atomcount += 1

    # sort by element followed by first coordinate (x) then second coordinate (y)
    sortedlines.sort(key=lambda x: (x[0], x[1], x[2]))

    sortedlabels = [line[0] for line in sortedlines]
    sortedcoords = [[line[1], line[2], line[3]] for line in sortedlines]
    sortedorder = [line[4] for line in sortedlines]

    NA = len(sortedlabels)
    return sortedlabels, sortedcoords, NA, sortedorder

def compare(in_a,in_b,simpleit=0,noHydrogens=0):
   try:
      a_labels = in_a.get_chemical_symbols()
      b_labels = in_b.get_chemical_symbols()
      a_coords = in_a.get_positions()
      b_coords = in_b.get_positions()
      NA_a = len(a_labels)
      NA_b = len(b_labels)
   except:
      from classes import Molecule
      a_labels, a_coords, NA_a = read_xyz_from_molecule(in_a, noHydrogens)
      b_labels, b_coords, NA_b = read_xyz_from_molecule(in_b, noHydrogens)

   b_init_labels = b_labels 
   b_init_coords = b_coords
   
   #Calculate the initial unsorted all-atom RMSD as a baseline
   A_all = np.array(a_coords)
   B_all = np.array(b_coords)

   #If the two molecules are of the same size, get 
   if NA_a == NA_b:
      InitRMSD_unsorted = kabsch(A_all,B_all)
   else:
      return("Error: unequal number of atoms. " + str(NA_a) + " is not equal to " + str(NA_b))
 

   """
   If the initial RMSD is zero (<0.001), then the structured are deemed identical already and 
   we don't need to do any reordering, swapping, or reflections
   """
   if InitRMSD_unsorted < 0.001:
      return(float(InitRMSD_unsorted))
   
   """
   Read in the original coordinates and labels of xyz1 and xyz2, 
   and sort them by atom labels so that atoms of the same label/name are grouped together

   Then, count how many types of atoms, and determine their numerical frequency
   """
   try:
      a_labels, a_coords, NA_a, order = sorted_xyz(in_a, noHydrogens)
      b_labels, b_coords, NA_b, junk = sorted_xyz(in_b, noHydrogens)
   except:
      a_labels, a_coords, NA_a, order = sorted_xyz_from_molecule(in_a, noHydrogens)
      b_labels, b_coords, NA_b, junk = sorted_xyz_from_molecule(in_b, noHydrogens)

   Uniq_a = list(set(a_labels))
   list.sort(Uniq_a)
   N_uniq_a = len(Uniq_a)
   Atom_freq_a = dict(Counter(a_labels))

   Uniq_b = list(set(b_labels))
   list.sort(Uniq_b)
   N_uniq_b = len(Uniq_b)
   Atom_freq_b = dict(Counter(b_labels))

   """
   If the number and type of atoms in the two structures are not equal, exit with 
   an error message
   """
   if (NA_a == NA_b) & (Uniq_a == Uniq_b) & (Atom_freq_a == Atom_freq_b) :
      num_atoms = NA_a
      num_uniq = N_uniq_a
      Uniq = Uniq_a  
      Atom_freq = Atom_freq_a
      Sorted_Atom_freq = sorted(Atom_freq.items(), key=operator.itemgetter(1), reverse=True)
      """
      Atom = sorted(Uniq, key=operator.itemgetter(0), reverse=True)
      print Atom
      print num_uniq
      """
   else:
      return("Unequal number or type of atoms. Exiting ... ")

   A_all = np.array(a_coords)
   A_all = A_all - sum(A_all) / len(A_all)
   B_all = np.array(b_coords)
   B_all = B_all - sum(B_all) / len(B_all)
   InitRMSD_sorted = kabsch(A_all,B_all)
   
   """
   Dynamically generate hashes of coordinates and atom indices for every atom type
   """
   a_Coords = {}
   a_Indices = {}
   b_Coords = {}
   b_Indices = {}
   Perm = {}
   for i in range(len(Uniq)):
      a_Coords[Uniq[i]]  = 'a_' + Uniq[i] + 'coords' 
      b_Coords[Uniq[i]]  = 'b_' + Uniq[i] + 'coords' 
      a_Indices[Uniq[i]] = 'a_' + Uniq[i] + 'indices' 
      b_Indices[Uniq[i]] = 'b_' + Uniq[i] + 'indices' 
      Perm[Uniq[i]]      = 'perm_' + Uniq[i]
      vars()[Perm[Uniq[i]]] = []
      vars()[a_Coords[Uniq[i]]]  = parse_for_atom(a_labels, a_coords, str(Uniq[i]))
      vars()[a_Indices[Uniq[i]]] = get_atom_indices(a_labels, str(Uniq[i]))
      vars()[b_Coords[Uniq[i]]]  = parse_for_atom(b_labels, b_coords, str(Uniq[i]))
      vars()[b_Indices[Uniq[i]]] = get_atom_indices(b_labels, str(Uniq[i]))

   l = 0 
   A = np.array(vars()[a_Coords[Uniq[l]]])
   A = A - sum(A) / len(A)
   B = np.array(vars()[b_Coords[Uniq[l]]])
   B = B - sum(B) / len(B)

   '''
   For each atom type, we can do a Kuhn-Munkres assignment in the initial 
   coordinates or the many swaps and reflections thereof

   If a single Kuhn-Munkres assignment is requested with a -s or --simple flag,
   no swaps and reflections are considered. Otherwise, the default is to perform 
   a combination of 6 axes swaps and 8 reflections and do Kuhn-Munkres assignment 
   on all 48 combinations. 
   '''

   swaps = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
   reflects = [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1),
               (-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1)]
   B_t = []
   for i in swaps:
      for j in reflects:
         B_t.append([transform_coords(B, i, j), i, j])

   rmsds = []
   # Performs the munkres algorithm on each set of transformed coordinates
   for i in range(len(B_t)):
      l = 0
      cost_matrix = np.array([[np.linalg.norm(a - b)
                                 for b in B_t[i][0]] for a in A])
      q1,q2,xx = lapjv.lapjv(cost_matrix)
      LAP = (q1,q2)
      vars()[Perm[Uniq[l]]] = []
      for j in range(len(LAP[0])):
         vars()[Perm[Uniq[l]]] += [(j,LAP[0][j])]
      
      vars()[Perm[Uniq[l]]] = sorted( vars()[Perm[Uniq[l]]], key = lambda x: x[0])
      vars()[Perm[Uniq[l]]] = [x[1] for x in vars()[Perm[Uniq[l]]]]
      # If there's more than one atom type, loop through each unique atom type 
      if num_uniq == 1:
         b_perm = permute_atoms(b_coords, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
         b_final = transform_coords(b_perm, B_t[i][1], B_t[i][2])
         rmsds.append([kabsch(a_coords, b_final), B_t[i][1], B_t[i][2], b_final, vars()[Perm[Uniq[l]]]])
         rmsds = sorted(rmsds, key = lambda x: x[0])
      else: 
         b_perm = permute_atoms(b_coords, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
         b_trans = transform_coords(b_perm, B_t[i][1], B_t[i][2])
         while l < num_uniq:
            if l > 0:
               vars()[b_Coords[Uniq[l]]] = parse_for_atom(b_labels, b_final, Uniq[l])
            else:
               vars()[b_Coords[Uniq[l]]] = parse_for_atom(b_labels, b_trans, Uniq[l])
            lol1 = np.array(vars()[b_Coords[Uniq[l]]])
            lol2 = np.array(vars()[a_Coords[Uniq[l]]])
            cost_matrix = np.array([[np.linalg.norm(a-b) for b in lol1] for a in lol2])
            q1,q2,xx = lapjv.lapjv(cost_matrix)
            LAP = (q1,q2) #New one gives the indices fucked flipped Q_Q
            vars()[Perm[Uniq[l]]] = []
            for k in range(len(LAP[0])):
               vars()[Perm[Uniq[l]]] += [(k,LAP[0][k])]
            vars()[Perm[Uniq[l]]] = sorted( vars()[Perm[Uniq[l]]], key = lambda x: x[0])
            vars()[Perm[Uniq[l]]] = [x[1] for x in vars()[Perm[Uniq[l]]]]
            b_final = permute_atoms(b_trans, vars()[Perm[Uniq[l]]], vars()[b_Indices[Uniq[l]]])
            b_trans = b_final
            l += 1
            q = l - 1 
            rmsds.append([kabsch(a_coords, b_final), B_t[i][1], B_t[i][2], b_final])
            rmsds = sorted(rmsds, key = lambda x: x[0])
   FinalRMSD = float(rmsds[0][0])
   if FinalRMSD < float(InitRMSD_unsorted): 
      return(float(rmsds[0][0]))
   else:
      return(float(InitRMSD_unsorted))
