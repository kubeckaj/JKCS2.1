
import os, sys
import numpy as np
import pandas as pd
from cluster_analysis import *

def printhelp():
  print('   Use this file like this:')
  print('           ClusterFilter.py [FILE] {OPTIONS}')
  print(' OPTIONS:')
  print('  -i/-itol 0.3     filter out if any intermolecular bond length has changed more than 30% [default=0.2]')
  print('  -d/-dmol 2.5     filter out if smallest distance between molecules greater than 2.5 [default=2.7]')
  print('  -h 2             filter clusters based on number of hydrogen bonds')
  print('  -reject sa,am    filter out subclusters consisting only specified molecules (separated with ",")')
  print('  -subcl 2         minimum size of returned subclusters (1 returns monomers) [default=2]')
  print('  ')

#######################################################################

internal_limit = 0.2    # largest allowed relative error
max_dist = 2.7          # maximum distance for H-bonding between molecules
min_dist = 1.4          # minimum allowed distance between molecules (otherwise a reaction has occured)
print_minimums = False
count_minimums = 1
subcl_size = 2
Hbonding = None

# Look for inputs
file_in = ""
last = ""
rejected_types = [] # e.g. sa,a,so4
n_arg = len(sys.argv)
for i in range(1,n_arg):
    arg=sys.argv[i]
    if str(arg) == "--help" or str(arg) == "-help":
        printhelp()
        exit()
    if arg == "-itol" or arg == '-i':
        last = "-itol"
        continue
    if last == "-itol":
        last = ""
        internal_limit = float(arg)
        continue
    if arg == "-dmol" or arg == '-d':
        last = "-dmol"
        continue
    if last == "-dmol":
        last = ""
        max_dist = float(arg)
        continue
    if arg == "-h":
        last = "-h"
        continue
    if last == "-h":
        last = ""
        Hbonding = int(arg)
        continue
    if arg == "-reject":
        last = "-reject"
        continue
    if last == "-reject":
        last = ""
        rejected_types = arg.split(',')
        continue
    if arg == "-subcl":
        last = "-subcl"
        continue
    if arg == "-subcl":
        last = ""
        subcl_size = int(arg)
        continue
    if last == "-minimums":
        last = ""
        print_minimums = False
        count_minimums = int(arg)
        continue
    file_in = arg

if ( file_in == "" ):
  print('ClusterFilter.py: Missing input file. EXITING')
  exit()

# check amount of rows
count = 0
for line in open(file_in): 
    count += 1
    if ( count > 1 ):
        break 
if ( count == 0 ):
    print('ClusterFilter.py: Empty input file. EXITING')
    exit()
if ( count == 1 ):
    print('ClusterFilter.py: Only one structure in the input file. EXITING')
    exit()

mol_file = "parameters.txt"
if not(os.path.isfile(mol_file)): 
    mol_file = "input.txt"
    if not(os.path.isfile(mol_file)): 
        print("input.txt or parameters.txt not found. Make sure you are in the correct folder!");exit()

d = bonding_paths(mol_file) # looks for molecules in parameters.txt or input.txt
mol_df = read_molecule_data(d)
n_mol = len(mol_df)

print("Found", n_mol, "molecules in parameters.txt.")

###############################################################################
print('ClusterFilter.py: Start.')

from time import time
t = time()

# output file to be printed  
file_out=file_in[:-4]+'_FILTERED.dat'

input_pkl = []
lines = open(file_in).readlines()
paths = pd.DataFrame([l.split()[0].split('/:EXTRACT:/') for l in lines], columns=['file', 'cluster'])

# Read pickle file(s)
clusters_df = read_pickled_data([file_in])

# Remove rows with missing structures (nan)
has_atoms = [type(x) != float for x in clusters_df[('xyz','structure')].values]
clusters_df = clusters_df[has_atoms]

# leave out clusters containing only (apart from 1) molecules in rejectd_type, e.g. "sa"
def rejected_cluster_type(molecules, counts):
    if len(rejected_types) == 0:
        return False
    
    total = np.sum(counts)
    n = 0
    for i in range(len(molecules)):   
        if molecules[i] in rejected_types:
            n += counts[i]

    if n >= total-1:
        return True
    
    return False

# Check that all monomers are in mol_df
monomers = np.unique(np.concatenate(clusters_df[("info", "components")].values))
for mol in monomers:
    if mol not in mol_df.index:
        print("Error: Molecule '"+mol+"' was not found in "+mol_file);exit()

# Separate clusters per type into subsets
if 'cluster_type' in clusters_df['info'].columns:
    unique_cluster_types = np.unique(clusters_df[("info", "cluster_type")].values)
    cluster_subsets = []
    for unique_cluster_type in unique_cluster_types:
        indexes = clusters_df[clusters_df[("info", "cluster_type")]==unique_cluster_type].index.values
        cluster_subsets.append(indexes)
else:
    print("Warning: No cluster types found")
    cluster_subsets = []
    indexes = clusters_df.index.values
    cluster_subsets.append(indexes)

"""
import matplotlib.pyplot as plt
d = distance_matrix(clusters_df)
data = []
for x in d:
    data.append(x.flatten())
data = np.concatenate(data)
data = data[data > 0.0]
print(len(data))

plt.hist(data, bins=1000)
plt.show(); exit()
"""

# Calculate distance matrices and bonding graphs
clusters_df[("xyz", "distances")] = distance_matrix(clusters_df)

subcl_df = pd.DataFrame()
k = 0
# Filter through each cluster type
for k0 in range(len(cluster_subsets)):
    unique_cl = clusters_df[clusters_df['info', 'cluster_type']==unique_cluster_types[k0]].iloc[0]
    
    components = []
    for i, c in enumerate(unique_cl['info', 'components']):
        ratio = unique_cl['info', 'component_ratio'][i]
        components += [c]*ratio


    if print_minimums:
        subset = clusters_df.loc[cluster_subsets[k0]]
        names = subset[("info", "file_basename")]
        clusters = subset[("xyz", "structure")]
        distances = subset[("xyz", "distances")]
        energies = subset[('log', 'electronic_energy')]    
        x = []
        y = []
        for cl in subset.index.values:
            n = len(components)
            s1 = 0
            min_distances=[]
            for m1 in range(n):
                e1 = s1 + mol_df.at[components[m1], 'size']
                s2 = e1

                for m2 in range(m1+1, n):
                    e2 = s2 + mol_df.at[components[m2], 'size']
                    AB = distances[cl][s1:e1,s2:e2]
                    i1,i2 = np.unravel_index(np.argmin(AB, axis=None), AB.shape)
                    t1 = mol_df.at[components[m1], 'xyz']['atom'][i1+1]
                    t2 = mol_df.at[components[m2], 'xyz']['atom'][i2+1]
                    bond_type = t1+'-'+t2
                    # May not work if mol size > 2
                    if bond_type != 'O-H' and bond_type != 'H-O':
                        print(names[cl], AB[i1,i2], energies[cl], t1+'-'+t2)
                    s2 = e2
                s1 = e1
        continue

    subset = clusters_df.loc[cluster_subsets[k0]]
    if Hbonding == None:
        # Filter out incorrectly bonded structures (where a chemical reaction has occured)
        filter = test_internal_bonds(subset, mol_df, components, internal_limit)
        cluster_subsets[k0] = cluster_subsets[k0][filter]
        # Filter out clusters where molecules are too far apart
        filter, sc = test_clustering(subset, mol_df, components, [min_dist, max_dist], subcl_size)
    else:
        mol_df = get_acceptors_and_donors(mol_df, all_oxygens=True)
        # Filter out clusters with too few H-bonds (doesn't select subclusters)
        filter, counts = test_Hbonds(subset, mol_df, components, internal_limit, num=Hbonding)
        cluster_subsets[k0] = cluster_subsets[k0][filter]
        #sc = subclusters(counts, components)
        continue

    ######### 
    # SUBCLUSTERS: rewrite cluster_type, components, ...
    sc = sc[sc['subsets'].isna() == False]
    df = clusters_df.loc[sc.index.values]
    for idx in sc.index.values:
        
        n = 0
        for set in sc.at[idx, 'subsets']:
            mol = components[set]
            components = np.unique(mol)
            component_ratio = [np.count_nonzero(mol == x) for x in components]
            cluster_type = ""
            for c in range(len(components)):
                cluster_type += str(component_ratio[c])+components[c]
            # Filter out unwanted cluster types
            if rejected_cluster_type(components, component_ratio):
                continue

            subcl_df = pd.concat([subcl_df, df.loc[[idx]]], reject_index=True)
            subcl_df.at[k, ('info', 'cluster_type')] = cluster_type
            
            subcl_df.at[k, ('info', 'components')] = components
            subcl_df.at[k, ('info', 'component_ratio')] = component_ratio

            basename = cluster_type+'-'+df.loc[idx][('info', 'file_basename')].split('-')[1]+'_'+str(n)
            subcl_df.at[k, ('info', 'file_basename')] = basename
            subcl_df.at[k, ('info', 'folder_path')] = df.loc[idx][('info','folder_path')]
            
            s = 0
            atoms = None
            for c in range(set[-1]+1):
                e = s + mol_df.at[components[c], 'size']
                if c in set:
                    if atoms == None:
                        atoms = df.at[idx, ('xyz', 'structure')][s:e]
                    else:
                        atoms = atoms + df.at[idx, ('xyz', 'structure')][s:e]
                s = e
            subcl_df.at[k, ('xyz', 'structure')] = atoms

            n += 1
            k += 1
    #########
"""
plt.ylim(ymax=180)
plt.ylabel('Bond angle (deg)')
plt.xlabel('Bond distance (Å)')
plt.legend()
plt.show()
    """
if print_minimums:
    exit()

t = time() - t

# Saving data
f = open(file_out, 'w')

len_df = 0
for k0 in range(len(cluster_subsets)):
    if len(cluster_subsets[k0]) == 0:
        continue
    len_df += len(cluster_subsets[k0])
    df = clusters_df.loc[cluster_subsets[k0]]
    ct = df[('info', 'file_basename')].values[0].split('-')[0]
    for i in range(len(lines)):
        if paths.at[i, 'cluster'].split('-')[0] == ct:
            if paths.at[i, 'cluster'] in df[('info', 'file_basename')].values:
                f.write(lines[i])
f.close()

print("Filtered", len_df, "clusters")
print("Time taken:", t, "seconds")

# Saving found subclusters to a pickle file
subcl_df.index = [str(j) for j in range(len(subcl_df))]

if len(subcl_df) > 0:
    file_subcl = file_in[:-4]+'_SUBCL'
    pd.to_pickle(subcl_df, file_subcl+'.pkl')
    print("Saved", len(subcl_df), "subclusters to "+file_subcl+'.pkl')
