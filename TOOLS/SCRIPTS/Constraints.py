
import os, sys
import numpy as np
import pandas as pd
from copy import deepcopy
from random import uniform, sample

def printhelp():
  print('   Use this file like this:')
  print('           Constraints.py calc.inp {OPTIONS}')
  print(' OPTIONS:')
  print('  -all       include all oxygen acceptors for hydrogen bonding')
  print('  -nosymm    do not remove symmetric cases')
  print('  -nfiles    maximum number of files to create [default=100]')
  print('  ')

#######################################################################

bond_limit = 0.1    # largest allowed relative error
#mol_dist = 2.7      # largest H-X bond distance between molecules
print_minimums = False
use_all = False
remove_symmetric = True
max_files = 100

# Look for inputs
file_in = sys.argv[1]
last = ""

n_arg = len(sys.argv)
for i in range(2,n_arg):
    arg=sys.argv[i]
    if str(arg) == "--help":
        printhelp()
        exit()
    if str(arg) == "-help":
        printhelp()
        exit()
    if arg == "-all":
        use_all = True
        last = ""
        continue
    if arg == "-nosymm":
        remove_symmetric = False
        last = ""
        continue
    if arg == "-nfiles":
        last = "-nfiles"
        continue
    if last == "-nfiles":
        max_files = int(arg)
        last = ""
        continue
    print('Constraints.py: Unknown argument "'+arg+'" EXITING'); exit()
  

if not(os.path.isfile(file_in)): print("File "+file_in+" not found. Make sure you are in the correct folder");exit()

###############################################################################
from cluster_analysis import *
print('Constraints.py: Start.')

d = bonding_paths(file_in)
if np.sum(d['n']) < 2:
    print('Not enough molecules in cluster (<2)!');exit()

mol_df = read_molecule_data(d)
cluster_size = mol_df['n'].sum()
n_mol = len(mol_df)
print("Found", n_mol, "molecule(s) in calc.inp.")

if use_all or cluster_size > 2:
    remove_symmetric = False
if remove_symmetric:
    # make sure rdkit is available
    try:
        import rdkit.Chem as Chem
        from rdkit.Chem import CanonicalRankAtoms
    except:
        remove_symmetric = False
if remove_symmetric:
    for i in mol_df.index.values:
        mol_df = get_SMILES(mol_df)
        smiles = mol_df.at[i, 'SMILES']
        xyz_df = mol_df.at[i, 'xyz']
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            ranks = list(CanonicalRankAtoms(mol, breakTies=False)) # symmetry class for each atom
            xyz_df['rank'] = ranks
        except:
            xyz_df['rank'] = np.arange(len(xyz_df))

###############################################################################

mol_df = get_acceptors_and_donors(mol_df, use_all)

components = []
for i in mol_df.index.values:
    for n in range(mol_df.at[i, 'n']):
        components.append(i)
cluster_size = len(components)
        
# Construct the cluster
cluster = construct_cluster(mol_df, components, remove_symmetric)
atoms = cluster['atoms']
bonds = cluster['bonds']

donors = [list(x.keys()) for x in cluster['donors']] # mostly hydrogens (sometimes N) bonded to O 
acceptors = [list(x.keys()) for x in cluster['acceptors']] # mostly oxygens (C-O-H, C-OO-H, C-O-C, C-O-O-C, C=O, NO2)
# {A:D} dictionary for preventing overlaps of acceptor-donor pairs 
AD = cluster['acceptors'][0] 
for x in cluster['acceptors'][1:]:
    AD.update(x)

############################################################
# removes overlapping OH--OH pairs
def remove_overlaps(pairs):
    k = AD.keys()
    to_remove = []

    for h1, o1 in deepcopy(pairs):
        if o1 not in k:
            continue

        for h2,o2 in deepcopy(pairs):
            if o2 not in k:
                continue
            if ((h1,o1) in to_remove) or ((h2,o2) in to_remove):
                continue
            
            if (h2 in AD[o1]) and (h1 in AD[o2]):
                # randomly select one pair
                if uniform(0,1) > 0.5: 
                    to_remove.append((h1,o1))
                else:
                    to_remove.append((h2,o2))

    for x in to_remove:
        pairs.remove(x)

    return pairs

nH = np.sum([len(x) for x in donors])
nO = np.sum([len(x) for x in acceptors])

combinations = []
while len(combinations) < 1e4:
    O = deepcopy(acceptors)
    H = deepcopy(donors)
    pair_list = []
    for i in range(nH):
        
        a,b = sample(range(cluster_size), 2) # pick 2 molecules a,b
        if (len(H[a]) == 0) or (len(O[b]) == 0):
            continue

        hi = sample(H[a],1)[0]; H[a].remove(hi)
        oj = sample(O[b],1)[0]; O[b].remove(oj)
        pair_list.append((hi,oj))
        
    pair_list = sorted(pair_list)
    combinations.append(remove_overlaps(pair_list))

combinations = list(np.unique(np.array(combinations, dtype=object)))

###########################################################

len_count = np.array([len(c) for c in combinations])
len_max = np.max(len_count)
# Remove combinations with low number of O--H bonds
if len(combinations) > 2*max_files:

    new_combinations = []
    for c in combinations:
        if len(c) == len_max:
            new_combinations.append(c)
    combinations = new_combinations

combinations = sample(combinations, min(2*max_files, len(combinations)))

# Pick only one constraint per symmetry group
if remove_symmetric:
    select = (atoms['atom'] == 'O') | ([a in np.concatenate(donors) for a in atoms.index.values])
    atoms = atoms[select]
    uniq, counts = np.unique(atoms['rank'], return_counts=True)
    
    g = []
    for u in uniq:
        g.append(atoms[atoms['rank'] == u].index.values)
    groups = pd.DataFrame(data={'count': counts, 'group': g}, index=uniq)
    
    def swap(i):
        r = atoms.at[i, 'rank']
        if groups.at[r, 'count'] == 1:
            return i
        for j in groups.at[r, 'group']:
            if j != i:
                return j
    
    new_combinations = combinations.copy()
    for c1 in combinations:
        n+=1
        if len(c1) < len_max:
            new_combinations.remove(c1)

        if c1 not in new_combinations:
            continue
        # swapping all H's in the same group returns an identical combination
        c2 = sorted([(swap(h),o) for h,o in c1])
        if c1 == c2:
            continue
        if c2 not in new_combinations:
            continue

        if uniform(0,1) > 0.5: 
            new_combinations.remove(c1)
        else:
            new_combinations.remove(c2)

    combinations = new_combinations

# save at most max_files combinations 
combinations = sample(combinations, min(max_files, len(combinations)))

# Create constraints[i].inp files
for i in range(len(combinations)):
    file_out = 'constraints'+str(i)+'.inp'
    with open(file_out, 'w') as f:
        f.write("$constrain\n")
        f.write("  force constant=0.001\n")
        for a, b in combinations[i]:
            #distances = {'H': 1.85, 'N': 2.9}
            if atoms.at[a, 'atom'] == 'H':
                f.write("  distance: "+str(a)+", "+str(b)+", 1.85\n")
            else:
                f.write("  distance: "+str(a)+", "+str(b)+", 2.9\n")
        f.write("$end\n")

print('Created '+str(len(combinations))+' constraint files')