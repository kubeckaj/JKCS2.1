import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # silences warning from pd.read_table

import numpy as np
import pandas as pd

# Look for paths to xyz-bonding files
def bonding_paths(file_in="parameters.txt"):
    
    if (file_in=="parameters.txt") or (file_in=="input.txt"):
        d = {'name': [], 'q': [], 'path': []}
        # Read parameters.txt
        with open(file_in) as f:
            while not f.readline().startswith("# name"):
                continue
            while True:
                line = f.readline()
                try:
                    name, q, path = line.split()
                    d['name'].append(name)
                    d['q'].append(int(q))
                    d['path'].append(path)
                except:
                    break
    elif file_in=='calc.inp':
        d = {'path': [], 'n': []}
        # Read from calc.inp
        with open(file_in) as f:
            while not f.readline().startswith("components"):
                continue
            line = f.readline()
            while not line.startswith("end"):
                path, n = line.split()
                d['path'].append(path)
                d['n'].append(int(n))

                for i in range(4):
                    line = f.readline()
    else: 
        print('Error: No valid input given!'); exit()

    if len(d['path']) < 1: 
        print('Error: No bonding files found in '+file_in);exit()

    return d

def get_SMILES(mol_df):
    n = len(mol_df)
    smiles = []
    for i in range(n):
        path = mol_df.iloc[i]['path']
        with open(path) as f:
            n_atoms = int(f.readline().strip())
            line = f.readline().split()
            smiles.append(line[-1])

    mol_df['SMILES'] = smiles
    return mol_df

def read_molecule_data(d):

    mol_df = pd.DataFrame(data=d)
    if 'name' in mol_df.columns:
        mol_df = mol_df.drop_duplicates(subset=['name'], keep='first').dropna()
        mol_df = mol_df.set_index(['name'])
    n = len(mol_df)

    sizes = np.zeros(n, dtype=int)
    xyzs = np.zeros(n, dtype=pd.DataFrame)
    bonds = np.zeros(n, dtype=pd.DataFrame)

    for i in range(n):
        path = mol_df.iloc[i]['path']
        if not(os.path.isfile(path)): print("File "+path+" not found.");exit()

        n_atoms = int(pd.read_table(path, nrows=0).columns[0])

        # read xyz data
        xyz_df = pd.read_table(path, skiprows=2, delim_whitespace=True,
                                names=['atom', 'x', 'y', 'z'], nrows=n_atoms)
        xyz_df.index += 1
        xyzs[i] = xyz_df
        sizes[i] = len(xyz_df)

        # read bonding data
        bond_df = pd.read_table(path, skiprows=2+n_atoms, index_col=0, header=None, names=list(range(9)),
                                delim_whitespace=True, nrows=n_atoms, usecols=[0,1,3,5,7], engine='python')
        bond_df = bond_df.fillna(0).astype('Int64')
        bond_types = pd.read_table(path, skiprows=2+n_atoms, index_col=0, header=None, names=list(range(9)),
                                    delim_whitespace=True, nrows=n_atoms, usecols=[0,2,4,6,8], engine='python')

        d = {'i1': [], 'i2': [], 'a1': [], 'a2': [], 'type': [], 'length': []}
        # compute bond distances
        for i1 in bond_df.index:
            for j in bond_df.columns:
                i2 = bond_df.at[i1,j]
                if (i2 > 0):
                    d['i1'].append(i1)
                    d['i2'].append(i2)
                    d['a1'].append(xyz_df.loc[i1]['atom'])
                    d['a2'].append(xyz_df.loc[i2]['atom'])
                    d['type'].append(int(bond_types.at[i1,j+1]))
                    r = xyz_df.loc[i1][['x','y','z']] - xyz_df.loc[i2][['x','y','z']]
                    d['length'].append(np.sqrt(r[0]**2 + r[1]**2 + r[2]**2))

        bonds[i] = pd.DataFrame(data=d)
    
    mol_df['size'] = sizes
    mol_df['xyz'] = xyzs
    mol_df['bonds'] = bonds

    return mol_df

def read_pickled_data(file_in, input_pkl=[]): 

    clusters_df = []
    # Read pickle file(s)
    for i in range(len(input_pkl)):
        if not(os.path.isfile(input_pkl[i])): print("File "+input_pkl[i]+" not found. Make sure you are in the correct folder");exit();
        
        newclusters_df = pd.read_pickle(input_pkl[i])
        newclusters_df = newclusters_df[newclusters_df['info']['file_basename'] != "No"]

        if i == 0: 
            newclusters_df.index = [str(j) for j in range(len(newclusters_df))]
            clusters_df = newclusters_df
        else:
            len_clusters_df = len(clusters_df)
            newclusters_df.index = [str(j+len_clusters_df) for j in range(len(newclusters_df))]
            from pandas import concat
            clusters_df = concat([clusters_df, newclusters_df.copy()], ignore_index=True)
            #clusters_df = clusters_df.append(newclusters_df)

    # Read data from .dat files
    lines = []
    clusters_df2 = []
    for f in file_in:
        input_pkl = []
        lines = np.append(lines, open(f).readlines())
        lines = np.unique([l.split()[0] for l in lines]) # Drop duplicate files
        
        paths = pd.DataFrame([l.split('/:EXTRACT:/') for l in lines], columns=['file', 'cluster'])
        input_pkl = pd.unique(paths['file'])
        
        # Read pickle file(s)
        if len(input_pkl) == 0:
            print("Error: no pickle files found", file_in)
        else:
            for i in range(len(input_pkl)):
                newclusters_df = pd.read_pickle(input_pkl[i])
                p = paths[paths["file"] == input_pkl[i]]
                filtered = [n in p['cluster'].values for n in newclusters_df['info']['file_basename'].values]
                newclusters_df = newclusters_df[filtered]

                if i == 0: 
                    newclusters_df.index = [str(j) for j in range(len(newclusters_df))]
                    clusters_df2 = newclusters_df
                else:
                    len_clusters_df = len(clusters_df2)
                    newclusters_df.index = [str(j+len_clusters_df) for j in range(len(newclusters_df))]
                    clusters_df2 = clusters_df2.append(newclusters_df)

    if len(clusters_df2) > 0:
        clusters_df2 = clusters_df2.drop_duplicates(subset=[("info", "file_basename"), 
                                                            ("info", "folder_path")], keep='first')

        len_clusters_df = len(clusters_df)
        clusters_df2.index = [str(j+len_clusters_df) for j in range(len(clusters_df2))]

        if len_clusters_df == 0:
            clusters_df = clusters_df2
        else:
            from pandas import concat
            clusters_df = concat([clusters_df, clusters_df2.copy()], ignore_index=True)
            #clusters_df = clusters_df.append(clusters_df2)
            
    return clusters_df

def distance_matrix(clusters_df):

    distances = []
    for atoms in clusters_df[("xyz", "structure")].values:
        # distance matrix
        try:
            d = atoms.get_all_distances()
            distances.append(d)
        except:
            distances.append(np.nan) 

    return distances


def get_acceptors_and_donors(mol_df, all_oxygens=False):
    bonds = mol_df['bonds']
    
    def add_to_dict(d, k, v):
        if k not in d.keys():
            if v == None:
                d[k] = []
                return
            d[k] = [v]
        elif (v != None) and (v not in d[k]):
            d[k].append(v)
        
    # Create lists of donors and acceptors to pair
    donors = [] 
    acceptors = []
    
    for mol in bonds:
        
        donors.append({})
        acceptors.append({})

        select_H = (mol['a1'] == 'H') | (mol['a2'] == 'H')
        select_O = (mol['a1'] == 'O') | (mol['a2'] == 'O')
        select_C = (mol['a1'] == 'C') | (mol['a2'] == 'C') 
        select_N = (mol['a1'] == 'N') | (mol['a2'] == 'N')
        select_S = (mol['a1'] == 'S') | (mol['a2'] == 'S')

        if not all_oxygens:
            # drop weaker oxygens
            select_C = select_C & (mol ['type'] == 2) # C=O
            select_N = select_N & (select_H | (mol ['type'] == 2)) # N-H or N=O
            
        # bonds with partial charges
        oh = mol[(select_O | select_N) & select_H] # O/N-H
        ox = mol[select_O & (select_C | select_S | select_N)] # O-C/N/S
        
        for row in oh.values:
            i1, i2, a1, a2, _, _ = row
            if a1 == 'H':
                add_to_dict(donors[-1], i1, i2)
                add_to_dict(acceptors[-1], i2, i1)
            else:
                add_to_dict(donors[-1], i2, i1)
                add_to_dict(acceptors[-1], i1, i2)

        for row in ox.values:
            i1, i2, a1, a2, _, _ = row

            if a1 == 'N' or a2 == 'N':
                if a1 == 'O':
                    add_to_dict(acceptors[-1], i1, i2)
                    add_to_dict(donors[-1], i2, i1)
                else:
                    add_to_dict(acceptors[-1], i2, i1)
                    add_to_dict(donors[-1], i1, i2)
            else:
                if a1 == 'O':
                    add_to_dict(acceptors[-1], i1, None)
                else:
                    add_to_dict(acceptors[-1], i2, None)

    mol_df['acceptors'] = acceptors
    mol_df['donors'] = donors
    
    return mol_df


bond_limits = pd.DataFrame(data={'D': ['H', 'H', 'N'], 'A': ['O', 'N', 'O'], 
                                'bond': [(1.4, 2.2), (1.4, 2.2), (2.5, 3.2)],
                                'angle': [(110, 180), (110, 180), (60, 120)]})
bond_limits = bond_limits.set_index(['D','A'])

# Returns a dictionary of general cluster properties using mol_df
def construct_cluster(mol_df, components, remove_symmetric=False):
    
    cluster_size = len(components)
    cluster = []
    bonds = []
    donors = [] # mostly hydrogens (sometimes N) bonded to O 
    acceptors = [] # mostly oxygens (C-O-H, C-OO-H, C-O-C, C-O-O-C, C=O, NO2)
    
    s = 0
    for i in range(len(components)):
        mol = mol_df.loc[components[i]]
        cluster.append(mol['xyz'].copy())
        cluster[-1]['mol'] = i
        bonds.append(mol['bonds'].copy())
        
        # shift indeces
        cluster[-1].index += s

        if remove_symmetric:
            cluster[-1]['rank'].index += s

        bonds[-1]['i1'] += s
        bonds[-1]['i2'] += s
        
        donors.append({ k+s: np.array(v)+s for k,v in mol['donors'].items()})
        acceptors.append({ k+s: np.array(v)+s for k,v in mol['acceptors'].items()})
        s += mol['size']
        
    if cluster_size > 1:
        cluster = pd.concat(cluster)
    else:
        cluster = cluster[0]

    cluster = {'atoms': cluster, 'bonds': bonds, 'donors': donors, 'acceptors': acceptors, 'limits': bond_limits}
    
    return cluster

def get_bonding(clusters_df, mol_df, components):
    if 'distances' in clusters_df['xyz'].keys():
        distances = clusters_df[("xyz", "distances")]
    else:
        distances = distance_matrix(clusters_df)

    n_cl = len(clusters_df)
    cluster_info = construct_cluster(mol_df, components)
    xyz = clusters_df[("xyz", "structure")]
    atoms = cluster_info['atoms'][['atom','mol']]

    don = np.concatenate([list(x.keys()) for x in cluster_info['donors']])
    acc = np.concatenate([list(x.keys()) for x in cluster_info['acceptors']])
    bond_limits = cluster_info['limits']
    
    DA = {}
    for x in cluster_info['donors']:
        DA.update(x)

    # list of donor-acceptor pairs to test
    pairs = []
    for i1 in don:
        for i2 in acc:
            if i2 in DA[i1]:
                continue
            pairs.append([i1,i2])
    pairs = np.array(pairs)
            
    internal = pd.concat(cluster_info['bonds'])
    i1 = internal['i1'].values-1
    i2 = internal['i2'].values-1
    # bonded indices
    used_idx = internal[['i1','i2']].to_numpy()-1
    used_idx = np.concatenate((used_idx, pairs), axis=0)
    used_idx = set(map(lambda x: tuple(x), np.sort(used_idx, axis=1)))

    # non-bonded indices
    other = set(zip(*np.triu_indices(len(atoms), 1)))
    other = np.array(list(other - used_idx))

    isH = (atoms['atom'].values=='H')[other]
    isH = np.any(isH, axis=1)
    otherH = other[isH]
    otherX = other[~isH]
    
    bond_lists = []
    bond_counts = []
    
    for cl in clusters_df.index.values:
        dist = distances[cl]
        bonds = internal.copy()
        bonds['d'] = dist[i1,i2]
        bonds['angle'] = np.nan
        
        n = len(components)
        counts=np.zeros((n,n)) # H-bond counts between each molecule
        
        for j1,j2 in pairs:
            a1,m1 = atoms.loc[j1]
            a2,m2 = atoms.loc[j2]
            
            d = dist[j1-1,j2-1]
            l,u = bond_limits.at[(a1,a2), 'bond']
            if d > u:
                continue
            if d < l:
                # unexpected covalent bond
                bonds.loc[len(bonds)] = j1,j2,a1,a2,np.nan,np.nan,d,np.nan
                break

            l,u = bond_limits.at[(a1,a2), 'angle']
            for j3 in DA[j1]:
                angle = xyz[cl].get_angle(j3-1, j1-1, j2-1)
                if angle < l or angle > u:
                    continue
            
            counts[m1,m2] += 1 # number of bonds between m1 and m2

            bonds.loc[len(bonds)] = j1,j2,a1,a2,'H',np.nan,d,angle

        # look for reactions (unexpected bonds) 
        for j1,j2 in otherX:
            d = dist[j1,j2]
            if d < 2.0:
                a1,m1 = atoms.loc[j1+1]
                a2,m2 = atoms.loc[j2+1]
                print(j1,j2,a1,a2,d)
                bonds.loc[len(bonds)] = j1+1,j2+1,a1,a2,'O',np.nan,d,np.nan
                break
        for j1,j2 in otherH:
            d = dist[j1,j2]
            if d < 1.4:
                a1,m1 = atoms.loc[j1+1]
                a2,m2 = atoms.loc[j2+1]
                print(j1,j2,a1,a2,d)
                bonds.loc[len(bonds)] = j1+1,j2+1,a1,a2,'O',np.nan,d,np.nan
                break
                
        bond_lists.append(bonds)
        bond_counts.append(counts)
        
    return bond_lists, bond_counts

def test_Hbonds(clusters_df, mol_df, components, rel_tol, num=None):

    n_cl = len(clusters_df)
    cl_size = len(components)
    cluster_info = construct_cluster(mol_df, components)
    distances = clusters_df[("xyz", "distances")]
    xyz = clusters_df[("xyz", "structure")]
    atoms = cluster_info['atoms']['atom']
    mols = cluster_info['atoms']['mol'].to_numpy()

    don = np.concatenate([list(x.keys()) for x in cluster_info['donors']])
    acc = np.concatenate([list(x.keys()) for x in cluster_info['acceptors']])
    
    DA = {}
    for x in cluster_info['donors']:
        DA.update(x)

    # list of donor-acceptor pairs to test
    p = {k: [] for k in bond_limits.index.values}
    used_idx = []
    for i1 in don:
        for i2 in acc:
            if i2 in DA[i1]:
                continue
            t1, t2 = (atoms[[i1, i2]])
            p[(t1,t2)].append([i1,i2])
            used_idx.append([i1,i2])
    
    pairs = {}
    triples = {}
    for k,v in p.items():
        if len(v) == 0:
            continue
        v = np.array(v).T
        pairs[k] = (v[0]-1, v[1]-1)
        triples[k] = np.array([DA[i]-1 for i in v[0]])

    internal = pd.concat(cluster_info['bonds'])
    i1 = internal['i1'].values-1
    i2 = internal['i2'].values-1
    int_length = internal['length'].values # reference covalent bond lengths

    # bonded indices
    used_idx = np.concatenate((used_idx, internal[['i1','i2']].to_numpy()), axis=0)-1
    used_idx = set(map(lambda x: tuple(x), np.sort(used_idx, axis=1)))
    
    # non-bonded indices
    other = set(zip(*np.triu_indices(len(atoms), 1)))
    other = np.array(list(other - used_idx))

    isH = np.any((atoms.values=='H')[other], axis=1) # which bonds contain H
    otherH = other[isH].T;  otherH = (otherH[0], otherH[1]) # hydrogens may be closer
    otherX = other[~isH].T; otherX = (otherX[0], otherX[1])

    # H-bond counts between each molecule
    counts = np.zeros((n_cl,cl_size,cl_size), dtype=int)
    
    for cl, dist in enumerate(distances):
        
        # test relative errors in intermolecular distances
        if not np.allclose(dist[i1,i2], int_length, rtol=rel_tol):
            continue

        # look for reactions (unexpected bonds)
        if np.any(dist[otherH] < 1.4):
            continue
        if np.any(dist[otherX] < 2.0):
            continue
        
        for a1,a2 in pairs.keys():
            l,u = bond_limits.at[(a1,a2), 'bond']
            j1,j2 = pairs[(a1,a2)]
            d = dist[j1,j2] # H-bonded distances
            
            if np.any(d < l):
                break

            select = d < u
            if not np.any(select):
                continue
            
            j1 = j1[select]; j2 = j2[select]
            j3s = triples[(a1,a2)][select]
            
            l,u = bond_limits.at[(a1,a2), 'angle']
            
            n = len(j1)
            angles = np.array([xyz[cl].get_angles([[j3,j1[i],j2[i]] for j3 in j3s[i]]) 
                                for i in range(n)])
            select = np.all((angles > l) & (angles < u), axis=1)
            
            j1 = j1[select]; j2 = j2[select]
            m1,m2 = mols[[j1,j2]]
            counts[cl,m1,m2] += 1 # number of bonds between m1 and m2
            
    # total number of H-bonds per cluster
    totals = np.sum(counts,axis=(1,2))
    
    if num == None:
        passed = np.array(totals) > 0
    else:
        # select highest counts
        min_count = max(np.max(totals) - num, 0)
        passed = np.array(totals) > min_count
        
    return passed, counts

# Test that molecules have not reacted internally
def test_internal_bonds(clusters_df, mol_df, components, rel_tol):
    
    if 'distances' in clusters_df['xyz'].keys():
        distances = clusters_df[("xyz", "distances")]
    else:
        distances = distance_matrix(clusters_df)
    
    passed = np.ones(len(distances), dtype=bool)
    for cl in range(len(distances)):
        # read cluster xyz data
        
        s = 0
        for mol in components:
            n = mol_df.at[mol, 'size']
            # reference bond lengths
            bond_df = mol_df.at[mol, 'bonds']
            
            a1 = bond_df['i1'] + s-1
            a2 = bond_df['i2'] + s-1
            # precalculated internal distances
            dist =  distances[cl][a1,a2] 

            # compute relative errors in bond lengths
            passed[cl] = np.allclose(dist, bond_df['length'], rel_tol)
            s += n

    return passed

# Test which molecules are properly clustered.
# Also returns any found subclusters.
def test_clustering(clusters_df, mol_df, components, limits=None, subcl_size=2):
    
    structures = clusters_df[("xyz", "structure")]
    if 'distances' in clusters_df['xyz'].keys():
        distances = clusters_df[("xyz", "distances")]
    else:
        distances = distance_matrix(clusters_df)
    
    if limits == None:
        min_dist, max_dist = 1.4, 3.0 # H s
    else:
        min_dist, max_dist = limits

    n_cl = len(clusters_df)
    passed = np.zeros(n_cl, dtype=bool)
    splitted = pd.DataFrame(columns=['subsets'], index=clusters_df.index.values)

    for cl in range(n_cl):
        n = len(components)
        s1 = 0
        
        # arrays for sub-clusters
        subcl = np.split(np.arange(n), n) # [[mol1], [mol2], ..., [molN]]
        
        for m1 in range(n):
            e1 = s1 + mol_df.at[components[m1], 'size']
            s2 = e1

            idxa = np.where([m1 in sc for sc in subcl])[0][0]
            for m2 in range(m1+1, n):
                e2 = s2 + mol_df.at[components[m2], 'size']
                AB = distances[cl][s1:e1,s2:e2]
                d = np.min(AB)

                if d < min_dist: # reaction has occured!
                    subcl = None
                    break
                
                idxb = np.where([m2 in sc for sc in subcl])[0][0]
                if idxa == idxb:
                    s2 = e2
                    continue
                
                if d < max_dist: 
                    if idxa > idxb:
                        idxa, idxb = idxb, idxa
                    # merge subclusters
                    subcl[idxa] = np.append(subcl[idxa], subcl.pop(idxb))
                s2 = e2
            s1 = e1

            if subcl == None:
                break
            if len(subcl) == 1:
                passed[cl] = True
                break

        if passed[cl] == False and subcl != None:
            idx = structures.index.values[cl]
            # take subclusters containing >=subcl_size molecules
            i = 0
            while i < len(subcl):
                if len(subcl[i]) < subcl_size:
                    subcl.pop(i)
                else:
                    i += 1
            if len(subcl) > 0:
                for i in range(len(subcl)):
                    subcl[i] = np.sort(subcl[i])
                splitted.at[idx,'subsets'] = subcl
    
    return passed, splitted[passed==False]
