import os
import numpy as np
import pandas as pd

# Look for paths to xyz-bonding files
def bonding_paths(file_in="parameters.txt"):
    
    if file_in=='calc.inp':
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
    elif file_in.endswith('.txt'):
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
    else: 
        print('Error: No valid input given!'); exit()

    if len(d['path']) < 1: 
        print('Error: No bonding files found in '+file_in);exit()

    return d

# Look for SMILES strings
def get_SMILES(mol_df):
    n = len(mol_df)
    smiles = []
    for i in range(n):
        path = mol_df.iloc[i]['path']
        with open(path) as f:
            if not path.startswith('.sdf'):
                n_atoms = int(f.readline().strip())
            line = f.readline().split()
            smiles.append(line[-1])
        
        if smiles[-1].endswith(('.mol', '.sdf')):
            if not os.path.isfile(smiles[-1]):
                samepath = os.path.abspath(os.path.dirname(path))+'/'+smiles[-1]
                if os.path.isfile(samepath):
                    smiles[-1] = samepath
                else:
                    print('Error: could not find file'+smiles[-1]+'from'+path);exit()
                    
    mol_df['SMILES'] = smiles
    return mol_df

def read_bonding_data(path, n_atoms, xyz_df, smiles=None):
    d = {'i1': [], 'i2': [], 'a1': [], 'a2': [], 'type': [], 'length': []}

    if smiles.endswith(('.sdf', '.mol')):
        path = smiles
        smiles = None

    if path.endswith(('.sdf', '.mol')):
        skip=0
        with open(path) as f:
            while "M" not in f.readline():
                continue
            skip = 1+len(f.readlines())
        bond_df = pd.read_table(path, skiprows=4+n_atoms, index_col=None, 
                                header=None, names=list(range(7)),sep='\s+', 
                                skipfooter=skip, engine='python')
    else:
        bond_df = pd.read_table(path, skiprows=2+n_atoms, index_col=None, 
                                header=None, names=list(range(9)),
                                sep='\s+', nrows=n_atoms, engine='python')
        
    if (len(bond_df) == 0) and smiles != None:
        from rdkit.Chem import MolFromSmiles, AddHs
        mol = AddHs(MolFromSmiles(smiles))
        for b in mol.GetBonds():
            i1 = b.GetBeginAtomIdx()+1
            i2 = b.GetEndAtomIdx()+1
            d['i1'].append(i1)
            d['i2'].append(i2)
            d['a1'].append(b.GetBeginAtom().GetSymbol())
            d['a2'].append(b.GetEndAtom().GetSymbol())
            d['type'].append(int(b.GetBondType()))
            r = xyz_df.loc[i1][['x','y','z']] - xyz_df.loc[i2][['x','y','z']]
            d['length'].append(np.sqrt(r.iloc[0]**2 + r.iloc[1]**2 + r.iloc[2]**2))
    else:
        if path.endswith(('.sdf','.mol')):
            bond_types = bond_df[[2]]
            bond_df = bond_df[[0,1]]
        else:
            bond_df = bond_df.fillna(0).astype('Int64')
            bond_types = bond_df[[2,4,6,8]]
            bond_df = bond_df[[0,1,3,5,7]]

        # compute bond distances
        for idx in bond_df.index.values:
            i1 = bond_df.at[idx,0]
            for j in bond_df.columns[1:]:
                i2 = bond_df.at[idx,j]
                if (i2 > 0):
                    d['i1'].append(i1)
                    d['i2'].append(i2)
                    d['a1'].append(xyz_df.loc[i1]['atom'])
                    d['a2'].append(xyz_df.loc[i2]['atom'])
                    d['type'].append(int(bond_types.at[idx,j+1]))
                    r = xyz_df.loc[i1][['x','y','z']] - xyz_df.loc[i2][['x','y','z']]
                    d['length'].append(np.sqrt(r.iloc[0]**2 + r.iloc[1]**2 + r.iloc[2]**2))
    
    return pd.DataFrame(data=d)

# read xyz and bonding data from files
def read_molecule_data(mol_file):
    if mol_file.endswith(('.xyz','.sdf')):
        d={'name': [mol_file.split('.')[0]], 'q':[0], 'path': mol_file}
    else:
        d = bonding_paths(mol_file)
    
    mol_df = pd.DataFrame(data=d)
    if 'name' in mol_df.columns:
        mol_df = mol_df.drop_duplicates(subset=['name'], keep='first').dropna()
        mol_df = mol_df.set_index(['name'])
    n = len(mol_df)
    
    try:
        get_SMILES(mol_df)
    except:
        print('Warning! Failed to read SMILES')
        
    sizes = np.zeros(n, dtype=int)
    xyzs = np.zeros(n, dtype=pd.DataFrame)
    bonds = np.zeros(n, dtype=pd.DataFrame)

    for i in range(n):
        path = mol_df.iloc[i]['path']
        if not(os.path.isfile(path)): print("File "+path+" not found."); continue

        if path.endswith('.sdf'):
            from ase.io import read
            atoms = read(path)
            n_atoms = len(atoms)
            x,y,z = atoms.get_positions().T
            xyz_df = pd.DataFrame(data={'atom': np.array(atoms.symbols), 'x': x, 'y':y, 'z':z})
            xyz_df.index += 1
            xyzs[i] = xyz_df
            sizes[i] = len(xyz_df)
        else:
            n_atoms = int(pd.read_table(path, nrows=0).columns[0])

            # read xyz data
            xyz_df = pd.read_table(path, skiprows=2, sep='\s+', names=['atom', 'x', 'y', 'z'], nrows=n_atoms)
            xyz_df.index += 1
            xyzs[i] = xyz_df
            sizes[i] = len(xyz_df)

        smi = mol_df.iloc[i]['SMILES'] if 'SMILES' in mol_df.columns else None
        # read bonding data
        bonds[i] = read_bonding_data(path, n_atoms, xyz_df, smi)
    
    mol_df['size'] = sizes
    mol_df['xyz'] = xyzs
    mol_df['bonds'] = bonds

    return mol_df

def read_xyz_data(xyz_files, noname=False) -> pd.DataFrame:
    from ase.io import read
    from re import split
    
    basenames = [f.split('/')[-1].split('.')[0] for f in xyz_files]
    cluster_types = [b.split('-')[0] for b in basenames]
    
    if noname:
        comp_ratio = np.ones(len(cluster_types))
        components = np.array(cluster_types)
    else:
        comp = [split(r"(\d+)", ct) for ct in cluster_types]
        comp_ratio = [list(np.array(c[1::2]).astype(int)) for c in comp]
        components = [c[2::2] for c in comp]
        
    atoms = [read(f, format='extxyz') for f in xyz_files]
    d = {('info', 'file_basename'): basenames,
            ('info', 'cluster_type'): cluster_types,
            ('info', 'components'): components, 
            ('info', 'component_ratio'): comp_ratio,
            ('xyz', 'structure'): atoms}
        
    clusters_df = pd.DataFrame(data=d)
    
    return clusters_df

def read_pickled_data(file_in, return_lines=False) -> pd.DataFrame: 
    if type(file_in) == str:
        file_in = [file_in]
        
    input_pkl=[f for f in file_in if f.endswith('.pkl')]
    file_in=[f for f in file_in if not f.endswith('.pkl')]
    
    clusters_df = []
    # Read pickle file(s)
    for i in range(len(input_pkl)):
        if not(os.path.isfile(input_pkl[i])): print("File "+input_pkl[i]+" not found. Make sure you are in the correct folder");exit();
        
        newclusters_df = pd.read_pickle(input_pkl[i])
        newclusters_df = newclusters_df[newclusters_df[('info', 'file_basename')] != "No"]

        if i == 0: 
            newclusters_df.index = [j for j in range(len(newclusters_df))]
            clusters_df = newclusters_df
        else:
            len_clusters_df = len(clusters_df)
            newclusters_df.index = [j+len_clusters_df for j in range(len(newclusters_df))]
            clusters_df = pd.concat([clusters_df, newclusters_df])

    if len(file_in) == 0:
        return clusters_df

    # Read data from .dat/.txt files
    lines = np.concatenate([open(f).readlines() for f in file_in])
    clusters_df2 = []
    input_pkl = []
    xyz_files = np.unique([l.split()[0] for l in lines if '.xyz' in l])
    lines = np.unique([l.split()[0] for l in lines if '/:EXTRACT:/' in l]) # Drop duplicates
    
    paths = pd.DataFrame([l.split('/:EXTRACT:/') for l in lines], columns=['file', 'cluster'])
    input_pkl = pd.unique(paths['file'])
    
    # Read pickle file(s)
    if len(input_pkl) == 0 and len(xyz_files) == 0:
        print("Error: no pickle files found", file_in)
    else:
        for i in range(len(input_pkl)):
            newclusters_df = pd.read_pickle(input_pkl[i])
            p = paths[paths["file"] == input_pkl[i]]
            filtered = newclusters_df[('info','file_basename')].isin(p['cluster'].values)
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
            clusters_df = clusters_df.append(clusters_df2)
    
    if len(xyz_files) > 0:
        newclusters_df = read_xyz_data(xyz_files)
        
        len_clusters_df = len(clusters_df)
        newclusters_df.index = [str(j+len_clusters_df) for j in range(len(newclusters_df))]
        if len_clusters_df == 0:
            clusters_df = newclusters_df
        else:
            clusters_df = clusters_df.append(newclusters_df)


    if return_lines: 
        return clusters_df, lines
    
    return clusters_df
