import os
import numpy as np
import pandas as pd
from pandas import DataFrame

from dataclasses import dataclass
from typing import Iterable, Union

# Look for paths to xyz-bonding files
def bonding_paths(file_in="parameters.txt"):
    
    if file_in=='calc.inp':
        d = {'name': [], 'path': [], 'n': []}
        i = 0
        # Read from calc.inp
        with open(file_in) as f:
            while not f.readline().startswith("components"):
                continue
            line = f.readline()
            while not line.startswith("end"):
                d['name'].append(i); i +=1
                path, n = line.split()
                d['path'].append(path)
                d['n'].append(int(n))

                for j in range(4):
                    line = f.readline()

    elif file_in.endswith('.txt'):
        d = {'name': [], 'q': [], 'path': []}
        # Read parameters.txt
        with open(file_in) as f:
            while not f.readline().startswith("# name"):
                continue
            for line in f:
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
def get_SMILES(path):
    with open(path) as f:
        if not path.endswith(('.mol', '.sdf')):
            n_atoms = int(f.readline().strip())
        line = f.readline().split()
        smiles = line[-1]
    
    if smiles.endswith(('.mol', '.sdf')):
        if not os.path.isfile(smiles[-1]):
            folder = os.path.abspath(os.path.dirname(path))
            samepath = folder+'/'+smiles
            if os.path.isfile(samepath):
                smiles = samepath
            else:
                print('Error: could not find file '+smiles+' from '+folder);exit()
                
    return smiles

def read_bonding_data(path, xyz_df, smiles=None):
    n_atoms = len(xyz_df)
    d = {'i1': [], 'i2': [], 'a1': [], 'a2': [], 'type': [], 'length': []}
    
    if smiles and smiles.endswith(('.sdf', '.mol')):
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
    
    return DataFrame(data=d)

@dataclass
class MolInfo:
    name: str
    q: int # charge, unless using Constraints
    path: str
    size: int = 0
    _smiles: Union[str, None] = None
    _xyz: Union[pd.DataFrame, None] = None
    _bonds: Union[pd.DataFrame, None] = None
    _donors: Union[dict, None] = None
    _acceptors: Union[dict, None] = None
    isomers: Union[Iterable[str], None] = None

    @property
    def xyz(self) -> DataFrame:
        if self._xyz is None:
            if not(os.path.isfile(self.path)): print("File "+self.path+" not found."); return

            if self.path.endswith(('.sdf', '.mol')):
                from ase.io import read
                atoms = read(self.path)
                n_atoms = len(atoms)
                x,y,z = atoms.get_positions().T
                self._xyz = DataFrame(data={'atom': np.array(atoms.symbols), 'x': x, 'y': y, 'z': z})
            else:
                n_atoms = int(pd.read_table(self.path, nrows=0).columns[0])

                # read xyz data
                self._xyz = pd.read_table(self.path, skiprows=2, sep='\s+', names=['atom', 'x', 'y', 'z'], nrows=n_atoms)
            
            self._xyz.index += 1
            self.size = len(self.xyz)

        return self._xyz

    @property
    def bonds(self) -> DataFrame:
        if self._bonds is None:
            self._bonds = read_bonding_data(self.path, self.xyz, self.smiles)
        return self._bonds
    
    @property
    def smiles(self) -> str:
        if self._smiles is None:
            try:
                self._smiles = get_SMILES(self.path)
            except: 
                print('Warning! Failed to read SMILES')

        return self._smiles
    
    @property
    def donors(self):
        if self._donors == None:
            self.get_acceptors_and_donors()

        return self._donors
    
    @property
    def acceptors(self):
        if self._acceptors == None:
            self.get_acceptors_and_donors()

        return self._acceptors
    
    def get_acceptors_and_donors(self, all_oxygens: bool=True):
        
        def add_to_dict(d, k, v):
            if k not in d.keys():
                if v == None:
                    d[k] = []
                    return
                d[k] = [v]
            elif (v != None) and (v not in d[k]):
                d[k].append(v)
            
        # Create lists of donors and acceptors to pair
        mol = self.bonds
        donors = {}
        acceptors = {}

        select_H = (mol['a1'] == 'H') | (mol['a2'] == 'H')
        select_O = (mol['a1'] == 'O') | (mol['a2'] == 'O')
        select_C = (mol['a1'] == 'C') | (mol['a2'] == 'C') 
        select_N = (mol['a1'] == 'N') | (mol['a2'] == 'N')
        select_S = (mol['a1'] == 'S') | (mol['a2'] == 'S')

        if not all_oxygens:
            # drop weakly charged oxygens
            select_C = select_C & (mol ['type'] == 2) # C=O
            select_N = select_N & (select_H | (mol ['type'] == 2)) # N-H or N=O
            
        # bonds with partial charges
        oh = mol[(select_O | select_N) & select_H] # O/N-H
        ox = mol[select_O & (select_C | select_S | select_N)] # O-C/N/S
        
        for row in oh.values:
            i1, i2, a1, a2, _, _ = row
            if a1 == 'H':
                add_to_dict(donors, i1, i2)
                add_to_dict(acceptors, i2, i1)
            else:
                add_to_dict(donors, i2, i1)
                add_to_dict(acceptors, i1, i2)

        for row in ox.values:
            i1, i2, a1, a2, _, _ = row

            if a1 == 'N' or a2 == 'N':
                if a1 == 'O':
                    add_to_dict(acceptors, i1, i2)
                    add_to_dict(donors, i2, i1)
                else:
                    add_to_dict(acceptors, i2, i1)
                    add_to_dict(donors, i1, i2)
            else:
                if a1 == 'O':
                    add_to_dict(acceptors, i1, None)
                else:
                    add_to_dict(acceptors, i2, None)

        self._acceptors = acceptors
        self._donors = donors

        return 

# read xyz and bonding data from files
def read_molecule_data(mol_file: str) -> dict[str, MolInfo]:
    if mol_file.endswith(('.xyz','.sdf','.mol')):
        d={'name': [mol_file.split('/')[-1].split('.')[0]], 'q': [0], 'path': mol_file}
    else:
        # input.txt or parameter.txt 
        d = bonding_paths(mol_file)
    
    mol_df = DataFrame(data=d)
    if 'name' in mol_df.columns:
        mol_df = mol_df.drop_duplicates(subset=['name'], keep='first').dropna()
        mol_df = mol_df.set_index(['name'])
    
    mol_info = {}
    for k in mol_df.index:
        q = mol_df.loc[k]['q'] if 'q' in mol_df.columns else mol_df.loc[k]['n'] 
        path = mol_df.loc[k]['path']
        mol_info[k] = MolInfo(k, q, path)

    return mol_info

def read_xyz_data(xyz_files, noname=False) -> DataFrame:
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
        
    clusters_df = DataFrame(data=d)
    
    return clusters_df

def read_pickled_data(file_in, return_lines=False) -> DataFrame: 
    if isinstance(file_in, str):
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
    
    paths = DataFrame([l.split('/:EXTRACT:/') for l in lines], columns=['file', 'cluster'])
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
                clusters_df2 = pd.concat([clusters_df2, newclusters_df])

    if len(clusters_df2) > 0:
        clusters_df2 = clusters_df2.drop_duplicates(subset=[("info", "file_basename"), 
                                                            ("info", "folder_path")], keep='first')

        len_clusters_df = len(clusters_df)
        clusters_df2.index = [str(j+len_clusters_df) for j in range(len(clusters_df2))]

        if len_clusters_df == 0:
            clusters_df = clusters_df2
        else:
            clusters_df = pd.concat([clusters_df, clusters_df2])
    
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
