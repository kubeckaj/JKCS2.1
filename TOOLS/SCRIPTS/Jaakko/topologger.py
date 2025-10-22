import os,sys
import pandas as pd
import numpy as np
from typing import Iterable

def get_ranks(mol_df):
    # make sure rdkit is available
    try:
        from rdkit.Chem import AddHs, MolFromSmiles, MolFromMolFile, CanonicalRankAtoms
    except:
        return False

    for i in mol_df.index.values:
        smiles = mol_df.at[i, 'SMILES']
        xyz_df = mol_df.at[i, 'xyz']
        try:
            if smiles.endswith('.mol'):
                mol = MolFromMolFile(smiles, removeHs=False)
            else:
                smiles = smiles.replace('/','').replace('\\','').replace('[C@H]','C')
                mol = MolFromSmiles(smiles)
                mol = AddHs(mol)
            #print(Chem.MolToMolBlock(mol))
            ranks = list(CanonicalRankAtoms(mol, breakTies=False)) # symmetry class for each atom
            xyz_df['rank'] = ranks
        except:
            xyz_df['rank'] = np.arange(len(xyz_df))
            
    return True

def clusters_to_smiles(clusters_df, components, mol_df, return_isomers=False, return_sorted=True):
    # make sure rdkit is available
    try:
        from rdkit.Chem import AddHs, RemoveHs,MolFromMolFile, MolFromSmiles, MolToSmiles, CanonSmiles
        from rdkit.Chem.AllChem import ETKDGv3, EmbedMolecule
        from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D
        from rdkit.Geometry.rdGeometry import Point3D
    except:
        print('Error! rdkit not found. Halting...')
        exit()

    n = len(components)
    n_atoms = mol_df.loc[components]['size'].values
    # SMILES patterns corresponding to the molecules
    smiles = mol_df.loc[components]['SMILES'].values
    
    # create initial Mol objects
    mols = []
    for smi in smiles:
        if smi.endswith('.mol'):
            mol = MolFromMolFile(smi,removeHs=False)
        elif smi.endswith('.sdf'):
            from rdkit.Chem import PandasTools
            mol = PandasTools.LoadSDF(smi, removeHs=False)['ROMol'].values[0]
        else: 
            AddHs(MolFromSmiles(smi))
        mols.append(mol)

    params = ETKDGv3()
    for m in mols:
        EmbedMolecule(m, params)
        
    xyz = clusters_df[('xyz', 'structure')].values
    names = clusters_df[('info', 'file_basename')].values
    # split clusters into monomers = [[monA_1,...], [monB_1,...], ...]
    monomers = [[]]*n
    s = 0
    for i in range(n):
        e = s+n_atoms[i]
        monomers[i] = [atoms[s:e] for atoms in xyz]
        s = e
    
    stereo_smiles = []
    # determine stereochemistry of each monomer
    for i in range(n):
        m = mols[i]
        for atoms in monomers[i]:
            for idx, pos in enumerate(atoms.get_positions()):
                point = Point3D(pos[0],pos[1],pos[2])
                m.GetConformer().SetAtomPosition(idx, point)

            AssignStereochemistryFrom3D(m)
            stereo_smiles.append(MolToSmiles(RemoveHs(m)))
    #for s in np.unique(stereo_smiles, return_counts=True):
    #for a in mol.GetAtoms():
    #    a.InvertChirality()
    
    # we can assume that chirality has no effect on binding energy
    #stereo_smiles = [CanonSmiles(smi.replace('@@','@').replace('[C@H]','C')).replace('[C@]','C') for smi in stereo_smiles]

    if return_isomers:
        # create new dataframe
        names = list(names)*n
        components = np.repeat(components, len(clusters_df))
        xyz = monomers[0]
        for m in monomers[1:]:
            xyz += m

        df = pd.DataFrame(data={'parent': names, 'component': components, 'xyz': xyz, 'SMILES': stereo_smiles})
        return df
    else: 
        if n == 1:
            return stereo_smiles
        s = int(len(stereo_smiles)/n)
        if return_sorted:
            # smiles are sorted so that A.B == B.A
            cluster_smiles = ['.'.join(np.sort(smi)) for smi in np.reshape(stereo_smiles, (n, s)).T]
        else:
            cluster_smiles = ['.'.join(smi) for smi in np.reshape(stereo_smiles, (n, s)).T]

        return cluster_smiles
    
def generate_rdkit_cluster(smiles):
    from rdkit.Chem import AddHs, MolFromMolFile, MolFromSmiles, CombineMols
    from functools import reduce

    # redefining nitro groups so that Os are interchangeable
    smiles = [smi.replace('N(=O)(=O)', 'N([O])([O])') for smi in smiles]
    # generating molecules separately and adding hydrogens

    mols = []
    for smi in smiles:
        if smi.endswith('.mol'):
            mol = MolFromMolFile(smi,removeHs=False)
        elif smi.endswith('.sdf'):
            from rdkit.Chem import PandasTools
            mol = PandasTools.LoadSDF(smi, removeHs=False)['ROMol'].values[0]
        else: 
            AddHs(MolFromSmiles(smi))
        mols.append(mol)

    # merging molecules to a cluster
    cluster = reduce(CombineMols, mols)

    return cluster

def get_cluster_fingerprint(fpgen, cluster, don, acc, xyz=None):
    from rdkit.Chem import RWMol, SanitizeMol, BondType
    cluster = RWMol(cluster)

    # Add a hydrogen bond between each D-A pair
    for d, a in zip(don, acc):
        # AddBond doesn't work without int()-conversion
        cluster.AddBond(int(d),int(a),BondType.HYDROGEN)
        
    SanitizeMol(cluster)
    
    """try:
    #    SanitizeMol(cluster)
    #except:
    if acc[0] == 52:
        print('Error: SanitizeMol failed.')
        print(np.stack((don, acc), axis=-1))
        from rdkit.Chem import Draw
        d = Draw.rdMolDraw2D.MolDraw2DSVG(600, 400) # or MolDraw2DCairo to get PNGs
        d.drawOptions().annotationFontScale = 0.7
        d.drawOptions().addAtomIndices = True
        d.DrawMolecule(cluster)

        #from ase.io import write
        #write('fail.xyz', xyz)
        d.FinishDrawing()
        #d.WriteDrawingText('cluster_labelled.png')   
        svg = d.GetDrawingText().encode()
        with open('cluster_error.svg', 'wb') as file:
            file.write(svg)
        exit()"""
    
    # return fingerprint in binary format
    return str(fpgen.GetCountFingerprint(cluster).ToBinary())

# filter clusters with the same isomers given in cluster_info
def filter_isomers(clusters_df, cluster_info):
    if len(clusters_df) == 0:
        return []
    if 'isomers' not in cluster_info.keys():
        return
    
    log = clusters_df['temp'][['SMILES']]
    log.index = np.arange(len(log))
    passed = np.repeat(True, len(clusters_df))

    isomer_smiles = []
    for smi in  cluster_info['isomers']:
        isomer_smiles = np.append(isomer_smiles, smi)
    
    for i, smiles in enumerate(log['SMILES'].values):
        passed[i] &= np.all([smi in isomer_smiles for smi in smiles.split('.')])
        
    return passed        

def filter_isomorphs(clusters_df, cluster_info, used_DA_pairs=None, el=('log','electronic_energy')):
    if len(clusters_df) == 0:
        return []
        
    if el not in clusters_df.columns:
        clusters_df[el] = np.nan
        
    log, el = el
    log = clusters_df[log][[el]]
    log['SMILES'] = clusters_df['temp']['SMILES']
    log.index = np.arange(len(log))

    # Using Morgan/circular fingerprint
    from rdkit.Chem.AllChem import GetMorganGenerator
    fpgen = GetMorganGenerator(radius=20, fpSize=2048)

    smiles = cluster_info['SMILES']
    cluster = generate_rdkit_cluster(smiles) # rdkit Mol object
    DA_pairs = clusters_df[('temp', 'Hbond_pairs')].values
    fps = [get_cluster_fingerprint(fpgen, cluster, don, acc) for (don, acc) in DA_pairs]
    
    log['TopFP'] = fps

    log = log.sort_values(by=[el])
    log = log.drop_duplicates(subset=['TopFP','SMILES'], keep='first')
    
    passed = np.repeat(False, len(clusters_df))
    passed[log.index.values] = True

    if isinstance(used_DA_pairs, Iterable) and len(used_DA_pairs) > 0:
        # TODO: should also check SMILES
        used_fps = [get_cluster_fingerprint(fpgen, cluster, don, acc)  for (don, acc,_,_) in used_DA_pairs]
        passed[log.index.values] &= ~log['TopFP'].isin(used_fps)
    
    return passed


"""def filter_isomorphs_old(clusters_df, mol_df, cluster_info, DA_pairs, used_DA_pairs=None):
    log = clusters_df['log'][['electronic_energy']]
    stereo_smiles = clusters_to_smiles(clusters_df, cluster_info['components'], mol_df)
    log['SMILES'] = stereo_smiles
    log.index = np.arange(len(log))
    
    # Select only clusters with given isomers
    if 'isomers' in cluster_info.keys():
        isomer_smiles = np.concatenate(cluster_info['isomers'])
        for i, smiles in enumerate(log['SMILES'].values):
            passed[i] &= np.all([smi in isomer_smiles for smi in smiles.split('.')])
        
    # Using Morgan/circular fingerprint
    from rdkit.Chem.AllChem import GetMorganGenerator
    fpgen = GetMorganGenerator(radius=20, fpSize=2048)

    smiles = cluster_info['SMILES']
    cluster = generate_rdkit_cluster(smiles)
    fps = [get_cluster_fingerprint(fpgen, cluster, don, acc) for (don, acc) in DA_pairs]
    log['TopFP'] = fps

    log = log.sort_values(by=['electronic_energy'])
    log = log.drop_duplicates(subset=['TopFP','SMILES'], keep='first')

    passed = np.repeat(False, len(clusters_df))
    passed[log.index.values] = True

    if isinstance(used_DA_pairs, Iterable) and len(used_DA_pairs) > 0:
        # TODO: also check SMILES
        used_fps = [get_cluster_fingerprint(fpgen, smiles, don, acc)  for (don, acc) in used_DA_pairs]
        for fp in used_fps:
            passed[log.index.values] &= (log['TopFP'].values != fp)

    return passed
"""