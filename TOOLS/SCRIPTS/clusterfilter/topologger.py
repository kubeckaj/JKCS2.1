import numpy as np
from pandas import DataFrame
from data_readers import MolInfo
from cluster_analysis import ClusterInfo
from typing import Iterable, Union

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

def clusters_to_smiles(clusters_df: DataFrame, 
                       cluster_info: ClusterInfo,
                       return_sorted: bool = False
    ) -> Union[list[str], DataFrame]:
    # make sure rdkit is available
    try:
        from rdkit.Chem import AddHs, RemoveHs, MolFromSmiles, MolToSmiles, CanonSmiles, SanitizeMol
        from rdkit.Chem.AllChem import ETKDGv3, EmbedMolecule
        from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D
        from rdkit.Geometry.rdGeometry import Point3D
    except:
        print('Error! rdkit not found. Halting...')
        exit()

    n = cluster_info.size
    n_atoms = cluster_info.component_sizes
    # SMILES patterns corresponding to the molecules
    smiles = cluster_info.smiles #[mol_df[k].smiles for k in components]
    
    # create initial Mol objects
    mols = []
    for smi in smiles:
        if smi.endswith('.mol'):
            from rdkit.Chem import MolFromMolFile, BondStereo, FindPotentialStereo
            mol = MolFromMolFile(smi, removeHs=False)

            for element in FindPotentialStereo(mol):
                if str(element.type) == "Bond_Double":
                    bond = mol.GetBondWithIdx(element.centeredOn)
                    bond.SetStereo(BondStereo.STEREOE)
                    bond.SetStereoAtoms(element.controllingAtoms[0],element.controllingAtoms[-1])
                    #mol.GetBondWithIdx(element.centeredOn).Specified()
                #print( f'chebi:57872 stereo from MOL: Type={element.type}, Which={element.centeredOn}, Specified={element.specified}, Descriptor={element.descriptor}, BondStereo={mol.GetBondWithIdx(element.centeredOn).GetStereo()}, Atoms={tuple(element.controllingAtoms)}' )
            
        elif smi.endswith('.sdf'):
            from rdkit.Chem import PandasTools
            mol = PandasTools.LoadSDF(smi, removeHs=False)['ROMol'].values[0]
        else: 
            mol = AddHs(MolFromSmiles(smi))
        mols.append(mol)

    params = ETKDGv3()
    for m in mols:
        EmbedMolecule(m, params)
        
    xyz = clusters_df[('xyz', 'structure')].values
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
            # we can assume that chirality has no effect on binding energy
            stereo_smiles1 = MolToSmiles(RemoveHs(m), canonical=False)
            #stereo_smiles1 = CanonSmiles(stereo_smiles1)
            for a in m.GetAtoms():
                a.InvertChirality()
            stereo_smiles2 = MolToSmiles(RemoveHs(m), canonical=False)
            #stereo_smiles2 = CanonSmiles(stereo_smiles2)
            
            stereo_smiles.append(sorted([stereo_smiles1, stereo_smiles2])[0])
            #print(stereo_smiles[-1])
            #stereo_smiles.append(stereo_smiles1)
            
    #print(np.unique(stereo_smiles, return_counts=True))
    # we can assume that chirality has no effect on binding energy
    #stereo_smiles = [CanonSmiles(smi.replace('@@','@').replace('[C@H]','C')).replace('[C@]','C') for smi in stereo_smiles]
    
    if n == 1:
        return stereo_smiles
    s = int(len(stereo_smiles)/n)
    if return_sorted:
        # smiles are sorted so that A.B == B.A
        cluster_smiles = ['.'.join(np.sort(smi)) for smi in np.reshape(stereo_smiles, (n, s)).T]
    else:
        cluster_smiles = ['.'.join(smi) for smi in np.reshape(stereo_smiles, (n, s)).T]

    return cluster_smiles
    
def tag_atoms(mol, offset=0):
    for i, atom in enumerate(mol.GetAtoms()):
        # Assign a unique ID based on atom position + an offset
        atom.SetAtomMapNum(i + offset)
    return mol

def generate_rdkit_cluster(smiles, tag_atoms: bool=False, return_parts: bool=False):
    from rdkit.Chem import AddHs, MolFromMolFile, MolFromSmiles, CombineMols
    from functools import reduce

    # redefining nitro groups so that Os are interchangeable
    smiles = [smi.replace('N(=O)(=O)', 'N([O])([O])')
                 .replace('[N+]([O-])=O', 'N([O])([O])') for smi in smiles]
    
    # generating molecules separately and adding hydrogens
    mols = []
    for smi in smiles:
        if smi.endswith('.mol'):
            mol = MolFromMolFile(smi, removeHs=False)
        elif smi.endswith('.sdf'):
            from rdkit.Chem import PandasTools
            mol = PandasTools.LoadSDF(smi, removeHs=False)['ROMol'].values[0]
        else: 
            mol = AddHs(MolFromSmiles(smi))

        mols.append(mol)

    if tag_atoms:
        offset = 0
        for mol in mols:
            # mark atoms in the molecule with IDs to map the atomic features
            mol = tag_atoms(mol, offset)
            offset += mol.GetNumAtoms() # update offset for next molecule

    # merging molecules to a cluster
    cluster = reduce(CombineMols, mols)

    if return_parts:
        return cluster, mols
    
    return cluster

def generate_rdkit_cluster_with_ids(smiles, return_parts: bool = False):
    return generate_rdkit_cluster(smiles, True, return_parts)


def add_Hbonds(cluster, don, acc, keepH=True, addId=False):
    from rdkit.Chem import RWMol, SanitizeMol, BondType, RemoveHs
    cl = RWMol(cluster) # make sure to copy the cluster

    i = 1
    # Add a hydrogen bond between each D-A pair
    for d, a in zip(don, acc):
        # AddBond doesn't work without int()-conversion
        cl.AddBond(int(d), int(a), BondType.HYDROGEN)
        if addId:
            cl.GetBondBetweenAtoms(int(d),int(a)).SetProp('bondNote', f'HB{i}')
            i += 1

    SanitizeMol(cl)
    return cl if keepH else RemoveHs(cl) # slow, removes non-Hbonded hydrogens 

"""def add_Hbonds_with_Hremoval(cluster, don, acc):
    from rdkit.Chem import RemoveHs
    cl = add_Hbonds(cluster, don, acc)
    cl = RemoveHs(cl)  # removes non-Hbonded hydrogens
    return cl
"""

def get_cluster_fingerprint(fpgen, cluster, don, acc):
    cluster = add_Hbonds(cluster, don, acc)
    
    # return fingerprint in binary format
    return str(fpgen.GetCountFingerprint(cluster).ToBinary())

# filter clusters with the same isomers given in cluster_info
def filter_isomers(clusters_df: DataFrame, cluster_info: ClusterInfo):
    if len(clusters_df) == 0:
        return []
    if cluster_info.isomers == None:
        return
    
    log = clusters_df['temp'][['SMILES']]
    log.index = np.arange(len(log))
    passed = np.repeat(True, len(clusters_df))
    for i, smiles in enumerate(log['SMILES'].values):
        for isomer_smiles, smi in zip(cluster_info.isomers, smiles.split('.')):
            if isomer_smiles is not None:
                passed[i] &= smi in isomer_smiles
        
    return passed        

def filter_isomorphs(clusters_df: DataFrame, cluster_info: ClusterInfo, 
                     used_DA_pairs: Union[Iterable, None ]= None, el=('log','electronic_energy')):
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
    
    cluster = generate_rdkit_cluster(cluster_info.smiles) # rdkit Mol object
    DA_pairs = clusters_df[('temp', 'Hbond_pairs')].values
    fps = [get_cluster_fingerprint(fpgen, cluster, don, acc) for (don, acc) in DA_pairs]
    
    log['TopFP'] = fps
    
    if isinstance(el, tuple):
        el = [el]
    
    log = log.sort_values(by=el)
    log = log.drop_duplicates(subset=['TopFP','SMILES'], keep='first')
    passed = np.repeat(False, len(clusters_df))
    passed[log.index.values] = True

    if isinstance(used_DA_pairs, Iterable) and len(used_DA_pairs) > 0:
        # TODO: should also check SMILES
        used_fps = [get_cluster_fingerprint(fpgen, cluster, don, acc)  for (don, acc,_,_) in used_DA_pairs]
        passed[log.index.values] &= ~log['TopFP'].isin(used_fps)
    
    return passed


def add_Hbonds_clusters(cf, clusters_df, keep_h=True):
    for i, ct in enumerate(cf.cluster_types):
        cluster_info = cf._cluster_info[i]

        cluster = generate_rdkit_cluster(cluster_info.smiles)
        subset = clusters_df[clusters_df[('info', 'cluster_type')] == ct]
        for idx in subset.index.values:
            don, acc = clusters_df.loc[idx, ('temp', 'Hbond_pairs')]
            mol = add_Hbonds(cluster, don, acc, keep_h)
            clusters_df.at[idx, ('temp', 'Mol')] = mol
            clusters_df.at[idx, ('temp', 'SMILES')] = cluster_info.smiles

    return clusters_df