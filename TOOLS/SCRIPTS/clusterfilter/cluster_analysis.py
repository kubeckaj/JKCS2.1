import numpy as np
#import pandas as pd
from pandas import DataFrame, concat
from data_readers import MolInfo

def distance_matrix(clusters_df):
    from ase import Atoms
    distances = []
    for atoms in clusters_df[("xyz", "structure")].values:
        if isinstance(atoms, Atoms):
            distances.append(atoms.get_all_distances())
        else:
            distances.append(None)
    return distances

def select_lowest(clusters_df, n, by=('log','electronic_energy')):
    #smiles = clusters_df[('temp', 'SMILES')]

    if len(clusters_df) > n:
        passed = np.array([False]*len(clusters_df))
        idx = np.argsort(clusters_df[by].values)[:n]
        passed[idx] = True
    else: 
        passed = np.array([True]*len(clusters_df))
        
    return passed

def test_convergence(clusters_df, fmax, invert=False):
    fcol = ('extra', 'forces')
    passed = np.array([True]*len(clusters_df))

    if fcol not in clusters_df.columns:
        print('Error! Forces not found in dataframe. Skipping convergence test.')
        return passed
    
    from ase.units import Ha
    for i, forces in enumerate(clusters_df[fcol].values):
        try:
            passed[i] = (fmax > np.linalg.norm(forces, axis=1).max()*Ha)
        except:
            passed[i] = False
            
    if invert: 
        passed = ~passed

    return passed

def cut_relative(clusters_df, cutoff, by=('log','electronic_energy')):

    min_value = np.min(clusters_df[by])
    passed = clusters_df[by] <= min_value + cutoff
    return passed

def default_Hbond_limits() -> DataFrame:
    df = DataFrame(data={'D': ['H', 'H', 'N'], 'A': ['O', 'N', 'O'], 
                                    'length': [(1.4, 2.2), (1.4, 2.2), (2.5, 3.2)],
                                    'angle': [(120, 180), (120, 180), (80, 100)]})
    df = df.set_index(['D','A'])
    return df

#donor = "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]"
#acceptor = """[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"""
#OH_pair = "O-[H]"

from dataclasses import dataclass
from typing import Iterable, Union

@dataclass
class ClusterInfo:
    size: int
    charge: int
    components: Iterable[str] # 5sa = [sa-1,sa-1,sa,sa,sa]
    component_sizes: Iterable[int]
    atoms: DataFrame
    bonds: DataFrame
    donors: Iterable[dict]
    acceptors: Iterable[dict]
    smiles: Union[Iterable[str], None] = None
    isomers: Union[Iterable[str], None] = None
    
# Returns a cluster info dictionary of general cluster properties using mol_df
def construct_cluster(mol_df: dict[str, MolInfo], 
                      components: Iterable[str], 
    ) -> ClusterInfo:

    cluster_size = len(components)
    cluster = []
    bonds = []
    donors = [] # mostly hydrogens (sometimes N) bonded to O 
    acceptors = [] # mostly oxygens (C-O-H, C-OO-H, C-O-C, C-O-O-C, C=O, NO2)
    smiles = []

    start = 0
    total_charge = 0
    for i, k in enumerate(components):
        mol_info = mol_df[k]
        total_charge += mol_info.q
        cluster.append(mol_info.xyz.copy())
        cluster[-1]['mol'] = i
        bonds.append(mol_info.bonds.copy())
        
        # shift indeces
        cluster[-1].index += start
        if 'rank' in cluster[-1].columns:
            cluster[-1]['rank'] += start

        bonds[-1]['i1'] += start
        bonds[-1]['i2'] += start
        
        donors.append({ k+start: np.array(v)+start for k,v in mol_info.donors.items()})
        acceptors.append({ k+start: np.array(v)+start for k,v in mol_info.acceptors.items()})
        start += mol_info.size

        smiles.append(mol_info.smiles)
    
    cluster = concat(cluster) if cluster_size > 1 else cluster[0]
    
    for smi in smiles:
        if smi==None:
            smiles = None
            break

    cluster_info = ClusterInfo(len(components), total_charge, components, [mol_df[k].size for k in components],
                               cluster, bonds, donors,  acceptors, smiles)
    return cluster_info

def select_fraction(cluster_info: ClusterInfo, counts, reacted, options, frac=0.0):
    n_mol = cluster_info.size

    # filtering options
    if isinstance(options, str):
        if options[-1].isdigit():
            num = int(options)
        else:
            options = [int(options[:-1]), options[-1]]
    if isinstance(options, list):
        num, filter_type = options
        if isinstance(filter_type, str):
            filter_type = filter_type.lower()
            
    # total number of H-bonds per cluster
    totals = np.sum(counts,axis=(1,2))
    
    if n_mol == 1:
        # monomers
        passed = np.ones(len(totals), dtype=bool)
    else:
        di1, di2 = np.diag_indices(n_mol)
        # the number of inter-molecular H-bonds per cluster
        inner_counts = np.sum(counts[:, di1, di2], axis=1)
        # test if H-bonds are only inter-molecular
        passed = (totals > inner_counts)
        # TODO: might leave out separated monomers!
    
    if filter_type == 'x':
        min_count = 0
        passed = ~reacted
    elif num == None:
        passed &= np.array(totals) >= 1
    else:
        from functools import reduce
        DA = reduce(lambda x, y: x|y, cluster_info.donors)
        don = list(DA.keys())
        if filter_type == 'a':
            acc = np.concatenate([list(x.keys()) for x in cluster_info.acceptors])
            min_count = len(acc) - num
        elif filter_type == 'd':
            min_count = len(don) - num
        elif filter_type == 'h':
            atoms = cluster_info.atoms['atom']
            nH = np.sum(atoms[don]=='H')
            min_count = nH - num
        else:
            min_count = np.max(totals) - num
        min_count = max(1, min_count) # cluster should have at least 1 H-bond
        
        n_passable = sum(passed)
        filtered = passed & (np.array(totals) >= min_count)
        while min_count > 1:
            # select clusters with highest H-bond counts
            if sum(filtered)/n_passable > frac:
                break    
            min_count -= 1
            filtered = passed & (np.array(totals) >= min_count)
            
        passed = filtered

    return passed, totals, min_count

def test_Hbonds(clusters_df: DataFrame, cluster_info: ClusterInfo, rel_tol: float=0.2, 
                bond_limits=None, options=None, frac: float=0.0, return_stats=False):
    
    n_cl = len(clusters_df)
    n_mol = cluster_info.size
    
    distances = clusters_df[("temp", "distances")].values
    xyz = clusters_df[("xyz", "structure")].values
    atoms = cluster_info.atoms['atom']
    mols = cluster_info.atoms['mol'].to_numpy()
    
    from functools import reduce
    DA = reduce(lambda x, y: x|y, cluster_info.donors)
    don = list(DA.keys())
    acc = np.concatenate([list(x.keys()) for x in cluster_info.acceptors])
    
    # list of donor-acceptor pairs to test
    pairs = {k: ([], []) for k in bond_limits.index.values}
    used_idx = []
    for i1 in don:
        for i2 in acc:
            if i2 in DA[i1]:
                continue
            t1, t2 = (atoms[[i1, i2]])
            pairs[(t1,t2)][0].append(i1)
            pairs[(t1,t2)][1].append(i2)
            used_idx.append([i1,i2])
    
    # start indeces by -1 to match distance matrices
    pairs = {k: (np.array(v[0])-1, np.array(v[1])-1) for k,v in pairs.items() if len(v[0]) > 0 }
    # atoms each donor is bonded to (for angle calculations)
    triples = {k: np.array([DA[i+1]-1 for i in v[0]],dtype=object) for k,v in pairs.items()}
    
    internal = concat(cluster_info.bonds)
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
    counts = np.zeros((n_cl,n_mol,n_mol), dtype=int)
    reacted = np.repeat(False, n_cl)
    Hbonds = []

    for cl, dist in enumerate(distances):
        Hbonds.append([[],[],[],[]])
        
        # test relative errors in intermolecular distances
        if not np.allclose(dist[i1,i2], int_length, rtol=rel_tol):
            reacted[cl] = True
            continue
        
        # look for reactions (unexpected bonds)
        if np.any(dist[otherH] < 1.4) or np.any(dist[otherX] < 2.0):
            reacted[cl] = True
            continue
        
        for a1,a2 in pairs.keys():
            l,u = bond_limits.at[(a1,a2), 'length']
            j1,j2 = pairs[(a1,a2)]
            d = dist[j1,j2] # H-bonded distances
            
            if np.any(d < l):
                break

            select = d < u
            if not np.any(select):
                continue
            
            j1 = j1[select]; j2 = j2[select]
            j3s = triples[(a1,a2)][select]
            d = d[select]
            
            l,u = bond_limits.at[(a1,a2), 'angle']
            angles = [xyz[cl].get_angles([[j3,j1[i],j2[i]] for j3 in j3s[i]]) 
                        for i in range(len(j1))]
            select = [np.all((a > l) & (a < u)) for a in angles]
            if not np.any(select):
                continue

            j1 = j1[select]; j2 = j2[select]
            m1, m2 = mols[[j1,j2]]
            np.add.at(counts[cl], (m1, m2), 1) # number of bonds between m1 and m2
            
            Hbonds[-1][0] += list(j1)
            Hbonds[-1][1] += list(j2)
            Hbonds[-1][2] += list(d[select])
            Hbonds[-1][3] += [angles[i] for i in range(len(select)) if select[i]]
            
            
    passed, totals, min_count = select_fraction(cluster_info, counts, reacted, options, frac)

    if return_stats:
        uq, uc = np.unique(totals, return_counts=True)
        stat_counts = np.zeros(np.max(uq)+1,dtype=int)
        for i,u in enumerate(uq):
            stat_counts[u] += uc[i]

        strings = [clusters_df[('info','cluster_type')].values[0], stat_counts, min_count, '{:2.2%}'.format(np.sum(passed) / len(passed))]
    else:
        strings = None
    
    return passed, counts, Hbonds, strings

# Test that molecules have not reacted internally
def test_internal_bonds(clusters_df, components,  mol_df, rel_tol):
    
    n = len(clusters_df)
    distances = clusters_df[("temp", "distances")].values
    passed = np.array([True]*n)
    for cl in range(n):
        s = 0
        for mol in components:
            # reference bond lengths
            bond_df = mol_df[mol].bonds
            
            a1 = bond_df['i1'] + s-1
            a2 = bond_df['i2'] + s-1
            # precalculated internal distances
            dist =  distances[cl][a1,a2] 

            # compute relative errors in bond lengths
            passed[cl] = np.allclose(dist, bond_df['length'], rel_tol)
            s += mol_df[mol].size

    return passed

def subclusters():

    return

# Test which molecules are properly clustered.
# Also returns any found subclusters.
def test_clustering(clusters_df, components, mol_df , limits=None, subn_mol=2):
    n_cl = len(clusters_df)

    if len(components) < 2:
        return [True]*n_cl

    structures = clusters_df[("xyz", "structure")]
    distances = clusters_df[("temp", "distances")].values
    
    if limits == None:
        min_dist, max_dist = 1.4, 3.0 # H s
    elif len(limits) == 2:
        min_dist, max_dist = limits
    elif len(limits) == 3:
        min_dist, minO_dist, max_dist = limits

    passed = [False]*n_cl
    splitted = DataFrame(columns=['subsets'], index=clusters_df.index.values)

    for cl in range(n_cl):
        n = len(components)
        s1 = 0
        
        # arrays for sub-clusters
        subcl = np.split(np.arange(n), n) # [[mol1], [mol2], ..., [molN]]
        
        for m1 in range(n):
            e1 = s1 + mol_df[components[m1]].size
            s2 = e1

            idxa = np.where([m1 in sc for sc in subcl])[0][0]
            for m2 in range(m1+1, n):
                e2 = s2 + mol_df[components[m1]].size
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
            # take subclusters containing >=subn_mol molecules
            i = 0
            while i < len(subcl):
                if len(subcl[i]) < subn_mol:
                    subcl.pop(i)
                else:
                    i += 1
            if len(subcl) > 0:
                for i in range(len(subcl)):
                    subcl[i] = np.sort(subcl[i])
                splitted.at[idx,'subsets'] = subcl
    
    #return passed, splitted[passed==False]
    return passed
