import numpy as np
#import pandas as pd
from pandas import DataFrame, concat

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

def get_acceptors_and_donors(mol_df, all_oxygens=False) -> DataFrame:
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
            # drop weakly charged oxygens
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

def default_Hbond_limits() -> DataFrame:
    df = DataFrame(data={'D': ['H', 'H', 'N'], 'A': ['O', 'N', 'O'], 
                                    'length': [(1.4, 2.2), (1.4, 2.2), (2.5, 3.2)],
                                    'angle': [(120, 180), (120, 180), (80, 100)]})
    df = df.set_index(['D','A'])
    return df

#donor = "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]"
#acceptor = """[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"""
#OH_pair = "O-[H]"

# Returns a cluster info dictionary of general cluster properties using mol_df
def construct_cluster(mol_df, components, bond_limits=None) -> dict:
    
    cluster_size = len(components)
    cluster = []
    bonds = []
    donors = [] # mostly hydrogens (sometimes N) bonded to O 
    acceptors = [] # mostly oxygens (C-O-H, C-OO-H, C-O-C, C-O-O-C, C=O, NO2)
    if not isinstance(bond_limits, DataFrame):
        bond_limits = default_Hbond_limits()

    s = 0
    for i in range(len(components)):
        mol = mol_df.loc[components[i]]
        cluster.append(mol['xyz'].copy())
        cluster[-1]['mol'] = i
        bonds.append(mol['bonds'].copy())
        
        # shift indeces
        cluster[-1].index += s
        if 'rank' in cluster[-1].columns:
            cluster[-1]['rank'] += s

        bonds[-1]['i1'] += s
        bonds[-1]['i2'] += s
        
        donors.append({ k+s: np.array(v)+s for k,v in mol['donors'].items()})
        acceptors.append({ k+s: np.array(v)+s for k,v in mol['acceptors'].items()})
        s += mol['size']

    if cluster_size > 1:
        cluster = concat(cluster)
    else:
        cluster = cluster[0]
    
    cluster = {'size': len(components), 'components': components, 'atoms': cluster, 'bonds': bonds, 
               'donors': donors, 'acceptors': acceptors, 'limits': bond_limits}
    
    if 'SMILES' in mol_df.columns:
        smiles = mol_df.loc[components]['SMILES'].values
        cluster.update({'SMILES': smiles})

    return cluster

def test_Hbonds(clusters_df: DataFrame, cluster_info: dict, rel_tol=0.2, options=None, return_stats=False):
    n_cl = len(clusters_df)
    n_mol = cluster_info['size']
    distances = clusters_df[("temp", "distances")].values
    xyz = clusters_df[("xyz", "structure")].values
    atoms = cluster_info['atoms']['atom']
    mols = cluster_info['atoms']['mol'].to_numpy()
    bond_limits = cluster_info['limits']
    
    from functools import reduce
    DA = reduce(lambda x, y: x|y, cluster_info['donors'])
    don = list(DA.keys())
    acc = np.concatenate([list(x.keys()) for x in cluster_info['acceptors']])

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
    
    # shift indeces by -1 to match distance matrices
    pairs = {k: (np.array(v[0])-1, np.array(v[1])-1) for k,v in pairs.items() if len(v[0]) > 0 }
    # atoms each donor is bonded to (for angle calculations)
    triples = {k: np.array([DA[i+1]-1 for i in v[0]],dtype=object) for k,v in pairs.items()}
    
    internal = concat(cluster_info['bonds'])
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
        # the number inter-molecular H-bonds per cluster
        inner_counts = np.sum(counts[:, di1, di2], axis=1)
        # test if H-bonds are only inter-molecular
        passed = (totals > inner_counts)
        
    if return_stats:
        uq, uc = np.unique(totals, return_counts=True)
        stat_counts = np.zeros(np.max(uq)+1,dtype=int)
        for i,u in enumerate(uq):
            stat_counts[u] += uc[i]

    if filter_type == 'x':
        min_count = 0
        passed = ~reacted
    elif num == None:
        passed &= np.array(totals) >= 1
    else:
        if filter_type == 'a':
            min_count = len(acc) - num
        elif filter_type == 'd':
            min_count = len(don) - num
        elif filter_type == 'h':
            nH = np.sum(atoms[don]=='H')
            min_count = nH - num
        else:
            min_count = np.max(totals) - num
        min_count = max(1, min_count) # cluster should have at least 1 H-bond

        # select clusters with highest H-bond counts
        passed &= np.array(totals) >= min_count
        
    if return_stats:
        strings = [clusters_df[('info','cluster_type')].values[0], stat_counts, min_count, '{:2.2%}'.format(np.sum(passed) / len(passed))]
        return passed, counts, Hbonds, strings
    
    return passed, counts, Hbonds, None

"""
import matplotlib.pyplot as plt
def test_Hbonding_old(clusters_df, mol_df, components, num, subn_mol=2):

    n_cl = len(clusters_df)
    bonds, counts = get_bonding(clusters_df, mol_df, components)
    
    passed = np.zeros(n_cl, dtype=bool)
    
    di = np.diag_indices(len(components))
    for cl in range(n_cl):
        if bonds[cl]['type'].values[-1] == 'O':
            continue
        c = counts[cl]
        c[di] = 0
        s = np.sum(c)
        passed[cl] = s > 0

    return passed, 0

    totals = np.zeros(n_cl)
    totals_w = np.zeros(n_cl)
    totals_w2 = np.zeros(n_cl)
    for cl in range(n_cl):
        n = len(components)
        counts=np.zeros((n,n)) # H-bond counts between each molecule
        
        for i1,i2 in pairs:
            if i2 in DA[i1]:
                continue
            a1,m1 = atoms.loc[i1]
            a2,m2 = atoms.loc[i2]
            
            d = distances[cl][i1-1,i2-1]
            closest = min(d, closest)
            l,u = bond_limits.at[(a1,a2), 'length']
            if d < l or d > u:
                continue

            i3 = DA[i1][0]
            angle = xyz[cl].get_angle(i3-1, i1-1, i2-1)
            
            l,u = bond_limits.at[(a1,a2), 'angle']
            if angle < l or angle > u:
                continue
            
            counts[m1,m2] += 1
            totals[cl] += 1
            q1 = charges[cl][i1-1] 
            q2 = charges[cl][i2-1]
            q3 = charges[cl][i3-1]

            r.append(d); theta.append(angle)

            #print(q1,q2,q3)

            w = -10*q1*q2*4/d**2
            alphas.append(w * (angle/150)**2)
            totals_w[cl] += w * (angle/150)**2.
            if a1 == 'N':
                angle = 180 - abs(angle-90)

            d2 = distances[cl][i3-1,i2-1] # O--O
            totals_w2[cl] += -20*(q1*q2*4/d**2. + q2*q3*4/d2**2.)
            
    totals_w /= np.mean(totals_w)/4
    totals_w2 /= np.mean(totals_w2)/4
    m = []
    for u in np.unique(totals):
        y = el[totals == u]
        m.append(np.mean(y))
        x = [u] * len(y)
        plt.plot(x,y,'bo', alpha=0.4)
    plt.plot(np.unique(totals),m,'r-')
    plt.plot(totals_w,el,'go', alpha=0.4)
    plt.plot(totals_w2,el,'ro', alpha=0.4)
    plt.show();exit()
    #plt.hist(r, bins=100, label=clusters_df['info']['cluster_type'][0])
    alphas -= np.min(alphas)
    alphas /= np.max(alphas)
    plt.scatter(r,theta, alpha=alphas, label=clusters_df['info']['cluster_type'][0])
    #plt.show()
    
    return
"""

# Test that molecules have not reacted internally
def test_internal_bonds(clusters_df, mol_df, components, rel_tol):
    
    n = len(clusters_df)
    distances = clusters_df[("temp", "distances")].values

    passed = np.array([True]*n)
    for cl in range(n):
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

def subclusters():

    return

# Test which molecules are properly clustered.
# Also returns any found subclusters.
def test_clustering(clusters_df, mol_df, components, limits=None, subn_mol=2):
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
