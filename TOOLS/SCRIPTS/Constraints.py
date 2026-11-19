import os, sys
from random import uniform, sample
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import chain, cycle

from time import time
import numpy as np
import argparse
from copy import deepcopy

from clusterfilter.data_readers import read_molecule_data
from clusterfilter.cluster_analysis import construct_cluster
from clusterfilter.topologger import generate_rdkit_cluster, add_Hbonds
from constraints_utils import get_run_config, real_path

def get_fingerprints(combinations, cluster, keep_h):
    clusters_to_test = []
    for pairs in combinations:
        don, acc = np.array(pairs).T - 1
        mol = add_Hbonds(cluster, don, acc, keepH=keep_h)
        clusters_to_test.append(mol)

    # Uniqueness filtering using Morgan/circular fingerprint
    from rdkit.Chem.AllChem import GetMorganGenerator
    fpgen = GetMorganGenerator(radius=20, fpSize=128)
    fps = [fpgen.GetCountFingerprint(cluster).ToBinary() for cluster in clusters_to_test]

    return fps, clusters_to_test


def cluster_sampling(donors, acceptors, cluster_size, nH, max_bonds, AD, DA, n_iter, seed):
    # np.random.default_rng(seed)
    np.random.seed(seed)

    combinations = []

    # Cluster combination sampling
    mol1, mol2 = np.meshgrid(range(cluster_size),range(cluster_size))
    mol_pairs = np.array([[m1,m2] for m1,m2 in zip(mol1.flatten(), mol2.flatten()) if m1!=m2])

    """def mol_generator():
        x = np.concatenate()
        while True:
            yield from x"""
    # generator for random pair permutations
    perm = [np.random.permutation(mol_pairs) for i in range(nH*5+1)]
    gen = cycle(chain(*perm))

    for i in range(n_iter):
        O = [np.random.permutation(a).tolist() for a in acceptors]
        H = [np.random.permutation(d).tolist() for d in donors]
        pair_list = []

        j = 0
        n = min(nH, max_bonds)
        while n > 0 and j < 100:
            j += 1
            a, b = next(gen) # pick two molecules
            if (len(H[a]) == 0) or (len(O[b]) == 0):
                continue

            # pick donor and acceptor
            pair_list.append([H[a].pop(), O[b].pop()])
            n -= 1

        pair_list = np.array(pair_list)
        combinations.append(remove_overlaps(pair_list, AD, DA))
    return combinations


def remove_overlaps(pairs, AD, DA):
    h1,o2 = pairs.T
    h2 = AD[o2]
    o1 = DA[h1]
    pairs2 = np.array([h2,o1]).T
    grid = np.array([(p==pairs).all(axis=1) for p in pairs2])
    """
    Pick lower triangle:
    [[-    False  True False]
    [False   -   False  True]
    [ True False   -   False]
    [False  True False   -]]
    """
    passed = (~np.tril(grid)).all(axis=0)
    return pairs[passed]

def get_additional_data(file_in, keep_h, n_procs=1):
    cluster_type = ''
    with open('rename.sh') as f:
        for line in f:
            if line.startswith('  mv $i'):
                cluster_type = line.split()[-1].split('-')[0]
                break

    if not real_path(file_in) or cluster_type == '':
        print('Pickled data not found in parent folder:', os.path.abspath('../..'))
        return None, None, None
    
    from clusterfilter.filter import ClusterFilter
    cf = ClusterFilter(file_in, mol_file=real_path('parameters.txt'))
    cf.extract_clusters(cluster_type)
    cf.set_Hbond('H','O', angle=(110, 180))
    cf.Hbonded(100)
    
    df = cf.get_filtered_data(return_temp=True)
    if len(df) == 0:
        return None, None, None
    
    lowest = np.sort(df[('log','electronic_energy')].values)[0]*627.5
    pairs = df[('temp', 'Hbond_pairs')].values
    combinations = np.array([np.array([a,b]).T + 1 for a,b in pairs], dtype=object)
    max_bonds = max(len(a) for a,b in pairs)
    cluster_info = cf._cluster_info[cluster_type]

    # parts = monomers
    cluster = generate_rdkit_cluster(cluster_info.smiles, add_ids=True)
    get_fingerprint_partial = partial(get_fingerprints, cluster=cluster, keep_h=keep_h)
    
    # Remove isomorphically identical duplicates
    with Pool(processes=n_procs) as pool:
        result = pool.map(get_fingerprint_partial, np.array_split(combinations, n_procs))
        fingerprints = np.concatenate([x[0] for x in result])

    return fingerprints, lowest, max_bonds


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate constraints files for a given cluster types")
    parser.add_argument('filename')
    parser.add_argument("-i", "-input", type=str, default=None, help="yaml file containing additional arguments")

    args, _ = parser.parse_known_args()

    ## check that input file is provided
    file_in = args.filename  # should be calc.inp
    if not real_path(file_in): print(
        "File " + file_in + " not found. Make sure you are in the correct folder");exit()
    
    ## check if config file is provided with default arguments
    cfg = get_run_config(args.i)

    ## use config file values as default values
    parser.add_argument("-nfiles", type=int, default=cfg.get("nfiles", 100), help="Maximum number of files to create")
    parser.add_argument("-all", nargs="?", default=cfg.get("all", False), const=True,
                        help="Include all oxygen as available for hydrogen bonding acceptors")
    parser.add_argument("-max", type=int, default=cfg.get("max", 100), help="Maximum number of H-bonds to create between molecules")
    parser.add_argument("-folder", type=str, default=cfg.get("folder", '.'), help="Folder to save files to")
    ## options used with model predictions
    parser.add_argument("-cutr", type=float, default=cfg.get('cutr', 10.0), help="Energy cutoff (kcal/mol) when using model predictions")
    parser.add_argument("-pkl", type=str, default=cfg.get('pkl', ''), help="Pickle database of sampled clusters")
    parser.add_argument("-min_files", type=int, default=100, help="Minimum number of predictions to save")
    parser.add_argument("-config_root", type=str, default=cfg.get("config_root", 'config'), help="Path to the config root directory")
    parser.add_argument("-model", type=str, default=cfg.get("model", 'gen_config'), help="Config file name to create the correct model object")
    parser.add_argument("-model_ckpt_folder", type=str, default=cfg.get("model_ckpt_folder", "ckpt/cluster/gcn_emb32_seed_42_lr_0.001_bs_64_time_20260113-111758"), help="Ckpt directory")
    parser.add_argument("-use_molgraphconvfeat", nargs="?", default=cfg.get("use_molgraphconvfeat", False), const=True,
                            help="A flag indicating if deepchem's MolGraphConvFeaturizer is used. Default is assumed to be a OGB featurizer (using AtomEncoder and BondEncoder)")
    parser.add_argument("-enc_monomers", nargs="?", default=cfg.get("enc_monomers", False), const=True,
                        help="Encoding individual monomers inside the cluster dataset instead of the full clusters")
    parser.add_argument("-dataset", type=str, default=cfg.get("dataset", "cluster"),
                        choices=["cluster", "flexibility"], help="Dataset to use (cluster or flexibility)")
    parser.add_argument("-data_path", type=str, default=cfg.get("data_path", "../data"),
                        help="Path to the root directory with data")
    #parser.add_argument("-max_flexibility", nargs="?", type=float, default=cfg.get("max_flexibility", None),
    #                    help="Maximum flexibility value for samples in the flexibility dataset, used in normalization")
    parser.add_argument("-include_mono", nargs="?", default=cfg.get("include_mono", False), const=True,
                        help="A boolean flag indicating if monomers should be included in the training data along with clusters (only for cluster dataset)")
    parser.add_argument("-keep_h", nargs="?", default=cfg.get("keep_h", False), const=True,
                        help="A boolean flag indicating if hydrogen atoms should be kept in the molecular graphs during pre-processing")
    parser.add_argument("-include_col", nargs="+", default=cfg.get("include_col", []),
                        help="Optional list of collections for clusters dataset separated by spaces")
    parser.add_argument("-only_col", nargs="?", default=cfg.get("only_col", False), const=True,
                        help="If set, only collections from include_col are used for clusters dataset")
    parser.add_argument("-eshift_by_monomers", nargs="?", default=cfg.get("eshift_by_monomers", False), const=True,
                        help="If set, the training targets for  the cluster dataset are shifted by the sum of the minimal energies of the individual monomers in the cluster, "
                             "instead of per atom energy shifts calculated on the train dataset.")
    parser.add_argument("-validate_model", nargs="?", default=cfg.get("validate_model", False), const=True,
                        help="A bool flag indicating if model validation should be run before evaluation to make sure that the ckpt is loaded correctly")

    ## parse arguments, CLI provided arguments override provided yaml values
    args = parser.parse_args()    
    print("Current arguments are", args)
    
    #args.config_root = real_path(args.model_ckpt_folder)
    use_all = args.all
    max_files = args.nfiles

    n_procs = min(cpu_count(), 8)
    print(f"Using {n_procs} processes")

    ###############################################################################
    print('Constraints.py: Start.')

    mol_df = read_molecule_data(file_in)
    n_mol = len(mol_df)
    print("Found", n_mol, "molecule(s) in calc.inp.")

    # read pickle data if available
    pickled_file = real_path(args.pkl)
    old_fps, lowest_el, max_bonds = get_additional_data(pickled_file, args.keep_h, n_procs)
    if max_bonds == None:
        max_bonds = args.max

    ###############################################################################
    
    components = []
    for name, mol in mol_df.items():
        components += [name] * mol.q
    cluster_size = len(components)

    # Construct the cluster
    cluster_info = construct_cluster(mol_df, components)
    atoms = cluster_info.atoms
    bonds = cluster_info.bonds

    donors = [list(x.keys()) for x in cluster_info.donors] # mostly hydrogens (sometimes N) bonded to O
    acceptors = [list(x.keys()) for x in cluster_info.acceptors] # mostly oxygens (C-O-H, C-OO-H, C-O-C, C-O-O-C, C=O, NO2)
    # {A:D} dictionary for preventing overlaps of acceptor-donor pairs
    AD = {}
    for x in cluster_info.acceptors:
        AD.update({k: v for k,v in x.items() if len(v) > 0})

    ############################################################
    # removes overlapping OH--OH pairs
    k = list(AD.keys())
    v = np.array([d[0] for d in AD.values()])
    AD = np.zeros(np.max([k,v])+1, dtype=int)
    AD[k] = v
    DA = np.zeros(np.max([k,v])+1, dtype=int)
    DA[v] = k

    t0 = time()

    nH = np.sum([len(x) for x in donors])
    nO = np.sum([len(x) for x in acceptors])
    
    combinations = []
    if cluster_size == 1:
        # Conformer sampling
        donors = donors[0]
        acceptors = acceptors[0]
        for i in range(int(1.5 * max_files)):
            pair_list = []
            O = deepcopy(acceptors)
            H = deepcopy(donors)

            for i in range(nH):
                # pick donor and acceptor
                hi = sample(H,1)[0]; H.remove(hi)
                oj = sample(O,1)[0]; #O.remove(oj)
                pair_list.append([hi,oj])

            pair_list = np.array(pair_list)
            combinations.append(remove_overlaps(pair_list, AD, DA))
    else:
        n = int(1.5 * max_files/n_procs)
        with Pool(processes=n_procs) as pool:
            result = pool.starmap(cluster_sampling,
                                  [(donors, acceptors, cluster_size, nH, max_bonds, AD, DA, n, i)
                                   for i in range(n_procs)])
            combinations = sum(result, [])

    combinations = np.array(combinations, dtype=object)
    print('Generation time:', time()-t0)

    ###########################################################
    lens = np.array([len(c) for c in combinations])
    ulen, counts = np.unique(lens, return_counts=True)
    if counts[-1] > 2*max_files:
        min_len = ulen[-1]
    else:
        ulen = np.flip(ulen)
        cumul = np.cumsum(np.flip(counts))
        min_len = ulen[cumul < 2*max_files][-1]-1

    # Select combinations with the highest number of bonds
    combinations = combinations[lens >= min_len]

    # parts = monomers
    try:
        cluster, parts = generate_rdkit_cluster(cluster_info.smiles, add_ids=True, return_parts=True)
        get_fingerprint_partial = partial(get_fingerprints, cluster=cluster, keep_h=args.keep_h)
    except:
        print('Error! Failed to generate fingerprints. Are SMILES missing?')
        cluster_size = 1

    if cluster_size > 1:
        # Remove isomorphically identical duplicates
        with Pool(processes=n_procs) as pool:
            result = pool.map(get_fingerprint_partial, np.array_split(combinations, n_procs))
            fingerprints = np.concatenate([x[0] for x in result])
            clusters_to_test = np.concatenate([x[1] for x in result])

        unique_fps, idx = np.unique(fingerprints, return_index=True)
        if isinstance(old_fps, np.ndarray):
            # remove already existing clusters
            idx = idx[~np.isin(fingerprints[idx], old_fps)]

        combinations = combinations[idx]
        selected_clusters = fingerprints[idx]
        clusters_to_test = clusters_to_test[idx]
        print('Filtered', len(selected_clusters), 'unique combinations')
        print('Filtering time:', time()-t0)

        # Predict energies with a trained neural network model
        args.model_ckpt_folder = '/'.join(args.model_ckpt_folder.split('/')[:-1])
        ckpt_folder = real_path(args.model_ckpt_folder)
        if args.model != '' and ckpt_folder and len(clusters_to_test) > 0:
            from utils.eval_utils import run_evaluation
            print("Calculate energies with a trained model from " + ckpt_folder)
            t = time()
            energies, eshift_values = run_evaluation(clusters_to_test, parts, args, ckpt_folder,
                                                     validate=args.validate_model, verbose=True)
            print('Evaluation time:', time()-t0)
            eshift = eshift_values[[a.GetAtomicNum() for a in cluster.GetAtoms()]].sum()
            energies = (energies + eshift)*627.5
            
            if not isinstance(lowest_el, float):
                lowest_el = min(energies)

            # select clusters with E <= E_min + 10 kcal/mol
            cutoff = lowest_el + args.cutr
            idx = np.argsort(energies)
            energies = energies[idx]
            if (energies < cutoff).sum() > args.min_files:
                combinations = combinations[idx][energies < cutoff]
            else:
                combinations = combinations[idx[:args.min_files]]

    ############################################################
    # sort and save at most max_files combinations
    combinations = sample(list(combinations), min(max_files, len(combinations)))
    combinations = [sorted([(a,b) for a,b in c]) for c in combinations]

    # Create output folder
    if not os.path.exists(args.folder):
        os.makedirs(args.folder)

    # Create constraints[i].inp files
    for i in range(len(combinations)):
        file_out = 'constraints'+str(i)+'.inp'
        with open(args.folder+'/'+file_out, 'w') as f:
            f.write("$constrain\n")
            f.write("  force constant=0.001\n")
            for a, b in combinations[i]:
                #distances = {'H': 1.85, 'C': 2.9, 'N': 2.9}
                if atoms.at[a, 'atom'] == 'H':
                    f.write("  distance: "+str(a)+", "+str(b)+", 1.85\n")
                else:
                    f.write("  distance: "+str(a)+", "+str(b)+", 2.9\n")
            f.write("$end\n")

    print('Created '+str(len(combinations))+' constraint files')