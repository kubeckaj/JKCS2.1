import os 
import numpy as np
from pandas import DataFrame
from typing import Callable, Iterable, Optional, Union
from data_readers import MolInfo, read_molecule_data, read_pickled_data
from cluster_analysis import ClusterInfo

class ClusterData():
    _unique_cluster_types = Iterable[str]

    def __init__(self, 
        file_in: Union[Optional[str], Iterable[str]] = None, 
        *,
        mol_file: Optional[str] = None,
        mol_df: Optional[MolInfo] = None,
        total_charge: int = 0
    ) -> None:
        self.file_in = file_in
        self.mol_file = mol_file
        self.clusters_df = None
        self._unique_cluster_types = []
        self.charge = total_charge

        if not isinstance(mol_df, MolInfo):
            self.read_molecule_data(mol_file)
        else:
            self.mol_df = mol_df
        self.read_cluster_data(file_in)

    def __len__(self):
        if isinstance(self.clusters_df, DataFrame):
            return len(self.clusters_df)
        else:
            return 0

    @property
    def cluster_types(self):
        return self._unique_cluster_types
    
    def read_molecule_data(self, 
        mol_file: Optional[str] = None
    ) -> DataFrame:
        
        from os.path import isfile
        # look for parameters.txt or input.txt
        if not isinstance(mol_file, str):
            mol_file = "parameters.txt"
            if not isfile(mol_file): 
                mol_file = "input.txt"
                if not(isfile(mol_file)): 
                    print("input.txt or parameters.txt not found. Make sure you are in the correct folder!")
                    self.mol_df = None
                    return
                
        elif not isfile(mol_file): 
            print('Warning! '+mol_file+" not found. Make sure --mol_file path is correct")
            self.mol_df = None
            return
            
        self.mol_file = mol_file
        self.mol_df = read_molecule_data(mol_file)
        self.n_mol = len(self.mol_df)
        
        return self.mol_df

    def read_cluster_data(self, 
        file_in: Union[Optional[str], Iterable[str]] = None,
    ) -> DataFrame:
            
        if isinstance(file_in, str):
            file_in = [file_in]
            
        if isinstance(file_in, Iterable):
            for f in file_in:
                if not f.endswith(('.txt', '.dat', '.pkl')):
                    print('Error! Cannot read data from: '+ f)
                    exit()
                
            self.clusters_df = read_pickled_data(file_in)

        if ("info", "cluster_type") in self.clusters_df.columns:
            cluster_types = self.clusters_df[("info", "cluster_type")].values
            self._unique_cluster_types = np.unique(cluster_types)
        else:
            self._unique_cluster_types = []

        return self.clusters_df
    
    # returns a copy of the loaded dataframe
    def get_data(self, return_temp=False):
        if (not return_temp) and ('temp' in self.clusters_df.columns):
            df = self.clusters_df.drop('temp',axis=1)
            return df
        else:
            return self.clusters_df.copy(deep=True)

    def save_pickle(self, file_out: str):
        if file_out.endswith('.pkl'):
            self.clusters_df.to_pickle(file_out)
        else:
            raise Exception('Output file name must end with .pkl')
        return 
    
    
class ClusterFilter(ClusterData):
    
    clusters_df: DataFrame
    mol_df: dict[str, MolInfo]
    _components: dict[str, Iterable[str]]
    _cluster_info: dict[str, ClusterInfo]

    def __init__(self, 
        file_in: Union[str, Iterable[str]] = None, 
        **kwargs,
    ) -> None:
        
        super().__init__(file_in, **kwargs)
        if isinstance(file_in, str) and file_in.endswith(('.dat','.txt')):
            self._lines = np.array([l for l in open(file_in).readlines() if '/:EXTRACT:/' in l])

        self.reset()
        self.reset_Hbond_limits()    

        # Check that all monomers are in mol_df
        monomers = np.unique(np.concatenate(list(self._components.values())))
        not_found = [mol for mol in monomers if mol not in self.mol_df.keys()]

        if len(not_found) == 1:
            print("Warning: Molecule "+not_found[0]+" was not found in "+self.mol_file)
        elif len(not_found) > 5:
            print("Warning: Molecules "+not_found[0]+',...,'+not_found[-1]+" were not found in "+self.mol_file)
        elif len(not_found) > 1:
            print("Warning: Molecules "+','.join(not_found)+" were not found in "+self.mol_file)
        

    def read_cluster_data(self, file_in = None) -> DataFrame:
        super().read_cluster_data(file_in)

        from cluster_analysis import distance_matrix
        # precalculating distance matrices
        self.clusters_df[("temp", "distances")] = distance_matrix(self.clusters_df)
        return self.clusters_df

    # reset filter (all values are set as passed)
    def reset(self):
        self._cluster_info = {}
        self.history = []

        # Separate clusters per type into subsets
        if ("info", "cluster_type") not in self.clusters_df.columns:
            raise Exception("Pickled data has no cluster types.")
        
        cluster_types = self.clusters_df[("info", "cluster_type")].values
        self._cluster_subsets = {}
        self._components = {}
        for ct in self.cluster_types:
            # Get indeces to each cluster type
            indexes = self.clusters_df[cluster_types==ct].index.values
            self._cluster_subsets[ct] = indexes
            
            # Get lists of components
            components = []
            unique_cl = self.clusters_df.loc[indexes[0]]
            for i, c in enumerate(unique_cl['info', 'components']):
                ratio = unique_cl['info', 'component_ratio'][i]
                components += [c]*ratio
            self._components[ct] = components
        
        # filter if atoms are missing
        self._filter(lambda df: [isinstance(x, Iterable) for x in df[("temp", "distances")].values])
        return
    
    @property
    def components(self):
        return self._components

    def cluster_info(self, cluster_type: str) -> ClusterInfo:
        if cluster_type not in self.components.keys():
            return None
        
        if cluster_type not in self._cluster_info.keys():
            from cluster_analysis import construct_cluster
            components = self.components[cluster_type]
            self._cluster_info[cluster_type] = construct_cluster(self.mol_df, components)
            
        return self._cluster_info[cluster_type]


    def _filter(self, func: Callable, *args, **constants):
        for ct in self.cluster_types:
            subset = self._cluster_subsets[ct]
            if len(subset) == 0:
                self._cluster_subsets[ct] = np.array([])
                continue
                
            # *args inputs are iterated over
            iterable_args = [arg[ct] for arg in args]
            passed = func(self.clusters_df.loc[subset], *iterable_args, **constants)
            if len(passed) == 0:
                self._cluster_subsets[ct] = np.array([])
            else:
                self._cluster_subsets[ct] = subset[passed]
        
        return

    # select n lowest clusters
    def select(self, n: int, sort_by: tuple=('log','electronic_energy')):
        if n > 0:
            from cluster_analysis import select_lowest
            self._filter(select_lowest, n=n, by=sort_by)
        return 
    
    # test clustering (molecular distance) and reactions between molecules
    def distance(self, minH=1.4, minO=2.0, max=3):
        self.history.append('distance')
        from cluster_analysis import test_clustering

        # Filter out clusters where molecules are too far apart (or too close)
        self._filter(test_clustering, self._components, mol_df=self.mol_df, limits=[minH, max])
        
        return 
    
    # test internal reactions
    def reacted(self, itol: float = 0.2):
        self.history.append('reacted')
        from cluster_analysis import test_internal_bonds

        # Filter out incorrectly bonded structures (where a chemical reaction has occured)
        self._filter(test_internal_bonds, self._components, mol_df=self.mol_df, rel_tol=itol)

        return 
    
    def converged(self, fmax: float = 0.01, invert: bool = False):
        self.history.append('converged')
        from cluster_analysis import test_convergence

        # Filter out (un)converged structures as defined by fmax (eV/Ang)
        self._filter(test_convergence, fmax=fmax, invert=invert)

        return 
    
    def extract_clusters(self, cluster_types: Union[str, Iterable[str]]):
        if isinstance(cluster_types, str):
            cluster_types = [cluster_types]

        func = lambda df, ct: [df[('info', 'cluster_type')].values[0] in ct]*len(df)
        self._filter(func, ct=cluster_types)
        return
    
    def except_clusters(self, cluster_types: Union[str, Iterable[str]]):
        if isinstance(cluster_types, str):
            cluster_types = [cluster_types]

        func = lambda df, ct: [df[('info', 'cluster_type')].values[0] not in ct]*len(df)
        self._filter(func, ct=cluster_types)
        return

    def monomers(self, extract=True):
        m = [k for k,v in self.components.items() if len(v)==1]
        if len(m) == 0:
            return
        
        if extract:
            self.extract_clusters(m)
        else:
            self.except_clusters(m)
        return

    @property
    def Hbond_limits(self):
        return self._Hbond_limits 

    def reset_Hbond_limits(self):
        from cluster_analysis import default_Hbond_limits
        self._Hbond_limits = default_Hbond_limits()
        
    def set_Hbond(self, 
        donor: str, 
        acceptor: str, 
        length: Optional[tuple[float]] = None, 
        angle: Optional[tuple[float]] = None
    ):
        from numpy import nan
        if (donor,acceptor) not in self._Hbond_limits.index:
            self._Hbond_limits.loc[(donor,acceptor), :] = nan
        if isinstance(length, tuple):
            self._Hbond_limits.at[(donor,acceptor), 'length'] = length
        if isinstance(angle, tuple):
            self._Hbond_limits.at[(donor,acceptor), 'angle'] = angle
        return
    
    def Hbonded(self, 
        min_Hbonds: int = 0,
        filter_type: str = 'D',
        rel_tol: float = 0.2,
        *,
        frac: float = 0.0,
        return_stats: bool = False,
    ):
        # add proton donors and acceptors to mol_df
        from cluster_analysis import test_Hbonds

        self.history.append('Hbonded')

        self._rel_tol = rel_tol
        self._Hbond_options=[min_Hbonds, filter_type]
        self.clusters_df[('temp', 'N_Hbonds')] = object()
        self.clusters_df[('temp', 'Hbond_pairs')] = object()
        self.clusters_df[('temp', 'Hbond_lengths')] = object()
        self.clusters_df[('temp', 'Hbond_angles')] = object()
        self._Hbond_counts = {}
        self._Hbond_stats = {}

        def func(clusters_df, **kwarks):
            ct = clusters_df[('info', 'cluster_type')].values[0]
            
            # dict which contains info needed by filter about the cluster type
            if len(clusters_df) > 0:
                cluster_info = self.cluster_info(ct)
                passed, Hbond_counts, Hbond_pairs, stats = test_Hbonds(clusters_df, cluster_info, **kwarks)

                # save these for later use
                for i, idx in enumerate(clusters_df.index.values):
                    self.clusters_df.at[idx, ('temp', 'N_Hbonds')] = np.sum(Hbond_counts[i])
                    self.clusters_df.at[idx, ('temp', 'Hbond_pairs')] = tuple(Hbond_pairs[i][:2])
                    self.clusters_df.at[idx, ('temp', 'Hbond_lengths')] = Hbond_pairs[i][2]
                    self.clusters_df.at[idx, ('temp', 'Hbond_angles')] = Hbond_pairs[i][3]
                
                self._Hbond_counts[ct] = Hbond_counts
                self._Hbond_stats[ct] = stats
                return passed
            else:
                self._Hbond_counts[ct] = None
                return []
        
        self._filter(func, rel_tol=self._rel_tol, options=self._Hbond_options, 
                     bond_limits=self._Hbond_limits, frac=frac, return_stats=return_stats)
        
        if return_stats:
            n = max([len(l[1]) for l in self._Hbond_stats.values()])
            
            from pandas import set_option
            set_option('display.max_columns', None)
            set_option('display.max_rows', None)
            set_option('display.width', None)
            
            count_table = np.array([np.pad(self._Hbond_stats[ct][1], 
                                           (0, n-len(self._Hbond_stats[ct][1])), 
                                           'constant', constant_values=(0,0)) 
                                    for ct in self._Hbond_stats.keys()]).astype(str)
            for i, ct in enumerate(self._Hbond_stats.keys()):
                l = self._Hbond_stats[ct]
                if l[2] >= count_table.shape[1]:
                    count_table[i, -1] += '|'
                else:
                    count_table[i, l[2]] ='|'+count_table[i, l[2]] 

            df = DataFrame(data={'cluster': [self._Hbond_stats[ct][0] for ct in self._Hbond_stats.keys()], 
                                 'passed': [self._Hbond_stats[ct][-1] for ct in self._Hbond_stats.keys()]})
            for i in range(n):
                df[i] = count_table[:,i]
            return df

        return 
    
    # uniqueness filtering (e.g. rg,el,dip)
    def uniq(self, *, rg=None, el=None, g=None, dip=None, dup=None):
        
        return 
    
    def cutr(self, cutoff: float, by: Union[str, tuple]):
        if by not in self.clusters_df:
            raise Exception('column',by,'not found in cluster data')
            
        from cluster_analysis import cut_relative
        if cutoff >= 0:
            self._filter(cut_relative, cutoff=cutoff, by=by)

        return 
    
    def topology(self, 
        n: int = 1,  # Selecting multiples not implemented yet
        selected_isomers_file: Optional[str] = None,
        excepted_topologies_file: Optional[str] = None,
        sort_by: Union[str, tuple] = ('log','electronic_energy')
    ) -> None:
        if 'Hbonded' not in self.history:
            raise Exception('H-bond filtering must be run before filtering topologies.')

        from topologger import clusters_to_smiles, filter_isomorphs
        self.clusters_df[('temp', 'SMILES')] = object()
        for ct in self.cluster_types:
            subset = self._cluster_subsets[ct]
            if len(subset) == 0:
                continue
            self.clusters_df.loc[subset, ('temp', 'SMILES')] =\
                clusters_to_smiles(self.clusters_df.loc[subset], self.cluster_info(ct))
            
        self.file_iso = selected_isomers_file
        if isinstance(self.file_iso, str):
            if not os.path.isfile(self.file_in):
                print('Error! Isomer file '+self.file_in+' not found.');exit()

            # Temporary filter to extract monomers isomers
            isomers_cf = ClusterFilter(self.file_iso, mol_file=self.mol_file)
            isomers_df = isomers_cf.get_unique_isomers(self.components)

            # Add smiles to MolInfos
            for ct, (c, ) in isomers_cf.components.items():
                smiles = isomers_df[isomers_df[('info', 'cluster_type')]==ct][('temp', 'SMILES')].values
                if len(smiles) > 0:
                    self.mol_df[c].isomers = smiles

            # Add smiles to ClusterInfos
            for cluster_info in self._cluster_info.values():
                iso_smiles = [self.mol_df[k].isomers for k in cluster_info.components]
                cluster_info.isomers = iso_smiles

            from topologger import filter_isomers
            self._filter(filter_isomers, self._cluster_info)
        else:
            self._isomers_df = None
        if n > 1:
            return

        self.file_not = excepted_topologies_file
        if isinstance(self.file_not, str):
            from cluster_analysis import distance_matrix, test_Hbonds

            # extract already found topologies 
            used_topo_df = read_pickled_data(self.file_not)
            used_topo_df[("temp", "distances")] = distance_matrix(used_topo_df)
            self._used_pairs = {}
            for ct in self.cluster_types:
                subset = self._cluster_subsets[ct]
                if len(subset) == 0:
                    continue
                used = used_topo_df[used_topo_df['info', 'cluster_type']==ct]
                if len(used) > 0:
                    f, _, used_topologies, _ = test_Hbonds(used, self._cluster_info[ct], 
                                                           self._rel_tol, bond_limits=self._Hbond_limits,
                                                           options=self._Hbond_options)
                    self._used_pairs[ct] = np.array(used_topologies, dtype=tuple)[f]

            # filter unique topologies and remove already found topologies 
            self._filter(filter_isomorphs, self._cluster_info, self._used_pairs, el=sort_by)
        else:
            # filter unique topologies 
            self._filter(filter_isomorphs, self._cluster_info, el=sort_by)
            
        return 
    
    def get_filtered_length(self) -> int:
        return sum([len(s) for s in self._cluster_subsets.values()])
        
    # return filtered dataframe
    def get_filtered_data(self, return_temp=False) -> DataFrame:
        
        idx = np.concatenate([v for v in self._cluster_subsets.values()])
        df = self.clusters_df.loc[idx].sort_index(axis=1)
        if (not return_temp) and ('temp' in self.clusters_df.columns):
            df = df.drop('temp',axis=1)
        return df
    
    def save_to(self, file_out: str):
        filtered_clusters = self.get_filtered_data()[('info', 'file_basename')].values
        paths = DataFrame([l.split()[0].split('/:EXTRACT:/') for l in self._lines], columns=['file', 'cluster'])        
        passed = paths['cluster'].isin(filtered_clusters)
        
        with open(file_out, 'w') as f:
            f.writelines(self._lines[passed])

        return 
    
    def get_unique_isomers(self, components: dict=None) -> DataFrame:
        self.reset()
        self.monomers()

        # Make sure we only select monomers
        if components is not None:
            extract = '1'+np.unique(np.concatenate(list(components.values())))
            self.extract_clusters(extract)
        else:
            self.monomers()
        
        self.Hbonded(100, 'X')
        self.topology()

        isomers_df = self.get_filtered_data(return_temp=True)
        if ('log', 'gibbs_free_energy') in isomers_df.columns:
            isomers_df = isomers_df.sort_values(by=[('log', 'gibbs_free_energy'), ('log', 'electronic_energy')])
        if ('log', 'electronic_energy') in isomers_df.columns:
            isomers_df = isomers_df.sort_values(by=('log', 'electronic_energy'))

        return isomers_df.drop_duplicates(subset=[('temp', 'SMILES')])
    
    def get_binding_energies(self, high_df: str=None):
        # 1. Extract monomers
        monomers = self.get_unique_isomers()

        # 2. Extract clusters
        self.reset()
        self.monomers(False)
        self.Hbonded(100, rel_tol=0.3)
        self.topology(sort_by=('log', 'gibbs_free_energy')) 

        columns = ['Cluster', 'Delta-E', 'Delta-G']
        #if isinstance(high, str):
        #    columns += ['High DE', 'Corrected DG']
        columns += ['SMILES']

        results = DataFrame(columns=columns)
        for ct in self.cluster_types:
            components = self._components[ct]
            if len(components) < 2:
                continue
            if len(self._cluster_subsets[ct])==0:
                results.loc[len(results), :] = [ct, np.nan, np.nan, np.nan]
                continue

            subset = self.clusters_df.loc[self._cluster_subsets[ct]]

            subset = subset.sort_values(by=[('log', 'gibbs_free_energy'), ('log', 'electronic_energy')])
            subset = subset.drop_duplicates(subset=[('temp', 'SMILES')])

            for i in range(len(subset)):
                dE, dG, smiles = subset[[('log', 'electronic_energy'), ('log', 'gibbs_free_energy'), ('temp', 'SMILES')]].values[i]
                for smi in smiles.split('.'):
                    m = monomers[monomers[('temp', 'SMILES')]==smi]
                    if len(m) == 0:
                        dE = np.nan
                        dG = np.nan
                        break
                    
                    dE -= m[('log', 'electronic_energy')].values[0]
                    dG -= m[('log', 'gibbs_free_energy')].values[0]
                
                results.loc[len(results), :] = ct, dE, dG, smiles

        results = results.sort_values(by='Cluster')
        #results = results.convert_dtypes()
        results[results.columns[1:-1]] = 627.5 * results[results.columns[1:-1]].astype(float)

        return results

        
        

if __name__ == '__main__':
    cf = ClusterFilter(['lowest.pkl', 'DADBsa_UMA.pkl'])
    table = cf.get_binding_energies()
    table.set_index('Cluster').to_csv('binding_energies.txt', float_format="%.5f")
    print(table)