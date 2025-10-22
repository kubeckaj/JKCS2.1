import numpy as np
from pandas import DataFrame
from typing import Callable, Iterable, Optional, Union

class ClusterData():
    def __init__(self, 
        file_in: Union[Optional[str], Iterable[str]] = None, 
        *,
        mol_file: Optional[str] = None,
    ) -> None:
        self.file_in = file_in
        self.mol_file = mol_file
        self.clusters_df = None
        
        self.read_molecule_data(mol_file)
        self.read_cluster_data(file_in)

    def __len__(self):
        if isinstance(self.clusters_df, DataFrame):
            return len(self.clusters_df)
        else:
            return 0

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
            
        from data_readers import read_molecule_data
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
                if not f.endswith(('.txt','.dat', '.pkl')):
                    print('Error! Cannot read data from: '+ f)
                    exit()
                
            from data_readers import read_pickled_data
            self.clusters_df = read_pickled_data(file_in)

        return self.clusters_df
    
    # returns a copy of the loaded dataframe
    def get_data(self):
        if 'temp' in self.clusters_df.columns:
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
    
    def binding_energies():

        return
    

class ClusterFilter(ClusterData):
    def __init__(self, 
        file_in: str, 
        **kwargs,
    ) -> None:
        
        super().__init__(file_in, **kwargs)
        if file_in.endswith(('.dat','.txt')):
            self._lines = np.array([l for l in open(file_in).readlines() if '/:EXTRACT:/' in l])
        self.reset()

        # Check that all monomers are in mol_df
        monomers = np.unique(np.concatenate(self._components))
        for mol in monomers:
            if mol not in self.mol_df.index:
                print("Warning: Molecule '"+mol+"' was not found in "+self.mol_file)

        from cluster_analysis import distance_matrix
        # precalculating distance matrices
        self.clusters_df[("temp", "distances")] = distance_matrix(self.clusters_df)
        self._filter(lambda df: [isinstance(x, Iterable) for x in df[("temp", "distances")].values])

    # reset filter (all values are set as passed)
    def reset(self):
        self._cluster_info = None
        self._Hbond_limits = None
        self.history = []

        # Separate clusters per type into subsets
        if ("info", "cluster_type") not in self.clusters_df.columns:
            raise Exception("Pickled data has no cluster types.")
        
        cluster_types = self.clusters_df[("info", "cluster_type")].values
        self._unique_cluster_types = np.unique(cluster_types)
        self._cluster_subsets = []
        self._components = []
        for ct in self._unique_cluster_types:
            # Get indeces to each cluster type
            indexes = self.clusters_df[cluster_types==ct].index.values
            self._cluster_subsets.append(indexes)
            
            # Get lists of components
            components = []
            unique_cl = self.clusters_df.loc[indexes[0]]
            for i, c in enumerate(unique_cl['info', 'components']):
                ratio = unique_cl['info', 'component_ratio'][i]
                components += [c]*ratio
            self._components.append(components)
        
        return
    
    @property
    def cluster_types(self):
        return self._unique_cluster_types
    
    def _filter(self, func: Callable, iter_names={}, **kwarks):
        for k, ct in enumerate(self.cluster_types):
            # adding variable names to **kwarks needed during iteration
            for key, var in iter_names.items():
                kwarks.update({key: self.__dict__['_'+var][k]})
                
            subset = self._cluster_subsets[k]
            passed = func(self.clusters_df.loc[subset], **kwarks)
            if len(passed) == 0:
                self._cluster_subsets[k] = np.array([])
            else:
                self._cluster_subsets[k] = subset[passed]
        
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
        self._filter(test_clustering, {'components':'components'}, mol_df=self.mol_df, limits=[minH, max])
        
        return 
    
    # test internal reactions
    def reacted(self, itol: float = 0.2):
        self.history.append('reacted')
        from cluster_analysis import test_internal_bonds

        # Filter out incorrectly bonded structures (where a chemical reaction has occured)
        self._filter(test_internal_bonds, {'components':'components'}, mol_df=self.mol_df, rel_tol=itol)

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

        func = lambda df, c, ct: [c in ct]*len(df)
        self._filter(func, {'c': 'unique_cluster_types'}, ct=cluster_types)
        return
    
    def except_clusters(self, cluster_types: Union[str, Iterable[str]]):
        if isinstance(cluster_types, str):
            cluster_types = [cluster_types]

        func = lambda df, c, ct: [c not in ct]*len(df)
        self._filter(func, {'c': 'unique_cluster_types'}, ct=cluster_types)
        return

    def set_Hbond(self, 
        donor: Optional[str] = None, 
        acceptor: Optional[str] = None, 
        length: Optional[tuple[float]] = None, 
        angle: Optional[tuple[float]] = None
    ):
        if not isinstance(self._Hbond_limits, DataFrame):
            from cluster_analysis import default_Hbond_limits
            self._Hbond_limits = default_Hbond_limits()
        
        if isinstance(donor, str) and isinstance(acceptor, str):
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
        return_stats: bool = False,
    ):
        self.history.append('Hbonded')
        if 'reacted' in self.history:
            return
        
        self.set_Hbond() # loads default H-bond limits if not done already
        
        # add proton donors and acceptors to mol_df
        from cluster_analysis import get_acceptors_and_donors, construct_cluster, test_Hbonds
        self.mol_df = get_acceptors_and_donors(self.mol_df, all_oxygens=True)

        self._rel_tol = rel_tol
        self._Hbond_options=[min_Hbonds, filter_type]
        self._cluster_info = []
        self.clusters_df[('temp', 'Hbond_pairs')] = object()
        self.clusters_df[('temp', 'Hbond_lengths')] = object()
        self.clusters_df[('temp', 'Hbond_angles')] = object()
        self._Hbond_counts = []
        self._Hbond_stats = []

        def func(clusters_df, components, **kwarks):
            # dict which contains info needed by filter about the cluster type
            if len(clusters_df) > 0:
                cluster_info = construct_cluster(self.mol_df, components, self._Hbond_limits)
                passed, Hbond_counts, Hbond_pairs, stats = test_Hbonds(clusters_df, cluster_info, **kwarks)
                
                # save these for later
                for i, idx in enumerate(clusters_df.index.values):
                    self.clusters_df.at[idx, ('temp', 'Hbond_pairs')] = tuple(Hbond_pairs[i][:2])
                    self.clusters_df.at[idx, ('temp', 'Hbond_lengths')] = Hbond_pairs[i][2]
                    self.clusters_df.at[idx, ('temp', 'Hbond_angles')] = Hbond_pairs[i][3]
                
                self._cluster_info.append(cluster_info)
                self._Hbond_counts.append(Hbond_counts)
                self._Hbond_stats.append(stats)
                return passed
            else:
                self._cluster_info.append([])
                self._Hbond_counts.append([])
                return []
        
        self._filter(func, {'components': 'components'}, rel_tol=self._rel_tol, 
                     options=self._Hbond_options, return_stats=return_stats)
        
        if return_stats:
            n = max([len(l[1]) for l in self._Hbond_stats])
            #m = max([max(l[1]) for l in self._Hbond_stats])+1
            from pandas import set_option
            set_option('display.max_columns', None)
            set_option('display.max_rows', None)
            set_option('display.width', None)
            
            count_table = np.array([np.pad(l[1], (0, n-len(l[1])), 'constant', constant_values=(0,0)) for l in self._Hbond_stats]).astype(str)
            for i, l in enumerate(self._Hbond_stats):
                if l[2] >= count_table.shape[1]:
                    count_table[i, -1] += '|'
                else:
                    count_table[i, l[2]] ='|'+count_table[i, l[2]] 

            df = DataFrame(data={'cluster': [l[0] for l in self._Hbond_stats], 'passed': [l[-1] for l in self._Hbond_stats]})
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
        n: int = 1,
        selected_isomers_file: Optional[str] = None,
        excepted_topologies_file: Optional[str] = None,
        sort_by: Union[str, tuple] = ('log','electronic_energy')
    ) -> None:
        if 'Hbonded' not in self.history:
            raise Exception('H-bond filtering must be run before filtering topologies.')

        from topologger import clusters_to_smiles, filter_isomorphs
        self.clusters_df[('temp', 'SMILES')] = object()
        for k, ct in enumerate(self.cluster_types):
            subset = self._cluster_subsets[k]
            if len(subset) == 0:
                continue
            self.clusters_df.loc[subset, ('temp', 'SMILES')] =\
                clusters_to_smiles(self.clusters_df.loc[subset], self._components[k], self.mol_df)
            
        self.file_iso = selected_isomers_file
        if isinstance(self.file_iso, str):
            from data_readers import read_pickled_data
            from topologger import filter_isomers
            
            self._isomers_df = read_pickled_data(self.file_iso)
            for mol in self.mol_df.index.values:
                subset = self._isomers_df[self._isomers_df[('info','cluster_type')]=='1'+mol]
                if len(subset) > 0:
                    self.mol_df.at[mol, 'isomers'] = clusters_to_smiles(subset, [mol], self.mol_df)
            
            for cluster_info in self._cluster_info:
                iso_smiles = self.mol_df.loc[cluster_info['components']]['isomers'].values
                cluster_info.update({'isomers': iso_smiles})

            self._filter(filter_isomers, {'cluster_info': 'cluster_info'})
        else:
            self._isomers_df = None
        if n>1:
            return

        self.file_not = excepted_topologies_file
        if isinstance(self.file_not, str):
            from data_readers import read_pickled_data
            from cluster_analysis import distance_matrix, test_Hbonds

            # extract already found topologies 
            used_topo_df = read_pickled_data(self.file_not)
            used_topo_df[("temp", "distances")] = distance_matrix(used_topo_df)
            self._used_pairs = [[]]*len(self.cluster_types)
            for k, ct in enumerate(self.cluster_types):
                subset = self._cluster_subsets[k]
                if len(subset) == 0:
                    continue
                used = used_topo_df[used_topo_df['info', 'cluster_type']==ct]
                if len(used) > 0:
                    f, _, used_topologies, _ = test_Hbonds(used, self._cluster_info[k], 
                                                           self._rel_tol, options=self._Hbond_options)
                    self._used_pairs[k] = np.array(used_topologies, dtype=tuple)[f]

            # filter unique topologies and remove already found topologies 
            self._filter(filter_isomorphs, {'cluster_info': 'cluster_info', 'used_DA_pairs': 'used_pairs'}, el=sort_by)
        else:
            # filter unique topologies 
            self._filter(filter_isomorphs, {'cluster_info': 'cluster_info'}, el=sort_by)
            
        return 
    
    def get_filtered_length(self) -> int:
        return sum([len(s) for s in self._cluster_subsets])
        
    # return filtered dataframe
    def get_filtered_data(self) -> DataFrame:
        idx = np.concatenate(self._cluster_subsets)
        df = self.clusters_df.loc[idx].sort_index(axis=1)
        if 'temp' in self.clusters_df.columns:
            df = df.drop('temp',axis=1)
        return df
    
    def save_to(self, file_out: str):
        filtered_clusters = self.get_filtered_data()[('info', 'file_basename')].values
        paths = DataFrame([l.split()[0].split('/:EXTRACT:/') for l in self._lines], columns=['file', 'cluster'])        
        passed = paths['cluster'].isin(filtered_clusters)
        
        with open(file_out, 'w') as f:
            f.writelines(self._lines[passed])

        return 
    
    def print_statistics():

        return
