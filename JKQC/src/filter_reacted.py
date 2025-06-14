def bonded(x,xA,xB, bonddistancethreshold = 2.0):
  #Are two atoms bonded?
  if xA == "C" and xB == "N" and x < 1.75:
    return 1
  elif xA == "N" and xB == "C" and x < 1.75:
    return 1
  elif xA == "N" and xB == "N" and x < 1.5:
    return 1
  elif xA == "S" and xB == "O" and x < 1.9:
    return 1
  elif xA == "O" and xB == "S" and x < 1.9:
    return 1
  elif xA == "O" and xB == "N" and x < 1.9:
    return 1
  elif xA == "N" and xB == "O" and x < 1.9:
    return 1
  elif xA == "H" or xB == "H":
    if xA == "H" and xB == "H" and x < 0.8:
      return 1
    else:
      return 0
  elif x < bonddistancethreshold:
    return 1
  else:
    return 0

def most_frequent(List):
  return max(set(List), key = List.count)

### REACTED ###
def filter_reacted(clusters_df, Qclustername = 0, Qreacted = 1, bonddistancethreshold = 2.0, Qout = 1):
  """for removing reacting structures
  clusters_df = JKQC pandas dataframe
  """
  from pandas import isna
  from numpy import sort, array, unique, dtype, sqrt, sum, asarray
  missing = float("nan")
  
  dt = dtype(object)
  #Are there some cluster types which I should distinguish?
  if Qclustername != 0:
    unique_cluster_types = unique(clusters_df.loc[:,("info","cluster_type")])
    cluster_subsets = []
    for unique_cluster_type in unique_cluster_types:
      indexes = clusters_df[clusters_df.loc[:,("info","cluster_type")]==unique_cluster_type].index.values
      cluster_subsets.append(indexes)
  else:
    cluster_subsets = []
    indexes = clusters_df.index.values
    cluster_subsets.append(indexes)   
  
  #Loop over all unique cluster subsets
  dist = lambda p1, p2: sqrt(sum(((p1-p2)**2)))
  for k0 in range(len(cluster_subsets)):
    all_molecules = []
    #Loop over the configuration of k0 subset
    for k in range(len(cluster_subsets[k0])):
      b = clusters_df.loc[cluster_subsets[k0][k],("xyz","structure")]
      if isna(b):
        all_molecules.append(str(missing))
      else:
        p = b.positions
        symb = array(b.get_chemical_symbols())
        ind = [i != 'test' for i in symb] 
   
        #distance matrix and bond matrix 
        dm = asarray([[dist(p1, p2) for p2 in p[ind]] for p1 in p[ind]])
        bm = array([[bonded(dm[i][j],symb[ind][i],symb[ind][j],bonddistancethreshold) for j in range(len(dm[i]))] for i in range(len(dm))])

        test = 0
        choosing_list = range(len(bm))
        molecules=[]
        if len(choosing_list) == 0:
          test = 1
        while test == 0:
          selected = [choosing_list[0]]
          test_chosen = 0
          j = -1
          while test_chosen == 0:
            j += 1
            chosen = selected[j]
            for i in choosing_list:
              if bm[chosen][i] == 1 and chosen != i and (not (i in selected)):
                selected.append(i)
            if len(selected)-1 == j:
              test_chosen = 1
          molecules.append("".join(sort([symb[ind][i] for i in selected])))
          choosing_list = [i for i in choosing_list if i not in selected]
          if len(choosing_list) == 0:
            test = 1
        all_molecules.append(str(sort(array(sort(molecules),dtype = dt))))
  
    
    mf = most_frequent(all_molecules)
    
    if Qreacted != 2:
      nind = [ i == mf for i in all_molecules ]
    else:
      nind = [ i != mf for i in all_molecules ]
    if k0 == 0:
      selclusters_df = clusters_df.loc[cluster_subsets[0][nind]].copy()
    else:
      from pandas import concat
      selclusters_df = concat([selclusters_df, clusters_df.loc[cluster_subsets[k0][nind]].copy()], ignore_index=True)
      #selclusters_df = selclusters_df.append(clusters_df.loc[cluster_subsets[k0][nind]].copy())
  if Qout >= 1:
    original_length = len(clusters_df)
  clusters_df = selclusters_df.copy()
  if Qout >= 1:
    print("Removing reacted: "+str(original_length)+" --> "+str(len(clusters_df)))
  return clusters_df
