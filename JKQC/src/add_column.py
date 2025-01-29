def add_column(clusters_df,addcolumn):
  from numpy import dtype,array
  
  dtype = dtype([("label", "U100"), ("value", object)])

  for i in range(len(addcolumn)):
    if addcolumn[i][0] == "SP" and addcolumn[i][1][-3:] == "pkl":
      from pandas import read_pickle
      loadedclustersdf = read_pickle(addcolumn[i][1])
      loadedmatrix = array([loadedclustersdf.loc[:,("info","file_basename")],loadedclustersdf.loc[:,("log","electronic_energy")]])
      loadeddictionary = {loadedmatrix[0][idx] : loadedmatrix[1][idx] for idx in range(len(loadedmatrix[1]))}
      tobeadded = clusters_df.loc[:,("info",'file_basename')].map(loadeddictionary)
      clusters_df.loc[tobeadded.index,("out","electronic_energy")] = tobeadded.values
    else:
      from numpy import loadtxt
      loadedmatrix = loadtxt(addcolumn[i][1], usecols=(0, 1), unpack=True, dtype=dtype)
      try:
        loadeddictionary = {loadedmatrix[0][idx] : loadedmatrix[1][idx] for idx in range(len(loadedmatrix[1]))}
      except:
        loadeddictionary = {loadedmatrix[0].item() : loadedmatrix[1].item()}
      tobeadded = clusters_df.loc[:,("info",'file_basename')].map(loadeddictionary)
      clusters_df.loc[tobeadded.index,("extra",addcolumn[i][0])] = tobeadded.values
 
  return clusters_df
