def take_glob(output, clusters_df, Qglob):
  from numpy import unique,array,min
  from pandas import isna
  missing = float("nan")

  uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")])
  indexes = []
  for i in uniqueclusters:
    indexes_i = clusters_df.loc[:,("info","cluster_type")] == i
    if Qglob == 1:
      GFE = clusters_df.loc[indexes_i,("log","gibbs_free_energy")].values
    elif Qglob == 2:
      GFE = clusters_df.loc[indexes_i,("log","gibbs_free_energy")] + clusters_df.loc[indexes_i,("out","electronic_energy")] - clusters_df.loc[indexes_i,("log","electronic_energy")]
    else:
      print("Qglob error. [EXITING]")
      exit()
    GFE = array([missing if str(x) == "NaN" else x for x in GFE])
    if len(GFE[~isna(GFE)]) != 0:
      globindex = clusters_df.index.values[clusters_df.loc[:,("info","cluster_type")] == i][GFE == min(GFE[~isna(GFE)])][0]
      indexes.append(int(array(range(len(clusters_df)))[clusters_df.index.values == globindex][0]))

  newoutput = []
  for j in range(output.shape[0]):
    toappend=[]
    for i in indexes:
      toappend.append(output[j][i])
    newoutput.append(array(toappend,dtype=object))
  output = array(newoutput)
  return output
