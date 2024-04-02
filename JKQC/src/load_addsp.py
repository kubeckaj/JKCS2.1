def load_addsp(clusters_df,input_pkl_sp,Qout):
  from pandas import read_pickle  
  from numpy import array

  for i in range(len(input_pkl_sp)):
    newclusters_df_sp = read_pickle(input_pkl_sp[i])
    if Qout >= 1:
      print("Number of files in "+input_pkl_sp[i]+": "+str(len(newclusters_df_sp)))
    if "log" in newclusters_df_sp:
      newclusters_df_sp = newclusters_df_sp.rename(columns={"log": "out"})
    if "out" in clusters_df:
      clusters_df = clusters_df.drop(["out"], axis=1, level=0)
    index_old = []
    index_new = []
    basenames_old = list(clusters_df.loc[:,('info', 'file_basename')].values)
    basenames_new = list(newclusters_df_sp.loc[:,('info', 'file_basename')].values)
    if len(basenames_old) != len(list(set(basenames_old))):
      print("The -addSP database has multiple occurances")
      exit()
    if len(basenames_new) != len(list(set(basenames_new))):
      print("The original database(s) has multiple occurances")
      exit()
    for ii in range(len(basenames_old)):
      try:
        index_new.append(basenames_new.index(basenames_old[ii]))
        index_old.append(ii)
      except:
        continue
    for ii in newclusters_df_sp.loc[:,"out"].columns.values:
      toadd = array(newclusters_df_sp.loc[index_new,("out",ii)], dtype=object)
      clusters_df.at[index_old,("out",ii)]=toadd

  return clusters_df
