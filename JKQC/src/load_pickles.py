def load_pickles(input_pkl,Qout,Qid):
  from pandas import DataFrame
  if len(input_pkl) == 0:
    clusters_df = DataFrame()
  else:
    from pandas import read_pickle
    import gc
    for i in range(len(input_pkl)):
      newclusters_df = read_pickle(input_pkl[i])
      if not isinstance(newclusters_df, DataFrame):
        print("File "+input_pkl[i]+" is not JKQC-compatible Pandas.DataFrame. Try to use JKTS.")
        exit()
      if Qout >= 1 and len(input_pkl) > 1:
        print("Number of files in "+input_pkl[i]+": "+str(len(newclusters_df)))
      if i == 0:
        clusters_df = newclusters_df
      else:
        from pandas import concat
        clusters_df = concat([clusters_df,newclusters_df.copy()], ignore_index=True)
        #clusters_df = clusters_df.append(newclusters_df, ignore_index=True)
      gc.collect()
    if Qout >= 2:
      print("Pickles collected, resetting index...")
    clusters_df = clusters_df.reset_index(drop=True)
    if Qout >= 2:
      print("Checking xyz_id1...")

  if Qid == 1 and not ("xyz","id1") in clusters_df.columns and ("xyz","structure") in clusters_df.columns:
    if Qout >= 2:
      print("Adding id1...") 
    from functions import df_add_iter
    from read_xyz import identify_1
    variables = [ identify_1(ase_ID) for ase_ID in clusters_df.loc[:,("xyz","structure")] ]
    clusters_df = df_add_iter(clusters_df,"xyz","id1", range(len(clusters_df)),variables) 
  if Qout >= 2:
    print("Sending back to the main file.")
  return clusters_df
