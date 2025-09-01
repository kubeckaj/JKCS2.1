def filter_select(clusters_df,Qselect,Qclustername,Qout):
  from numpy import unique

  if len(clusters_df) == 0:
    print("Empty database.")
    exit()
  if Qclustername != 0:
    uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")].values)
  else:
    uniqueclusters = "1"
  newclusters_df = []
  for i in uniqueclusters:
    if Qclustername != 0:
      selected_df = clusters_df.loc[clusters_df.loc[:,("info","cluster_type")] == i][0:Qselect]
    else:
      selected_df = clusters_df.iloc[0:Qselect]
    if len(newclusters_df) == 0:
      newclusters_df = selected_df
    else:
      from pandas import concat
      #newclusters_df = concat([newclusters_df,selected_df.copy()], ignore_index=True)
      newclusters_df = concat([newclusters_df,selected_df.copy()])
      #newclusters_df = newclusters_df.append(selected_df)

  if Qout >= 1:
    print("Selecting/Sampling: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
  return newclusters_df
