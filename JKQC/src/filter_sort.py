def filter_sort(clusters_df,Qsort):
  if Qsort == "g":
    Qsort = "gibbs_free_energy"
  if Qsort == "el":
    Qsort = "electronic_energy"
  if Qsort == "gout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = (clusters_df.loc[:,("log","gibbs_free_energy")]-clusters_df.loc[:,("log","electronic_energy")]+clusters_df.loc[:,("out","electronic_energy")]).sort_values().index
      clusters_df = clusters_df.loc[sorted_indices]
    else:
      print("Cannot sort this")
      exit()
  elif Qsort == "elout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = clusters_df.loc[:("out","electronic_energy")].sort_values().index
      clusters_df = clusters_df.loc[sorted_indices]
    else:
      print("Cannot sort this")
      exit()
  elif str(Qsort) == "no":
    clusters_df = clusters_df
  elif Qsort == "b":
    clusters_df = clusters_df.sort_values([('info','file_basename')])
  else:
    clusters_df = clusters_df.sort_values([("log",Qsort)])
  return clusters_df
