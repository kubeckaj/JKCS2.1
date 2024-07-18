def filter_sort(clusters_df,Qsort,Qreverse):
  if Qsort == "g":
    Qsort = "gibbs_free_energy"
  if Qsort == "el":
    Qsort = "electronic_energy"
  if Qsort == "gout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = (clusters_df.loc[:,("log","gibbs_free_energy")]-clusters_df.loc[:,("log","electronic_energy")]+clusters_df.loc[:,("out","electronic_energy")]).sort_values(ascending = Qreverse).index
      clusters_df = clusters_df.loc[sorted_indices]
    else:
      print("Cannot sort this")
      exit()
  elif Qsort == "elout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = clusters_df.loc[:,("out","electronic_energy")].sort_values(ascending = Qreverse).index
      clusters_df = clusters_df.loc[sorted_indices]
    else:
      print("Cannot sort this")
      exit()
  elif str(Qsort) == "errpa":
    err = []
    for ind in clusters_df.index:
      try:
        aseCL = clusters_df.loc[ind,("xyz","structure")]
        symb = aseCL.get_chemical_symbols()
        atoms = len(symb)
        err.append(float(atoms))
      except:
        err.append(float("nan"))
    sorted_indices = (clusters_df.loc[:,("extra","error")]/err).sort_values(ascending = Qreverse).index
    clusters_df = clusters_df.loc[sorted_indices]
  elif str(Qsort) == "no":
    clusters_df = clusters_df
  elif Qsort == "b":
    clusters_df = clusters_df.sort_values([('info','file_basename')], ascending = Qreverse)
  elif len(Qsort.split(",")) == 2:
    clusters_df = clusters_df.sort_values([(Qsort.split(",")[0],Qsort.split(",")[1])], ascending = Qreverse)
  else:
    clusters_df = clusters_df.sort_values([("log",Qsort)], ascending = Qreverse)
  return clusters_df
