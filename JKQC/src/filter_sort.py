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
  elif str(Qsort) == "rg":
    from numpy import array, tile, errstate, unique, sum
    rg = []
    for aseCL in clusters_df.loc[:,("xyz","structure")].values:
      try:
        rg.append((sum(sum((aseCL.positions-tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/sum(aseCL.get_masses()))**0.5)
      except:
        rg.append(float("nan"))
    #sorted_indices = rg.sort_values(ascending = Qreverse).index
    #sorted_indices = sorted(range(len(rg)), key=lambda x: rg[x])
    #HA HA THIS IS SO STUPID CODE
    #print(rg)
    from pandas import DataFrame
    df = DataFrame({"rg":rg},index=range(len(rg)))
    df = df.loc[:,("rg")].sort_values(ascending = Qreverse)
    sorted_indices = df.index
    #print(DataFrame({"rg":rg},index=range(len(rg))).loc[:,"rg"])
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
