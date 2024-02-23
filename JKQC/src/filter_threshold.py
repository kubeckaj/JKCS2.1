def filter_threshold(clusters_df,Qcut,Qclustername,Qout):
  from numpy import array, tile, errstate, unique, sum
  from pandas import isna
  original_length = len(clusters_df)


  if Qclustername != 0:
    uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")].values)
  else:
    uniqueclusters = "1"
  newclusters_df = []

  myNaN = lambda x : missing if x == "NaN" else x
  for i in uniqueclusters:
    if Qclustername != 0:
      preselected_df = clusters_df.loc[clusters_df.loc[:,("info","cluster_type")] == i]
    else:
      preselected_df = clusters_df

    for i in range(len(Qcut)):
      if Qcut[i][1] == "relative":
        multiplier = 627.503
      else:
        multiplier = 1
      if Qcut[i][2] == "el":
        what = multiplier*preselected_df.loc[:,("log","electronic_energy")].values
      elif Qcut[i][2] == "g":
        what = multiplier*preselected_df.loc[:,("log","gibbs_free_energy")].values
      elif Qcut[i][2] == "elout":
        if not ("out","electronic_energy") in preselected_df.columns:
          print("The out-el.energy column is not present.")
          exit()
        what = multiplier*preselected_df.loc[:,("out","electronic_energy")].values
      elif Qcut[i][2] == "gout":
        if not ("out","electronic_energy") in preselected_df.columns:
          print("The out-el.energy column is not present.")
          exit()
        what = multiplier*(preselected_df.loc[:,("log","gibbs_free_energy")].values-preselected_df.loc[:,("log","electronic_energy")].values+preselected_df.loc[:,("out","electronic_energy")].values)
      elif Qcut[i][2] == "lf":
        what = array([array(ii) if isna([ii]).any() else array(ii[0]) for ii in preselected_df.loc[:,("log","vibrational_frequencies")].values], dtype=object)
      elif Qcut[i][2] == "rg":
        rg = []
        for aseCL in preselected_df.loc[:,("xyz","structure")].values:
          try:
            rg.append((sum(sum((aseCL.positions-tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/sum(aseCL.get_masses()))**0.5)
          except:
            rg.append(missing)
        what = array(rg)
      else:
        what = preselected_df.loc[:,("log",Qcut[i][2])].values
      if Qcut[i][1] == "relative":
        minimum = min(what)
      else:
        minimum = 0
      if Qcut[i][3] == "nan" or Qcut[i][3] == "NA" or Qcut[i][3] == "na" or Qcut[i][3] == "NaN":
        if Qcut[i][0] == ">":
          preselected_df = preselected_df.loc[isna(what-minimum)]
        else:
          preselected_df = preselected_df.drop(index=preselected_df[isna(what-minimum)].index)
      else:
        with errstate(invalid='ignore'):
          if Qcut[i][0] == ">":
            preselected_df = preselected_df.loc[what-minimum > float(Qcut[i][3])]
          else:
            preselected_df = preselected_df.loc[what-minimum <= float(Qcut[i][3])]
      print(preselected_df)
      #APPEND SELECTED
      if len(newclusters_df) == 0:
        newclusters_df = preselected_df.copy()
      else:
        newclusters_df = newclusters_df.append(preselected_df.copy())

  if Qout >= 1:
    print("Threshold filtering: "+str(original_length)+" --> "+str(len(newclusters_df)))
  return newclusters_df
