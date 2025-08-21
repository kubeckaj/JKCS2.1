def filter_threshold(clusters_df,Qcut,Qclustername,Qout):
  from numpy import array, tile, errstate, unique, sum
  from pandas import isna
  missing = float("nan")
  original_length = len(clusters_df)


  if Qclustername != 0:
    uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")].values)
  else:
    uniqueclusters = "1"
  newclusters_df = []

  bonded_incr = 0
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
      elif Qcut[i][2] == "errpa":
        errpa = []
        for ind in clusters_df.index:
          try:
            aseCL = clusters_df.loc[ind,("xyz","structure")]
            symb = aseCL.get_chemical_symbols()
            err =  clusters_df.loc[ind,("extra","error")] 
            atoms = len(symb)
            errpa.append(float(err/atoms))
          except:
            errpa.append(float("nan"))
        what = array(errpa)
      elif Qcut[i][2] == "lf":
        what = array([array(ii) if (isna([ii]).any() or type(ii)==type(array([]))) else array(ii[0]) for ii in preselected_df.loc[:,("log","vibrational_frequencies")].values], dtype=object)
        #what = array([array(ii) if isna([ii]).any() else array(ii[0]) for ii in preselected_df.loc[:,("log","vibrational_frequencies")].values], dtype=object)
      elif Qcut[i][2] == "rg":
        rg = []
        for aseCL in preselected_df.loc[:,("xyz","structure")].values:
          try:
            rg.append((sum(sum((aseCL.positions-tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/sum(aseCL.get_masses()))**0.5)
          except:
            rg.append(missing)
        what = array(rg)
      elif Qcut[i][2] == "bonded":
        bonded = []
        for ind in clusters_df.index:
          try:
            aseCL = clusters_df.loc[ind,("xyz","structure")]
            positions = aseCL.positions
            symb = aseCL.get_chemical_symbols()
            dist = lambda p1, p2: sum((p1-p2)**2)**0.5
            symb_ind = array(aseCL.get_chemical_symbols())
            mask1 = symb_ind == Qcut[i][3][1]
            mask2 = symb_ind == Qcut[i][3][2]
            dm = [dist(p1, p2) for p2 in positions[mask1] for p1 in positions[mask2]]
            bonds = sum(test_i <= float(Qcut[i][3][0]) for test_i in dm)
            if str(Qcut[i][3][1]) == str(Qcut[i][3][2]):
              bonds = (bonds-sum(mask1))/2
            bonded.append(int(bonds))
          except:
            bonded.append(float("nan"))
        what = array(bonded)       
        bonded_incr = 1
      elif len(Qcut[i][2].split(",")) == 2:
        what = preselected_df.loc[:,(Qcut[i][2].split(",")[0],Qcut[i][2].split(",")[1])].values
      else:
        what = preselected_df.loc[:,("log",Qcut[i][2])].values
      if Qcut[i][1] == "relative":
        minimum = min(what)
      else:
        minimum = 0
      if Qcut[i][3+bonded_incr] == "nan" or Qcut[i][3+bonded_incr] == "NA" or Qcut[i][3+bonded_incr] == "na" or Qcut[i][3+bonded_incr] == "NaN":
        if Qcut[i][0] == "==" or Qcut[i][0] == ">=" or Qcut[i][0] == "<=":
          #mask = preselected_df[what].apply(lambda x: pd.isna(x) or (hasattr(x, "__len__") and len(x) == 0))
          #preselected_df = preselected_df.loc[mask]
          preselected_df = preselected_df.loc[isna(what-minimum)]
        else:
          preselected_df = preselected_df.drop(index=preselected_df[isna(what-minimum)].index)
          #mask = preselected_df[array(what)].apply(lambda x: pd.isna(x) or (hasattr(x, "__len__") and len(x) == 0))
          #preselected_df = preselected_df.loc[~mask]
      else:
        with errstate(invalid='ignore'):
          if Qcut[i][0] == ">":
            preselected_df = preselected_df.loc[what-minimum > float(Qcut[i][3+bonded_incr])]
          elif Qcut[i][0] == ">=":
            preselected_df = preselected_df.loc[what-minimum >= float(Qcut[i][3+bonded_incr])]
          elif Qcut[i][0] == "==":
            preselected_df = preselected_df.loc[what-minimum == float(Qcut[i][3+bonded_incr])]
          elif Qcut[i][0] == "!=":
            preselected_df = preselected_df.loc[what-minimum != float(Qcut[i][3+bonded_incr])]
          elif Qcut[i][0] == "<":
            preselected_df = preselected_df.loc[what-minimum < float(Qcut[i][3+bonded_incr])] 
          elif Qcut[i][0] == "<=":
            preselected_df = preselected_df.loc[what-minimum <= float(Qcut[i][3+bonded_incr])]
          else:
            print("Something weird with the filtering")
            exit()
    #APPEND SELECTED
    if len(newclusters_df) == 0:
      newclusters_df = preselected_df.copy()
    else:
      from pandas import concat
      #newclusters_df = concat([newclusters_df,preselected_df.copy()], ignore_index=True)
      newclusters_df = concat([newclusters_df,preselected_df.copy()])
      #newclusters_df = newclusters_df.append(preselected_df.copy())

  if Qout >= 1:
    print("Threshold filtering: "+str(original_length)+" --> "+str(len(newclusters_df)))
  return newclusters_df
