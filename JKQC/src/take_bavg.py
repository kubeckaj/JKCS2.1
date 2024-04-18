def is_averagable(input_array): # :-D
  try:
    [float(i) for i in input_array]
    test = 0
    for i in range(len(input_array)):
      if input_array[i] != input_array[0]:
        test = 1
        break
    if test == 0:
      return False
    else:
      return True
  except ValueError:
    return False

def is_the_same(input_array):
  test = 0
  for i in range(len(input_array)):
    if input_array[i] != input_array[0]:
      test = 1
      break
  if test == 0:
    return input_array[0]
  else:
    return missing

def myif(cond,opt1,opt2,opt3,Pout):
  if Pout[cond] == "-g" or Pout[cond] == "-gout":
    return opt2
  elif Pout[cond] == "-s":
    return opt3
  else:
    return opt1

def take_bavg(output, clusters_df, Qbavg, Qt, QUenergy, Pout, Qclustername):
  from numpy import unique,min,exp,log,sum,array
  from pandas import isna
  missing = float("nan")
   
  k = 1.380662*10**-23 # [J/K]
  uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")].values)
  portions = []
  freeenergies = []
  entropies = []
  indexes = []
  for i in uniqueclusters:
    if Qbavg == 1:
      GFE = clusters_df.loc[clusters_df.loc[:,("info","cluster_type")] == i,("log","gibbs_free_energy")].values
    elif Qbavg == 2:
      GFE = clusters_df.loc[:,("log","gibbs_free_energy")] + clusters_df.loc[:,("out","electronic_energy")] - clusters_df.loc[:,("log","electronic_energy")]
      GFE = GFE.loc[clusters_df.loc[:,("info","cluster_type")] == i].values
    else:
      print("Qglob error. [EXITING]")
      exit()
    nonans = ~isna(GFE)
    GFE = GFE[nonans]
    try:
      minimum = min(GFE)
    except:
      minimum = missing
    GFE = GFE-minimum

    preportions = [exp(-GFE[i]*43.60*10**-19/k/Qt) for i in range(GFE.shape[0])]
    toaddfreeenergies = []
    try:
      addfreeenergy = QUenergy*(minimum - 1/43.60/10**-19*k*Qt*log(sum([exp(-GFE[i]*43.60*10**-19/k/Qt) for i in range(GFE.shape[0])])))
    except:
      addfreeenergy = missing
    freeenergies.append(addfreeenergy)
    sumpreportions = sum(preportions)
    portions.append(preportions/sumpreportions)
    indexes.append(array(range(len(clusters_df)))[clusters_df.loc[:,("info","cluster_type")] == i][nonans])

  newoutput = []
  skip = []
  for l in range(output.shape[0]):
    toappend = []
    for i in range(len(portions)):
      if i in skip:
        continue
      try:
        appendable = myif(l,sum([float(portions[i][j])*float(output[l,indexes[i][j]]) for j in range(len(portions[i]))]),freeenergies[i],missing,Pout) if is_averagable(output[l][indexes[i]]) else is_the_same(output[l,indexes[i]])
      except:
        appendable = missing
      if l == 0:
        if isna(appendable):
          skip.append(i)
          continue
      toappend.append(appendable)
    newoutput.append(array(toappend, dtype=object))
  output = array(newoutput)
  return output
