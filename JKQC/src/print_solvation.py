def myFunc(e):
  try:
    numero = sum([int(i) for i in seperate_string_number(dash_comment(e[0])[0])[0::2]])
  except:
    numero = 0
  return numero+len(e[0])/100.

def seperate_string_number(string):
    # "1sa_215" -> ["1","sa","_","215"]
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i
        previous_character = i
        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups

def dash_comment(input_array):
  #removes comment
  #['1','sa','24','b','_','25',"-","26"] -> ['1','sa','24','b']
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "-" or input_array[element] == "_":
      output_array = []
      partbefore = input_array[0:element]
      output_array.append(partbefore)
      break
  return output_array

def print_solvation(output,Qsolvation,Qt,Qrh,QPsolvent,Pout,QUenergy,Qout):
  from numpy import apply_along_axis,transpose,array,unique,sum,exp,min
  from pandas import isna

  if Qout >= 1:
    print("#####################################",flush = True)
    print("##########  SOLVATION  ##############",flush = True)
    print("#####################################",flush = True)
  if len(Pout)>1:
    if Pout[1] != "-g" and Pout[1] != "-gout":
      print("Please use -ct -g or -ct -gout as the first two arguments.")
      exit()

  #SORT OUTPUT AND CLUSTERS
  output = transpose(transpose(output)[apply_along_axis(myFunc, axis=1, arr=transpose(output)).argsort()])
  cluster_types = [dash_comment(seperate_string_number(i))[0] for i in output[0]]
  #UNIQUE CLUSTERS
  no_solvent_cluster_types = []
  solvent_content = []
  for i in cluster_types:
    solvent_content_i = 0
    no_solvent_cluster_types_i = []
    for j in range(int(len(i)/2)):
      if i[2*j+1] == Qsolvation:
        solvent_content_i = i[2*j]
      else:
        no_solvent_cluster_types_i.append(i[2*j]) 
        no_solvent_cluster_types_i.append(i[2*j+1]) 
    no_solvent_cluster_types.append(no_solvent_cluster_types_i)
    solvent_content.append(solvent_content_i)
  cluster_types = ["".join(i) for i in cluster_types]
  no_solvent_cluster_types = ["".join(i) for i in no_solvent_cluster_types]
  unique_clusters = [x for x in array(unique(no_solvent_cluster_types), dtype=object) if x != ""]
  free_energies = output[1]
 
  #FIND SOLVENT
  index_of_solvent = -1
  for i in range(len(cluster_types)):
    if cluster_types[i] == "1"+Qsolvation:
      index_of_solvent = i
  if index_of_solvent == -1:
    print("Missing solvent")
    exit()

  #SET CONDITIONS
  if not isna(Qt):
    Temp = Qt
  else:
    Temp = 298.15
  if Qout >= 1:
    print(f"Temp = %.2f K; "%(Temp), end = "")
  if not isna(QPsolvent):
    p_solvent = QPsolvent
  else:
    if not isna(Qrh):
      rh = Qrh
    else:
      rh = 1.0
    if Qout >= 1:
      print(f"RH = %.2f %%; "%(rh*100), end = "") 
    #Actually not the Antoine equation - see the link
    #https://www.omnicalculator.com/chemistry/vapour-pressure-of-water
    p_solvent = rh*100*6.1094*exp(17.625*(Temp-273.15)/(Temp-273.15+243.04))
  if Qout >= 1:
    print(f"p_solvent = %.2f Pa "%(p_solvent))
  p_ref = 101325
  R = 1.987 #cal/mol/K #=8.31441

  #PRINT SOLVATION
  unique_clusters.append('')
  for i in unique_clusters:
    indexes = []
    for j in range(len(no_solvent_cluster_types)):
      if i == no_solvent_cluster_types[j]:
        indexes.append(j)
    tot_conc = 0
    nominators = []
    free_energies_i = []
    for j in indexes:
      free_energies_i.append(output[1][j]-float(solvent_content[j])*output[1][index_of_solvent])
    minimum = min(free_energies_i)
    free_energies_i = free_energies_i - minimum
    for j in range(len(indexes)):
      nominator = (p_solvent/p_ref)**float(solvent_content[indexes[j]])*exp(-627.503/QUenergy*free_energies_i[j]/R*1000/Temp)
      nominators.append(nominator)
    denominator = sum(nominators)
    
    for j in range(len(indexes)):
      print(f"%8s "%(cluster_types[indexes[j]]), end="")
    print("")
    for j in range(len(indexes)):
      print(f"%8.1f "%(nominators[j]/denominator*100),end="")
    print("")
    if Qout >= 1:
      print("----------------------------------------")
    else:
      print("")
