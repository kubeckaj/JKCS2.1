def seperate_string_number(string):
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i == "_" or previous_character == "_":
            newword += i
        elif i.isalpha() and previous_character.isalpha():
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
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "-" or input_array[element] == "_":
      output_array = []
      partbefore = input_array[0:element]
      output_array.append(partbefore)
      break
  return output_array

def adjust_monomers(chosen_cluster_type,old_monomers,new_monomers):
  #test whether the monomers are really present (with positive numbers)
  for i in range(1,len(old_monomers),2):
    if old_monomers[i] in chosen_cluster_type[1::2]:
      if int(chosen_cluster_type[chosen_cluster_type.index(old_monomers[i])-1])-int(old_monomers[i-1]) >= 0:
        continue
      else:
        return chosen_cluster_type
    else:
      return chosen_cluster_type
  new_cluster_type = []
  for j in range(1,len(chosen_cluster_type),2):
    if chosen_cluster_type[j] in old_monomers[1::2]:
      subsctracted = int(chosen_cluster_type[j-1]) - int(old_monomers[old_monomers.index(chosen_cluster_type[j])-1])
      if chosen_cluster_type[j] in new_monomers[1::2]:
        subsctracted += int(new_monomers[new_monomers.index(chosen_cluster_type[j])-1])
        del new_monomers[new_monomers.index(chosen_cluster_type[j])-1:new_monomers.index(chosen_cluster_type[j])+1]
      if subsctracted != 0:
        new_cluster_type.append(str(subsctracted))
        new_cluster_type.append(chosen_cluster_type[j])
    elif chosen_cluster_type[j] in new_monomers[1::2]:
      added = int(chosen_cluster_type[j-1]) + int(new_monomers[new_monomers.index(chosen_cluster_type[j])-1])
      del new_monomers[new_monomers.index(chosen_cluster_type[j])-1:new_monomers.index(chosen_cluster_type[j])+1]
      if added != 0:
        new_cluster_type.append(str(added))
        new_cluster_type.append(chosen_cluster_type[j])
    else:
      new_cluster_type.append(chosen_cluster_type[j-1])
      new_cluster_type.append(chosen_cluster_type[j])
  for j in range(1,len(new_monomers),2): 
    new_cluster_type.append(new_monomers[j-1])
    new_cluster_type.append(new_monomers[j])
  return adjust_monomers(new_cluster_type,old_monomers,new_monomers)

def myFunc(e):
  from numpy import sum
  try:
    numero = sum([int(i) for i in seperate_string_number(dash_comment(e[0])[0])[0::2]])
  except:
    numero = 0
  return numero+len(e[0])/100.

##########################################################################################
##########################################################################################
##########################################################################################

def print_formation(output, Qout=1, Qt = 298.15, Qp = 101325.0, Qconc = 0, conc = [], CNTfactor = 0, QUenergy = 1):
  """print formation/binding properties
  output = np.array of output
  """
  from numpy import transpose,apply_along_axis,array,sum,dtype
  missing = float("nan") 
  if Qconc > 0:
    from numpy import log
    R = 1.987 #cal/mol/K #=8.31441
 
  if Qout >= 1:
    print("#####################################",flush = True)
    print("##########  FORMATION  ##############",flush = True)
    print("#####################################",flush = True)
  #SORT OUTPUT
  output = transpose(transpose(output)[apply_along_axis(myFunc, axis=1, arr=transpose(output)).argsort()])
  cluster_types = [dash_comment(seperate_string_number(i))[0] for i in output[0]]
  cluster_names = ["".join(dash_comment(seperate_string_number(i))[0]) for i in output[0]]

  ########################        
  ### SOLVING PROTONATION
  ########################        
  monomer_types = list(set([item for sublist in cluster_types for item in sublist[1::2]]))
  all_positive = []
  possible_positive = ["gd","eda","tma","dma","am","gly","buta","dbma","dea","dhexa","dmea","dpenta","dpropa","ibuta","IIebuta","ipropa","ipropea","mea","nona","propa","sbuta","tbuta","tea","tibuta","tpropa","w"]
  for possible_positive_i in range(len(possible_positive)):
    if possible_positive[possible_positive_i] in monomer_types:
      if "1"+possible_positive[possible_positive_i]+"1p" in cluster_names or "1p1"+possible_positive[possible_positive_i] in cluster_names:
        all_positive.append(possible_positive[possible_positive_i])
  if len(all_positive) > 0:
    for j in range(0,len(all_positive)):
      for i in range(len(cluster_types)):
        if j == 0:
          cluster_types[i] = adjust_monomers(cluster_types[i],["1",all_positive[j],"1","p"],["1","protonated_"+all_positive[0]])
        else:
          cluster_types[i] = adjust_monomers(cluster_types[i],["1",all_positive[j],"1","p"],["1","protonated_"+all_positive[0],"-1",all_positive[0],"1",all_positive[j]])
  ########################        
  ## SOLVING DEPROTONATION
  ########################        
  monomer_types = list(set([item for sublist in cluster_types for item in sublist[1::2]]))
  all_negative = []
  all_negatives_neutral = []
  possible_negative = ["b","mb","glyt","acc","bzc","cc","fc","mbtcc","oc","pc","suc","tbc","tpc","ttc"]
  possible_negatives_neutral = ["sa","msa","gly","aca","bza","ca","fa","mbtca","oa","pa","sua","tba","tpa","tta"]
  for possible_negative_i in range(len(possible_negative)):
    if possible_negative[possible_negative_i] in monomer_types:
      all_negative.append(possible_negative[possible_negative_i])
      all_negatives_neutral.append(possible_negatives_neutral[possible_negative_i ]) 
  if len(all_negative) > 1:
    for j in range(1,len(all_negative)):
      for i in range(len(cluster_types)):
        cluster_types[i] = adjust_monomers(cluster_types[i],["1",all_negative[j]],["1",all_negative[0],"-1",all_negatives_neutral[0],"1",all_negatives_neutral[j]])
  ########################        
  ########################        
  ########################        

  cluster_types_sorted = [sorted([extract_i[i:i + 2] for i in range(0, len(extract_i), 2)],key=lambda x: x[1]) for extract_i in cluster_types]
  monomers = array([sum([abs(int(j)) for j in array(i)[:,0]]) for i in cluster_types_sorted]) == 1
  dt = dtype(object)
  monomer_nonunique_types = [i[1] for i in array(cluster_types,dtype=dt)[monomers]]
  monomer_types = list(set(monomer_nonunique_types))
  for i in monomer_types:
    found = 0
    for j in range(len(monomer_nonunique_types)):
      if i == monomer_nonunique_types[j]:
        found += 1
        if found > 1:
          monomers[j] = False
  if Qout >= 1:
    print("ANCHOR MONOMERS: " + " ".join(monomer_types),flush = True)
  new_output = []
  #TODO slow for large number of files
  for i in range(len(output[0])):
    line = array(output)[:,i]
    cluster_total_number = sum([int(sel[0]) for sel in cluster_types_sorted[i]])
    for j in range(len(cluster_types_sorted[i])):
      cluster_molecule = cluster_types_sorted[i][j][1]
      cluster_molecule_number = cluster_types_sorted[i][j][0]
      test_monomer = 0
      #for k in range(len(array(output)[:,monomers][0])):
      for k in range(sum(monomers)):
        selected_monomer = array(cluster_types_sorted,dtype=dt)[monomers][k][0][1]
        if cluster_molecule == selected_monomer:
          for line_i in range(1,len(line)):
            if type(line[line_i]) != type("str"):
              try:
                line[line_i] = float(line[line_i]) - float(cluster_molecule_number) * float(array(output)[:,monomers][line_i,k])
                if Qconc > 0:
                  for conc_j in range(len(conc)):
                    if conc[conc_j][0] == selected_monomer:         
                      conc_mon=float(eval(conc[conc_j][1].replace("ppt","*10**-12*"+str(Qp)).replace("ppb","*10**-9*"+str(Qp)).replace("^","**").replace("cmmc","*10**6*1.380662*10**-23*"+str(Qt)) ))
                      line[line_i] = float(line[line_i]) - QUenergy*(float(cluster_molecule_number) - CNTfactor*float(cluster_molecule_number)/cluster_total_number) * R/1000/627.503 * Qt * log( conc_mon / Qp)
              except:
                line[line_i] = missing
          test_monomer = 1
      if test_monomer == 0:
        line[1:] = missing 
    new_output.append(line)
  new_output = transpose(array(new_output))
  toprint = list(zip(*new_output)) #[row for row in list(zip(*out))]
  if len(toprint) > 0:
    column_widths = [max(len(str(row[i])) for row in toprint) for i in range(len(toprint[0]))]
    for row in toprint:
      formatted_row = [str(row[i]).ljust(column_widths[i]) for i in range(len(row))]
      print(" ".join(formatted_row),flush = True)

