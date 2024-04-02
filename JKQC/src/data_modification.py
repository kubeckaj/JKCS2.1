def is_nameable(input_array):
  from numpy import mod

  nameable_test = True
  if mod(len(input_array),2) == 0:
    for input_array_i in input_array[0::2]:
      if not input_array_i.isnumeric():
        nameable_test = False
        break
    for input_array_i in input_array[1::2]:
      if input_array_i.isnumeric():
        nameable_test = False
        break
  else:
    nameable_test = False
  return nameable_test

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

def zeros(input_array):
  output_string = ""
  skip = 0
  for i in range(len(input_array)):
    if skip == 1:
      skip = 0
      continue
    if input_array[i] == "0":
      skip = 1
      continue
    output_string += input_array[i]
  return output_string

def replace_first_occurrence(string, old_char, new_char):
  index = string.find(old_char)  # Find the index of the first occurrence
  if index != -1:  # If the character is found
   string = string[:index] + new_char + string[index+1:]  # Replace the character
  return string

def data_modification(clusters_df, Qunderscore, Qrename, Qclustername, QrenameWHAT, Qiamifo, Qrebasename, Qdrop, Qout2log,Qpresplit,Qindex):

  #SPLIT THE DATABASE AT START  
  if Qpresplit > 0:
    clusters_df = clusters_df.sample(n=Qpresplit, random_state=42)
  
  #SELECT INDEX
  if Qindex != "-1":
    try:
      clusters_df = eval(f"clusters_df[{Qindex}]")
    except:
      try:
        Qindex = Qindex+":"+str(int(Qindex)+1)
        clusters_df = eval(f"clusters_df[{Qindex}]")
      except:
        print("Error selecting given index")
        exit()

  if Qrename == 1:
    if Qclustername == 1:
      import re
      from numpy import array
      l1 = []
      l2 = []
      l3 = []
      l4 = []
      for cluster_id in clusters_df.index:
        file_i_BASE = clusters_df.loc[cluster_id]['info']['file_basename']
        cluster_type_array = seperate_string_number(file_i_BASE.split("-")[0].split("_")[0])
        if is_nameable(cluster_type_array):
          for QrenameONE in QrenameWHAT:
            cluster_type_array = [w.replace(QrenameONE[0],QrenameONE[1]) if QrenameONE[0] == w else w for w in cluster_type_array]
          l4.append("".join(cluster_type_array)+file_i_BASE.replace(file_i_BASE.split("-")[0].split("_")[0],""))
          cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
          cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
          cluster_type_string = zeros(cluster_type_array_sorted)
          components = cluster_type_array[1::2]
          component_ratio = [int(i) for i in re.split('(\d+)', file_i_BASE.split("-")[0])[1:][0::2]]
          l1.append(cluster_type_string)
          l2.append(components)
          l3.append(component_ratio)
      clusters_df.loc[:,("info", "cluster_type")] =  array(l1)
      clusters_df.loc[:,("info", "components")] = array(l2, dtype=object)
      clusters_df.loc[:,("info", "component_ratio")] = array(l3, dtype=object)
      clusters_df.loc[:,("info", "file_basename")] = array(l4)

  if Qunderscore == 1:
    for cluster_id in clusters_df.index:
      clusters_df.loc[cluster_id,('info','file_basename')] = replace_first_occurrence(clusters_df.loc[cluster_id,('info','file_basename')],"_","-")

  if Qiamifo > 0:
    predashstrings = [i.split("-")[0] for i in clusters_df.loc[:,("info","file_basename")].values]
    clusters_df.loc[:,("info", "cluster_type")] = predashstrings

  if Qrebasename == 1:
    from numpy import delete
    values=clusters_df.loc[:,("info","file_basename")].values
    for i in range(len(clusters_df)):
      if values[i] in delete(values, i, axis=0):
        version=1
        values=clusters_df.loc[:,("info","file_basename")].values
        while values[i]+"-v"+str(version) in values:
          version += 1
        clusters_df.at[i,("info","file_basename")] = values[i]+"-v"+str(version)

  if Qdrop != "0":
    if Qdrop in clusters_df:
      clusters_df = clusters_df.drop([Qdrop], axis=1, level=0)
  
  if Qout2log == 1:
    clusters_df = clusters_df.rename({'out': 'log'}, axis='columns')

  return clusters_df
