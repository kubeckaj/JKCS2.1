
# EXTRACT
def listToString(s,spaces):
  # initialize an empty string
  str1 = ""
  # traverse in the string
  for ele in s:
      str1 += ele
      str1 += spaces
  # return string
  return str1

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

def is_nameable(input_array):
  nameable_test = True
  if len(input_array) % 2 == 0:
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

def dash(input_array):
  #['1','sa',"23",'-','24','b','_','25',"-","26"]
  # -- > [['1','sa','23','b','_','25','-','26'],['1','sa','24','b','_','25','-','26']]
  output_array = [input_array]
  Qnounderscore = 0
  for element in range(len(input_array)):
    if input_array[element] == "_":
      Qnounderscore = 1
    if input_array[element] == "-" and Qnounderscore == 0:
      partbefore = input_array[0:element-1]
      partafter = input_array[element+2:]
      output_array_1 = []
      try:
        num1=int(input_array[element-1])
        num2=int(input_array[element+1])+1
        output_array = []
      except:
        break
      for thenumber in range(num1,num2):
        output_array_1.append(partbefore+[str(thenumber)]+partafter)
      for i in output_array_1:
        for j in dash(i): 
          output_array.append(j)
      break
  return output_array

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

def comma(input_array):
  # 1sa,2sa-02 -> ["1sa","2sa-02"]
  output_array = []
  output_string = ""
  bracket_val = 0
  for element in range(len(input_array)):
    if input_array[element] == "(":
      bracket_val += 1
    if input_array[element] == ")":
      bracket_val -= 1
    if input_array[element] == ",":
      if bracket_val == 0:
        output_array.append(output_string)
        output_string = ""
        continue 
    output_string += input_array[element]
  output_array.append(output_string)
  return output_array

def comma2(input_array):
  # (1,3)... -> [[1,...],[3,...]]
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "(":
      elementL = element
    if input_array[element] == ")":
      output_array = []
      elementR = element
      partbefore = input_array[0:elementL]
      partafter = input_array[elementR+1:]
      thenumbers = listToString(input_array[elementL+1:elementR],"").split(",")
      for thenumber in thenumbers:
        for j in comma2(partbefore+[str(thenumber)]+partafter):
          output_array.append(j)
      break
  return output_array

def unique(list1):
  # initialize a null list
  unique_list = []
  # traverse for all elements
  for x in list1:
    # check if exists in unique_list or not
    if x not in unique_list:
      unique_list.append(x)
  # return list
  return unique_list

def zeros(input_array):
  #remnoves zeros from name
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

def my_special_compare_with_asterix(string1,string2_array,howmany,what):
  string1_array = seperate_string_number(string1)
  tothowmany=howmany
  for indx in range(int(len(string1_array)/2)):
    if string1_array[1::2][indx] in string2_array[1::2]:
      found_position = string2_array[1::2].index(string1_array[1::2][indx])
      if string1_array[::2][indx] != string2_array[::2][found_position]:
        return False
      else:
        string2_array=string2_array[:2*found_position]+string2_array[2*found_position+2:]
        continue
    else:
      if what[0] == "whatever":
        tothowmany-=int(string1_array[::2][indx])
        if howmany>0 and tothowmany<0:
          return False
        else:
          continue
      elif string1_array[1::2][indx] in what:
        tothowmany-=int(string1_array[::2][indx])
        continue
      else:
        return False
  if tothowmany == 0 or howmany<0:
    if len(string2_array) > 0:
      return False
    else:
      return True
  else:
    return False

###########################################################
###########################################################
###########################################################

def extract_clusters(clusters_df,Qextract,Pextract,Qclustername,Qout):
  from pandas import DataFrame,concat
  from re import split

  #COMMA SEPARATED 
  # 1sa,2sa-02 -> ["1sa","2sa-02"]
  Pextract_comma = []
  for extract_i in Pextract:
    #EXTRACT IF STRING IS PICKLE FILE:
    if extract_i[-4:] == ".pkl":
      from pandas import read_pickle
      for_extracting = read_pickle(extract_i)
      comma_separated = for_extracting.loc[:,("info","file_basename")].values
    else:
      comma_separated = comma(extract_i)
    for separated_i in range(len(comma_separated)):
      Pextract_comma.append(comma_separated[separated_i])
  
  #DEALING WITH THE EXTRACTION INPUT  
  if Qclustername == 0:
    Pextract_ultimate = unique(Pextract_comma)
  else:
    #CONVERT DASHES TO RANGE
    Pextract_dash = []
    for extract_i in Pextract_comma:
      separated = seperate_string_number(extract_i)
      dash_separated = dash(separated)
      for separated_i in range(len(dash_separated)):
        Pextract_dash.append(dash_separated[separated_i])
    #REMOVE COMMENTS
    Pextract_dash = [ dash_comment(i)[0] for i in Pextract_dash ]
    #CONVERTS MULTI-CHOICE COMMA
    Pextract_comma2 = []
    for extract_i in Pextract_dash:
      comma2_separated = comma2(extract_i)
      for separated_i in range(len(comma2_separated)):
        Pextract_comma2.append(comma2_separated[separated_i])
    #SORTING NAMES
    Pextract_sorted = []
    for extract_i in Pextract_comma2:
      if len(extract_i) > 1:
        array_sorted = sorted([extract_i[i:i + 2] for i in range(0, len(extract_i), 2)],key=lambda x: x[1])
        Pextract_sorted.append([item for sublist in array_sorted for item in sublist])
      else:
        Pextract_sorted.append(extract_i)
        Qclustername = 0
    #ADDING TO UNIQUE ULTIMATE
    Pextract_final = []
    for extract_i in Pextract_sorted:
      corrected = zeros(extract_i)
      if len(corrected) > 0:
        Pextract_final.append(corrected)
    Pextract_ultimate = unique(Pextract_final)

  #JK note: I got once fucked up that this was not sorted so I am sorting it here:
  for index in clusters_df.index:
    file_basename_split = clusters_df.loc[index,("info","cluster_type")]
    split_numbers_letters = split('(\d+)',file_basename_split)[1:]
    cluster_type_array = seperate_string_number(file_basename_split)
    if is_nameable(cluster_type_array):
      cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
      cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
      cluster_type = zeros(cluster_type_array_sorted)
      clusters_df.loc[index,("info","cluster_type")] = cluster_type    

  #SEARCHING FOR CLUSTERS IN THE DATABASE
  newclusters_df = DataFrame()
  for extract_i in Pextract_ultimate: 
    if Qclustername == 0:
      extracted_df = clusters_df[clusters_df["info"]["file_basename"].values == extract_i].copy()
    else:
      if "*" not in extract_i:
        extracted_df = clusters_df[clusters_df["info"]["cluster_type"].values == extract_i].copy()
      else:
        try:
          asterix_position=seperate_string_number(extract_i)[1::2].index('*')
          new_extract_i=seperate_string_number(extract_i)[:2*asterix_position]+seperate_string_number(extract_i)[asterix_position*2+1+1:]
          try:
            howmany=int(extract_i[asterix_position*2])
          except:
            howmany=-1
          what=["whatever"]
        except:
          what=[]
          while 1==1:
            try:
              asterix_position=seperate_string_number(extract_i)[::2].index('*')
              try:
                if len(asterix_position)>1:
                  asterix_position=asterix_position[0]
              except:
                howmany=-1
              new_extract_i=seperate_string_number(extract_i)[:2*asterix_position]+seperate_string_number(extract_i)[asterix_position*2+1+1:]
              howmany=-1
              what.append(seperate_string_number(extract_i)[asterix_position*2+1])
              extract_i=new_extract_i
            except:
              break
        extracted_df = clusters_df[[my_special_compare_with_asterix(value_i,new_extract_i,howmany,what) for value_i in clusters_df["info"]["cluster_type"].values]].copy()
   
    if len(extracted_df) > 0:
      if len(newclusters_df) == 0:
        newclusters_df = extracted_df.copy()
      else:
        newclusters_df = concat([newclusters_df, extracted_df.copy()], axis=0, ignore_index=False)
  
  #RETURNING CLUSTERS_DF
  prelen=len(clusters_df)
  if Qextract == 2:
    from pandas import Index
    a1 = Index(clusters_df.index)
    a2 = Index(newclusters_df.index)
    clusters_df = clusters_df.iloc[a1.difference(a2, sort=False)]
  else:
    clusters_df = newclusters_df.copy()
  if Qout >= 1:
    print("Extracting: "+str(prelen)+" --> "+str(len(clusters_df)))
 
  return clusters_df
