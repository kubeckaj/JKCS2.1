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

def combine_same(input_array):
  output_array = []
  for i in range(0,len(input_array)-1,2):
    test = 0 
    for j in range(0,len(output_array)-1,2):
      if input_array[i+1] == output_array[j+1]:
        test = 1
        output_array[j] = str(int(output_array[j]) + int(input_array[i]))
    if test == 0:
      output_array.append(input_array[i])
      output_array.append(input_array[i+1])
  return output_array

def replace_first_occurrence(string, old_char, new_char):
  index = string.find(old_char)  # Find the index of the first occurrence
  if index != -1:  # If the character is found
   string = string[:index] + new_char + string[index+1:]  # Replace the character
  return string

def rearrange_formula(formula):
    import re
    # Define a regular expression pattern to capture elements and their counts
    pattern = r'([A-Z][a-z]?)(\d*)'
    
    # Use re.findall to extract element symbols and their counts from the formula
    elements = re.findall(pattern, formula)
    
    # Initialize lists to store rearranged parts of the formula
    counts = []
    symbols = []
    
    # Iterate through the extracted elements
    for element, count in elements:
        if count == '':  # If no count is specified, set count to 1
            count = '1'
        counts.append(count)
        symbols.append(element)
    
    # Rearrange the elements and counts to the desired format 'count1symbol1count2symbol2...'
    rearranged_formula = ''.join(counts[i] + symbols[i] for i in range(len(elements)))
    
    return rearranged_formula,counts,symbols

def fitPlaneSVD(XYZ, res):
    import numpy as np
    XYZ = np.array([np.array(i) for i in XYZ])
    res = np.array(res) 
    rows, cols = XYZ.shape
    A = np.hstack([XYZ, 0*np.ones((rows, 1))])
    A_transpose_A_inv = np.linalg.pinv(np.dot(A.T, A))
    refined_coeffs = np.dot(A_transpose_A_inv, np.dot(A.T, res))
    #print(refined_coeffs)
    refined_plane_coeffs = refined_coeffs[:-1]
    
    return refined_plane_coeffs

def mergeDictionary(dict_1, dict_2):
  if len(dict_1) == 0:
    return dict_2
  dict_3 = {**dict_1, **dict_2}
  for key, value in dict_3.items():
    if key in dict_1 and key in dict_2:
      dict_3[key] = value + dict_1[key]
    elif key in dict_1:
      dict_3[key] = [float("nan")] + value
    else:
      ml = max([len(x) for x in dict_1.values()])
      dict_3[key] = value + ml*[float("nan")]
  return dict_3

def data_modification(clusters_df, Qunderscore, Qrename, Qclustername, QrenameWHAT, Qiamifo, Qrebasename, Qdrop, Qout2log,Qpresplit,Qindex, seed,Qatomize):

  #SPLIT THE DATABASE AT START  
  if Qpresplit != 0:
    if Qpresplit > 0:
      clusters_df = clusters_df.sample(n=Qpresplit, random_state=seed)
    else:
      clusters_df = clusters_df.sample(n=len(clusters_df), random_state=seed)
  
  #SELECT INDEX
  if Qindex != "no":
    try:
      clusters_df = eval(f"clusters_df[{Qindex}]")
    except:
      try:
        if Qindex == "-1":
          Qindex = str(len(clusters_df)-1)
        Qindex = Qindex+":"+str(int(Qindex)+1)
        clusters_df = eval(f"clusters_df[{Qindex}]")
      except:
        print("Error selecting given index")
        exit()
    clusters_df = clusters_df.reset_index(drop=True)

  if Qrename == 1:
    if Qclustername == 1:
      from read_files import seperate_string_number,is_nameable
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
          cluster_type_array = combine_same(cluster_type_array)
          l4.append("".join(cluster_type_array)+file_i_BASE.replace(file_i_BASE.split("-")[0].split("_")[0],""))
          cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
          cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
          cluster_type_string = zeros(cluster_type_array_sorted)
          components = cluster_type_array[1::2]
          component_ratio = [int(i) for i in re.split('(\d+)', file_i_BASE.split("-")[0])[1:][0::2]]
          l1.append(cluster_type_string)
          l2.append(components)
          l3.append(component_ratio)
      from pandas import options
      options.mode.chained_assignment = None
      for i in range(len(clusters_df)):
        clusters_df.at[i,("info", "cluster_type")] =  l1[i]
        clusters_df.at[i,("info", "file_basename")] = l4[i]
        clusters_df.at[i,("info", "components")] = array(l2[i], dtype=object)
        clusters_df.at[i,("info", "component_ratio")] = array(l3[i], dtype=object)
      #clusters_df.loc[:,("info", "components")] = array(l2, dtype=object)
      #clusters_df.loc[:,("info", "component_ratio")] = array(l3, dtype=object)

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

  if Qatomize == 1:
    from numpy import array
    overall_formulas = []
    overall_counts = []
    overall_symbols = []
    overall_properties = []
    #clusters_df.loc[:,("info", "component_ratio")].astype(object)
    #clusters_df.loc[:,("info", "components")].astype(object)
    #clusters_df.loc[:,("info", "cluster_type")].astype(object)
    from pandas import option_context
    with option_context('mode.chained_assignment', None):
      if ("info", "component_ratio") in clusters_df:
        clusters_df.loc[:,("info", "component_ratio")] = clusters_df.loc[:,("info", "component_ratio")].astype(object)
      else:
        clusters_df.loc[:,("info", "component_ratio")] = array(len(clusters_df)*[float('nan')]).astype(object)
      if ("info", "components") in clusters_df:
        clusters_df.loc[:,("info", "components")] = clusters_df.loc[:,("info", "components")].astype(object)
      else:
        clusters_df.loc[:,("info", "components")] = array(len(clusters_df)*[float('nan')]).astype(object)
      if ("info", "cluster_type") in clusters_df:
        clusters_df.loc[:,("info", "cluster_type")] = clusters_df.loc[:,("info", "cluster_type")].astype(object)
      else:
        clusters_df.loc[:,("info", "cluster_type")] = array(len(clusters_df)*[float('nan')]).astype(object)
    for cluster_id in clusters_df.index:
      structure = clusters_df.loc[cluster_id,('xyz','structure')]
      output_formula,output_counts,output_symbols = rearrange_formula(structure.get_chemical_formula())
      overall_formulas.append(output_formula)
      overall_counts.append(output_counts)
      overall_symbols.append(output_symbols)
      if ('log','electronic_energy') in clusters_df:
        overall_properties.append(clusters_df.loc[cluster_id,('log','electronic_energy')])
      clusters_df.at[cluster_id,("info", "component_ratio")] = [int(i) for i in output_counts]
      clusters_df.at[cluster_id,("info", "components")] = output_symbols
      clusters_df.at[cluster_id,("info", "cluster_type")] = output_formula
    overall_counts_new = []
    overall_symbols_new = list(set([item for sublist in overall_symbols for item in sublist]))
    for i in range(len(overall_counts)):
      toappend = []
      for j in range(len(overall_symbols_new)):
        if overall_symbols_new[j] in overall_symbols[i]:
          ind = overall_symbols[i].index(overall_symbols_new[j])
          toappend.append(int(overall_counts[i][ind]))
        else:
          toappend.append(int(0))
      overall_counts_new.append(toappend)
    if ('log','electronic_energy') in clusters_df:
      fitted = fitPlaneSVD(overall_counts_new,overall_properties)
    mons_dic = {}
    from numpy import array
    from ase import Atoms
    for i in range(len(overall_symbols_new)):
      columns = ["folder_path","file_basename","cluster_type","components","component_ratio"]
      folder_path = float("nan")
      file_basename = "1"+overall_symbols_new[i]+"-artificial"
      cluster_type = "1"+overall_symbols_new[i]
      components = [overall_symbols_new[i]]
      component_ratio = [int(1)]
      all_locals = locals()
      dic = {("info",column):[all_locals.get(column)] for column in columns}
      if ('log','electronic_energy') in clusters_df:
        dic.update({("log","electronic_energy"):[fitted[i]]})
      dic.update({("xyz","structure"):[Atoms(overall_symbols_new[i], positions = [[0,0,0]])]})
      mons_dic = mergeDictionary(mons_dic,dic)
    from pandas import DataFrame
    mons_df = DataFrame(mons_dic,index=range(len(overall_symbols_new)))
    mons_df.to_pickle("atoms.pkl")
  elif Qatomize == 2:
    for cluster_id in clusters_df.index:
      from read_files import seperate_string_number,is_nameable
      from re import split
      file_basename = clusters_df.loc[cluster_id,("info","file_basename")]
      file_basename_split = file_basename.split("-")[0].split("_")[0] 
      split_numbers_letters = split('(\d+)',file_basename_split)[1:]
      cluster_type_array = seperate_string_number(file_basename_split)
      if is_nameable(cluster_type_array):
        cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
        cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
        cluster_type = zeros(cluster_type_array_sorted)
        components = split_numbers_letters[1::2]
        component_ratio = [int(i) for i in split_numbers_letters[0::2]]
      else:
        cluster_type = float("nan")
        components = float("nan")
        component_ratio = float("nan")
      clusters_df.at[cluster_id,("info","cluster_type")] = cluster_type
      clusters_df.at[cluster_id,("info","components")] = components
      clusters_df.at[cluster_id,("info","component_ratio")] = component_ratio

  if Qdrop != "0":
    if Qdrop in clusters_df:
      clusters_df = clusters_df.drop([Qdrop], axis=1, level=0)
  
  if Qout2log == 1:
    clusters_df = clusters_df.rename({'out': 'log'}, axis='columns')
  if Qout2log == -1:
    clusters_df = clusters_df.rename({'log': 'out'}, axis='columns')
  
  return clusters_df
