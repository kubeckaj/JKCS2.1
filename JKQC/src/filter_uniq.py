def seperate_string_number2(string):
    ### rg3,el2.4,g -> [rg,3,el,2.4,g]
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i == "_" or previous_character == "_":
            newword += i
        elif i != "-" and i.isalpha() and previous_character.isalpha():
            newword += i
        elif (i.isnumeric() or i == ".") and (previous_character.isnumeric() or previous_character == "." or previous_character == "-"):
            newword += i
        else:
            if previous_character != ",":
              if previous_character.isnumeric() or previous_character == ".":
                groups[-1]=[groups[-1],newword]
              else:
                groups.append(newword)
            if i != ",":
              newword = i
        previous_character = i
        if x == len(string) - 2:
            if previous_character.isnumeric() or previous_character == ".":
              groups[-1]=[groups[-1],newword]
            else:
              groups.append(newword)
            newword = ''
    return groups

def filter_uniq(clusters_df,Quniq,Qclustername,Qsample,Qout):
  from numpy import unique, transpose, sum, tile, floor
  missing = float("nan")

  if Quniq == "dup":
    newclusters_df = clusters_df.copy()
    newclusters_df = newclusters_df.drop_duplicates(subset=[("info","file_basename")])
  else:
    if Qclustername != 0:
      uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")].values)
    else:
      uniqueclusters = "1"
    newclusters_df = []

    myNaN = lambda x : missing if x == "NaN" else x
    for i in uniqueclusters:
       if Qclustername != 0:
         preselected_df = clusters_df[clusters_df.loc[:,("info","cluster_type")] == i]
       else:
         preselected_df = clusters_df
       separated_inputs = seperate_string_number2(str(Quniq))
       compare_list = []
       compare_list_num = []
       for separated_input in separated_inputs:
         if isinstance(separated_input,list):
           if separated_input[0] == "rg":
             import warnings
             warnings.filterwarnings("ignore", category=RuntimeWarning)
             compare_list.append("rg")
           elif separated_input[0] == "el":
             compare_list.append("electronic_energy")
           elif separated_input[0] == "g":
             compare_list.append("gibbs_free_energy")
           elif separated_input[0] == "d" or separated_input[0] == "dip":
             compare_list.append("dipole_moment")
           else:
             compare_list.append(separated_input[0])
           compare_list_num.append(float(separated_input[1]))
         else:
           #compare_list.append(separated_input)
           if separated_input == "rg":
             import warnings
             warnings.filterwarnings("ignore", category=RuntimeWarning)
             compare_list.append("rg")
             compare_list_num.append(2)
           elif separated_input == "mass":
             compare_list.append("mass")
             compare_list_num.append(2)
           elif separated_input == "el":
             compare_list.append("electronic_energy")
             compare_list_num.append(3)
           elif separated_input == "g":
             compare_list.append("gibbs_free_energy")
             compare_list_num.append(3)
           elif separated_input == "d" or separated_input == "dip":
             compare_list.append("dipole_moment")
             compare_list_num.append(1)
           else:
             compare_list.append(separated_input)
             compare_list_num.append(1)
       scale_test=1
       scale=0
       scale_scale = 0.1
       scale_iter = 0
       while scale_test == 1:
         scale_iter = scale_iter + 1
         tocompare = []
         for compare_element_i in range(len(compare_list)):
           j = compare_list[compare_element_i]
           jj = scale+compare_list_num[compare_element_i]
           if j == "rg":
             rg = []
             for aseCL in preselected_df["xyz"]["structure"]:
               try:
                 rg.append((sum(sum((aseCL.positions-tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/sum(aseCL.get_masses()))**0.5)
               except:
                 rg.append(missing)
             values = [floor(myNaN(o)*10**jj) for o in rg]
           elif j == "mass":
             mass = []
             for aseCL in preselected_df["xyz"]["structure"]:
               try:
                 mass.append(sum(aseCL.get_masses()))
               except:
                 mass.append(missing)
             values = [floor(myNaN(o)*10**jj) for o in mass]
           elif j == "gout":
             gout = []
             for Gouti in range(len(preselected_df)):
               try:
                 gout.append(preselected_df.loc[Gouti,("log","gibbs_free_energy")].values-preselected_df.loc[Gouti,("log","electronic_energy")].values+preselected_df.loc[Gouti,("out","electronic_energy")].values)
               except:
                 gout.append(missing)
             values = [floor(myNaN(o)*10**jj) for o in gout]
           else:
             values = [floor(float(myNaN(o))*10**jj) for o in preselected_df.loc[:,("log",j)].values]
           tocompare.append(values)
         tocompare = transpose(tocompare)
         uniqueindexes = unique(tocompare,axis = 0,return_index=True)[1]
         if Qsample > 0 and len(uniqueindexes) > Qsample and scale_iter < 100:
           if scale_scale > 0:
             scale_scale = -0.9*scale_scale
           else:
             scale_scale = 1.1*scale_scale
           scale = scale - abs(scale_scale)
         elif Qsample > 0 and len(uniqueindexes) < Qsample and len(preselected_df) >= Qsample and scale_iter < 100:
           if scale_scale < 0:
             scale_scale = -0.9*scale_scale
           else:
             scale_scale = 1.1*scale_scale       
           scale = scale + abs(scale_scale)
         else:
           scale_test = 0
       selected_df = preselected_df.iloc[uniqueindexes]
       if len(newclusters_df) == 0:
         newclusters_df = selected_df
       else:
         from pandas import concat
         newclusters_df = concat([newclusters_df,selected_df.copy()], ignore_index=True)
         #newclusters_df = newclusters_df.append(selected_df)
  if Qout >= 1:
    if Qsample > 0:
      print("Sampled: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
    else:
      print("Uniqueness: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
  return newclusters_df
