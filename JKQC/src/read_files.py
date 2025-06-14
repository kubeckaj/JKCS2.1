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

def read_files(clusters_df, files, orcaextname = "out", orcaext = "out", turbomoleext = "log", Qclustername = 1, Qforces = 0, Qanharm = 0, Qdisp_electronic_energy = 0, Qdisp_forces = 0):
  from os import path
  from re import split
  from pandas import DataFrame
  from mmap import mmap,ACCESS_READ

  Q_ORCA_used = 0
  Q_G16_used = 0
  Q_XTB_used = 0
  Q_ABC_CREST_used = 0
  Q_TURBOMOLE_used = 0
  Q_MRCC_used = 0

  clusters_dict = {}
  rem_orcaextname = orcaextname
  for file_i in files:
    #############################
    ### FILE AND FOLDER NAMES ###
    #############################
    folder_path       = path.abspath(file_i)[::-1].split("/",1)[1][::-1]+"/"
    file_basename     = file_i[:-4][::-1].split("/",1)[0][::-1]
    file_i_ABC_CREST  = folder_path+file_basename+".log"
    file_i_XTB        = folder_path+file_basename+".log"
    file_i_engrad     = folder_path+file_basename+".engrad"
    file_i_G16        = folder_path+file_basename+".log"
    file_i_XYZ        = folder_path+file_basename+".xyz"
    #ORCA && MRCC
    file_i_MRCC       = folder_path+file_basename+".out"
    orcaextname = rem_orcaextname
    file_i_ORCA2 = folder_path+file_basename+"."+"bullshit"
    if not path.exists(folder_path+file_basename+".log"):
      file_i_ORCA = folder_path+file_basename+"."+orcaext
      orcaextname = "log"
      mrccextname = "log"
    else:
      if orcaext == "out":
        if not path.exists(folder_path+file_basename+"."+orcaext):
          file_i_ORCA = folder_path+file_basename+".log"
          orcaextname = "log"
        else:
          file_i_ORCA = folder_path+file_basename+"."+orcaext
      else:
        file_i_ORCA = folder_path+file_basename+"."+orcaext
        file_i_ORCA2 = folder_path+file_basename+"."+"out"
        orcaextname = "log"
        orcaextname2 = "out"
    ##
    file_i_TURBOMOLE = folder_path+file_basename+"."+turbomoleext
    file_i_INFO      = folder_path+"info.txt"

    ###############
    ### INFO ######
    ###############
    columns = ["folder_path","file_basename"]
    if Qclustername == 1:
      columns = columns + ["cluster_type","components","component_ratio"]
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
    all_locals = locals()
    dic = {("info",column):[all_locals.get(column)] for column in columns}
  
    ### EXTRA INFO FILE ###
    if path.exists(file_i_INFO):
      file = open(file_i_INFO, "r")
      for line in file:
        #TODO not sure whether this still works
        splitted_line = line.split(" ",1)
        dic.update({("info",str(splitted_line[0])):[splitted_line[-1].strip()]})
      file.close()    

    ################
    #### XYZ #######
    ################
    if path.exists(file_i_XYZ):
      from read_xyz import read_xyz,identify_1
      out = read_xyz(file_i_XYZ)
      dic.update({("xyz","structure"):[out]})
      out = identify_1(out)
      dic.update({("xyz","id1"):[out]})

    for file_test_ext in list(set(["log","out",turbomoleext,orcaext])):
      file_test = folder_path+file_basename+"."+file_test_ext
      if path.exists(file_test):
        with open(file_test, "r", encoding="utf-8") as f:
          mm = mmap(f.fileno(), 0, access=ACCESS_READ)

          ###############
          #### G16 ######
          ###############
          if file_test == file_i_G16:
            testG16 = mm.find(rb'Gaussian(R)')+1
            if testG16 > 0:
              if Q_G16_used == 0:
                from read_g16 import read_g16,read_g16_init
                read_g16_init(Qforces = Qforces, Qanharm = Qanharm)
                Q_G16_used = 1
              dic_g16 = read_g16(mm, Qforces = Qforces, Qanharm = Qanharm)
              dic.update(dic_g16)
              continue

          ###############
          ### ORCA ######
          ###############
          if file_test == file_i_ORCA:
            testORCA = mm.find(rb'O   R   C   A')+mm.find(rb'ORCA')+mm.find(rb'SHARK')+3
            if testORCA > 0:
              if Q_ORCA_used == 0:
                from read_orca import read_orca,read_orca_init
                read_orca_init(Qforces = Qforces, Qanharm = Qanharm, Qdisp_forces = Qdisp_forces)
                Q_ORCA_used = 1
              dic_orca = read_orca(mm, orcaextname, Qforces = Qforces, Qanharm = Qanharm, Qdisp_electronic_energy = Qdisp_electronic_energy, Qdisp_forces = Qdisp_forces)
              dic.update(dic_orca)
              continue
          ###############
          ### ORCA ######
          ###############
          if file_test == file_i_ORCA2:
            testORCA = mm.find(rb'O   R   C   A')+mm.find(rb'ORCA')+mm.find(rb'SHARK')+3
            if testORCA > 0:
              if Q_ORCA_used == 0:
                from read_orca import read_orca,read_orca_init
                read_orca_init(Qforces = Qforces, Qanharm = Qanharm, Qdisp_forces = Qdisp_forces)
                Q_ORCA_used = 1
              dic_orca = read_orca(mm, orcaextname2, Qforces = Qforces, Qanharm = Qanharm, Qdisp_electronic_energy = Qdisp_electronic_energy, Qdisp_forces = Qdisp_forces)
              dic.update(dic_orca)
              continue

          ###############
          ### XTB ######
          ###############
          if file_test == file_i_XTB:
            testXTB = mm.find(rb'|                           x T B                           |')+1
            if testXTB > 0:
              if Q_XTB_used == 0:
                from read_xtb import read_xtb,read_xtb_init
                read_xtb_init()
                Q_xtb_used = 1
              dic_xtb = read_xtb(mm)
              dic.update(dic_xtb)
              continue
 
          ######################
          #### ABC/CREST #######
          ######################
          if file_test == file_i_ABC_CREST:
            testABC_CREST = mm.find(rb'ABC')+mm.find(rb'JXYZ')+2
            if testABC_CREST > 0:
              if Q_ABC_CREST_used == 0:
                from read_abc_crest import read_abc_crest,read_abc_crest_init
                read_abc_crest_init()
                Q_ABC_CREST_used = 1
              dic_abc_crest = read_abc_crest(mm)
              dic.update(dic_abc_crest)
              continue
  
          ######################
          ####### MRCC #########
          ######################
          if file_test == file_i_MRCC:
            testMRCC = mm.find(rb'MRCC program system')+1
            if testMRCC > 0:
              if Q_MRCC_used == 0:
                from read_mrcc import read_mrcc,read_mrcc_init
                read_mrcc_init()
                Q_MRCC_used = 1
              dic_mrcc = read_mrcc(mm)
              dic.update(dic_mrcc)
              continue         
  
    ###############
    ### ENGRAD ####
    ###############
    ### This part is not needed and is turned off
    if path.exists(file_i_engrad) and Qforces == 1:
      from numpy import array
      if path.exists(file_i_engrad):
        file = open(file_i_engrad, "r")
        for gradi in range(3):
          file.readline()
        out_NAtoms = int(file.readline())
        for gradi in range(7):
          file.readline()
        try:
          save_forces = []
          for gradi in range(3*out_NAtoms):
            save_forces.append(-float(file.readline())/0.529177)
          out_forces = [array([save_forces[i],save_forces[i+1],save_forces[i+2]]) for i in range(0,len(save_forces),3)]
        except:
          out_forces = float("nan")
        dic.update({("extra","forces"):[out_forces]})
        file.close()

    #combine into large dictionary
    clusters_dict = mergeDictionary(clusters_dict,dic)

  newclusters_df = DataFrame(clusters_dict,index=range(len(clusters_df),len(clusters_df)+len(files)))
  if len(clusters_df) > 0:
    from pandas import concat
    clusters_df = concat([clusters_df,newclusters_df.copy()], ignore_index=True)
    #clusters_df = clusters_df.append(newclusters_df)
  else:
    clusters_df = newclusters_df

  return clusters_df

  ################
  #### TURBOMOLE #
  ################
  #if path.exists(file_i_TURBOMOLE):
  #  file = open(file_i_TURBOMOLE, "r")
  #  testTURBOMOLE = 0
  #  for i in range(5):
  #    if re.search("TURBOMOLE", file.readline()):
  #      testTURBOMOLE = 1
  #      break
  #  if testTURBOMOLE == 1:
  #    out = missing
  #    for line in file:
  #      if re.search("Final CCSD\(F12\*\)\(T\) energy", line):
  #        out = float(line.split()[5])
  #      if re.search("Final MP2 energy", line):
  #        out = float(line.split()[5])
  #    file.close()
  #    clusters_df = df_add_iter(clusters_df, turbomoleextname, "electronic_energy", [str(cluster_id)], [out])
