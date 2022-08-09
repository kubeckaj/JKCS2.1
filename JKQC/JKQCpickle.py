####################################################################################################
# Reading input

import os
import sys
import numpy as np

def listToString(s,spaces): 
    # initialize an empty string
    str1 = ""
    # traverse in the string  
    for ele in s: 
        str1 += ele  
        str1 += spaces 
    # return string  
    return str1

# Checking if correct environmnent is loaded
if str(os.environ.get("VIRTUAL_ENV").split("/")[-1]) != "JKCS":
  print("Trying to load JKCS environment by myself since you were lazy ass and did not do it!.")
  from subprocess import call
  if not call(['/bin/bash', '-i', '-c', "JKpython; python "+listToString(sys.argv," ")]):
    print("This time, JKCS environment was loaded. Please, load it next time by yourself (JKpython).")
  else:
    print("You did not load JKCS environment and I was not able to do it for you.")
  exit()

def print_help():
  print("##################################################")
  print("### python JKQC.py [FILES] [ARGUMENTS] [PRINT] ###")
  print("\nFILES:")
  print(" any .xyz .log .out files")
  print(" any .pkl (you can also use -in/-out)")
  print("\nARGUMENTS:")
  print(" -help      print this help")
  print(" -in X.pkl  read data from a database")
  print(" -out X.pkl save data to XX.pkl (-noex = do not print example)")
  print(" -folder X  takes in all X/*.log files")
  print(" -noname    the file names are not analysed (e.g. 1000-10_1.xyz)")
  print(" -extract X prints only selected clusters (e.g. 1sa1w,1sa3-4w or 1sa1w-1_0)")
  print(" -except X  prints only non-selected clusters (e.g. 1sa1w,1sa3-4w or 1sa1w-1_0)")
  print(" -reacted   removes (the minority of) reacted structures")
  print(" -noexample,-noex does not print an example")
  print("\nPRINT:")
  print(" -ct                cluster type (e.g. 1sa,3sa2am)")
  print(" -b                 basename (e.g. 1sa-3, 3sa2am-4_508)")
  print(" -nXYZ,-nLOG,-nOUT  name with its extension")
  print(" -pXYZ,-pLOG,-pOUT  full_path/name with its extension")
  print(" -ePKL              pkl_full_filepath/:EXTRACT:/basename")
  print(" -el       electronic energy from .log [Eh]   -elSP        single-point el.en. (1. in .log) [Eh]")
  print(" -elout    el. energy from .out [Eh]          -rot         RMSD of rotational constants [GHz]") 
  print(" -elc      elc=elout-el [Eh]                  -rots        rotational constants [GHz]")
  print(" -zpec     ZPE correction [Eh]                -mult        muliplicity [-]")
  print(" -zpe      ZPE=el+zpec [Eh]                   -char,-chrg  charge [-electron_charge]")
  print(" -zpeout   ZPEout=elout+zpec [Eh]             -dip         dipole moment [Debye]")
  print(" -uc       energy thermal correction [Eh]     -dips        dipole moments [Debye]")
  print(" -u        internal energy U=el+uc [Eh]       -pol         polarizability [Angstrom^3]")
  print(" -uout     Uout=elout+uc [Eh]                 -templog     temperature used in .log [K]")
  print(" -hc       enthalpy th. corr. hc=uc+kT [Eh]   -preslog     pressure used in .log [atm]")
  print(" -h        enthalpy energy H=el+hc [Eh]       -lf          the lowest vib. freq. [1/cm]")
  print(" -hout     Hout=elout+hc [Eh]                 -f           array of vibration freq. [1/cm]")
  print(" -s        entropy [Eh/K]                     -rsn         rotational symmetry number [-]")
  print(" -gc       Gibbs free energy th. corr. [Eh]   -t,-time     total computational (elapsed) time [mins]")
  print(" -g        Gibbs free energy [Eh]             -termination count lines with \"Normal termination\" status")
  print(" -gout     G with el.en.corr.:Gout=G+elc [Eh] -mi          moments of inertia")
  print(" -mull     Mulliken charges [-el.charge]      -ami         average moment of inertia")
  print(" -xyz      save all XYZ files                 -rg          radius of gyration [Angstrom]")
  print(" -movie    save all XYZs to movie.xyz         -radius      approx. radius of cluster size [Angstrom]")
  print("                                              -radius0.5   radius with +0.5 Angstrom correction")
  print("\nPOST-CALCULATIONS:")
  print(" -fc [value in cm^-1] frequency cut-off for low-vibrational frequencies")
  print(" -temp [value in K]   recalculate for different temperature")
  print(" -v,-as [value]       anharmonicity scaling factor")
  print(" -unit                converts units [Eh] -> [kcal/mol] (for entropy: [Eh/K] -> [cal/mol/K])")
  print("\nFILTERING:")
  print(" -sort <str>          sort by: g,gout,el")
  print(" -select <int>        selects <int> best structures from each cluster")
  print(" -uniq,-unique <str>  selects only unique based on, e.g.: rg,el or rg,g or rg,el,dip")
  print("\nFORMATION PROPERTIES:")
  print(" -glob OR -globout       prints only values for clusters with the lowest -g OR -gout")
  print(" -bavg OR -bavgout       prints a value that is Boltzmann average over each cluster using -g OR -gout")
  print("                         NOTE: -g/-gout is treated correctly + -s not treated; use (G - H)/T")
  print(" -formation              print values as formation ")
  print(" <input_file> -formation print formations for the input file (no averaging though)")
  print(" -conc sa 0.00001        dG at given conc. (use -cnt for self-consistent dG)")

folder = "./"	
files = []  
input_pkl = []
output_pkl = "mydatabase.pkl"

Qclustername = 1 #Analyse file names for cluster definition?
Qextract = 0 #Do I want to extarct only some cluster_type(s)?
Pextract = []
Qreacted = 0 #Remove reacted structures

Qout = 0 #Do I want to save output.pkl? 0=NO,1=YES,2=YES but do not print example
Pout = []

missing = np.nan

QUenergy = 1 #if you want to energy 
QUentropy = 1 #if you want to energy 

Qqha = 0 #Run the QHA
Qt = missing
Qp = missing
Qfc = 0 #Run QHA with vib. frequency cutoff
Qanh = 1

Qsort = 0 # 0=no sorting, otherwise string
Qselect = 0 # 0=nothing, otherwise the number of selected structures
Quniq = 0 # uniqie based on given arguments
formation_input_file = ""

Qglob = 0 # 1=values for lowest -g, 2=values for lowest -gout
Qbavg = 0 # 1=Boltzmann avg over -g, 2=Boltzmann avg over -gout
Qformation = 0
Qconc = 0
conc = []
CNTfactor = 0 #see Wyslouzil

orcaext = "out"
orcaextname = "out"
turbomoleext = "out"
turbomoleextname = "out"

last = ""
for i in sys.argv[1:]:
  #HELP
  if i == "-help" or i == "--help":
    print_help()
    exit()
  #FOLDER
  if i == "-folder":
    last = "-folder"
    continue
  if last == "-folder":
    last = ""
    if os.path.exists(i):
      folder = i
      continue
    else:
      print("Folder "+i+" does not exist. [EXITING]")
      exit()
  #INPKL
  if i == "-in" or i == "-i":
    last = "-in"
    continue
  if last == "-in":
    last = ""
    if os.path.exists(i):    
      input_pkl.append(i)
      continue
    else:
      print("File "+i+" does not exist. Sorry [EXITING]")
      exit()
  #OUTPKL
  if i == "-out" or i == "-o":
    last = "-out"
    continue
  if last == "-out":
    last = ""
    if output_pkl != "mydatabase.pkl":
      print("Hey are you trying to use two different outputs? "+output_pkl+"/"+i+" [EXITING]")
      exit()
    else:
      output_pkl = i
      if Qout == 0:
        Qout = 1
      continue
  #ORCA EXTENSION
  if i == "-orcaext" or i == "-orca":
    last = "-orcaext"
    continue
  if last == "-orcaext":
    last = ""
    orcaext = str(i)
    orcaextname = "log"
    continue
  #TURBOMOLE EXTENSION
  if i == "-turbomoleext" or i == "-turbomole":
    last = "-turbomoleext"
    continue
  if last == "-turbomoleext":
    last = ""
    turbomoleext = str(i)
    turbomoleextname = "log"
    continue
  #NONAME
  if i == "-noname":
    Qclustername = 0
    continue
  #REACTED
  if i == "-reacted":
    Qreacted = 1
    continue
  #NOEXAMPLE
  if i == "-noexample" or i == "-noex":
    Qout = 2
    continue
  #EXTRACT
  if i == "-extract":
    last = "-extract"
    continue
  if last == "-extract":
    last = ""
    Qextract = 1
    Pextract.append(i)
    continue
  #EXCEPT
  if i == "-except":
    last = "-except"
    continue
  if last == "-except":
    last = ""
    Qextract = 2
    Pextract.append(i)
    continue
  #INPUT FILES
  if len(i) > 3:
    ext = i[-4:]
    if ext == ".xyz" or ext == ".log" or ext == ".out":
      if os.path.exists(i):
        files.append(i)
        continue
      else:
        print("File "+i+" does not exist. Sorry [EXITING]")
        exit()
    if ext == ".pkl":
      if os.path.exists(i):
        input_pkl.append(i)
        continue
      else:
        if output_pkl != "mydatabase.pkl":
          print("Hey either one input database does not exist or you are trying to use two different outputs? "+output_pkl+"/"+i+" [EXITING]")
          exit()
        else:
          output_pkl = i
          continue
  ########
  # SORT 
  if i == "-sort" or i == "--sort":
    last = "-sort"
    continue
  if last == "-sort":
    last = ""
    Qsort = str(i)
    continue
  # SELECT
  if i == "-select" or i == "--select":
    last = "-select"
    continue
  if last == "-select":
    last = ""
    Qselect = int(i)
    continue  
  # UNIQUE
  if i == "-unique" or i == "--unique" or i == "-uniq" or i == "--uniq":
    last = "-uniq"
    continue
  if last == "-uniq":
    last = ""
    Quniq = str(i)
    continue
  ########
  # INFO
  if i == "-info" or i == "--info":
    Pout.append("-info")
    continue
  # CITE
  if i == "-cite" or i == "--cite":
    Pout.append("-cite")
    continue
  # XYZ
  if i == "-xyz" or i == "--xyz" or i == "-XYZ" or i == "--XYZ":
    Pout.append("-xyz")
    continue
  # MOVIE
  if i == "-movie" or i == "--movie":
    Pout.append("-movie")
    continue
  # RG
  if i == "-rg" or i == "--rg" or i == "-Rg" or i == "--Rg":
    Pout.append("-rg")
    continue
  if i == "-radius" or i == "--radius":
    Pout.append("-radius")
    continue
  if i == "-radius0.5" or i == "--radius0.5":
    Pout.append("-radius0.5")
    continue
  # INFO
  if i == "-ct" or i == "-clustertype" or i == "--ct" or i == "--clustertype":
    Pout.append("-ct")
    continue
  if i == "-b" or i == "-basename" or i == "--b" or i == "--basename":
    Pout.append("-b")
    continue
  if i == "-nOUT" or i == "--nOUT":
    Pout.append("-nOUT")
    continue
  if i == "-nLOG" or i == "--nLOG":
    Pout.append("-nLOG")
    continue
  if i == "-nXYZ" or i == "--nXYZ":
    Pout.append("-nXYZ")
    continue
  if i == "-pOUT" or i == "--pOUT":
    Pout.append("-pOUT")
    continue
  if i == "-pLOG" or i == "--pLOG":
    Pout.append("-pLOG")
    continue
  if i == "-pXYZ" or i == "--pXYZ":
    Pout.append("-pXYZ")
    continue
  if i == "-ePKL" or i == "--ePKL":
    Pout.append("-ePKL")
    continue
  # LOG & OUT
  if i == "-elSP" or i == "-elsp" or i == "--elSP" or i == "--elsp":
    Pout.append("-elsp")
    continue
  if i == "-el" or i == "-elen" or i == "--el" or i == "--elen":
    Pout.append("-el")
    continue
  if i == "-elout" or i == "--elout":
    Pout.append("-elout")
    continue
  if i == "-elc" or i == "--elc":
    Pout.append("-elc")
    continue
  if i == "-g" or i == "-gibbs" or i == "--g" or i == "--gibbs":
    Pout.append("-g")
    continue
  if i == "-h" or i == "-enthalpy" or i == "--h" or i == "--enthalpy":
    Pout.append("-h")
    continue
  if i == "-s" or i == "-entropy" or i == "--s" or i == "--entropy":
    Pout.append("-s")
    continue
  if i == "-u" or i == "--u":
    Pout.append("-u")
    continue
  if i == "-uc" or i == "--uc":
    Pout.append("-uc")
    continue
  if i == "-uout" or i == "--uout":
    Pout.append("-uout")
    continue
  if i == "-hc" or i == "--hc":
    Pout.append("-hc")
    continue
  if i == "-hout" or i == "--hout":
    Pout.append("-hout")
    continue
  if i == "-gc" or i == "--gc":
    Pout.append("-gc")
    continue
  if i == "-g" or i == "--g" or i == "-gibbs" or i == "--gibbs" or i == "-G" or i == "--G":
    Pout.append("-g")
    continue
  if i == "-gout" or i == "--gout" or i == "-GDh" or i == "--GDh" or i == "-GD" or i == "--GD":
    Pout.append("-gout")
    continue
  if i == "-zpec" or i == "-zpecorr" or i == "--zpec" or i == "--zpecorr":
    Pout.append("-zpec")
    continue
  if i == "-zpe" or i == "-ZPE" or i == "--zpe" or i == "--ZPE":
    Pout.append("-zpe")
    continue
  if i == "-zpeout" or i == "-ZPEout" or i == "--zpeout" or i == "--ZPEout":
    Pout.append("-zpeout")
    continue
  if i == "-lf" or i == "-lfreq" or i == "--lf" or i == "--lfreq":
    Pout.append("-lf")
    continue
  if i == "-f" or i == "-freq" or i == "--f" or i == "--freq":
    Pout.append("-f")
    continue
  if i == "-rot" or i == "--rot":
    Pout.append("-rot")
    continue
  if i == "-rots" or i == "--rots":
    Pout.append("-rots")
    continue
  if i == "-mult" or i == "--mult":
    Pout.append("-mult")
    continue
  if i == "-char" or i == "-chrg" or i == "--char" or i == "--chrg":
    Pout.append("-char")
    continue
  if i == "-mull" or i == "--mull":
    Pout.append("-mull")
    continue
  if i == "-dip" or i == "--dip":
    Pout.append("-dip")
    continue
  if i == "-dips" or i == "--dips":
    Pout.append("-dips")
    continue
  if i == "-pol" or i == "--pol":
    Pout.append("-pol")
    continue
  if i == "-templog" or i == "--templog":
    Pout.append("-templog")
    continue
  if i == "-preslog" or i == "--preslog":
    Pout.append("-preslog")
    continue
  if i == "-mi" or i == "--mi":
    Pout.append("-mi")
    continue
  if i == "-ami" or i == "--ami":
    Pout.append("-ami")
    continue
  if i == "-rsn" or i == "--rsn":
    Pout.append("-rsn")
    continue
  if i == "-t" or i == "--t" or i == "-time" or i == "--time":
    Pout.append("-t")
    continue
  if i == "-termination" or i == "--termination":
    Pout.append("-termination")
    continue
  #PRE_EXTRACT_DATA MANIPULATION
  if i == "-fc" or i == "--fc":
    Qqha = 1
    last = "-fc"
    continue
  if last == "-fc":
    last = ""
    Qfc = float(i)
    continue
  if i == "-temp" or i == "--temp":
    Qqha = 1
    last = "-temp"
    continue
  if last == "-temp":
    last = ""
    Qt = float(i)
    continue
  if i == "-as" or i == "--as" or i == "-v" or i == "--v":
    Qqha = 1
    last = "-as"
    continue
  if last == "-as":
    last = ""
    Qanh = float(i)
    continue
  if i == "-unit" or i == "--unit":
    QUenergy = 627.503
    QUentropy = 627503.0
    continue
  #FORMATION
  if i == "-glob" or i == "--glob":
    Qglob = 1
    continue
  if i == "-globout" or i == "--globout":
    Qglob = 2
    continue
  if i == "-bavg" or i == "--bavg":
    Qbavg = 1
    continue
  if i == "-bavgout" or i == "--bavgout":
    Qbavg = 2
    continue
  if i == "-formation" or i == "--formation":
    Qformation = 1
    continue
  if i == "-conc" or i == "--conc":
    Qconc = 1
    last = "-conc"
    continue
  if last == "-conc":
    last = "-conc2"
    remember = str(i)
    continue
  if last == "-conc2":
    last = ""
    conc.append(np.array([remember, float(i)]))
    continue
  if os.path.exists(i):
    formation_input_file = i
    continue
  if i == "-cnt" or i == "--cnt":
    CNTfactor = 1
    continue
  #UNKNOWN ARGUMENT   
  print("I am sorry but I do not understand the argument: "+i+" [EXITING]")
  exit()  

####################################################################################################
# Checking for arguments

if len(files) == 0:
  if len(input_pkl) == 0 and len(formation_input_file) == 0:
    import glob
    files = glob.glob(folder+"/*.log")
    if len(files) == 0:
      print("No inputs. No *.log files. [EXITING]")
      exit()
  if len(input_pkl) > 1:
    if Qout == 0:
      Qout = 1
else:
  if len(Pout) == 0:
    if Qout==0:
      Qout = 1
    print("Number of files: "+str(len(files)))
  else:
    Qout = 2

import pandas as pd
import re
from ase.io import read,write


####################################################################################################
# THIS IS TAKEN FROM tool.py IN JKCS.py 
### Append to dataframe 
def df_add_append(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    newdataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
    dataframe = dataframe.append(newdataframe)
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )

  return dataframe

### Add to dataframe via iteration
def df_add_iter(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    df_l = dataframe.shape[0]
    var_l = len(variables)
    for i in range(var_l):
      dataframe[label,name][df_l-var_l+i] = variables[i]
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
    for i in range(len(variables)):
      dataframe[label,name][i] = variables[i]

  return dataframe
####################################################################################################

# Loading input pikles
if len(input_pkl) == 0:
  clusters_df = pd.DataFrame()
else:
  for i in range(len(input_pkl)):
    newclusters_df = pd.read_pickle(input_pkl[i])
    newclusters_df.index = [str(j) for j in range(len(newclusters_df))]
    if Qout == 1:
      print("Number of files in "+input_pkl[i]+": "+str(len(newclusters_df)))
    if i == 0: 
      clusters_df = newclusters_df
    else:
      len_clusters_df = len(clusters_df)
      newclusters_df.index = [str(j+len_clusters_df) for j in range(len(newclusters_df))]
      clusters_df = clusters_df.append(newclusters_df)

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

# Reading the new files in
for file_i in files:
  folder = os.path.abspath(file_i)[::-1].split("/",1)[1][::-1]+"/"
  file_i_BASE = file_i[:-4][::-1].split("/",1)[0][::-1]
  file_i_ABC  = folder+file_i_BASE+".log"
  file_i_XTB  = folder+file_i_BASE+".log"
  file_i_G16  = folder+file_i_BASE+".log"
  file_i_XYZ  = folder+file_i_BASE+".xyz"
  file_i_ORCA = folder+file_i_BASE+"."+orcaext
  file_i_TURBOMOLE = folder+file_i_BASE+"."+turbomoleext
  file_i_INFO = folder+"info.txt" 
 
  ###############
  ### INFO ######
  ###############
  cluster_id = len(clusters_df)
  clusters_df = df_add_append(clusters_df, "info", "folder_path", [str(cluster_id)], folder)
  clusters_df = df_add_iter(clusters_df, "info", "file_basename", [str(cluster_id)], [file_i_BASE])
  # sorted cluster type
  if Qclustername == 1:
    cluster_type_array = seperate_string_number(file_i_BASE.split("-")[0]) 
    if np.mod(len(cluster_type_array),2) == 0:
      cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
      cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
      cluster_type_string = zeros(cluster_type_array_sorted)
      clusters_df = df_add_iter(clusters_df, "info", "cluster_type", [str(cluster_id)], [cluster_type_string])
      #
      components = re.split('(\d+)', file_i_BASE.split("-")[0])[1:][1::2]
      clusters_df = df_add_iter(clusters_df, "info", "components", [str(cluster_id)], [components])
      component_ratio = [int(i) for i in re.split('(\d+)', file_i_BASE.split("-")[0])[1:][0::2]]
      clusters_df = df_add_iter(clusters_df, "info", "component_ratio", [str(cluster_id)], [component_ratio])
  if os.path.exists(file_i_INFO):
    file = open(file_i_INFO, "r")  
    for line in file:
      clusters_df = df_add_iter(clusters_df, "info", str(line.split(" ",1)[0]), [str(cluster_id)], [line.split(" ",1)[-1].strip()])    
    file.close()

  ###############
  ### ABC #######
  ###############
  if os.path.exists(file_i_ABC):
    file = open(file_i_ABC, "r")
    testABC = 0
    for i in range(1):
      if re.search("ABC", file.readline()):
        testABC = 1
        break
    if testABC == 1:
      out_electronic_energy = missing  #1
      for line in file:
        if re.search("ABC energy:", line): #1
          try:
            out_electronic_energy = float(line.split()[2]) 
          except:
            out_electronic_energy = missing
      clusters_df = df_add_iter(clusters_df, "log", "electronic_energy", [str(cluster_id)], [out_electronic_energy])
    file.close()
      

  ###############
  ### XTB #######
  ###############
  if os.path.exists(file_i_XTB):
    file = open(file_i_XTB, "r")
    testXTB = 0
    for i in range(5):
      if re.search("\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-", file.readline()):
        testXTB = 1
        break
    if testXTB == 1:
      out_dipole_moment = missing      #0
      out_electronic_energy = missing  #1
      out_enthalpy_energy = missing                      #T8
      out_gibbs_free_energy = missing  #2
      out_entropy = missing                              #E1
      for line in file:
        if re.search("Debye", line): #0
          try:
            out_dipole_moment = float(line.split()[5])
          except:
            out_dipole_moment = missing
          continue
        if re.search("TOTAL ENERGY", line): #1
          try:
            out_electronic_energy = float(line.split()[3])
          except:
            out_electronic_energy = missing
          continue
        if re.search("total E", line): #1
          try:
            out_electronic_energy = float(line.split()[3])
          except:
            out_electronic_energy = missing
          continue
        if re.search("^H\(T\)", line): #T8
          try:
            out_enthalpy_energy = out_electronic_energy + float(line.split()[1])
          except:
            out_enthalpy_energy = missing
        if re.search("^T\*S ", line): #E1
          try:
            out_entropy = float(line.split()[1])*1000*627.503/298.15
          except:
            out_entropy = missing
        #if re.search("^G\(T\)", line): 
        #  try:
        #    out_gibbs_free_energy = out_electronic_energy + float(line.split()[1])
        #  except:
        #    out_gibbs_free_energy = missing
        if re.search("TOTAL FREE ENERGY", line): #2
          try:
            out_gibbs_free_energy = float(line.split()[0])
          except:
            out_gibbs_free_energy = missing
          continue
      clusters_df = df_add_iter(clusters_df, "log", "dipole_moment", [str(cluster_id)], [out_dipole_moment])
      clusters_df = df_add_iter(clusters_df, "log", "electronic_energy", [str(cluster_id)], [out_electronic_energy])
      clusters_df = df_add_iter(clusters_df, "log", "gibbs_free_energy", [str(cluster_id)], [out_gibbs_free_energy])
      clusters_df = df_add_iter(clusters_df, "log", "enthalpy_energy", [str(cluster_id)], [out_enthalpy_energy])
      clusters_df = df_add_iter(clusters_df, "log", "entropy", [str(cluster_id)], [out_entropy])
    file.close()

  ###############
  ### G16 #######
  ###############
  if os.path.exists(file_i_G16):
    file = open(file_i_G16, "r")
    testG16 = 0
    for i in range(5):
      if re.search("Gaussian", file.readline()):
        testG16 = 1
        break
    if testG16 == 1:
      out_time = missing                                 #TIME
      out_termination = 0                                #TERMINATION
      out_charge = missing                               #I1
      out_multiplicity = missing                         #I2
      out_NAtoms = missing                               #I3
      out_rotational_constants = missing                 #O1
      out_rotational_constant = missing                  #O2
      out_sp_electronic_energy = missing                 #O3x
      out_electronic_energy = missing                    #O3
      out_mulliken_charges = missing                     #O4
      out_dipole_moment = missing                        #O5
      out_dipole_moments = missing                       #O6
      out_polarizability = missing                       #O7
      out_vibrational_frequencies = missing              #V1
      out_temperature = missing                          #V2
      out_pressure = missing                             #V3
      out_moments_of_inertia = missing                   #V4
      out_rotational_symmetry_number = missing           #T1
      out_zero_point_correction = missing                #T2
      out_energy_thermal_correction = missing            #T3
      out_enthalpy_thermal_correction = missing          #T4
      out_gibbs_free_energy_thermal_correction = missing #T5
      out_zero_point_energy = missing                    #T6
      out_internal_energy = missing                      #T7
      out_enthalpy_energy = missing                      #T8
      out_gibbs_free_energy = missing                    #T9
      out_entropy = missing                              #E1
      search=-1
      save_mulliken_charges=0
      save_something=""
      for line in file:
        #TIME
        if re.search("Elapsed time",line): 
          if np.isnan(out_time):
            out_time = 0
          try:
            out_time += float(line.split()[2])*24*60+float(line.split()[4])*60+float(line.split()[6])+float(line.split()[8])/60 
          except:
            out_time = missing
          continue
        #TERMINATION
        if re.search("Normal termination",line):
          out_termination += 1
          continue
        #INITIAL SEARCH
        if search==-1:
          if re.search(" Charge = ", line): #I1/I2
            try:
              out_charge = int(line.split()[2])
              out_multiplicity = int(line.split()[5])
            except:
              out_charge = missing
              out_multiplicity = missing
            search += 1
            continue
        #OPTIMIZATION SEARCH
        ## I will convert the strings into floats/ints later !!!!!!!
        if search==0 or search==-1: 
          if np.isnan(out_NAtoms): 
            if re.search(" NAtoms=", line): #I3
              try:
                out_NAtoms = int(line.split()[1])
              except:
                out_NAtoms = missing
              continue
          if re.search("Rotational constants ", line): #O1/O2
            try:
              out_rotational_constants = [str(line.split()[3]),str(line.split()[4]),str(line.split()[5])]
            except:
              out_rotational_constants = missing
            continue
          if re.search("SCF Done", line): #O3
            try:
              out_electronic_energy = float(line.split()[4])
            except:
              out_electronic_energy = missing
            if np.isnan(out_sp_electronic_energy):
              out_sp_electronic_energy = out_electronic_energy
            continue
          ## MULLIKEN
          if re.search(" Mulliken charges:", line): #O4
            save_mulliken_charges = out_NAtoms+1
            save_something = "mulliken_charges"
            try:
              out_mulliken_charges = ["0"]*out_NAtoms
            except:
              out_mulliken_charges = missing
            continue
          if save_something == "mulliken_charges":
            try:
              if save_mulliken_charges<=out_NAtoms:
                out_mulliken_charges[out_NAtoms-save_mulliken_charges] = str(line.split()[2])
            except:
              out_mulliken_charges = missing
            save_mulliken_charges-=1
            if save_mulliken_charges == 0:
              save_something = ""
            continue
          ## DIPOLE
          if re.search('Dipole moment \(field-independent basis, Debye\):', line): #O5/O6
            save_something = "dipole_moment"
            continue
          if save_something == "dipole_moment":
            try:
              out_dipole_moment = str(line.split()[7])
            except:
              out_dipole_moment = missing
            try:
              out_dipole_moments = [str(line.split()[1]), str(line.split()[3]), str(line.split()[5])]
            except:
              out_dipole_moments = missing
            save_something = ""
            continue
          # POLARIZABILITY 
          if re.search("Exact polarizability", line): #O7
            try:
              pol = [float(i) for i in line.split()[2:]]
              pol_mat = np.matrix([[pol[0],pol[1],pol[3]],[pol[1],pol[2],pol[4]],[pol[3],pol[4],pol[5]]])
              out_polarizability = 0.14818471147*sum(np.linalg.eigh(pol_mat)[0])/3.
            except:
              out_polarizability = missing
            search+=1
            #+unit conversion
            try:
              out_rotational_constants = [float(i) for i in out_rotational_constants]
            except:
              out_rotational_constants = missing
            out_rotational_constant = np.linalg.norm(out_rotational_constants)
            out_electronic_energy = float(out_electronic_energy)
            try:
              out_mulliken_charges = [float(i) for i in out_mulliken_charges]
            except:
              out_mulliken_charges = missing
            out_dipole_moment = float(out_dipole_moment)
            try:
              out_dipole_moments = [float(i) for i in out_dipole_moments]
            except:
              out_dipole_moments = missing
            continue       
        #VIBRATIONAL FREQUENCIES
        if search==1:
          if re.search("Frequencies", line): #V1
            if np.all(np.isnan(out_vibrational_frequencies)):
              out_vibrational_frequencies = []
            try:
              out_vibrational_frequencies += [float(i) for i in line.split()[2:]]
            except:
              out_vibrational_frequencies = missing
            continue
          if re.search("Kelvin.  Pressure", line): #V2/V3
            try:
              out_temperature = float(line.split()[1])
              out_pressure = float(line.split()[4])
            except:
              out_temperature = missing
              out_pressure = missing
            continue
          if re.search("Eigenvalues -- ", line): #V4
            if np.isnan(out_moments_of_inertia):
              out_moments_of_inertia = []
            try:
              out_moments_of_inertia += [float(i) for i in line.split()[2:]]
            except:
              out_moments_of_inertia = missing
            search+=1
            continue
          if re.search("Zero-point correction=", line):
            search+=1
        #THERMOCHEMISTRY
        if search==2:
          if re.search("Rotational symmetry number", line): #T1
            try:
              out_rotational_symmetry_number = float(line.split()[3])
            except:
              out_rotational_symmetry_number = missing
            continue
          if re.search("Zero-point correction=", line): #T2
            try:
              out_zero_point_correction = float(line.split()[2])
            except:
              out_zero_point_correction = missing
            continue
          if re.search("Thermal correction to Energy", line): #T3
            try:
              out_energy_thermal_correction = float(line.split()[4])
            except:
              out_energy_thermal_correction = missing
            continue
          if re.search("Thermal correction to Enthalpy", line): #T4
            try:
              out_enthalpy_thermal_correction = float(line.split()[4])
            except:
              out_enthalpy_thermal_correction = missing
            continue
          if re.search("Thermal correction to Gibbs Free Energy", line): #T5
            try:
              out_gibbs_free_energy_thermal_correction = float(line.split()[6])
            except:
              out_gibbs_free_energy_thermal_correction = missing
            continue
          if re.search("Sum of electronic and zero-point Energies", line): #T6
            try:
              out_zero_point_energy = float(line.split()[6])
            except:
              out_zero_point_energy = missing
            continue
          if re.search("Sum of electronic and thermal Energies", line): #T7
            try:
              out_internal_energy = float(line.split()[6])
            except:
              out_internal_energy = missing
            continue
          if re.search("Sum of electronic and thermal Enthalpies", line): #T8
            try:
              out_enthalpy_energy = float(line.split()[6])
            except:
              out_enthalpy_energy = missing
            continue
          if re.search("Sum of electronic and thermal Free Energies", line): #T9
            try:
              out_gibbs_free_energy = float(line.split()[7])
            except:
              out_gibbs_free_energy = missing
            search += 1
            continue
        #ENTROPY
        if search==3:
          if re.search("Total", line): #E1
            try:
              out_entropy = float(line.split()[3])
            except:
              out_entropy = missing
            search += 1 
            continue
      #SAVE
      clusters_df = df_add_iter(clusters_df, "log", "time", [str(cluster_id)], [out_time]) #TIME
      clusters_df = df_add_iter(clusters_df, "log", "termination", [str(cluster_id)], [out_termination]) #TERMINATION
      clusters_df = df_add_iter(clusters_df, "log", "charge", [str(cluster_id)], [out_charge]) #I1
      clusters_df = df_add_iter(clusters_df, "log", "multiplicity", [str(cluster_id)], [out_multiplicity]) #I2
      clusters_df = df_add_iter(clusters_df, "log", "NAtoms", [str(cluster_id)], [out_NAtoms]) #I3
      clusters_df = df_add_iter(clusters_df, "log", "rotational_constants", [str(cluster_id)], [out_rotational_constants]) #O1
      clusters_df = df_add_iter(clusters_df, "log", "rotational_constant", [str(cluster_id)], [out_rotational_constant]) #O2
      clusters_df = df_add_iter(clusters_df, "log", "sp_electronic_energy", [str(cluster_id)], [out_sp_electronic_energy]) #O3
      clusters_df = df_add_iter(clusters_df, "log", "electronic_energy", [str(cluster_id)], [out_electronic_energy]) #O3
      clusters_df = df_add_iter(clusters_df, "log", "mulliken_charges", [str(cluster_id)], [out_mulliken_charges]) #O4
      clusters_df = df_add_iter(clusters_df, "log", "dipole_moment", [str(cluster_id)], [out_dipole_moment]) #O5
      clusters_df = df_add_iter(clusters_df, "log", "dipole_moments", [str(cluster_id)], [out_dipole_moments]) #O6
      clusters_df = df_add_iter(clusters_df, "log", "polarizability", [str(cluster_id)], [out_polarizability]) #O7
      clusters_df = df_add_iter(clusters_df, "log", "vibrational_frequencies", [str(cluster_id)], [out_vibrational_frequencies]) #V1
      clusters_df = df_add_iter(clusters_df, "log", "temperature", [str(cluster_id)], [out_temperature]) #V2
      clusters_df = df_add_iter(clusters_df, "log", "pressure", [str(cluster_id)], [out_pressure]) #V3
      clusters_df = df_add_iter(clusters_df, "log", "moments_of_inertia", [str(cluster_id)], [out_moments_of_inertia]) #V4
      clusters_df = df_add_iter(clusters_df, "log", "rotational_symmetry_number", [str(cluster_id)], [out_rotational_symmetry_number]) #T1
      clusters_df = df_add_iter(clusters_df, "log", "zero_point_correction", [str(cluster_id)], [out_zero_point_correction]) #T2
      clusters_df = df_add_iter(clusters_df, "log", "energy_thermal_correction", [str(cluster_id)], [out_energy_thermal_correction]) #T3
      clusters_df = df_add_iter(clusters_df, "log", "enthalpy_thermal_correction", [str(cluster_id)], [out_enthalpy_thermal_correction]) #T4
      clusters_df = df_add_iter(clusters_df, "log", "gibbs_free_energy_thermal_correction", [str(cluster_id)], [out_gibbs_free_energy_thermal_correction]) #T5
      clusters_df = df_add_iter(clusters_df, "log", "zero_point_energy", [str(cluster_id)], [out_zero_point_energy]) #T6
      clusters_df = df_add_iter(clusters_df, "log", "internal_energy", [str(cluster_id)], [out_internal_energy]) #T7
      clusters_df = df_add_iter(clusters_df, "log", "enthalpy_energy", [str(cluster_id)], [out_enthalpy_energy]) #T8
      clusters_df = df_add_iter(clusters_df, "log", "gibbs_free_energy", [str(cluster_id)], [out_gibbs_free_energy]) #T9
      clusters_df = df_add_iter(clusters_df, "log", "entropy", [str(cluster_id)], [out_entropy]) #E1
    file.close()

  ###############
  ### XYZ #######
  ###############
  if os.path.exists(file_i_XYZ):
    try:
      out = read(file_i_XYZ)
    except:
      out = missing
  else:
    out = missing
  clusters_df = df_add_iter(clusters_df, "xyz", "structure", [str(cluster_id)], [out])

  ###############
  ### ORCA ######
  ###############
  if os.path.exists(file_i_ORCA):
    file = open(file_i_ORCA, "r")
    out = missing
    for line in file:
      if re.search("FINAL SINGLE POINT ENERGY", line):
        out = float(line.split()[4])
    file.close()
    clusters_df = df_add_iter(clusters_df, orcaextname, "electronic_energy", [str(cluster_id)], [out])

  ###############
  ### TURBOMOLE #
  ###############
  if os.path.exists(file_i_TURBOMOLE):
    file = open(file_i_TURBOMOLE, "r")
    testTURBOMOLE = 0
    for i in range(5):
      if re.search("TURBOMOLE", file.readline()):
        testTURBOMOLE = 1
        break
    if testTURBOMOLE == 1:
      out = missing
      for line in file:
        if re.search("Final CCSD\(F12\*\)\(T\) energy", line):
          out = float(line.split()[5])
        if re.search("Final MP2 energy", line):
          out = float(line.split()[5])
      file.close()
      clusters_df = df_add_iter(clusters_df, turbomoleextname, "electronic_energy", [str(cluster_id)], [out])

## EXTRACT
def dash(input_array):
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "-":
      partbefore = input_array[0:element-1]
      partafter = input_array[element+2:]
      output_array_1 = []
      try:
        num1=int(input_array[element-1])
        num2=int(input_array[element+1])+1
        output_array = []
      except:
        break
      for thenumber in range(int(input_array[element-1]),int(input_array[element+1])+1):
        output_array_1.append(partbefore+[str(thenumber)]+partafter)
      for i in output_array_1:
        for j in dash(i): 
          output_array.append(j)
      break
  return output_array

def dash_comment(input_array):
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "-":
      output_array = []
      partbefore = input_array[0:element]
      output_array.append(partbefore)
      break
  return output_array

def comma(input_array):
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

if Qextract > 0:
  #print(Pextract)
  Pextract_comma = []
  for extract_i in Pextract:
    comma_separated = comma(extract_i)
    for separated_i in range(len(comma_separated)):
      Pextract_comma.append(comma_separated[separated_i])
  #print(Pextract_comma)    
  if Qclustername == 0:
    Pextract_ultimate = Pextract_comma
  else:
    #print(Pextract_comma)    
    Pextract_dash = []
    for extract_i in Pextract_comma:
      separated = seperate_string_number(extract_i)
      dash_separated = dash(separated)
      for separated_i in range(len(dash_separated)):
        Pextract_dash.append(dash_separated[separated_i])
    #print([listToString(i,"") for i in Pextract_dash])
    Pextract_dash = [ dash_comment(i)[0] for i in Pextract_dash ]
    #print(Pextract_dash)
    #print([listToString(i,"") for i in Pextract_dash])
    Pextract_comma2 = []
    for extract_i in Pextract_dash:
      comma2_separated = comma2(extract_i)
      for separated_i in range(len(comma2_separated)):
        Pextract_comma2.append(comma2_separated[separated_i])
    #print(Pextract_comma2)
    Pextract_sorted = []
    for extract_i in Pextract_comma2:
      if len(extract_i) > 1:
        array_sorted = sorted([extract_i[i:i + 2] for i in range(0, len(extract_i), 2)],key=lambda x: x[1])
        Pextract_sorted.append([item for sublist in array_sorted for item in sublist])
      else:
        Pextract_sorted.append(extract_i)
        Qclustername = 0
    #print(Pextract_sorted)
    Pextract_final = []
    for extract_i in Pextract_sorted:
      corrected = zeros(extract_i)
      if len(corrected) > 0:
        Pextract_final.append(corrected)
    Pextract_ultimate = unique(Pextract_final)
    #print(Pextract_ultimate)
  newclusters_df = pd.DataFrame()
  for extract_i in Pextract_ultimate:
    if Qclustername == 0:
      extracted_df = clusters_df[clusters_df["info"]["file_basename"].values == extract_i]
    else:
      extracted_df = clusters_df[clusters_df["info"]["cluster_type"].values == extract_i]
    if len(extracted_df) > 0:
      if len(newclusters_df) == 0:
        newclusters_df = extracted_df
      else:
        #print(newclusters_df)
        newclusters_df = newclusters_df.append(extracted_df)
        #print(newclusters_df)
  prelen=len(clusters_df)
  if Qextract == 2:
    a1 = pd.Index(clusters_df.index)
    a2 = pd.Index(newclusters_df.index)
    clusters_df = clusters_df.iloc[a1.difference(a2, sort=False)]
  else:
    clusters_df = newclusters_df
  if Qout == 1:
    print("Extracting: "+str(prelen)+" --> "+str(len(clusters_df)))

if Qreacted > 0:
  all_molecules = []
  dt = np.dtype(object)
  a = clusters_df
  for k in range(len(clusters_df)):
    b = a["xyz"]["structure"][k]
    if pd.isna(b):
      all_molecules.append(str(missing))
    else:
      p = b.positions
      symb = np.array(b.get_chemical_symbols())
      ind = [i != 'H' for i in symb]
  
      dist = lambda p1, p2: np.sqrt(np.sum(((p1-p2)**2)))
      dm = np.asarray([[dist(p1, p2) for p2 in p[ind]] for p1 in p[ind]])
  
      def bonded(x):
        if x < 1.75:
          return 1
        else:
          return 0
  
      bm = np.array([[bonded(j) for j in i] for i in dm])
  
      test = 0
      choosing_list = range(len(bm))
      molecules=[]
      if len(choosing_list) == 0:
        test = 1
      while test == 0:
        selected = [choosing_list[0]]
        test_chosen = 0
        j = -1
        while test_chosen == 0:
          j += 1
          #print(str(j)+"/"+str(len(bm)))
          #print(selected)
          chosen = selected[j]
          for i in choosing_list:
            if bm[choosing_list[j]][i] == 1 and choosing_list[j] != i and not (i in selected):
              selected.append(i)
          if len(selected)-1 == j:
            test_chosen = 1
        molecules.append(list(np.sort([symb[ind][i] for i in selected])))
        choosing_list = [i for i in choosing_list if i not in selected]
        if len(choosing_list) == 0:
          test = 1
      all_molecules.append(str(np.sort(np.array(molecules,dtype = dt))))

  def most_frequent(List):
      return max(set(List), key = List.count)
  
  #print(most_frequent(all_molecules))
  mf = most_frequent(all_molecules)
  #print(mf)
  #print(all_molecules[0])
  nind = [ i == mf for i in all_molecules]
  clusters_df = clusters_df[nind]

## SAVE OUTPUT.pkl ##
if Qout > 0:
  original_clusters_df = clusters_df.copy()

## PRE_EXTRACT_DATA MANIPULATION ##
if Qqha == 1:
  h = 6.626176*10**-34 #m^2 kg s^-1
  R = 1.987 #cal/mol/K #=8.31441
  k = 1.380662*10**-23 #m^2 kg s^-2 K^-1
  if Qanh != 1:
    # VIBRATIONAL FREQ MODIFICATION e.g. anharmonicity (vib.freq.,ZPE,ZPEc,U,Uc,H,Hc,S // G,Gc)
    for i in range(len(clusters_df)):
      try:
        QtOLD = clusters_df["log","temperature"].values[i]
        if pd.isna(clusters_df["log"]["vibrational_frequencies"].values[i]).any():
          continue
      except:
        QtOLD = clusters_df["log","temperature"].values[i]
      try:
        lf = float(clusters_df["log"]["vibrational_frequencies"].values[i][0])
      except:
        lf = 0
      if lf <= 0:
        clusters_df["log","entropy"][i] = missing
        clusters_df["log","enthalpy_energy"][i] = missing
        clusters_df["log","enthalpy_thermal_correction"][i] = missing
        clusters_df["log","internal_energy"][i] = missing
        clusters_df["log","energy_thermal_correction"][i] = missing
        clusters_df["log","zero_point_correction"][i] = missing
        clusters_df["log","zero_point_energy"][i] = missing
        continue
      Sv_OLD = np.sum([R*h*vib*2.99793*10**10/k/QtOLD/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K
      Ev_OLD = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
      #
      clusters_df["log","vibrational_frequencies"][i] = [Qanh * j for j in clusters_df["log","vibrational_frequencies"].values[i]]
      #
      Sv = np.sum([R*h*vib*2.99793*10**10/k/QtOLD/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K  
      Ev = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
      ###
      clusters_df["log","zero_point_correction"][i] = np.sum([0.5*h*vib*2.99793*10**10 for vib in clusters_df["log","vibrational_frequencies"][i]])*0.00038088*6.022*10**23/1000
      clusters_df["log","zero_point_energy"][i] = clusters_df["log","electronic_energy"][i] + clusters_df["log","zero_point_correction"][i]  
      clusters_df["log","internal_energy"][i] += (Ev - Ev_OLD)/1000/627.503    
      clusters_df["log","energy_thermal_correction"][i] += (Ev - Ev_OLD)/1000/627.503    
      clusters_df["log","enthalpy_energy"][i] += (Ev - Ev_OLD)/1000/627.503
      clusters_df["log","enthalpy_thermal_correction"][i] += (Ev - Ev_OLD)/1000/627.503
      clusters_df["log","entropy"][i] += Sv - Sv_OLD      
      ###
    
  
  # NEW TEMPERATURE (T,S,H,Hc,U,Uc // G,Gc)
  if ~np.isnan(Qt):
    for i in range(len(clusters_df)):
      try:
        if pd.isna(clusters_df["log"]["vibrational_frequencies"].values[i]):
          clusters_df["log","temperature"][i] = Qt
          continue
      except:
        QtOLD = clusters_df["log","temperature"].values[i]
      if Qt != QtOLD:
        try:
          lf = float(clusters_df["log"]["vibrational_frequencies"].values[i][0])
        except:
          lf = 0
        if lf <= 0:
          clusters_df["log","temperature"][i] = Qt
          clusters_df["log","entropy"][i] = missing
          clusters_df["log","enthalpy_energy"][i] = missing
          clusters_df["log","enthalpy_thermal_correction"][i] = missing
          clusters_df["log","internal_energy"][i] = missing
          clusters_df["log","energy_thermal_correction"][i] = missing
          continue
        Sv_OLD = np.sum([R*h*vib*2.99793*10**10/k/QtOLD/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K
        Sv = np.sum([R*h*vib*2.99793*10**10/k/Qt/(np.exp(h*vib*2.99793*10**10/k/Qt)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/Qt)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K
        Ev_OLD = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
        Ev = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/Qt)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
        ###
        clusters_df["log","temperature"][i] = Qt
        clusters_df["log","entropy"][i] += Sv - Sv_OLD + 4*R*np.log(Qt/QtOLD)
        clusters_df["log","enthalpy_energy"][i] += (Ev - Ev_OLD + 4*R*(Qt-QtOLD))/1000/627.503
        clusters_df["log","enthalpy_thermal_correction"][i] += (Ev - Ev_OLD + 4*R*(Qt-QtOLD))/1000/627.503
        clusters_df["log","internal_energy"][i] += (Ev - Ev_OLD + 3*R*(Qt-QtOLD))/1000/627.503
        clusters_df["log","energy_thermal_correction"][i] += (Ev - Ev_OLD + 3*R*(Qt-QtOLD))/1000/627.503
        ###
  
  # LOW VIBRATIONAL FREQUNECY TREATMENT (S // G,Gc)
  if Qfc > 0:
    for i in range(len(clusters_df)):
      try:
        if pd.isna(clusters_df["log"]["vibrational_frequencies"].values[i]):
          continue
      except:
        lf = 0
      try:
        lf = float(clusters_df["log"]["vibrational_frequencies"].values[i][0])
      except:
        lf = 0
      if lf <= 0:
        clusters_df["log","entropy"][i] = missing
        continue
      Qt = clusters_df["log"]["temperature"].values[i]
      vibs = clusters_df["log"]["vibrational_frequencies"].values[i]
      structure = clusters_df["xyz"]["structure"].values[i]
      
      mu = [float(h/(8*np.pi**2*2.99793*10**10*vib)) for vib in vibs]
      try:
        mi = np.mean(structure.get_moments_of_inertia())
      except:
        mi = missing
      #print(mu)
      #print(mi)
      Sr = [R*(0.5+np.log((8*np.pi**2.99793*(mu[j]*mi/(mu[j]+mi))*k*Qt/h**2)**0.5)) for j in range(len(mu))]  #cal/mol/K
      #print(Sr,flush=True)
      #print(vibs)
      Sv = [R*h*vib*2.99793*10**10/k/Qt/(np.exp(h*vib*2.99793*10**10/k/Qt)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/Qt)) for vib in vibs] #cal/mol/K
      w = [1/(1+(Qfc/vib)**4) for vib in vibs]
      Sv_corr = np.sum([w[j]*Sv[j]+(1-w[j])*Sr[j] for j in range(len(w))])
      Sv_each = np.sum(Sv)  #cal/mol/K
      #St_each = [R*np.log((2*np.pi*0.001*np.sum(structure.get_masses())*k**2/8.31441*Qt/h**2)**(3/2)*k*Qt/101325)+5/2 for structure in clusters_df["xyz"]["structure"].values]
      ###
      clusters_df["log","entropy"][i] = clusters_df["log","entropy"][i]+(Sv_corr-Sv_each)
    ###

#  if Qfc > 0:
#    Qt = clusters_df["log"]["temperature"].values
#    vibs = clusters_df["log"]["vibrational_frequencies"].values
#    structures = clusters_df["xyz"]["structure"].values
#    mu = [[h/(8*np.pi**2*2.99793*10**10*vib) for vib in vibs[i]] for i in range(len(vibs))]
#    mi = [np.mean(structure.get_moments_of_inertia()) for structure in structures]
#    Sr = [[R*(0.5+np.log((8*np.pi**2.99793*(mu[i][j]*mi[i]/(mu[i][j]+mi[i]))*k*Qt[i]/h**2)**0.5)) for j in range(len(mu[i]))] for i in range(len(mu))] #cal/mol/K
#    Sv = [[R*h*vib*2.99793*10**10/k/Qt[i]/(np.exp(h*vib*2.99793*10**10/k/Qt[i])-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/Qt[i])) for vib in vibs[i]] for i in range(len(vibs))] #cal/mol/K
#    w = [[1/(1+(Qfc/vib)**4) for vib in vibs[i]] for i in range(len(vibs))]
#    Sv_corr = [np.sum([w[i][j]*Sv[i][j]+(1-w[i][j])*Sr[i][j] for j in range(len(w[i]))]) for i in range(len(w))]
#    Sv_each = [np.sum(i) for i in Sv] #cal/mol/K
#    #St_each = [R*np.log((2*np.pi*0.001*np.sum(structure.get_masses())*k**2/8.31441*Qt/h**2)**(3/2)*k*Qt/101325)+5/2 for structure in clusters_df["xyz"]["structure"].values]
#    ###
#    clusters_df["log","entropy"] = [clusters_df["log","entropy"][i]+(Sv_corr[i]-Sv_each[i]) for i in range(len(clusters_df))]
#    ###

  ## CORRECTIONS FOR GIBBS FREE ENERGY
  for i in range(len(clusters_df)):
    try:
      clusters_df["log","gibbs_free_energy"][i] = clusters_df["log","enthalpy_energy"][i] - clusters_df["log","entropy"][i]/1000/627.503 * clusters_df["log","temperature"][i]
    except:
      clusters_df["log","gibbs_free_energy"][i] = missing
    try:
      clusters_df["log","gibbs_free_energy_thermal_correction"][i] = clusters_df["log","gibbs_free_energy"][i] - clusters_df["log","electronic_energy"][i]
    except:
      clusters_df["log","gibbs_free_energy_thermal_correction"][i] = missing

  #  clusters_df["log","gibbs_free_energy"] = [clusters_df["log","enthalpy_energy"][i] - clusters_df["log","entropy"][i]/1000/627.503 * clusters_df["log","temperature"][i] for i in range(len(clusters_df))] 
   #  clusters_df["log","gibbs_free_energy_thermal_correction"] = [clusters_df["log","gibbs_free_energy"][i] - clusters_df["log","electronic_energy"][i] for i in range(len(clusters_df))] 

###### FILTERING ######
if ( str(Qselect) != "0" or str(Quniq) != "0" ) and str(Qsort) == "0":
  Qsort = "g"
if str(Qsort) != "0":
  if Qsort == "g":
    Qsort = "gibbs_free_energy"
  if Qsort == "el":
    Qsort = "electronic_energy"
  clusters_df = clusters_df.sort_values([("log",Qsort)])
if str(Quniq) != "0":
  uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
  newclusters_df = []
  myNaN = lambda x : missing if x == "NaN" else x
  for i in uniqueclusters:
     preselected_df = clusters_df[clusters_df["info"]["cluster_type"] == i]
     tocompare = []
     for j in ["rg","electronic_energy"]:#,"gibbs_free_energy"]:
       if j == "rg":
         rg = []
         for aseCL in preselected_df["xyz"]["structure"]:
           try:
             rg.append((np.sum(np.sum((aseCL.positions-np.tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/np.sum(aseCL.get_masses()))**0.5)
           except:
             rg.append(missing)
         values = [np.floor(myNaN(o)*1e2) for o in rg]
       else:  
         values = [np.floor(myNaN(o)*1e4) for o in preselected_df["log"][j].values]
       tocompare.append(values)
     tocompare = np.transpose(tocompare)
     uniqueindexes = np.unique(tocompare,axis = 0,return_index=True)[1]
     selected_df = preselected_df.iloc[uniqueindexes]
     if len(newclusters_df) == 0:
       newclusters_df = selected_df
     else:
       #print(newclusters_df)
       newclusters_df = newclusters_df.append(selected_df)
       #print(newclusters_df)
  if Qout == 1:
    print("Uniqueness: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
  clusters_df = newclusters_df
if str(Qsort) != "0":
  clusters_df = clusters_df.sort_values([("log",Qsort)])
if str(Qselect) != "0":
  uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
  newclusters_df = []
  for i in uniqueclusters:
     selected_df = clusters_df[clusters_df["info"]["cluster_type"] == i][0:Qselect] 
     if len(newclusters_df) == 0:
       newclusters_df = selected_df
     else:
       #print(newclusters_df)
       newclusters_df = newclusters_df.append(selected_df)
       #print(newclusters_df)
  if Qout == 1:
    print("Selecting/Sampling: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
  clusters_df = newclusters_df

## SAVE OUTPUT.pkl ##
if Qout > 0:
  tosave = original_clusters_df.loc[clusters_df.index]
  tosave.to_pickle(output_pkl)
  if Qout == 1:
    if len(tosave) == 0:
      print(tosave)
      print("No files in the input!")
    else:
      print("Example output:")
      print(tosave.iloc[0])
      print("Number of files in "+output_pkl+": "+str(len(tosave)))

## EXTRACT DATA ##
output = []
for i in Pout:
  if i == "-cite":
    try:
      output.append(clusters_df["info"]["citation"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  #XYZ
  if i == "-xyz":
    for ind in clusters_df.index:
      write(clusters_df["info"]["file_basename"][ind]+".xyz",clusters_df["xyz"]["structure"][ind])
    continue
  #MOVIE
  if i == "-movie":
    f = open("movie.xyz","w")
    for ind in clusters_df.index:
      write(".movie.xyz",clusters_df["xyz"]["structure"][ind])
      f2 = open(".movie.xyz", 'r')
      f.write(f2.read())
      f2.close()
    f.close()
    continue
  #Rg
  if i == "-rg":
    rg = []
    for aseCL in clusters_df["xyz"]["structure"]:
      try:
        rg.append((np.sum(np.sum((aseCL.positions-np.tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/np.sum(aseCL.get_masses()))**0.5)
      except:
        rg.append(missing)
    output.append(rg)
    continue
  if i == "-radius" or i == "-radius0.5":
    radius = []
    for aseCL in clusters_df["xyz"]["structure"]:
      try:
        dist = lambda p1, p2: np.sqrt(np.sum((p1-p2)**2))
        centered = aseCL.positions-aseCL.positions.mean(axis = 0)
        ratios = np.sqrt(np.linalg.eigvalsh(np.dot(centered.transpose(),centered)/len(centered)))
        maxdist = np.max(np.asarray([[dist(p1, p2) for p2 in aseCL.positions] for p1 in aseCL.positions]))
        if i == "-radius0.5":
          if ratios[0] < 0.5:
            ratios[0] = 0.5
          if ratios[1] < 0.5:
            ratios[1] = 0.5
          if ratios[2] < 0.5:
            ratios[2] = 0.5
          if maxdist < 0.5:
            maxdist = 0.5
        if np.max(ratios) > 0:
          ratios = ratios / np.max(ratios)
        else:
          ratios = missing
        radius.append((np.prod(maxdist*ratios))**(1/3))	
      except:
        radius.append(missing)
    output.append(radius)
    continue
  #INFO
  if i == "-ct":
    try: 
      output.append(clusters_df["info"]["cluster_type"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-b": 
    try:
      output.append(clusters_df["info"]["file_basename"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-nOUT":
    try:
      output.append(clusters_df["info"]["file_basename"].values+".out")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-nLOG":
    try:
      output.append(clusters_df["info"]["file_basename"].values+".log")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-nXYZ":
    try:
      output.append(clusters_df["info"]["file_basename"].values+".xyz")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pOUT":
    try:
      output.append(clusters_df["info"]["folder_path"].values+clusters_df["info"]["file_basename"].values+".out")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pLOG":
    try:
      output.append(clusters_df["info"]["folder_path"].values+clusters_df["info"]["file_basename"].values+".log")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pXYZ":
    try:
      output.append(clusters_df["info"]["folder_path"].values+clusters_df["info"]["file_basename"].values+".xyz")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-ePKL":
    if Qout > 0 or len(input_pkl) == 1:
      if Qout > 0:
        try:
          output.append(os.path.abspath(output_pkl)+"/:EXTRACT:/"+clusters_df["info"]["file_basename"].values)
        except:
          output.append([missing]*len(clusters_df))
      else:
        try:
          output.append(os.path.abspath(input_pkl[0])+"/:EXTRACT:/"+clusters_df["info"]["file_basename"].values)
        except:
          output.append([missing]*len(clusters_df))
    else:
      print("Sorry but it seems to me that you are taking this from more file or no file is formed (yet)")
      exit()
    continue
  if i == "-elsp":
    try:
      output.append(QUenergy*clusters_df["log"]["sp_electronic_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-el":
    try:
      output.append(QUenergy*clusters_df["log"]["electronic_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-elout":
    try:
      output.append(QUenergy*clusters_df["out"]["electronic_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-elc":
    try:
      output.append(QUenergy*(clusters_df["out"]["electronic_energy"].values-clusters_df["log"]["electronic_energy"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-uc":
    try:
      output.append(QUenergy*clusters_df["log"]["energy_thermal_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-u":
    try:
      output.append(QUenergy*(clusters_df["log"]["electronic_energy"].values+clusters_df["log"]["energy_thermal_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-uout":
    try:
      output.append(QUenergy*(clusters_df["out"]["electronic_energy"].values+clusters_df["log"]["energy_thermal_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-zpec": 
    try:
      output.append(QUenergy*clusters_df["log"]["zero_point_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-zpe":
    try:
      output.append(QUenergy*(clusters_df["log"]["electronic_energy"].values+clusters_df["log"]["zero_point_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-zpeout":
    try:
      output.append(QUenergy*(clusters_df["out"]["electronic_energy"].values+clusters_df["log"]["zero_point_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-g":
    try:
      output.append(QUenergy*clusters_df["log"]["gibbs_free_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-gc":
    try:
      output.append(QUenergy*clusters_df["log"]["gibbs_free_energy_thermal_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-gout":
    try:
      output.append(QUenergy*(clusters_df["log"]["gibbs_free_energy"].values+clusters_df["out"]["electronic_energy"].values-clusters_df["log"]["electronic_energy"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-h":
    try:
      entalpies = clusters_df["log"]["enthalpy_energy"].values
      entalpies[entalpies == "NaN"] = missing 
      output.append(QUenergy*entalpies)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-hc":
    try:
      output.append(QUenergy*clusters_df["log"]["enthalpy_thermal_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-hout":
    try:
      output.append(QUenergy*(clusters_df["log"]["enthalpy_energy"].values+clusters_df["out"]["electronic_energy"].values-clusters_df["log"]["electronic_energy"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-s":
    try:
      entropies = clusters_df["log"]["entropy"].values
      entropies[entropies == "NaN"] = missing
      output.append(QUentropy*entropies/1000/627.503)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-lf": 
    lowestfreq = []
    for aseCL in clusters_df["log"]["vibrational_frequencies"]:
      try:
        lowestfreq.append(aseCL[0])
      except:
        lowestfreq.append(missing)
    output.append(lowestfreq)
    continue
  if i == "-f":
    try:
      output.append(clusters_df["log"]["vibrational_frequencies"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-rot":
    try:
      output.append(clusters_df["log"]["rotational_constant"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-rots":
    try:
      output.append(clusters_df["log"]["rotational_constants"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-mult":
    try:
      output.append(clusters_df["log"]["multiplicity"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-char":
    try:
      output.append(clusters_df["log"]["charge"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-mull":
    try:
      output.append(clusters_df["log"]["mulliken_charges"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-dip":
    try:
      output.append(clusters_df["log"]["dipole_moment"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-dips":
    try:
      output.append(clusters_df["log"]["dipole_moments"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pol":
    try:
      output.append(clusters_df["log"]["polarizability"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-templog":
    try:
      output.append(clusters_df["log"]["temperature"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-preslog":
    try:
      output.append(clusters_df["log"]["pressure"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-mi":
    try:
      output.append([structure.get_moments_of_inertia() for structure in clusters_df["xyz"]["structure"].values])
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-ami":
    try:
      output.append([np.mean(structure.get_moments_of_inertia()) for structure in clusters_df["xyz"]["structure"].values])
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-rsn":
    try:
      output.append(clusters_df["log"]["rotational_symmetry_number"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-t":
    try:
      output.append(clusters_df["log"]["time"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-termination":
    try: 
      output.append(clusters_df["log"]["termination"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  output.append(["UNKNOWN_ARGUMENT"]*len(clusters_df))

## PRINT DATA ##
if not len(output) == 0:
  output = np.array(output)

  #TAKING GLOBAL MINIMA ONLY
  if Qglob == 1 or Qglob == 2:
    uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
    indexes = []
    for i in uniqueclusters:
      if Qglob == 1:
        GFE = clusters_df["log","gibbs_free_energy"].values[clusters_df["info"]["cluster_type"] == i]
      elif Qglob == 2:
        GFE = clusters_df["log","gibbs_free_energy"].values + clusters_df["out"]["electronic_energy"].values - clusters_df["log"]["electronic_energy"].values 
        GFE = np.array(GFE)[clusters_df["info"]["cluster_type"] == i]
      else:
        print("Qglob error. [EXITING]")
        exit()
      GFE[GFE == "NaN"] = missing
      if len(GFE[~pd.isna(GFE)]) != 0:
        globindex = clusters_df.index.values[clusters_df["info"]["cluster_type"] == i][GFE == np.min(GFE[~pd.isna(GFE)])][0]
        indexes.append(int(np.array(range(len(clusters_df)))[clusters_df.index.values == globindex][0]))
   
    output = [[output[j][i] for i in indexes] for j in range(output.shape[0])]
  
  #TAKING BOLTZMANN AVERAGE OVER ALL MINIMA 
  if Qbavg == 1 or Qbavg == 2:
    k = 1.380662*10**-23 # [J/K]
    uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
    portions = []
    freeenergies = [] 
    entropies = []
    indexes = []
    for i in uniqueclusters:
      if Qbavg == 1:
        GFE = clusters_df["log","gibbs_free_energy"].values[clusters_df["info"]["cluster_type"] == i]
        temps = clusters_df["log","temperature"].values[clusters_df["info"]["cluster_type"] == i]
      elif Qbavg == 2:
        GFE = clusters_df["log","gibbs_free_energy"].values + clusters_df["out"]["electronic_energy"].values - clusters_df["log"]["electronic_energy"].values
        GFE = np.array(GFE)[clusters_df["info"]["cluster_type"] == i]
        temps = clusters_df["log","temperature"].values[clusters_df["info"]["cluster_type"] == i]
      else:
        print("Qglob error. [EXITING]")
        exit()
      nonans = ~pd.isna(GFE)
      GFE = GFE[nonans]
      minimum = np.min(GFE)
      GFE = GFE-minimum
      preportions = [np.exp(-GFE[i]*43.60*10**-19/k/temps[nonans][i]) for i in range(GFE.shape[0])]
      freeenergies.append(QUenergy*(minimum - 1/43.60/10**-19*k*temps[nonans][0]*np.log(np.sum([np.exp(GFE[i]*43.60*10**-19/k/temps[nonans][i]) for i in range(GFE.shape[0])]))))
      sumpreportions = np.sum(preportions)
      portions.append(preportions/sumpreportions)
      indexes.append(np.array(range(len(clusters_df)))[clusters_df["info"]["cluster_type"] == i][nonans])
    
    def is_averagable(input_array): # :-D
      try:
        [float(i) for i in input_array]
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

    def myif(cond,opt1,opt2,opt3):
      if Pout[cond] == "-g" or Pout[cond] == "-gout":
        return opt2
      elif Pout[cond] == "-s":
        return opt3
      else:
        return opt1
    
    output = [[myif(l,np.sum([portions[i][j]*output[l,indexes[i][j]] for j in range(len(portions[i]))]),freeenergies[i],missing) if is_averagable(output[l][indexes[i]]) else is_the_same(output[l,indexes[i]]) for i in range(len(portions))] for l in range(output.shape[0])]

  fn = ".help"+str(np.random.randint(100000,size=1)[0])+".txt" 
  f = open(fn, "w")
  [f.write(" ".join(map(str,row))+"\n") for row in list(zip(*output))]
  f.close()
  #TODO can you make this working using only python?
  if Qout != 2 or Qformation == 0:
    os.system("cat "+fn+" | column -t")
    os.remove(fn)

def myFunc(e):
  try:
    numero = np.sum([int(i) for i in seperate_string_number(dash_comment(e[0])[0])[0::2]])
  except:
    numero = 0
  return numero+len(e[0])/100.
dt = np.dtype(object)

if Qformation == 1:
  if Qout != 2:
    print("#####################################",flush = True)
    print("##########  FORMATION  ##############",flush = True)
    print("#####################################",flush = True)
  if len(formation_input_file) > 0:
    f = open(formation_input_file, "r")
    output = []
    for line in f.readlines():
      output.append(line.split())
    output = np.array(output,dtype=dt).transpose()
    f.close()
  #SORT OUTPUT
  output = np.transpose(np.transpose(output)[np.apply_along_axis(myFunc, axis=1, arr=np.transpose(output)).argsort()])
  #
  cluster_types = [dash_comment(seperate_string_number(i))[0] for i in output[0]]
  ######## SOLVING PROTONATION
  for i in range(len(cluster_types)):
    chosen_cluster_type = cluster_types[i]
    if "p" in chosen_cluster_type:
      protonated_base = ""
      if "gd" in chosen_cluster_type:
        protonated_base = "gd"
      elif "eda" in chosen_cluster_type:
        protonated_base = "eda"
      elif "tma" in chosen_cluster_type:
        protonated_base = "tma"
      elif "dma" in chosen_cluster_type:
        protonated_base = "dma"
      elif "am" in chosen_cluster_type:
        protonated_base = "am"
      elif "w" in chosen_cluster_type:
        protonated_base = "w"
      new_cluster_type = []
      for j in range(1,len(chosen_cluster_type),2):
        if chosen_cluster_type[j] == protonated_base:
          subsctracted = int(chosen_cluster_type[j-1]) - 1
          if subsctracted > 0:
            new_cluster_type.append(str(subsctracted))
            new_cluster_type.append(chosen_cluster_type[j])
        elif chosen_cluster_type[j] == "p":
          if int(chosen_cluster_type[j-1]) > 1:
            new_cluster_type.append(str(int(chosen_cluster_type[j-1]) - 1))
            new_cluster_type.append(chosen_cluster_type[j])
          else:
            continue
        else:
          new_cluster_type.append(chosen_cluster_type[j-1])
          new_cluster_type.append(chosen_cluster_type[j])
      new_cluster_type.append("1")
      new_cluster_type.append("protonated_"+protonated_base)       
      cluster_types[i] = new_cluster_type 
  ########################        
  #print(cluster_types)
  cluster_types_sorted = [sorted([extract_i[i:i + 2] for i in range(0, len(extract_i), 2)],key=lambda x: x[1]) for extract_i in cluster_types]
  #print(cluster_types_sorted)
  monomers = np.array([np.sum([int(j) for j in np.array(i)[:,0]]) for i in cluster_types_sorted]) == 1
  dt = np.dtype(object)
  #print(np.array(cluster_types,dtype=dt)[monomers])
  monomer_types = [i[1] for i in np.array(cluster_types,dtype=dt)[monomers]]
  #print(monomer_types)
  if Qout != 2:
    print("MONOMERS: " + " ".join(monomer_types),flush = True)
  f = open(".help.txt", "w")
  for i in range(len(output[0])):
    line = np.array(output)[:,i]
    cluster_total_number = np.sum([int(sel[0]) for sel in cluster_types_sorted[i]])
    for j in range(len(cluster_types_sorted[i])):
      cluster_molecule = cluster_types_sorted[i][j][1]
      cluster_molecule_number = cluster_types_sorted[i][j][0]
      #print(np.array(output)[:,monomers])
      test_monomer = 0
      for k in range(len(np.array(output)[:,monomers][0])):
        #selected_monomer = seperate_string_number(np.array(output)[:,monomers][0,k])[1]
        selected_monomer = np.array(cluster_types_sorted,dtype=dt)[monomers][k][0][1]
        #print("--------------------------")
        #print(selected_monomer)
        #print(cluster_molecule)
        if cluster_molecule == selected_monomer:
          for line_i in range(1,len(line)):
            try:
              line[line_i] = float(line[line_i]) - float(cluster_molecule_number) * float(np.array(output)[:,monomers][line_i,k])
              if Qconc > 0:
                for conc_j in range(len(conc)):
                  if conc[conc_j][0] == selected_monomer:         
                    R = 1.987 #cal/mol/K #=8.31441
                    if np.isnan(Qt):
                      try:
                        Qt = clusters_df["log","temperature"].values[i]
                      except:
                        Qt = 298.15
                    if np.isnan(Qp):
                      try:
                        Qp = 101325.0*float(clusters_df["log","pressure"].values[i])
                      except:
                        Qp = 101325.0
                    line[line_i] = float(line[line_i]) - QUenergy*(float(cluster_molecule_number) - CNTfactor/cluster_total_number) * R/1000/627.503 * Qt * np.log( float(conc[conc_j][1]) / Qp)
            except:
              line[line_i] = missing
          test_monomer = 1
      if test_monomer == 0:
        line[1:] = missing 
    f.write(" ".join(map(str,line))+"\n")
  f.close()
  os.system("cat .help.txt | column -t")
  os.remove(".help.txt")
