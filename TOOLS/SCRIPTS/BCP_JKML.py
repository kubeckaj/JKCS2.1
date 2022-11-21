import sys
import os.path
import subprocess

def help():
  print("-method <str>                      direct|delta [default = delta]", flush = True)
  print("-size <int>                        random number of samples for training set", flush = True)
  print("-train <file_HIGH> [<file_LOW>]    train on given pikled files", flush = True)
  print("-trained <file_VARS-PKL>           take pre-trained ML", flush = True)
  print("-eval <file_STRS> [<file_LOW>]     evaluate energies of NEW structures in pickled file(s)", flush = True)
  print("-test <file_HIGH> [<file_LOW>]     validate ML on energies of structures in pickled file(s)", flush = True)
  print("-monomers <file_HIGH> [<file_LOW>] binding energies with respect to monomer(s) in in pickled file(s)", flush = True)
  print("    /or/  none                     training directly on el.energies (not good for mix of clusters)", flush = True)
  print("-sigma <X> -lambda <Y>             hyperparameters [default: sigma = 1.0 and lambda = 1e-4]", flush = True)
  print("OTHER: -split X, -startsplit X, -finishsplit X, -array, -optimize", flush = True)
  print("OUTPUTFILES: -out X.pkl [def = predicted_QML.pkl], -varsout X.pkl [def = vars.pkl]")
  
#PREDEFINED ARGUMENTS
method = "delta"
size = "full"
TRAIN_HIGH = ""
TRAIN_LOW = ""
TEST_DATABASE = ""
TEST_LOW = ""
TEST_HIGH = ""
Qtrain = 0 #0=nothing (fails), 1=train, 2=trained
Qsplit = 1 #1=no split, how many splits to do; ONLY FOR TRAINING
Qsplit_i = 1; Qsplit_j = 1; 
Qeval = 0 #0=nothing (possible), 1=validate, 2=eval
Qopt = 0 #0=nothing (possible), 1=optimize
Qmonomers = 0 #0=monomers taken from database, 1=monomers in separate files, 2=no monomer subtraction

#PREDEFINED QML
sigmas = [1.0]
lambdas = [1e-4]*len(sigmas)
#kernel_args = {
#            "cut_distance": 1e1,
#            "cut_start": 1.0,
#            "two_body_width": 1.0, #0.1
#            "two_body_scaling": 200.0,#2.0
#            "two_body_power": 6.0,#6.0
#            "three_body_width": 3.0,#3.0
#            "three_body_scaling": 20.0,
#            "three_body_power": 3.0,
#            "alchemy": "off", #"periodic-table"
#            "alchemy_period_width": 1.0,
#            "alchemy_group_width": 1.0,#1.0
#            "fourier_order": 3}  
kernel_args = {}

outfile="predicted_QML.pkl"
varsoutfile="vars.pkl"

#TREATING ARGUMENTS
last=""
for i in sys.argv[1:]:
  #HELP
  if i == "-help":
    help()
    exit()
  #VARS OUT FILE
  if i == "-varsout":
    last = "-varsout"
    continue
  if last == "-varsout": 
   varsoutfile = i
   varsoutfile = varsoutfile.split(".pkl")[0]+".pkl"
   last = ""
   continue
  #OUT FILE
  if i == "-out":
    last = "-out"
    continue
  if last == "-out":
   outfile = i
   outfile = outfile.split(".pkl")[0]+".pkl"
   last = ""
   continue
  #METHOD
  if i == "-method":
    last = "-method"
    continue
  if i == "direct" or last == "-method" or i == "delta":
    last = ""
    if i == "direct":
      method = "direct"
    if i == "delta":
      method = "delta"
    continue
  #HYPERPARAMETERS
  if i == "-sigma":
    last = "-sigma"
    continue
  if last == "-sigma":
    last = ""
    sigmas = [float(i)]
    continue
  if i == "-lambda":
    last = "-lambda"
    continue
  if last == "-lambda":
    last = ""
    lambdas = [float(i)]*len(sigmas)
    continue
  #TRAINING SIZE
  if i == "-size":
    last = "-size"
    continue
  if last == "-size":
    last = ""
    try:
      if i != "full":
        size = int(i)
    except:
      print("Wrong argument for size [EXITING]", flush = True)
      exit()
    continue
  #DATABASES
  if i == "-train":
    last = "-train"
    continue
  if i == "-test":
    last = "-test"
    continue
  if i == "-eval":
    last = "-eval"
    continue
  if i == "-opt":
    last = "-opt"
    continue
  if i == "-monomers":
    last = "-monomers"
    continue
  if i == "-trained":
    if Qtrain == 0:
      Qtrain = 2
    else:
      print("cannot take trained if training [EXITING]")
      exit
    last = "-trained"
    continue
  #SPLIT
  if i == "-finishsplit":
    Qsplit=-1
  if i == "-split" or i == "-startsplit" or i == "-finishsplit":
    last = "-split"
    continue
  if last == "-split":
    last = "-splitt"
    Qsplit = int(i)
    continue
  if last == "-splitt":
    last = "-splittt"
    Qsplit_i = int(i)-1
    continue
  if last == "-splittt":
    last = ""
    Qsplit_j = int(i)-1
    continue
  #TRAIN DATABASE(S)
  if last == "-trained":
    last = ""
    VARS_PKL = i
    continue
  if last == "-train":
    TRAIN_HIGH = i
    if Qtrain == 0:
      Qtrain = 1
    else:
      print("cannot train if taking trained [EXITING]")
      exit
    last = "-train2"
    continue
  if last == "-train2":
    TRAIN_LOW = i
    last = ""
    continue
  #TEST/EVAL/OPT DATABASE(S)
  if last == "-eval":
    TEST_HIGH = i
    Qeval = 2
    last = "-test2"
    continue
  if last == "-test":
    TEST_HIGH = i
    Qeval = 1
    last = "-test2"
    continue
  if last == "-opt":
    TEST_HIGH = i
    Qopt = 1
    last = "-test2"
    continue
  if last == "-test2":
    TEST_LOW = i
    last = ""
    continue
  #MONOMER DATABASE(S)
  if last == "-monomers":
    if i == "0" or i == "none" or i == "no":
      Qmonomers=2
      last = ""
    else:
      MONOMERS_HIGH = i
      Qmonomers = 1
    last = "-monomers2"
    continue
  if last == "-monomers2":
    MONOMERS_LOW = i
    last = ""
    continue
  #UNKNOWN ARGUMENT
  print("Sorry cannot understand this argument: "+i)
  exit()

## TEST IF ALL REQUIRED FILES EXIST
if Qtrain == 1:
  if not os.path.exists(TRAIN_HIGH):
    print("Error reading file. TRAIN_HIGH = "+TRAIN_HIGH)
    exit()
  if method == "delta":
    if not os.path.exists(TRAIN_LOW):
      print("Error reading file. TRAIN_LOW = "+TRAIN_LOW)
      exit()
if Qtrain == 2:
  if not os.path.exists(VARS_PKL):
    print("Error reading file. VARS_PKL = "+VARS_PKL)
    exit()
if Qeval == 1 or Qeval == 2:
  if not os.path.exists(TEST_HIGH):
    print("Error reading file. + TEST_HIGH = "+TEST_HIGH)
    exit()
  if method == "delta":
    if not os.path.exists(TEST_LOW):
      print("Error reading file. + TEST_LOW = "+TEST_LOW)
      exit()
if Qmonomers == 1:
  if not os.path.exists(MONOMERS_HIGH):
    print("Error reading file. MONOMERS_HIGH = "+MONOMERS_HIGH)
    exit()
  if method == "delta":
    if not os.path.exists(MONOMERS_LOW):
      print("Error reading file. MONOMERS_LOW = "+MONOMERS_LOW)
      exit()

#IMPORTING THE REST OF LIBRARIES
import pandas as pd
import numpy as np
from ase import Atoms
from ase.io import read, write
from sklearn.model_selection import train_test_split
import pickle
#import scipy

#from qml import fchl
#from qml import Compound
from qml.fchl import generate_representation
from qml.fchl import get_local_symmetric_kernels
from qml.fchl import get_local_kernels
#from qml.fchl import get_atomic_symmetric_kernels
#from qml.fchl import get_atomic_kernels
#from qml.fchl import get_global_symmetric_kernels
#from qml.fchl import get_global_kernels
#from qml.kernels import gaussian_kernel
from qml.math import cho_solve
#import matplotlib.pyplot as plt

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

######################
###### TRAINING ######
######################

if Qtrain == 0:
  print("training option is missing [EXITING]", flush = True)
  exit()

#TRAIN
if Qtrain == 1:
  ### DATABASE LOADING ###
  ## The high level of theory
  clusters_df = pd.read_pickle(TRAIN_HIGH).sort_values([('info','file_basename')])
  ens = (clusters_df["log"]["electronic_energy"]).values.astype("float")
  strs = clusters_df["xyz"]["structure"]
  ## The low level of theory
  if method == "delta":
    clusters_df2 = pd.read_pickle(TRAIN_LOW).sort_values([('info','file_basename')])
    ens2 = (clusters_df2["log"]["electronic_energy"]).values.astype("float")
    #str2 should be the same as str by princip
  
  ### REPRESENTATION CALCULATION ###
  repres_dataframe = pd.DataFrame(index = strs.index, columns = ["xyz"])
  max_atoms = max([len(strs[i].get_atomic_numbers()) for i in range(len(strs))])
  for i in range(len(repres_dataframe)):
    repres_dataframe["xyz"][i] = generate_representation(strs[i].get_positions(), strs[i].get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=10.0)
  fchl_representations = np.array([mol for mol in repres_dataframe["xyz"]])

  #some info about the full representation
  print(fchl_representations.shape, flush = True)
  print(ens.shape, flush = True)
  if method == "delta":
    print(ens2.shape, flush = True)
    if ens.shape != ens2.shape:
      print("The LOW and HIGH method train sizes do not match. [EXITING]")
      exit()

  ### BINDING ENERGIES CALCULATION ###
  if Qmonomers == 2: #no difference from monomers at all
    ens_correction = [0]*len(ens)
    if method == "delta":
      ens2_correction = [0]*len(ens2)
  else:
    #HIGH LEVEL
    if Qmonomers == 1:
      clusters_df0 = pd.read_pickle(MONOMERS_HIGH)
      ins_monomers = [np.sum(clusters_df0["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df0))]
      monomers_df = clusters_df0[ins_monomers]
    else:
      ins_monomers = [np.sum(clusters_df["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df))]
      monomers_df = clusters_df[ins_monomers]
    ens_correction = [0]*len(ens)
    n = 0
    for i in np.transpose([clusters_df["info"]["components"].values,clusters_df["info"]["component_ratio"].values]):
      nn = 0
      for j in i[0]:
        test = 0
        for k in range(len(monomers_df)):
          monomer_k = monomers_df.iloc[k]
          monomer_k_name = monomer_k["info"]["components"]
          if j == monomer_k_name[0]:
            ens_correction[n] += monomer_k["log"]["electronic_energy"]*i[1][nn]
            test = 1
            break
        if test == 0:
          print("OMG; monomer "+j+" was not found in:")
          print(monomers_df["info"]["components"].values)
          exit()
        nn += 1
      n += 1
    #LOW LEVEL
    if method == "delta":
      if Qmonomers == 1:
        clusters_df0 = pd.read_pickle(MONOMERS_LOW)
        ins_monomers = [np.sum(clusters_df0["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df0))]
        monomers_df = clusters_df0[ins_monomers]
      else:
        ins_monomers = [np.sum(clusters_df2["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df2))]
        monomers_df = clusters_df2[ins_monomers]
      ens2_correction = [0]*len(ens2)
      n = 0
      for i in np.transpose([clusters_df2["info"]["components"].values,clusters_df2["info"]["component_ratio"].values]):
        nn = 0
        for j in i[0]:
          test = 0
          for k in range(len(monomers_df)):
            monomer_k = monomers_df.iloc[k]
            monomer_k_name = monomer_k["info"]["components"]
            if j == monomer_k_name[0]:
              ens2_correction[n] += monomer_k["log"]["electronic_energy"]*i[1][nn]
              test = 1
              break
          if test == 0:
            print("OMG; monomer "+j+" was not found in:")
            print(monomers_df["info"]["components"].values)
            exit()
          nn += 1
        n += 1
  
  #The binding (formation) energy calculation
  form_ens = ens - ens_correction
  #print(form_ens, flush = True)
  if method == "delta":
    form_ens2 = ens2 - ens2_correction
    #print(form_ens2, flush = True)

  ### DEFINING THE TRAINING Xs and Ys:  Y = QML(X) 
  #the full set
  X_train0 = fchl_representations
  if method == "delta":
    Y_train0 = form_ens - form_ens2
  else: 
    Y_train0 = form_ens
  #Do we take only subset for training?
  if size != "full":
    X_train, X_trash, idx, idx_trash = train_test_split(X_train0, range(len(X_train0)), test_size=(len(X_train0)-size)/len(X_train0), random_state=1)
    Y_train = Y_train0[idx]
    size = "full" #IT IS BECAUSE I DO WANT TO MAKE MY TEST SET SMALLER
  else:
    X_train = X_train0
    Y_train = Y_train0
  
  ### QML TRAINING ###
  #TODO if splitting is required include this:
  #splits=int(sys.argv[3])
  #split_taken_1=int(sys.argv[4])-1
  #split_taken_2=int(sys.argv[5])-1
  #X_train_1 = np.array_split(X_train,splits)[split_taken_1]
  #X_train_2 = np.array_split(X_train,splits)[split_taken_2]
  #Y_train_1 = np.array_split(Y_train,splits)[split_taken_1]
  #Y_train_2 = np.array_split(Y_train,splits)[split_taken_2]
  #if split_taken_1 == split_taken_2:
  #  K = get_local_symmetric_kernels(X_train_2, sigmas, **kernel_args)
  #  K = [K[i] + lambdas[i]*np.eye(len(K[i])) for i in range(len(sigmas))]
  #  #alpha = [cho_solve(Ki, Y_train_2) for Ki in K]
  #else:
  #  K = get_local_kernels(X_train_1,X_train_2, sigmas, **kernel_args)
  
  #ONLY FOR JOINING ALL THE SPLITS AND CHOLESKY DECOMPOSITION
  if Qsplit == -1:
    splits = Qsplit_i+1
    K = [];
    for i in range(0,splits):
      Kx = []
      for j in range(0,splits):
        if i < j:
          s1 = j
          s2 = i
        else:
          s1 = i
          s2 = j
        f = open(varsoutfile.split(".pkl")[0]+"_"++str(splits)+"_"+str(s1)+"_"+str(s2)+".pkl","rb")
        Kcell, Y_train = pickle.load(f)
        if i > j:
          Kcell = np.transpose(Kcell[0])
        else:
          Kcell = Kcell[0]
        if len(Kx) == 0:
          Kx = Kcell
        else:
          Kx = np.concatenate((Kx,Kcell))
        f.close()
      if len(K) == 0:
        K = Kx
      else:
        K = np.concatenate((K,Kx),axis = 1)
    alpha = [cho_solve(K, Y_train)]
    #alpha = [cho_solve(Ki, Y_train_1) for Ki in K]
    f = open(varsoutfile,"wb")
    pickle.dump([X_train, sigmas, alpha],f)
    f.close()
    print("Training completed.", flush = True)
  elif Qsplit == 1:
    K = get_local_symmetric_kernels(X_train, sigmas, **kernel_args)       #calculates kernel
    K = [K[i] + lambdas[i]*np.eye(len(K[i])) for i in range(len(sigmas))] #corrects kernel
    alpha = [cho_solve(Ki, Y_train) for Ki in K]                          #calculates regression coeffitients

    #I will for now everytime save the trained QML
    f = open(varsoutfile,"wb")
    pickle.dump([X_train, sigmas, alpha],f)
    f.close()
    print("Training completed.", flush = True)
  else:
    X_train_i = np.array_split(X_train,Qsplit)[Qsplit_i]
    X_train_j = np.array_split(X_train,Qsplit)[Qsplit_j]
    if Qsplit_i == Qsplit_j:
      K = get_local_symmetric_kernels(X_train_i, sigmas, **kernel_args)       #calculates kernel
      K = [K[i] + lambdas[i]*np.eye(len(K[i])) for i in range(len(sigmas))] #corrects kernel
    else:
      K = get_local_kernels(X_train_i, X_train_j, sigmas, **kernel_args)
    f = open(varsoutfile.split(".pkl")[0]+"_"+str(Qsplit)+"_"+str(Qsplit_i)+"_"+str(Qsplit_j)+".pkl","wb")
    pickle.dump([K,Y_train],f)
    f.close()
    exit()

#LOAD TRAINING
#TODO collect splitting: /home/kubeckaj/ML_SA_B/ML/TRAIN/TEST/SEPARATE/cho_solve.pkl
if Qtrain == 2:
  f = open(VARS_PKL,"rb")
  X_train, sigmas, alpha = pickle.load(f)
  if len(alpha)!=1:
    alpha = [alpha]
  f.close()
  print("Training loaded.", flush = True)


######################
###### EVALUATE ######
######################

if Qeval == 1 or Qeval == 2:
  ### DATABASE LOADING ###
  ## The high level of theory
  clusters_df = pd.read_pickle(TEST_HIGH).sort_values([('info','file_basename')])
  try:
    ens = (clusters_df["log"]["electronic_energy"]).values.astype("float")
    Qeval = 2 #2=Compare ML prediction with QC
  except:
    Qeval = 1 #1=Only predicts
  strs = clusters_df["xyz"]["structure"]
  ## The low level of theory 
  if method == "delta":
    if Qeval == 1:
      clusters_df2 = clusters_df
    else:
      clusters_df2 = pd.read_pickle(TEST_LOW).sort_values([('info','file_basename')])
    ens2 = (clusters_df2["log"]["electronic_energy"]).values.astype("float")
    #str2 should be the same as str by princip

  ### REPRESENTATION CALCULATION ###
  repres_dataframe = pd.DataFrame(index = strs.index, columns = ["xyz"])
  max_atoms = max([len(strs[i].get_atomic_numbers()) for i in range(len(strs))])
  for i in range(len(repres_dataframe)):
    repres_dataframe["xyz"][i] = generate_representation(strs[i].get_positions(), strs[i].get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=10.0)
  fchl_representations = np.array([mol for mol in repres_dataframe["xyz"]])
 
  #some info about the full representation
  print(fchl_representations.shape, flush = True)
  if Qeval == 2:
    print(ens.shape, flush = True)
  if method == "delta":
    print(ens2.shape, flush = True)
    if Qeval == 2:
      if ens.shape != ens2.shape:
        print("The LOW and HIGH method test sizes do not match. [EXITING]", flush = True)
        exit()
 
  ### BINDING ENERGIES CALCULATION ###
  #TODO clustername given as argument
  #monomers_df = pd.read_pickle("/home/kubeckaj/ML_SA_B/DATABASES/database_XTB_monomers_at_DFT.pkl")
  #import re
  #clustername=np.array(re.split(r'(\d+)', sys.argv[3])[1:])
  #clustername=clustername.reshape(int(len(clustername)/2),2)
  #clustername=[clustername[:,(1,0)]]*len(clusters_df)
  #print(clustername)
  #print(clustername[0][0])
  #for i in clustername: #np.transpose([clusters_df["info"]["components"].values,clusters_df["info"]["component_ratio"].values]):
  if Qmonomers == 2: #no difference from monomers at all
    ens_correction = [0]*len(clusters_df)
    if method == "delta":
      ens2_correction = [0]*len(ens2)
  else:
    #HIGH LEVEL
    if Qmonomers == 1:
      clusters_df0 = pd.read_pickle(MONOMERS_HIGH)
      ins_monomers = [np.sum(clusters_df0["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df0))]
      monomers_df = clusters_df0[ins_monomers]
    else:
      ins_monomers = [np.sum(clusters_df["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df))]
      monomers_df = clusters_df[ins_monomers]
    ens_correction = [0]*len(strs)
    n = 0
    for i in np.transpose([clusters_df["info"]["components"].values,clusters_df["info"]["component_ratio"].values]):
      nn = 0
      for j in i[0]:
        test = 0
        for k in range(len(monomers_df)):
          monomer_k = monomers_df.iloc[k]
          monomer_k_name = monomer_k["info"]["components"]
          if j == monomer_k_name[0]:
            ens_correction[n] += monomer_k["log"]["electronic_energy"]*i[1][nn]
            test = 1
            break
        if test == 0:
          print("OMG; monomer "+j+" was not found in:", flush = True)
          print(monomers_df["info"]["components"].values, flush = True)
          exit()
        nn += 1
      n += 1
    #LOW LEVEL
    if method == "delta":
      if Qmonomers == 1:
        clusters_df0 = pd.read_pickle(MONOMERS_LOW)
        ins_monomers = [np.sum(clusters_df0["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df0))]
        monomers_df = clusters_df0[ins_monomers]
      else:
        ins_monomers = [np.sum(clusters_df2["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df2))]
        monomers_df = clusters_df2[ins_monomers]
      ens2_correction = [0]*len(ens2)
      n = 0
      for i in np.transpose([clusters_df2["info"]["components"].values,clusters_df2["info"]["component_ratio"].values]):
        nn = 0
        for j in i[0]:
          test = 0
          for k in range(len(monomers_df)):
            monomer_k = monomers_df.iloc[k]
            monomer_k_name = monomer_k["info"]["components"]
            if j == monomer_k_name[0]:
              ens2_correction[n] += monomer_k["log"]["electronic_energy"]*i[1][nn]
              test = 1
              break
          if test == 0:
            print("OMG; monomer "+j+" was not found in:", flush = True)
            print(monomers_df["info"]["components"].values, flush = True)
            exit()
          nn += 1
        n += 1
      
  #The binding (formation) energy calculation
  if Qeval == 2:
    form_ens = ens - ens_correction
    #print(form_ens, flush = True)
  if method == "delta":
    form_ens2 = ens2 - ens2_correction
    #print(form_ens2, flush = True) 

  ### DEFINING THE EVALUATION Xs:  Y = QML(X) 
  #the full set
  X_test0 = fchl_representations
  #TODO Do we take only subset for the testing?
  if size != "full":
    X_test, X_trash, idx, idx_trash = train_test_split(X_test0, range(len(X_test0)), test_size=(len(X_test0)-size)/len(X_test0))
    clusters_df = clusters_df[idx]
    ens_correction = ens_correction[idx]
    if Qeval == 2:
      form_ens = form_ens[idx]
    if method == "delta":
      form_ens2 = form_ens2[idx]
  else:
    X_test = X_test0

  ### CORRECTING THE FCHL MATRIX SIZES
  #X_train = fchl_representations0
  #IF YOU ARE EXTENDING THIS WILL MAKE THE MATRIXES OF THE SAME SIZE
  if X_train.shape[1] != X_test.shape[1]:
    if X_train.shape[1] > X_test.shape[1]:
      small = X_test
      large = X_train
    else:
      small = X_train
      large = X_test
    newmatrix = np.zeros([small.shape[0],large.shape[1],5,large.shape[3]])
    newmatrix[:,:,0,:] = 1e+100
    newmatrix[0:small.shape[0],0:small.shape[1],0:5,0:small.shape[3]] = small
    if X_train.shape[1] > X_test.shape[1]:
      X_test = newmatrix
    else:
      X_train = newmatrix
 
  ### THE EVALUATION
  Ks = get_local_kernels(X_test, X_train, sigmas, **kernel_args)
  Y_predicted = [np.dot(Ks[i], alpha[i]) for i in range(len(sigmas))]

  ### PRINTING THE RESULTS
  #print("Ypredicted = {" + ",".join([str(i) for i in Y_predicted[0]])+"};")
  if Qeval == 2:
    print("Ytest = {" + ",".join([str(i) for i in ens])+"};", flush = True)
  if method != "delta":
    print("Ypredicted = {" + ",".join([str(i) for i in Y_predicted[0]+ens_correction])+"};", flush = True)
  else:
    #print("TEST",flush=True)
    #print(Y_predicted[0],flush=True)
    #print("TEST",flush=True)
    #print(form_ens2,flush=True)
    #print("TEST",flush=True)
    #print(ens_correction,flush=True)
    #print("TEST",flush=True)
    print("Ypredicted = {" + ",".join([str(i) for i in Y_predicted[0]+form_ens2+ens_correction])+"};", flush = True)
  #print("formensval = {" + ",".join([str(i) for i in form_ens[idx]])+"};")
  #print("ens_correction = {" + ",".join([str(i) for i in ens_correction])+"};")
 
  # If possible print MAE and RMSE
  if Qeval == 2:
    if method != "delta":
      Y_validation = form_ens
    else:
      Y_validation = form_ens - form_ens2
    # Calculate mean-absolute-error (MAE):
    mae = [627.503*np.mean(np.abs(Y_predicted[i] - Y_validation))  for i in range(len(sigmas))]
    print("mae = " + ",".join([str(i) for i in mae])+" kcal/mol", flush = True)
    # Calculate root-mean-squared-error (RMSE):
    rmse = [627.503*np.sqrt(np.mean(np.abs(Y_predicted[i] - Y_validation)**2))  for i in range(len(sigmas))]
    print("rmse = " + ",".join([str(i) for i in rmse])+" kcal/mol", flush = True)

  ### PRINTING THE QML PICKLES
  clustersout_df = clusters_df.copy()
  for i in range(len(clustersout_df)):
    if method != "delta":
      clustersout_df.loc[clustersout_df.iloc[i].name,("log","electronic_energy")] = Y_predicted[0][i]+ens_correction[i]
    else:
      clustersout_df.loc[clustersout_df.iloc[i].name,("log","electronic_energy")] = Y_predicted[0][i]+form_ens2[i]+ens_correction[i]
  clustersout_df.to_pickle(outfile)
########

######################
###### OPTIMIZE ######
######################

if Qopt == 1:
  ### DATABASE LOADING ###
  ## The high level of theory
  clusters_df = pd.read_pickle(TEST_HIGH).sort_values([('info','file_basename')])
  strs = clusters_df["xyz"]["structure"]
 
  print("Preparing optimization", flush = True) 
  xyz=strs[0].get_positions()
  #F=np.zeros_like(xyz)
  #for i in range(len(xyz)):
  #      for j in range(3):
  #            F[i,j]=float(13.0)
  #print(F[0,0], flush = True)
  #print(xyz[0,0]+0.23, flush = True)
  #print(type(xyz), flush = True)
  maxsteps=8
  xyzdeviation=0.05
  shift=0.3  
  print("Starting optimization", flush = True)
  for step in range(maxsteps):
    ### GENERATE SHIFTED STRUCTURES
    
    R=[xyz]
    if step != maxsteps-1:
      for i in range(len(xyz)):
        for j in range(3):
          ch=np.zeros_like(xyz)
          ch[i,j]=+xyzdeviation #THIS IS THE SHIFT OF 0.05 Angstrom
          R.append(xyz+ch)
          #ch=-ch
          #R.append(xyz+ch)
    
    RR=pd.DataFrame(np.zeros(len(R)),index=range(len(R)))
    RR[0]=R
    #print(RR, flush = True)  

    if method == "delta":
      for RR_iter in range(len(RR[0])): 
        #print(RR_iter,flush = True)
        os.system("mkdir test;")
        tocalc = strs[0].copy()
        tocalc.set_positions(RR[0][RR_iter])
        write("test/test.xyz",tocalc)
        os.system("cd test;xtb test.xyz --sp --gfn 1 > test.log 2>&1 ;cd ..;JKQC -folder test -out JKMLtest.pkl -noex;rm -r test")
        os.system("JKQC JKMLtest.pkl -el > .en")
        with open(".en", "r") as ping:
          en=float(ping.read().rstrip())
        #print(en, flush=True)
        if RR_iter == 0:
          all_ens = [en]
        else:
          all_ens.append(en)
        #print(all_ens, flush = True)
      #print(all_ens, flush = True)

    #print(RR.values[0][0], flush = True)
    ### REPRESENTATION CALCULATION ###
    repres_dataframe = pd.DataFrame(index = RR.index, columns = ["xyz"])
    max_atoms = max([len(strs[i].get_atomic_numbers()) for i in range(len(strs))])
    for i in range(len(repres_dataframe)):#TODO strs[0] cannot be define like that for different molecules, i.e. I can optimize only 1 molecule
      repres_dataframe["xyz"][i] = generate_representation(RR.values[i][0], strs[0].get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=10.0)
    fchl_representations = np.array([mol for mol in repres_dataframe["xyz"]])
  
    #some info about the full representation
    #print(fchl_representations.shape, flush = True)
  
    ### DEFINING THE EVALUATION Xs:  Y = QML(X) 
    #the full set
    X_test = fchl_representations
  
    ### CORRECTING THE FCHL MATRIX SIZES
    #X_train = fchl_representations0
    #IF YOU ARE EXTENDING THIS WILL MAKE THE MATRIXES OF THE SAME SIZE
    if X_train.shape[1] != X_test.shape[1]:
      if X_train.shape[1] > X_test.shape[1]:
        small = X_test
        large = X_train
      else:
        small = X_train
        large = X_test
      newmatrix = np.zeros([small.shape[0],large.shape[1],5,large.shape[3]])
      newmatrix[:,:,0,:] = 1e+100
      newmatrix[0:small.shape[0],0:small.shape[1],0:5,0:small.shape[3]] = small
      if X_train.shape[1] > X_test.shape[1]:
        X_test = newmatrix
      else:
        X_train = newmatrix
  
    ### THE EVALUATION
    Ks = get_local_kernels(X_test, X_train, sigmas, **kernel_args)
    Y_predicted = [np.dot(Ks[i], alpha[i]) for i in range(len(sigmas))]

    if method == "delta":
      new_save_energy = Y_predicted[0][0]+all_ens[0]
    else:
      new_save_energy = Y_predicted[0][0]
    if step != 0:
      if new_save_energy > save_energy:     
        xyz = xyz + change
        shift = shift/2
        continue

    xyzold = xyz
    F=np.zeros_like(xyz)
    if step != maxsteps-1:
      Fp=np.zeros_like(xyz)
      #Fm=np.zeros_like(xyz)
      for i in range(len(xyz)):
        for j in range(3):
          if method == "delta":
            #print(Y_predicted[0][0])
            #print(all_ens[0])
            #print(wtf1)
            #print(F[i,j])
            F[i,j]=Y_predicted[0][0]+all_ens[0]
            #print(F[i,j])
            Fp[i,j]=Y_predicted[0][1+j+3*i]+all_ens[1+j+3*i]
          else:
            F[i,j]=Y_predicted[0][0]
            Fp[i,j]=Y_predicted[0][1+j+3*i]
          #Fp[i,j]=Y_predicted[0][1+2*j+6*i]
          #Fm[i,j]=Y_predicted[0][2+2*j+6*i]
  
      
      #print(np.linalg.inv((Fm+Fp-2*F)/(2*xyzdeviation)))
      #print(xyz+0.5*np.matmul(np.linalg.inv((Fm+Fp-2*F)/(2*xyzdeviation)),(Fp-Fm)/(2*xyzdeviation)), flush = True)
      change=(Fp-F)/xyzdeviation*shift
      xyz = xyz - change
      #print("TEST:",flush = True)
      #print(change, flush = True)
      #print(change*change, flush = True)
      #print(np.transpose(change*change),flush=True)
      #print(sum(np.transpose(change*change)),flush=True)
      #print(np.sqrt(sum(np.transpose(change*change))),flush = True)
      maxdev=max(np.sqrt(sum(np.transpose(change*change))))

    if step == 0:
      print("step \tenergy [Eh]        \tmax.shift [A]", flush = True)
      clustersout_df = pd.DataFrame()
    save_energy = new_save_energy
    cluster_id = len(clustersout_df)
    clustersout_df = df_add_append(clustersout_df, "info", "folder_path", [str(cluster_id)], os.path.abspath(TEST_HIGH)[::-1].split("/",1)[1][::-1]+"/")
    clustersout_df = df_add_iter(clustersout_df, "log", "electronic_energy", [str(cluster_id)], [save_energy])
    newxyz = strs[0].copy()
    newxyz.set_positions(xyzold)
    clustersout_df = clustersout_df = df_add_iter(clustersout_df, "xyz", "structure", [str(cluster_id)], [newxyz])
    print(str(step)+" \t"+str(save_energy)+" \t"+str(maxdev),flush = True)
    #print(xyz, flush = True)
    

  ### PRINTING THE RESULTS
  #print("Ypredicted = {" + ",".join([str(i) for i in save_energy])+"};", flush = True)

  ### PRINTING THE QML PICKLES
  #clustersout_df = clusters_df.copy()
  #for i in range(len(clustersout_df)):
  #  clustersout_df.loc[clustersout_df.iloc[i].name,("log","electronic_energy")] = Y_predicted[0][i]
  clustersout_df.to_pickle(outfile)
########


