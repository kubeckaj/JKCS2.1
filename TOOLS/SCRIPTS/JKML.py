import sys
import os.path
import subprocess

def help():
  print("#################################################################", flush = True)
  word = "JKML"
  colors = ["\033[34m", "\033[31m", "\033[32m", "\033[33m"]
  # 31red, 32green, 33yellow, 34blue, 35magenta, 36cyan
  def pJKML():
    for i, letter in enumerate(word): 
      color = colors[i % len(colors)]  # cycle through colors
      print(f"{color}{letter}\033[0m", end="", flush = True)
  pJKML();print(" HELP:", flush = True)
  print("This script interfaces ML model and pickled structures.", flush = True)
  print("", flush = True)
  print(f"  ",end="");pJKML();print(" [OPTION(s)]", flush = True)
  print("", flush = True)
  print("  HELP OPTIONS:", flush = True)
  print("    -help                            prints basic help", flush = True)
  print("    -help_nn                         prints help for neural network methods (e.g. PaiNN,SchNet)", flush = True)
  print("    -help_krr                        prints help for kernel ridge regression methods (e.g. FCHL)", flush = True)
  print("    -help_adv                        print some advanced features")
  print("", flush = True)
  print("  OPTIONS:", flush = True)
  print("    -qml                             use KRR with FCHL19 [set by default]", flush = True)
  print("    -mbdf                            use KRR with MBDF", flush = True)
  print("    -nn,-painn                       switch to NN = neural network with PaiNN", flush = True)
  print("    -schnet                          switch to NN = neural network with SchNet", flush = True)
  print("", flush = True)
  print("  INPUT FILES:", flush = True)
  print("    -train <HIGH.pkl> [<LOW.pkl>]    train on given pikled files", flush = True)
  print("    -trained <model.pkl>             take pre-trained ML model", flush = True)
  print("    -monomers <HIGH.pkl> [<LOW.pkl>] properties with respect to monomer(s) in pickled file(s)", flush = True)
  print("    -test <HIGH.pkl> [<LOW.pkl>]     validate ML on energies of structures in pickled file(s)", flush = True)
  print("    -eval <STRS.pkl> [<LOW.pkl>]     evaluate energies of NEW structures in pickled file(s)", flush = True)
  print("", flush = True)
  print("  OUTPUT FILES:", flush = True)
  print("    -out <file.pkl>                  name of file with predicted values [def = predicted.pkl]", flush = True)
  print("    -modelout <file.pkl>             name of file with saved model [def = model.pkl]", flush = True)
  print("", flush = True)
  print("  EXAMPLES:", flush = True)
  print(f"    ",end="");pJKML();print(" -loc -train collected_structures.pkl -modelout trained_example.pkl", flush = True)
  print(f"    ",end="");pJKML();print(" -par q64 -cpu 64 -train tr_high.pkl tr_low.pkl -test te_high.pkl te_low.pkl", flush = True)
  print(f"    ",end="");pJKML();print(" -par qgpu -cpu 2 -nn -epochs 10000 -train x.pkl -eval str.pkl -monomers monomers.pkl", flush = True)
  print(f"    ",end="");pJKML();print(" -loc -trained model_for_atomization_en.pkl -eval diff_molecules.pkl -monomers atoms.pkl ", flush = True)
  print("", flush = True)

  #SO FAR I WANT THE MONOMERS SEPARATELY: print("    /or/  none                     training directly on el.energies (not good for mix of clusters)", flush = True)
  #METHOD IS NOT NECESSARY: print("-method <str>                      direct|delta [default]", flush = True)

def help_adv():
  print("  ADVANCED OPTIONS:", flush = True)
  print("    -column <str> <str> selects a different column from the pickled files (e.g. log entropy)", flush = True)
  print("    -noforces           forces are not trained/tested even if it is possible")
  print("    -size <int>         randomly selects only portion of the trainin database (e.g. 200)", flush = True)
  print("    -seed <int>         seed for random number generators [def = 42]")
  print("    -categorize <int>   selects structures with similar MBTR (our similarity definition)", flush = True)
  print("    -similarity <int>   selects structures with similar FCHL (uses kernel)", flush = True)
  print("    -forcemonomers      adds (extra) monomers to sampleeach/selection", flush = True)
  print("    -printforces        print out all forces (this might be a lot of numbers)", flush = True)
  print("    -opt <STRS.pkl>     optimize structure based on model [NN]", flush = True)
  print("    -opt_maxs <float>   max step in Angstrom in optimization [def = 0.02]")
  print("    -md <STRS.pkl>      run md starting from provided structure(s) based on model [NN]", flush = True)
  print("    -md_temperature     temperature of the simulation [def = 300.0]", flush = True)
  print("", flush = True)
  print("  EXTRA ADVANCED OPTIONS:", flush = True)
  print("    -so3net             switch to NN = neural network with SO3net (from SchNetPack)", flush = True)
  print("    -split <int>        only with -krr/-fchl how many splits of KRR matrix do you do", flush = True)
  print("    -startsplit <int>   the same like above but only construct kernels", flush = True)
  print("    -finishsplit <int>  (see -split) combines the splitted kernel and creates model", flush = True)
  print("    -wolfram            prints {} instead of []", flush = True)
  print("", flush = True)

def help_krr():
  print("  OPTIONS FOR KERNEL RIDGE REGRESSION:", flush = True)
  print("    -sigma <int>        Gaussian width hyperparameter [def = 1.0]", flush = True)
  print("    -lambda <int>       numerical stability (for matrix inv.) hyperparameter [def = 1e-4]")
  print("    -laplacian          switch to Laplacian kernel (FCHL)")
  print("    -krr_cutoff <float> cutoff function (Angstrom) [def = 10.0]", flush = True)
  print("", flush = True)
  
def help_nn():
  print("  OPTIONS FOR NEURAL NETWORKS:", flush = True)
  print("    -epochs <int>              number of epochs [def = 1000]", flush = True)
  print("    -nn_train <float>          portion of training data (exlc. validation) [def = 0.9]", flush = True)
  print("    -nn_ESpatience <int>       Early-Stop patience of epochs for no improvement [def = 200]", flush = True)
  print("    -nn_energytradeoff <float> trade-off [energy, force] = [<float>, 1] [def = 0.01]")
  print("    -nn_lr                     learning rate [def = 1e-4]", flush = True)
  print("    -nw                        number of workers for database manipulation [def = 1]")
  print("", flush = True)
  print("  OPTIONS FOR REPRESENTATION:", flush = True)
  print("    -nn_ab <int>       number of atom basis/features/size of embeding vector [def = 256]", flush = True)
  print("    -nn_int <int>      number of interaction blocks [def = 5]", flush = True)
  print("    -nn_rb <int>       number of radial basis for exp. int. dist. [def = 20]", flush = True)
  print("    -nn_cutoff <float> cutoff function (Angstrom) [def = 5.0]", flush = True)
  print("", flush = True)

#PREDEFINED ARGUMENTS
method = "direct"
size = "full"
seed = 42
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
Qwolfram = 0 #1 means that {} will be printed insead of []
Qprintforces = 0 #print forces?
Qmonomers = 2 #0=monomers taken from database, 1=monomers in separate files, 2=no monomer subtraction
Qsampleeach = 0
Qforcemonomers = 0
column_name_1 = "log"
column_name_2 = "electronic_energy" 
Qifforces = 1 #IF forces exist use them in calculations

#krr : fchl
#nn : painn
#nn : 
Qmethod = "krr"
Qrepresentation = "fchl"

#PREDEFINED QML
Qkernel = "Gaussian"
sigmas = [1.0]
lambdas = [1e-4]*len(sigmas)
krr_cutoff = 10.0
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

#PREDEFINED NN - PaiNN
nn_rbf = 20
nn_tvv = 0.9
nn_cutoff = 5.0
nn_atom_basis = 256
nn_interactions = 5
nn_epochs = 1000
Qlearningrate = 1e-4
Qearlystop = 200
Qenergytradoff = 0.01 #if forces are trained on: [energy, force] = [X, 1]
nw = 1

#OPT/MD
md_temperature = 300.0
opt_maxstep = 0.02

#OUTPUT FILES
outfile="predicted.pkl"
varsoutfile="model.pkl"

#TREATING ARGUMENTS
last=""
for i in sys.argv[1:]:
  #HELP
  if i == "-help" or i == "-help_nn" or i == "-help_adv" or i == "-help_krr":
    help()
    if i == "-help_nn":
      help_nn()
    if i == "-help_adv":
      help_adv()
    if i == "-help_krr":
      help_krr()
    print("#################################################################", flush = True)
    exit()
  #VARS OUT FILE
  if i == "-varsout" or i == "-modelout":
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
    elif i == "delta":
      method = "delta"
    else:
      print("I do not understand your input:")
      print(i)
      exit()
    continue
  #FORCE MONOMERS
  if i == "-forcemonomers":
    Qforcemonomers = 1
    continue
  #SEED
  if i == "-seed":
    last = "-seed"
    continue
  if last == "-seed":
    last = ""
    seed = int(i)
    continue
  #WOLFRAM
  if i == "-wolfram":
    Qwolfram = 1
    continue
  #PRINT FORCES
  if i == "-printforces":
    Qprintforces = 1
    continue
  #COLUMN
  if i == "-column":
    last = "-column"
    continue
  if last == "-column":
    column_name_1 = i
    last = "-column2"
    continue  
  if last == "-column2":
    column_name_2 = i
    last = ""
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
  if i == "-opt" or i == "-optimize":
    Qopt = 1
    last = "-opt"
    continue
  if i == "-md":
    Qopt = 2
    last = "-opt"
    continue
  if i == "-monomers" or i == "-mon":
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
  #SAMPLEEACH 
  if i == "-sampleeach" or i == "-se" or i == "-categorize":
    last = "-sampleeach"
    continue
  if last == "-sampleeach":
    last = ""
    if Qtrain == 0:
      Qtrain = 1
    if Qeval == 0:
      Qeval = 1
    Qsampleeach = int(i)
    continue
  #NO FORCES
  if i == "-noforces":
    Qifforces = 0
    continue
  #SIMILARITY
  if i == "-similarity" or i == "-sim":
    last = "-similarity"
    continue
  if last == "-similarity":
    last = ""
    if Qtrain == 0:
      Qtrain = 1
    if Qeval == 0:
      Qeval = 1
    Qsampleeach = -int(i)
    continue
  #LAPLACIAN
  if i == "-laplacian":
    Qkernel = "Laplacian"
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
    if os.path.exists(i):
      method = "delta"
      TRAIN_LOW = i
      last = ""
      continue
  #TEST/EVAL/OPT DATABASE(S)
  if last == "-eval":
    TEST_HIGH = i
    Qeval = 1
    last = "-test2"
    continue
  if last == "-test":
    TEST_HIGH = i
    Qeval = 2
    last = "-test2"
    continue
  if last == "-opt":
    TEST_HIGH = i
    last = "-test2"
    continue
  if last == "-test2":
    if os.path.exists(i):
      method = "delta"
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
    if os.path.exists(i):
      method = "delta"
      MONOMERS_LOW = i
      last = ""
      continue
  #MODELS AND REPRESENTATIONS:
  if i == "-painn" or i == "-nn":
    Qmethod = "nn"
    Qrepresentation = "painn"
    continue
  if i == "-schnet":
    Qmethod = "nn"
    Qrepresentation = "schnet"
    continue
  if i == "-so3net":
    Qmethod = "nn"
    Qrepresentation = "so3net"
    continue
  #EPOCHS
  if i == "-nn_epochs" or i == "-epochs":
    last = "-nn_epochs"
    continue
  if last == "-nn_epochs":
    last = ""
    nn_epochs = int(i)
    continue
  #EPOCHS
  if i == "-nn_tvv" or i == "-nn_train":
    last = "-nn_tvv"
    continue
  if last == "-nn_tvv":
    last = ""
    nn_tvv = float(i)
    continue
  #RADIAL BASIS
  if i == "-nn_rbf" or i == "-nn_rb":
    last = "-nn_rbf"
    continue
  if last == "-nn_rbf":
    last = ""
    nn_rbf = int(i)
    continue
  #NN CUTOFF
  if i == "-nn_cutoff" or i == "-krr_cutoff" or i == "-cutoff":
    last = "-cutoff"
    continue
  if last == "-cutoff":
    last = ""
    nn_cutoff = float(i)
    krr_cutoff = float(i)
    continue
  #OPT MAXSTEP
  if i == "-opt_maxstep" or i == "-opt_maxs":
    last = "-opt_maxstep"
    continue
  if last == "-opt_maxstep":
    last = ""
    opt_maxstep = float(i)
    continue
  #ATOM BASIS
  if i == "-nn_ab" or i == "-nn_atom_basis":
    last = "-nn_ab"
    continue
  if last == "-nn_ab":
    last = ""
    nn_atom_basis = int(i)
    continue
  #NN INTERACTIONS
  if i == "-nn_int" or i == "-nn_interctions":
    last = "-nn_int"
    continue
  if last == "-nn_int":
    last = ""
    nn_interactions = int(i)
    continue
  #NN INTERACTIONS
  if i == "-nw":
    last = "-nw"
    continue
  if last == "-nw":
    last = ""
    nw = int(i)
    continue
  #EPOCHS
  if i == "-nn_epochs" or i == "-epochs":
    last = "-nn_epochs"
    continue
  if last == "-nn_epochs":
    last = ""
    nn_epochs = int(i)
    continue
  #ENERGY TRADEOFF
  if i == "-nn_energytradeoff":
    last = "-nn_energytradeoff"
    continue
  if last == "-nn_energytradeoff":
    last = ""
    nn_energytradeoff = float(i)
    continue
  #EARLY STOP
  if i == "-nn_ESpatience" or i == "-nn_espatience":
    last = "-nn_espatience"
    continue
  if last == "-nn_espatience":
    last = ""
    Qearlystop = int(i)
    continue
  #LEARNING RATE
  if i == "-nn_lr":
    last = "-nn_lr"
    continue
  if last == "-nn_lr":
    last = ""
    Qlearningrate = float(i)
    continue
  #MD temperature
  if i == "-md_temperature" or i == "-temperature":
    last = "-md_temperature"
    continue
  if last == "-md_temperature":
    last = ""
    md_temperature = float(i)
    continue
  #KRR (by default)
  if i == "-krr" or i == "-fchl" or i == "-qml":
    Qmethod = "krr"
    Qrepresentation = "fchl"
    continue
  if i == "-mbdf":
    Qmethod = "krr"
    Qrepresentation = "mbdf"
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
      if Qeval == 1:
        TEST_LOW = TEST_HIGH
      else:
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
import numpy as np
from ase import Atoms
from ase.io import read, write
from sklearn.model_selection import train_test_split
import pickle
import pandas as pd

#LOADING THE DATABASES
train_high_database = "none"
train_low_database = "none"
monomers_high_database = "none"
monomers_low_database = "none"
if Qtrain == 1:
  train_high_database = pd.read_pickle(TRAIN_HIGH).sort_values([('info','file_basename')])
  if method == "delta":
    train_low_database = pd.read_pickle(TRAIN_LOW).sort_values([('info','file_basename')])
if Qeval == 1 or Qeval == 2 or Qopt > 0:
  test_high_database = pd.read_pickle(TEST_HIGH).sort_values([('info','file_basename')])
  if method == "delta":
    test_low_database = pd.read_pickle(TEST_LOW).sort_values([('info','file_basename')])  
if Qmonomers == 1:
  monomers_high_database = pd.read_pickle(MONOMERS_HIGH).sort_values([('info','file_basename')])
  if method == "delta":
    monomers_low_database = pd.read_pickle(MONOMERS_LOW).sort_values([('info','file_basename')])

#LIBRARIES FOR GIVEN METHOD AND REPRESENTATION
if Qmethod == "krr":
  from qml.math import cho_solve
  if Qrepresentation == "fchl":
    from qml.fchl import generate_representation
    if Qkernel == "Gaussian":
      from qml.fchl import get_local_symmetric_kernels
      from qml.fchl import get_local_kernels
      JKML_sym_kernel = get_local_symmetric_kernels
      JKML_kernel = get_local_kernels
    else:
      from qml.kernels import laplacian_kernel
      from qml.kernels import laplacian_kernel_symmetric
      JKML_sym_kernel = laplacian_kernel_symmetric
      JKML_kernel = laplacian_kernel
  elif Qrepresentation == "mbdf":
    from MBDF import generate_mbdf
    generate_representation = generate_mbdf
    from qml.kernels import get_local_symmetric_kernel_mbdf, get_local_kernel_mbdf
    JKML_sym_kernel = get_local_symmetric_kernel_mbdf
    JKML_kernel = get_local_kernel_mbdf
  #from qml.fchl import get_atomic_symmetric_kernels
  #from qml.fchl import get_atomic_kernels
  #from qml.fchl import get_global_symmetric_kernels
  #from qml.fchl import get_global_kernels
  #from qml.kernels import gaussian_kernel
elif Qmethod == "nn" and Qrepresentation == "painn":
  if Qtrain > 0:
    from schnetpack.data import ASEAtomsData
    from schnetpack.data import AtomsDataModule
    import logging
    import schnetpack.transform as trn
    import torchmetrics
    import pytorch_lightning as pl
  if Qeval > 0 or Qopt > 0 or Qopt > 0:
    from ase.units import Ha, Bohr
    from schnetpack.interfaces import SpkCalculator
  import torch
  import schnetpack as spk
else:
  print("Wrong method or representation chosen.")
  exit()

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

def substract_monomers(the_clusters_df,the_monomers_df):
  if Qmonomers == 2: #no difference from monomers at all
    the_ens_correction = [0]*len(the_clusters_df)
  else:
    if Qmonomers == 1:
      clusters_df0 = the_monomers_df
      ins_monomers = [np.sum(clusters_df0["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df0))]
      monomers_df = clusters_df0[ins_monomers]
    else:
      ins_monomers = [np.sum(the_clusters_df["info"]["component_ratio"].values[i]) == 1 for i in range(len(the_clusters_df))]
      monomers_df = the_clusters_df[ins_monomers]
    the_ens_correction = [0]*len(the_clusters_df)
    n = 0
    for i in np.transpose([the_clusters_df["info"]["components"].values,the_clusters_df["info"]["component_ratio"].values]):
      nn = 0
      for j in i[0]:
        test = 0
        for k in range(len(monomers_df)):
          monomer_k = monomers_df.iloc[k]
          monomer_k_name = monomer_k["info"]["components"]
          if j == monomer_k_name[0]:
            the_ens_correction[n] += monomer_k[column_name_1][column_name_2]*i[1][nn]
            test = 1
            break
        if test == 0:
          print("OMG; monomer "+j+" was not found in:")
          print(monomers_df["info"]["components"].values)
          exit()
        nn += 1
      n += 1
  return the_ens_correction

###################################
###### SAMPLEEACH/SIMILARITY ######
###################################

if Qtrain == 0:
  print("training option is missing [EXITING]", flush = True)
  exit()

# SAMPLE EACH SECTION:
if Qsampleeach > 0:
  from dscribe.descriptors import MBTR
  def flatten(l):
      return [item for sublist in l for item in sublist]
  chemsyms_uniques = list(set(flatten([i.get_chemical_symbols() for i in train_high_database["xyz"]["structure"]])))
  #max_atoms=max([len(i.get_chemical_symbols()) for i in pd.read_pickle(TRAIN_HIGH).sort_values([('info','file_basename')])["xyz"]["structure"]])
  k2min, k2max, k2n = 0.7, 2.0, 100
  mbtr = MBTR(
        species=chemsyms_uniques,
        #ATOMIC NUMBER
        #k1={
        #    "geometry": {"function": "atomic_number"},
        #    "grid": {"min": 0, "max": 16, "n": 100, "sigma": 0.1},
        #},
        #DISTANCE OR INVERSE DISTANCE
        k2={
            "geometry": {"function": "distance"},
            "grid": {"min": k2min, "max": k2max, "n": k2n, "sigma": 0.000000001},
            "weighting": {"function": "exp", "scale": 0.5, "threshold": 3e-3},
            #"weighting": {"function": "unity"},#, "scale": 0.5, "threshold": 3e-3},
        },
        #ANGLE OR COSINE(ANGLE)
        #TODO to be tested as well
        k3={
            "geometry": {"function": "cosine"},
            "grid": {"min": -1, "max": 1, "n": 100, "sigma": 0.000000001},
            "weighting": {"function": "exp", "scale": 0.5, "threshold": 3e-3},
            #"weighting": {"function": "unity"},# "scale": 0.5, "threshold": 3e-3},
        },
        periodic=False,
        normalization="l2_each",
        flatten=False  # without this it's gonna be a dict and then should find out which part corresponds to which stuff
  )

  def compare_mbtr(structa, structb):
    '''
    Expects two xyz pickle structures
    Returns a float equal the sum of the average mean square deviation
    between each bond type for the twombtr inputs.
    Is 0 if equal, lower values equal similiarty.
    '''
    return np.sum(np.square(structa-structb))

  def cr_mbtr_k2(struct):
    '''
    Expects xyz pickle structure
    Returns the mbtr for bond lenghts. 
    '''
    return mbtr.create(struct)["k2"][:][:]

  print("MBTR must be calculated (just once) for both train and test datasets", flush = True)
  from joblib import Parallel, delayed
  import multiprocessing
  
  def task_train(arg):
    #mbtr_train[arg] = cr_mbtr_k2(train_database[arg])
    #return cr_mbtr_k2(train_database[arg])
    return cr_mbtr_k2(arg)
  
  def task_test(arg):
    #mbtr_test[arg] = cr_mbtr_k2(test_database[arg])
    #return cr_mbtr_k2(test_database[arg])
    return cr_mbtr_k2(arg)
  
  num_cores = multiprocessing.cpu_count()
  print("Trying to use "+str(num_cores)+" CPUs for MBTR. (If less are available I hope that nothing gets fucked up.)")
  
  mbtr_train = Parallel(n_jobs=num_cores)(delayed(task_train)(i) for i in train_high_database["xyz"]["structure"])
  mbtr_test = Parallel(n_jobs=num_cores)(delayed(task_train)(i) for i in test_high_database["xyz"]["structure"])
  print("MBTR done", flush = True)
  sampleeach_all = range(len(mbtr_test))
elif Qsampleeach < 0: 
  #SAMLEEACH = SIMILARITY based on FCHL 
  from joblib import Parallel, delayed
  import multiprocessing
  def task(arg):
    return generate_representation(arg.get_positions(),arg.get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=krr_cutoff)

  num_cores = multiprocessing.cpu_count()
  print("Trying to use "+str(num_cores)+" CPUs for FCHL. (If less are available I hope that nothing gets fucked up.)")
  max_atoms = max([len(i.get_atomic_numbers()) for i in train_high_database["xyz"]["structure"]])
  max_atoms2 = max([len(i.get_atomic_numbers()) for i in test_high_database["xyz"]["structure"]])
  if max_atoms2 > max_atoms:
    max_atoms = max_atoms2
  fchl_train = Parallel(n_jobs=num_cores)(delayed(task)(i) for i in train_high_database["xyz"]["structure"])
  fchl_test = Parallel(n_jobs=num_cores)(delayed(task)(i) for i in test_high_database["xyz"]["structure"])
  sampleeach_all = range(len(fchl_test))
  print("FCHL done", flush = True)
else:
  sampleeach_all = ["once"]


########################################################################################################################
########################################################################################################################
########################################################################################################################
### looping over all test structures / or / doing the whole process process for all

for sampleeach_i in sampleeach_all:
  ### SIMILARITY AND SAMPLE EACH
  if Qsampleeach > 0:
    dist = np.array([compare_mbtr(mbtr_train[i],mbtr_test[sampleeach_i]) for i in range(np.shape(mbtr_train)[0])])
    sampledist = dist.argsort()[:Qsampleeach] 
  elif Qsampleeach < 0:
    simil = JKML_kernel(np.array([m for m in [fchl_test[sampleeach_i]]]), np.array([m for m in fchl_train]), [0.001], **kernel_args)[0][0]
    simil = [ simil[i]/np.sqrt(JKML_kernel(np.array([m for m in [fchl_train[i]]]),np.array([m for m in [fchl_train[i]]]))[0][0][0]) for i in range(len(train_high_database))]
    dist = np.array([-m for m in simil])
    sampledist = dist.argsort()[:-Qsampleeach]

  ### TRAININING
  if Qtrain == 1:
    ### DATABASE LOADING ###
    ## The high level of theory
    clusters_df = train_high_database
    if Qsampleeach != 0:
      clusters_df = clusters_df.iloc[sampledist]
    if Qforcemonomers == 1:
      clusters_df0 = monomers_high_database
      clusters_df = clusters_df.append(clusters_df0, ignore_index=True)    
    ## The low level of theory
    if method == "delta":
      clusters_df2 = train_low_database
      if Qsampleeach != 0:
        clusters_df2 = clusters_df2.iloc[sampledist]
      if Qforcemonomers == 1:
        clusters_df0l = monomers_low_database
        clusters_df2 = clusters_df2.append(clusters_df0l, ignore_index=True)

    #Do we take only subset for training?
    if size != "full":
      if len(clusters_df) <= int(size):
        size = "full"
    if size != "full":
      clusters_df, clusters_df_trash, idx, idx_trash = train_test_split(clusters_df, range(len(clusters_df)), test_size=(len(clusters_df)-size)/len(clusters_df), random_state=seed)
      if method == "delta":
        clusters_df2 = clusters_df2.iloc[idx]
      size = "full" #IT IS BECAUSE I DO NOT WANT TO MAKE MY TEST SET SMALLER

    ### ENERGIES = VARIABLES / STRUCTURES
    ens = (clusters_df[column_name_1][column_name_2]).values.astype("float")
    strs = clusters_df["xyz"]["structure"]
    if method == "delta":
      ens2 = (clusters_df2[column_name_1][column_name_2]).values.astype("float")
      #str2 should be the same as str by principle

    ### FORCES 
    if ("extra","forces") in clusters_df.columns and Qifforces == 1:
      F_train = clusters_df["extra"]["forces"].values
      Qforces = 1
    else:
      Qforces = 0
    
    print(ens.shape, flush = True)
    if method == "delta":
      print(ens2.shape, flush = True)
      if ens.shape != ens2.shape:
        print("The LOW and HIGH method train sizes do not match. [EXITING]")
        exit()
  
    ### BINDING PROPERTIES CALCULATION (i.e. relative to monomers) ###
    #HIGH LEVEL
    ens_correction = substract_monomers(clusters_df,monomers_high_database)
    #LOW LEVEL
    if method == "delta":
      ens2_correction = substract_monomers(clusters_df2,monomers_low_database)
    
    #The binding (formation) energy calculation (or final property)
    form_ens = ens - ens_correction
    #print(form_ens, flush = True)
    if method == "delta":
      form_ens2 = ens2 - ens2_correction
      #print(form_ens2, flush = True)
    if method == "delta":
      Y_train = form_ens - form_ens2
    else:
      Y_train = form_ens
    
    if Qmethod == "krr":
      ### REPRESENTATION CALCULATION ###
      if Qrepresentation == "fchl":
        repres_dataframe = pd.DataFrame(index = strs.index, columns = ["xyz"])
        max_atoms = max([len(strs[i].get_atomic_numbers()) for i in range(len(strs))])
        for i in range(len(repres_dataframe)):
          repres_dataframe["xyz"][i] = generate_representation(strs[i].get_positions(), strs[i].get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=krr_cutoff)
        representations = np.array([mol for mol in repres_dataframe["xyz"]])
      elif Qrepresentation == "mbdf":
        X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
        representations = generate_representation(np.array([i.get_atomic_numbers() for i in strs]), np.array([i.get_positions() for i in strs]), cutoff_r = krr_cutoff, normalized = False)
      
      #some info about the full representation
      print(representations.shape, flush = True)
  
      ### DEFINING THE TRAINING X:  Y = QML(X) 
      X_train = representations
      
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
            f = open(varsoutfile.split(".pkl")[0]+"_"+str(splits)+"_"+str(s1)+"_"+str(s2)+".pkl","rb")
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
        f = open(varsoutfile,"wb")
        pickle.dump([X_train, sigmas, alpha],f)
        f.close()
        print("Training completed.", flush = True)
      elif Qsplit == 1:
        if Qrepresentation == "fchl":
          K = JKML_sym_kernel(X_train, sigmas, **kernel_args)       #calculates kernel
        elif Qrepresentation == "mbdf":
          K = [JKML_sym_kernel(X_train, X_atoms, sigmas[0], **kernel_args)]  #calculates kernel
        K = [K[i] + lambdas[i]*np.eye(len(K[i])) for i in range(len(sigmas))] #corrects kernel
        alpha = [cho_solve(Ki, Y_train) for Ki in K]                          #calculates regression coeffitients
  
        #I will for now everytime save the trained QML
        f = open(varsoutfile,"wb")
        if Qrepresentation == "fchl":
          pickle.dump([X_train, sigmas, alpha],f)
        elif Qrepresentation == "mbdf":
          pickle.dump([X_train, X_atoms, sigmas, alpha],f)
        f.close()
        print("Training completed.", flush = True)
      else:
        X_train_i = np.array_split(X_train,Qsplit)[Qsplit_i]
        X_train_j = np.array_split(X_train,Qsplit)[Qsplit_j]
        if Qsplit_i == Qsplit_j:
          K = JKML_sym_kernel(X_train_i, sigmas, **kernel_args)       #calculates kernel
          K = [K[i] + lambdas[i]*np.eye(len(K[i])) for i in range(len(sigmas))] #corrects kernel
        else:
          K = JKML_kernel(X_train_i, X_train_j, sigmas, **kernel_args)
        f = open(varsoutfile.split(".pkl")[0]+"_"+str(Qsplit)+"_"+str(Qsplit_i)+"_"+str(Qsplit_j)+".pkl","wb")
        pickle.dump([K,Y_train],f)
        f.close()
        exit()

    #####################################
    ### TRAINING NN #####################
    #####################################
    elif Qmethod == "nn":
      #PREPARING TRAINING DATABASE
      temperary_file_name = "training.db"
      if os.path.exists(temperary_file_name):
        os.remove(temperary_file_name)
      if Qforces == 0:
        new_dataset = ASEAtomsData.create(temperary_file_name,
          distance_unit='Ang',
          property_unit_dict={'energy':'eV', 'total_charge': 'e'},
          atomrefs = {'energy': [0]*100}
          )
        properties = [{'energy': np.array([i]), 'total_charge': np.array([0], dtype=np.float32)} for i in Y_train]
        target_properties = [spk.properties.energy]
        tradoffs = [1]
      else:
        new_dataset = ASEAtomsData.create(temperary_file_name,
          distance_unit='Ang',
          property_unit_dict={'energy':'eV', 'forces': 'eV/Ang', 'total_charge': 'e'},
          atomrefs = {'energy': [0]*100}
          )
        properties = [{'energy': 27.2107*np.array([Y_train[i]]), 'forces': 27.2114*np.array(F_train[i]), 'total_charge': np.array([0], dtype=np.float32)} for i in range(len(Y_train))]
        target_properties = [spk.properties.energy, spk.properties.forces]
        tradoffs = [Qenergytradoff, 1]
      new_dataset.add_systems(properties, strs)
     
      n_train = int(np.round(nn_tvv*len(strs)))
      n_val = len(strs) - n_train
      pl.seed_everything(seed)
      dataset = AtomsDataModule(temperary_file_name,
          batch_size=16,
          num_train=n_train,
          num_val=n_val,
          #num_test=n_test,
          transforms=[
              trn.ASENeighborList(cutoff=nn_cutoff),
              trn.RemoveOffsets(spk.properties.energy, remove_mean=True, remove_atomrefs=False),
              trn.CastTo32()
          ],
          num_workers=nw,#TODO how does this work?
          #split_file=split_file,
          data_workdir="./"
      )
      dataset.prepare_data()
      dataset.setup()
      logging.info(f'Number of train-val-test data: {dataset.num_train} - {dataset.num_val}')
      properties = dataset.dataset[0]
      logging.info('Loaded properties:')
      for k, v in properties.items():
        logging.info(f'     {k:20s} : {v.shape}')      

      # PainNN representation
      pairwise_distance = spk.atomistic.PairwiseDistances()
      radial_basis = spk.nn.GaussianRBF(n_rbf=nn_rbf, cutoff = nn_cutoff)
       
      if Qrepresentation == "painn": 
        X_train = spk.representation.PaiNN(
            n_atom_basis=nn_atom_basis,
            n_interactions=nn_interactions,
            radial_basis=radial_basis,
            cutoff_fn=spk.nn.CosineCutoff(nn_cutoff)
        )      
      elif Qrepresentation == "schnet":
        X_train = spk.representation.SchNet(
            n_atom_basis=nn_atom_basis,
            n_interactions=nn_interactions,
            radial_basis=radial_basis,
            cutoff_fn=spk.nn.CosineCutoff(nn_cutoff)
        )
      elif Qrepresentation == "so3net":
        X_train = spk.representation.SO3net(
            n_atom_basis=nn_atom_basis,
            n_interactions=nn_interactions,
            radial_basis=radial_basis,
            cutoff_fn=spk.nn.CosineCutoff(nn_cutoff)
        )
      else:
        print("You have probably enetered some weird molecular representation. [EXIT]")
        exit()

      output_modules = []
      output_losses = []
      for p, w in zip(target_properties, tradoffs):
          if p == spk.properties.energy:
              pred = spk.atomistic.Atomwise(n_in=nn_atom_basis, output_key=p)
          elif p == spk.properties.forces:
              pred = spk.atomistic.Forces(energy_key=spk.properties.energy, force_key=spk.properties.forces)
          elif p == spk.properties.dipole_moment:
              pred = spk.atomistic.DipoleMoment(n_in=nn_atom_basis, return_charges=True)
          else:
              raise NotImplementedError(f'{p} property does not exist')
      
          loss = spk.task.ModelOutput(
                  name=p,
                  loss_fn=torch.nn.MSELoss(),
                  loss_weight=w,
                  metrics={
                      "MAE": torchmetrics.MeanAbsoluteError(),
                      "RMSE": torchmetrics.MeanSquaredError(squared=False)
                  }
              )
          output_modules.append(pred)
          output_losses.append(loss)

      #MODEL (this could be for instance also Atomistic Model)
      nnpot = spk.model.NeuralNetworkPotential(
          representation=X_train,
          input_modules=[pairwise_distance],
          output_modules=output_modules,
          postprocessors=[
              trn.CastTo64(),
              trn.AddOffsets(spk.properties.energy, add_mean=True, add_atomrefs=False)
          ]
      )
 
      task = spk.task.AtomisticTask(
          model=nnpot,
          outputs=output_losses,
          optimizer_cls=torch.optim.AdamW,
          optimizer_args={"lr": Qlearningrate},
          scheduler_cls=spk.train.ReduceLROnPlateau,
          scheduler_args={'factor': 0.5, 'patience': 20, 'min_lr': 1e-7},
          scheduler_monitor = 'val_loss'
      )

      logger = pl.loggers.TensorBoardLogger(save_dir="./")
      #logger = pl.loggers.CSVLogger(save_dir=model_dir, flush_logs_every_n_steps=1)
      callbacks = [
          spk.train.ModelCheckpoint(
              model_path=os.path.join("./", varsoutfile),
              save_top_k=1,
              monitor="val_loss",
              save_last=True,
          ),
          pl.callbacks.EarlyStopping(
              monitor="val_loss",
              patience=Qearlystop,
          ),
          pl.callbacks.LearningRateMonitor(logging_interval='epoch')
      ]

      if torch.cuda.is_available():
          device = 'gpu'
          logging.info(torch.cuda.get_device_name(0))
      else:
          device = 'cpu'
      logging.info(f'Using device {device}')
      
      trainer = pl.Trainer(
          accelerator=device,
          #devices=1,
          callbacks=callbacks,
          logger=logger,
          default_root_dir="./",
          max_epochs=nn_epochs,
          #log_every_n_steps=1,
      )

      trainer.fit(task, datamodule=dataset)
    ######################
    #You should not get below this one to reach the error
    else:
      print("Wrong method or representation chosen.")
      exit()
  
  
  #LOAD TRAINING
  #TODO collect splitting: /home/kubeckaj/ML_SA_B/ML/TRAIN/TEST/SEPARATE/cho_solve.pkl
  if Qtrain == 2:
    if not os.path.exists(VARS_PKL):
      print("Error reading trained model. VARS_PKL = "+VARS_PKL)
      exit()
    if Qmethod == "krr":
      f = open(VARS_PKL,"rb")
      if Qrepresentation == "fchl":
        X_train, sigmas, alpha = pickle.load(f)
      elif Qrepresentation == "mbdf":
        X_train, X_atoms, sigmas, alpha = pickle.load(f)
      if len(alpha)!=1:
        alpha = [alpha]
      f.close()
      print("Trained model loaded.", flush = True)
    elif Qmethod == "nn":
      varsoutfile = VARS_PKL
      print("Trained model found.")
    else:
      print("Wrong method or representation chosen.")
      exit()
  
  ######################
  ###### EVALUATE ######
  ######################
  
  if Qeval == 1 or Qeval == 2:
    ### DATABASE LOADING ###
    ## The high level of theory
    clusters_df = test_high_database
    if Qsampleeach != 0:
      clusters_df = clusters_df.iloc[[sampleeach_i]]
    if method == "delta":
      clusters_df2 = test_low_database
      if Qsampleeach != 0:
        clusters_df2 = clusters_df2.iloc[[sampleeach_i]]

    #Do we take only subset for testing?
    if size != "full":
      if len(clusters_df) <= int(size):
        size = "full"
    if size != "full":
      clusters_df, clusters_df_trash, idx, idx_trash = train_test_split(clusters_df, range(len(clusters_df)), test_size=(len(clusters_df)-size)/len(clusters_df), random_state=seed)
      if method == "delta":
        clusters_df2 = clusters_df2.iloc[idx]
    clustersout_df = clusters_df.copy()

    if Qeval == 2:
      try:
        ens = (clusters_df[column_name_1][column_name_2]).values.astype("float")
        #Qeval = 2 #2=Compare ML prediction with QC
      except:
        Qeval = 1 #1=Only predicts
    if method == "delta":
      ens2 = (clusters_df2[column_name_1][column_name_2]).values.astype("float")
    strs = clusters_df["xyz"]["structure"]
      #str2 should be the same as str by princip
  
    ### FORCES 
    if ("extra","forces") in clusters_df.columns and Qifforces == 1:
      F_test = clusters_df["extra"]["forces"].values
      Qforces = 1
    elif Qprintforces == 1 and Qeval == 1 and Qifforces == 1:
      Qforces = 1
    else:
      Qforces = 0

    if Qeval == 2:
      print(ens.shape, flush = True)
    if method == "delta":
      print(ens2.shape, flush = True)
      if Qeval == 2:
        if ens.shape != ens2.shape:
          print("The LOW and HIGH method test sizes do not match. [EXITING]", flush = True)
          exit()
  
    ### BINDING PROPERTIES CALCULATION (i.e. relative to monomers) ###
    #HIGH LEVEL
    ens_correction = substract_monomers(clusters_df,monomers_high_database)
    #LOW LEVEL
    if method == "delta":
      ens2_correction = substract_monomers(clusters_df2,monomers_low_database)
    #TODO clustername given as argument
    #monomers_df = pd.read_pickle("/home/kubeckaj/ML_SA_B/DATABASES/database_XTB_monomers_at_DFT.pkl")
    #import re
    #clustername=np.array(re.split(r'(\d+)', sys.argv[3])[1:])
    #clustername=clustername.reshape(int(len(clustername)/2),2)
    #clustername=[clustername[:,(1,0)]]*len(clusters_df)
    #print(clustername)
    #print(clustername[0][0])
    #for i in clustername: #np.transpose([clusters_df["info"]["components"].values,clusters_df["info"]["component_ratio"].values]):
        
    #The binding (formation) energy calculation
    if Qeval == 2:
      form_ens = ens - ens_correction
      #print(form_ens, flush = True)
    if method == "delta":
      form_ens2 = ens2 - ens2_correction
      #print(form_ens2, flush = True) 

    ##################
    ### KRR + FCHL ###
    ################## 
    if Qmethod == "krr": 
      ### REPRESENTATION CALCULATION ###
      if Qrepresentation == "fchl":
        repres_dataframe = pd.DataFrame(index = strs.index, columns = ["xyz"])
        max_atoms = max([len(strs[i].get_atomic_numbers()) for i in range(len(strs))])
        for i in range(len(repres_dataframe)):
          repres_dataframe["xyz"][i] = generate_representation(strs[i].get_positions(), strs[i].get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=krr_cutoff)
        representations = np.array([mol for mol in repres_dataframe["xyz"]])
      elif Qrepresentation == "mbdf":
        X_test_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
        representations = generate_representation(np.array([i.get_atomic_numbers() for i in strs]), np.array([i.get_positions() for i in strs]), cutoff_r = krr_cutoff, normalized = False)

      #some info about the full representation
      print(representations.shape, flush = True)
      
      ### DEFINING THE EVALUATION Xs:  Y = QML(X) 
      #the full set
      X_test = representations
  
      ### CORRECTING THE FCHL MATRIX SIZES
      #IF YOU ARE EXTENDING THIS WILL MAKE THE MATRIXES OF THE SAME SIZE
      if Qrepresentation == "fchl":
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
      if Qrepresentation == "fchl":
        Ks = JKML_kernel(X_test, X_train, sigmas, **kernel_args)
      elif Qrepresentation == "mbdf":
        #Ks = [JKML_kernel(X_test, X_train, X_test_atoms, X_atoms, sigmas[0], **kernel_args)]
        Ks = [JKML_kernel(X_train, X_test, X_atoms, X_test_atoms, sigmas[0], **kernel_args)]
      Y_predicted = [np.dot(Ks[i], alpha[i]) for i in range(len(sigmas))]
    
    ##################
    ### NN + PaiNN ###
    ##################
    elif Qmethod == "nn":
      if torch.cuda.is_available():
        device = 'cuda'
      else:
        device = 'cpu'
      if Qforces == 0:
        spk_calc = SpkCalculator(
          model_file=varsoutfile,
          device=device,
          neighbor_list=spk.transform.ASENeighborList(cutoff=nn_cutoff),
          #neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
          #transforms=spk.transform.atomistic.SubtractCenterOfMass(),
          energy_key='energy',
          energy_unit="eV",#YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
          position_unit="Ang",
          )
      else:
        spk_calc = SpkCalculator(
          model_file=varsoutfile,
          device=device,
          neighbor_list=spk.transform.ASENeighborList(cutoff=nn_cutoff),
          #neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
          #transforms=spk.transform.atomistic.SubtractCenterOfMass(),
          energy_key='energy',
          force_key='forces',
          energy_unit="eV",#YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
          force_unit="eV/Ang",#YEAH I have no idea what the output is :-D
          position_unit="Ang",
          )
      Y_predicted = []
      F_predicted = []
      for i in range(len(clusters_df)):
        atoms = clusters_df["xyz"]["structure"].values[i].copy()
        atoms.calc = spk_calc
        Y_predicted.append(atoms.get_potential_energy())
        if Qforces == 1:
          F_predicted.append(atoms.get_forces())
      if Qforces == 1:
        Y_predicted = [0.0367493*np.array(Y_predicted)] #Hartree
        F_predicted = [0.0367493*np.array(F_predicted)] #Hartree/Ang
      else:
        Y_predicted = [np.array(Y_predicted)]

    else:
      print("Wrong method or representation chosen.")
      exit()
  
    ### PRINTING THE RESULTS
    if Qwolfram == 0:
      lb = "["
      rb = "]"
    else:
      lb = "{"
      rb = "}"
    #print("Ypredicted = {" + ",".join([str(i) for i in Y_predicted[0]])+"};")
    if Qeval == 2:
      print("Ytest = " + lb + ",".join([str(i) for i in ens]) + rb + ";", flush = True)
      if Qforces == 1:
        if Qprintforces == 0:
          print("Ftest = I will not print the forces becuase it would be too many numbers. (Use -printforces)", flush = True)
        else:
          print("Ftest = " + lb + ",".join([lb+",".join([lb+",".join([str(k) for k in j])+rb for j in i])+rb for i in F_test])+rb+";", flush = True)
    if method != "delta":
      print("Ypredicted = "+ lb + ",".join([str(i) for i in Y_predicted[0]+ens_correction])+rb+";", flush = True)
    else:
      #print("TEST",flush=True)
      #print(Y_predicted[0],flush=True)
      #print("TEST",flush=True)
      #print(form_ens2,flush=True)
      #print("TEST",flush=True)
      #print(ens_correction,flush=True)
      #print("TEST",flush=True)
      print("Ypredicted = "+ lb + ",".join([str(i) for i in Y_predicted[0]+form_ens2+ens_correction])+rb+";", flush = True)
    #print("formensval = {" + ",".join([str(i) for i in form_ens[idx]])+"};")
    #print("ens_correction = {" + ",".join([str(i) for i in ens_correction])+"};")
    if Qforces == 1:
      if Qprintforces == 0:
        print("Fpredicted = I will not print the forces becuase it would be too many numbers. (Use -printforces)", flush = True)
      else:
        print("Fpredicted = "+lb + ",".join([lb+",".join([lb+",".join([str(k) for k in j])+rb for j in i])+rb for i in F_predicted[0]])+rb+";", flush = True)
   
    # If possible print MAE and RMSE
    if Qeval == 2:
      if method != "delta":
        Y_validation = form_ens
      else:
        Y_validation = form_ens - form_ens2
      ### Calculate mean-absolute-error (MAE):
      if column_name_1 == "log" and column_name_2 == "electronic_energy":
        multiplier = 627.503
        units = " kcal/mol"
      else:
        multiplier = 1.0
        units = " [?]"
      print("", flush = True)
      print("Results:", flush = True)
      mae = [multiplier*np.mean(np.abs(Y_predicted[i] - Y_validation))  for i in range(len(sigmas))]
      print("Mean Absolute Error:", flush = True)
      print("mae = " + ",".join([str(i) for i in mae])+units, flush = True)
      ### Calculate root-mean-squared-error (RMSE):
      rmse = [multiplier*np.sqrt(np.mean(np.abs(Y_predicted[i] - Y_validation)**2))  for i in range(len(sigmas))]
      print("Root Mean Squared Error", flush = True)
      print("rmse = " + ",".join([str(i) for i in rmse])+units, flush = True)
      ### Calculate mean-absolute-relative-error (MArE):
      diff = [np.mean(Y_predicted[i]) - np.mean(Y_validation) for i in range(len(sigmas))]
      mare = [multiplier*np.mean(np.abs(Y_predicted[i] - Y_validation - diff[i]))  for i in range(len(sigmas))]
      print("Mean Absolute (mean-)Relative Error:", flush = True)
      print("mare = " + ",".join([str(i) for i in mare])+units, flush = True)
      ### Calculate root-mean-squared-relative-error (RMSrE):
      rmsre = [multiplier*np.sqrt(np.mean(np.abs(Y_predicted[i] - Y_validation - diff[i])**2))  for i in range(len(sigmas))]
      print("Root Mean Squared (mean-)Relative Error", flush = True)
      print("rmsre = " + ",".join([str(i) for i in rmsre])+units, flush = True)
      diff = [np.median(Y_predicted[i]) - np.median(Y_validation) for i in range(len(sigmas))]
      mare = [multiplier*np.mean(np.abs(Y_predicted[i] - Y_validation - diff[i]))  for i in range(len(sigmas))]
      print("Mean Absolute (median-)Relative Error:", flush = True)
      print("mare = " + ",".join([str(i) for i in mare])+units, flush = True)
      ### Calculate root-mean-squared-relative-error (RMSrE):
      rmsre = [multiplier*np.sqrt(np.mean(np.abs(Y_predicted[i] - Y_validation - diff[i])**2))  for i in range(len(sigmas))]
      print("Root Mean Squared (median-)Relative Error", flush = True)
      print("rmsre = " + ",".join([str(i) for i in rmsre])+units, flush = True)
      if Qforces == 1:
        print("", flush = True)
        print("Results for forces:", flush = True)
        mae = [np.mean(np.abs(np.array([np.array(j) for j in F_predicted[i]]).flatten() - np.array([np.array(j) for j in F_test]).flatten()))  for i in range(len(sigmas))]
        print("MAE = " + ",".join([str(i) for i in mae])+" [Eh/Angstrom]", flush = True)
        rmse = [np.sqrt(np.mean(np.abs(np.array([np.array(j) for j in F_predicted[i]]).flatten() - np.array([np.array(j) for j in F_test]).flatten())**2))  for i in range(len(sigmas))]
        print("RMSE = " + ",".join([str(i) for i in rmse])+" [Eh/Angstrom]", flush = True)
  
    ### PRINTING THE QML PICKLES
    #print(type(clustersout_df["xyz"]["structure"].values[0]))
    #print(type(Y_predicted[0][0]))
    #print(type(ens_correction[0]))
    for i in range(len(clustersout_df)):
      if method != "delta":
        clustersout_df.loc[clustersout_df.iloc[i].name,(column_name_1,column_name_2)] = Y_predicted[0][i]+ens_correction[i]
      else:
        clustersout_df.loc[clustersout_df.iloc[i].name,(column_name_1,column_name_2)] = Y_predicted[0][i]+form_ens2[i]+ens_correction[i]
      if Qforces == 1:
        clustersout_df = df_add_iter(clustersout_df, "extra", "forces", [clustersout_df.iloc[i].index], [F_predicted[i]])
    clustersout_df.to_pickle(outfile)
    if Qsampleeach > 0:
      if sampleeach_i == 0:
        os.system("JKQC "+outfile+" -out predicted_QML_FULL.pkl -noex")
      else:
        os.system("JKQC "+outfile+" predicted_QML_FULL.pkl -out predicted_QML_FULL.pkl -noex")
  ########
  
  ######################
  ###### OPTIMIZE ######
  ######################
  
  if Qopt > 0:
    ### DATABASE LOADING ###
    ## The high level of theory
    clusters_df = test_high_database
    strs = clusters_df["xyz"]["structure"]
    
    ##################
    ### KRR + FCHL ###
    ################## 
    if Qmethod == "krr" and Qrepresentation == "fchl":
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
          repres_dataframe["xyz"][i] = generate_representation(RR.values[i][0], strs[0].get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=krr_cutoff)
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
        Ks = JKML_kernel(X_test, X_train, sigmas, **kernel_args)
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
        clustersout_df = df_add_iter(clustersout_df, column_name_1, column_name_2, [str(cluster_id)], [save_energy])
        newxyz = strs[0].copy()
        newxyz.set_positions(xyzold)
        clustersout_df = clustersout_df = df_add_iter(clustersout_df, "xyz", "structure", [str(cluster_id)], [newxyz])
        print(str(step)+" \t"+str(save_energy)+" \t"+str(maxdev),flush = True)
        #print(xyz, flush = True)

    ##################
    ### NN + PaiNN ###
    ##################
    elif Qmethod == "nn":
      if torch.cuda.is_available():
        device = 'cuda'
      else:
        device = 'cpu'
      spk_calc = SpkCalculator(
        model_file=varsoutfile,
        device=device,
        neighbor_list=spk.transform.ASENeighborList(cutoff=nn_cutoff),
        #neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
        #transforms=spk.transform.atomistic.SubtractCenterOfMass(),
        energy_key='energy',
        force_key='forces',
        energy_unit="eV",#YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
        force_unit="eV/Ang",#YEAH I have no idea what the output is :-D
        position_unit="Ang",
        )
      Y_predicted = []
      F_predicted = []
      from ase.optimize import BFGS
      from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
      from ase.md.verlet import VelocityVerlet
      from ase.md.langevin import Langevin
      from ase import units
      from ase.io import read,write
      for i in range(len(clusters_df)):
        atoms = clusters_df["xyz"]["structure"].values[i].copy()
        atoms.calc = spk_calc
        if Qopt == 1:
          dyn = BFGS(atoms,maxstep=opt_maxstep)
          def printenergy(a=atoms):
            write("opt.xyz", a, append = True)
          dyn.attach(printenergy, interval=1)
          dyn.run(fmax=1e-6)
        else: 
          # Set the momenta corresponding to T
          MaxwellBoltzmannDistribution(atoms, temperature_K=md_temperature)
          # We want to run MD with constant energy using the VelocityVerlet algorithm.
          #dyn = VelocityVerlet(atoms, 1 * units.fs)  # 5 fs time step.
          dyn = Langevin(atoms, 0.1*units.fs, md_temperature*units.kB, 0.002) #friction coeffitient 0.002
          def printenergy(a=atoms):  # store a reference to atoms in the definition.
            """Function to print the potential, kinetic and total energy."""
            epot = a.get_potential_energy() / len(a)
            ekin = a.get_kinetic_energy() / len(a)
            write("traj.xyz", a, append = True)
            print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
                  'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
          # Now run the dynamics
          dyn.attach(printenergy, interval=1)
          printenergy()
          dyn.run(1000000)
 
      #  Y_predicted.append(atoms.get_potential_energy())
      #  if Qforces == 1:
      #    F_predicted.append(atoms.get_forces())
      #Y_predicted = [np.array(Y_predicted)]
      #if Qforces == 1:
      #  F_predicted = [np.array(F_predicted)]

    else:
      print("Wrong method or representation chosen.")
      exit()  
  
    ### PRINTING THE RESULTS
    #print("Ypredicted = {" + ",".join([str(i) for i in save_energy])+"};", flush = True)
  
    ### PRINTING THE QML PICKLES
    #clustersout_df = clusters_df.copy()
    #for i in range(len(clustersout_df)):
    #  clustersout_df.loc[clustersout_df.iloc[i].name,(column_name_1,column_name_2)] = Y_predicted[0][i]
    #clustersout_df.to_pickle(outfile)
  ########


