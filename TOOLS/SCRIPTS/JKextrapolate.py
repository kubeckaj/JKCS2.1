from ase.io import read
from sys import argv
from os import path
from pandas import DataFrame, to_pickle, read_pickle
from read_files import mergeDictionary
import numpy as np

def two_point_extrapolation(LOWEST_CARDINALNUMBER,lowPKL,highPKL,alpha,beta):
    E_SCF_X=lowPKL.loc[:,("log","scf_energy")].values
    E_SCF_Y=highPKL.loc[:,("log","scf_energy")].values
    if ("log","correlation_energy") in lowPKL.columns:
      E_CORR_X=lowPKL.loc[:,("log","correlation_energy")].values
      E_CORR_Y=highPKL.loc[:,("log","correlation_energy")].values
    else:
      E_CORR_X=lowPKL.loc[:,("log","electronic_energy")].values-E_SCF_X
      E_CORR_Y=highPKL.loc[:,("log","electronic_energy")].values-E_SCF_Y
    ep2_eX = np.exp(-alpha*np.sqrt(LOWEST_CARDINALNUMBER))
    ep2_eY = np.exp(-alpha*np.sqrt(LOWEST_CARDINALNUMBER+1))
    ep2_CBS_SCF = (E_SCF_X*ep2_eY - E_SCF_Y*ep2_eX)/(ep2_eY-ep2_eX)
    ep2_CBS_CORR = (LOWEST_CARDINALNUMBER**beta*E_CORR_X - (LOWEST_CARDINALNUMBER+1)**beta*E_CORR_Y)/(LOWEST_CARDINALNUMBER**beta-(LOWEST_CARDINALNUMBER+1)**beta)
    return (ep2_CBS_CORR+ep2_CBS_SCF)

def cps_laf_extrapolation(E_low,E_high,F):
    return E_low+F*(E_high-E_low)


def readAndSortPickles(files):
  '''
  Check if list of files exist and that they have matching file_basenames.
  Return the read pickles (dataframes) with coloumns stored by the file_basename. 
  '''
  for file in files:
    if not path.exists(file):
      print(f"The file '{file}' does not exist.", flush=True)
      exit()
      
    pickles = list((f.sort_values(by=[('info','file_basename')]) for f in [read_pickle(f) for f in files])) #sort according to basename
    for idx in range(len(pickles) - 1):
      if not (pickles[idx].loc[:,('info','file_basename')].values == pickles[idx+1].loc[:,('info','file_basename')].values).any():
        print('Pickles have mismatched entries')
        exit()
  return pickles



OUTPUT = []
HEADER = ['basename']
####HANDLE INPUT#########
method = argv[1]
if method == '2p' or method == '2pf':
  doLNOCPS = False
  files = [*argv[2].split(','),*argv[3].split(',')]  #Always make a list of files even if input is given as 
  lowestCardinalNumber = int(argv[4])
  alpha = float(argv[5])
  beta = float(argv[6])
  if len(files) == 4:
      doLNOCPS = True
      F_val = float(argv[7])
  
  pickles = readAndSortPickles(files) 
  OUTPUT.append(pickles[0].loc[:,('info','file_basename')].values)

  if doLNOCPS == True:
    OUTPUT.append(two_point_extrapolation(lowestCardinalNumber,pickles[0],pickles[2],alpha,beta))
    HEADER.append(f'E_CBS_LowerSettings')
    OUTPUT.append(two_point_extrapolation(lowestCardinalNumber,pickles[1],pickles[2],alpha,beta))
    HEADER.append('E_CBS_HigherSettings')
    OUTPUT.append(cps_laf_extrapolation(OUTPUT[-2],OUTPUT[-1],F_val))
    HEADER.append('E_CBS_CPS/LAF')
  else:
    OUTPUT.append(two_point_extrapolation(lowestCardinalNumber,pickles[0],pickles[1],alpha,beta))
    HEADER.append('E(CBS)')
      
elif method == '2s' or method == '2sf':
  files = [*argv[2].split(','),*argv[3].split(',')]
  F_val = float(argv[4])
  pickles = readAndSortPickles(files) 
  OUTPUT.append(pickles[0].loc[:,('info','file_basename')].values)
  E_X=pickles[0].loc[:,("log","electronic_energy")].values
  E_Y=pickles[1].loc[:,("log","electronic_energy")].values
  OUTPUT.append(cps_laf_extrapolation(E_X,E_Y,F_val))
  HEADER.append('E_CPS/LAF')
  
else:
  print(f"Unknown method '{method}'")
  exit()

file_name = "0"
if method == "2pf":
    if doLNOCPS == True:
        file_name = argv[8]
    else:
        file_name = argv[7]
if method == "2sf":
    file_name = argv[5]
if file_name != "0":
    df = read_pickle(file_name)
    print(f'Will write results to {file_name} under the "out" "electronic_energy" column.')    

print(*HEADER)
for x in zip(*OUTPUT):
    if file_name != "0":
        mask = df[('info', 'file_basename')] == x[0]
        df.loc[mask, ('out', 'electronic_energy')] = x[-1]
    print(*x)

if file_name != "0":
    df.to_pickle(file_name)
    
  
  


