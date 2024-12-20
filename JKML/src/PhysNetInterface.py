# import sys
import os.path

## Lib 1
import numpy as np
import time

########################################################################################################################
########################################################################################################################
########################################################################################################################
def training_nn(Qforces,
                Y_train,
                F_train,
                D_dipole,
                Q_charge,
                strs,
                nn_atom_basis, nn_rbf, nn_interactions, nn_cutoff, nn_tvv, Qbatch_size, seed, nn_epochs, Qlearningrate, varsoutfile):
#                Qenergytradoff,
#                nw,
#                Qrepresentation,
#                parentdir,
#                Qearlystop,
#                Qcheckpoint,
#                Qtime):
   
  # PREPARING TRAINING DATABASE 
  Nmax = np.max([len(i) for i in strs])
  num = len(strs)
  N = np.zeros([num], dtype=int) #Number of atoms
  E = np.zeros([num], dtype=float) #Potential energy with respect to atoms
  Q = np.zeros([num], dtype=float) #Total charge
  D = np.zeros([num, 3], dtype=float) #Dipole 
  Z = np.zeros([num, Nmax], dtype=int) #Nuclear charges/atomic numbers of nuclei
  R = np.zeros([num, Nmax, 3], dtype=float) #Cartesian coordinates of nuclei
  F = np.zeros([num, Nmax, 3], dtype=float) #Forces

  np.savez("database.npz", R=np.array([i.get_positions() for i in strs]), Q=Q_charge, D=0.2081943*np.array([np.array(i) for i in D_dipole]), E=27.2107*Y_train, F=27.2107*np.array([np.array(i) for i in F_train]), Z=np.array([i.get_atomic_numbers() for i in strs]), N=np.array([len(i) for i in strs]))
 
  import os,sys
 
  #print('sys.argv[0] =', sys.argv[0])             
  pathname = os.path.dirname(sys.argv[0])        
  #print('path =', pathname)
  #print('full path =', os.path.abspath(pathname)) 
  command  = "export PYTHONUNBUFFERED=1; python "
  command += os.path.abspath(pathname)+"/src/PhysNet_DER/run_train.py @input.inp "

  if not os.path.isfile("input.inp"): 
    from torch import cuda
    if cuda.is_available():
      device = 'cuda'
      print(cuda.get_device_name(0))
    else:
      device = 'cpu'

    with open('input.inp', 'w') as file:
      # Print directly to the file
      print("--restart=No", file=file)
      print("--num_features="+str(nn_atom_basis), file=file)
      print("--num_basis="+str(nn_rbf), file=file)
      print("--num_blocks="+str(nn_interactions), file=file)
      print("--num_residual_atomic=2", file=file)
      print("--num_residual_interaction=3", file=file)
      print("--num_residual_output=1", file=file)
      print("--cutoff="+str(nn_cutoff), file=file)
      print("--use_electrostatic=1", file=file)
      print("--use_dispersion=1", file=file)
      print("--grimme_s6=0.5", file=file)
      print("--grimme_s8=0.2130", file=file)
      print("--grimme_a1=0.0", file=file)
      print("--grimme_a2=6.0519", file=file)
      print("--dataset=database.npz", file=file)
      print("--num_train="+str(round(nn_tvv*len(Y_train))), file=file)
      print("--num_valid="+str(round((1-nn_tvv)*len(Y_train))), file=file)
      print("--batch_size="+str(Qbatch_size), file=file)
      print("--valid_batch_size="+str(Qbatch_size), file=file)
      print("--seed="+str(seed), file=file)
      print("--max_steps="+str(nn_epochs), file=file)
      print("--learning_rate="+str(Qlearningrate), file=file)
      print("--decay_steps="+str(nn_epochs), file=file)
      print("--decay_rate=0.1", file=file)
      print("--lambda_conf=0.4", file=file)
      print("--max_norm=1000.0", file=file)
      print("--ema_decay=0.999", file=file)
      print("--rate=0.0", file=file)
      print("--summary_interval=5", file=file)
      print("--validation_interval=5", file=file)
      print("--show_progress=False", file=file)
      print("--save_interval=5", file=file)
      print("--record_run_metadata=0", file=file)
      print("--device="+device, file=file)
  else:
    print("JKML: THE OLD input.inp FILE IS USED INSTEAD OF ANY NEW COMMANDS!!!")
    from torch import cuda 
    if cuda.is_available():
      os.system("sed -i 's/--device=cpu/--device=cuda/g' input.inp")
    else:
      os.system("sed -i 's/--device=cuda/--device=cpu/g' input.inp")
     
  os.system(command)
  os.system("cp */best/best_model.pt "+varsoutfile)

########################################################################################################################
########################################################################################################################
########################################################################################################################

def evaluating_nn(varsoutfile,
                  clusters_df,
                  method,
                  Qmin):

  import os,sys
  pathname = os.path.dirname(sys.argv[0])
  sys.path.append(pathname + "/src/PhysNet_DER/")
  from PNcalculator import PhysNetCalculator

  from torch import cuda 
  if cuda.is_available():
    os.system("sed -i 's/--device=cpu/--device=cuda/g' input.inp")
  else:
    os.system("sed -i 's/--device=cuda/--device=cpu/g' input.inp")

  Y_predicted = []
  F_predicted = []
  for i in range(len(clusters_df)):
    atoms = clusters_df["xyz"]["structure"].values[i].copy()
    if ("log","charge") in clusters_df.columns:
      charge = clusters_df["log"]["charge"].values[i]
    else:
      charge = 0

    calc = PhysNetCalculator(
      checkpoint=varsoutfile,
      atoms=atoms,
      charge=charge,
      config='input.inp')

    atoms.set_calculator(calc)
    Y_predicted.append(0.0367493 * atoms.get_potential_energy())
    F_predicted.append(0.0367493 *atoms.get_forces())  # Hartree/Ang

  Y_predicted = [np.array(Y_predicted)]  # Hartree 0.0367493 *
  F_predicted = [F_predicted]

  if method == "min":
    Y_predicted[0] += Qmin

  return Y_predicted, F_predicted
