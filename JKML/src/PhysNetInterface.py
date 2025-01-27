####################################################################################################
####################################################################################################

def training(Qforces,Y_train,F_train,D_dipole,Q_charge,Q_charges,strs,nn_atom_basis,nn_rbf,nn_interactions,nn_cutoff,nn_tvv,Qbatch_size,seed,nn_epochs,Qlearningrate,varsoutfile):
  #nw, Qrepresentation,Qearlystop,Qcheckpoint,Qtime):

  import os.path
  import numpy as np
  import os
  from torch import cuda
  import sys
   
  # PREPARING TRAINING DATABASE 
  Nmax = np.max([len(i) for i in strs])
  num = len(strs)
  N = np.zeros([num], dtype=int) #Number of atoms
  E = np.zeros([num], dtype=float) #Potential energy with respect to atoms
  Q = np.zeros([num], dtype=float) #Total charge
  Qa= np.zeros([num, Nmax], dtype=float) #TPartial charges
  D = np.zeros([num, 3], dtype=float) #Dipole 
  Z = np.zeros([num, Nmax], dtype=int) #Nuclear charges/atomic numbers of nuclei
  R = np.zeros([num, Nmax, 3], dtype=float) #Cartesian coordinates of nuclei
  F = np.zeros([num, Nmax, 3], dtype=float) #Forces
  for i in range(num):
    N[i] = len(strs.values[i])
    E[i] = 27.2107*Y_train[i]
    Q[i] = Q_charge[i]
    Qa[i,:N[i]] = Q_charges[i]
    D[i] = 0.2081943*np.array(D_dipole[i])
    Z[i,:N[i]] = np.array(strs.values[i].get_atomic_numbers())
    R[i,:N[i],:] = np.array(strs.values[i].get_positions())
    F[i,:N[i],:] = 27.2107*np.array(F_train[i])
  
  np.savez("database.npz", R=R, Q=Q, Qa=Qa, D=D, E=E, F=F, Z=Z, N=N)
 
  pathname = os.path.dirname(sys.argv[0])        

  if not os.path.isfile("input.inp"): 
    if cuda.is_available():
      device = 'cuda'
      print("JKML(PhysNetInterface): Using GPU: "+str(cuda.get_device_name(0)))
    else:
      device = 'cpu'
      print("JKML(PhysNetInterface): Using CPU")

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
      print("--decay_steps="+str(1000), file=file)
      print("--decay_rate=0.1", file=file)
      print("--lambda_conf=0.4", file=file)
      print("--max_norm=1000.0", file=file)
      print("--ema_decay=0.999", file=file)
      print("--rate=1.0", file=file)
      print("--summary_interval=5", file=file)
      print("--validation_interval=5", file=file)
      print("--show_progress=False", file=file)
      print("--save_interval=5", file=file)
      print("--record_run_metadata=0", file=file)
      print("--device="+device, file=file)
  else:
    print("JKML(PhysNetInterface): THE OLD input.inp FILE IS USED INSTEAD OF ANY NEW COMMANDS!!!")
    if cuda.is_available():
      os.system("sed -i 's/--device=cpu/--device=cuda/g' input.inp")
    else:
      os.system("sed -i 's/--device=cuda/--device=cpu/g' input.inp")
  
  command  = "export PYTHONUNBUFFERED=1; python "
  command += os.path.abspath(pathname)+"/src/PhysNet_DER/run_train.py @input.inp "
  os.system(command)
  os.system("cp */best/best_model.pt "+varsoutfile)

####################################################################################################
####################################################################################################

def evaluate(varsoutfile,clusters_df,method,Qmin):

  import os,sys
  from PNcalculator import PhysNetCalculator
  from torch import cuda 

  pathname = os.path.dirname(sys.argv[0])
  sys.path.append(pathname + "/src/PhysNet_DER/")

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

    #NOT APPRECIATED: atoms.set_calculator(calc)
    atoms.calc = calc
    Y_predicted.append(0.0367493 * atoms.get_potential_energy()) #Hartree
    F_predicted.append(0.0367493 * atoms.get_forces())  # Hartree/Ang

  Y_predicted = [np.array(Y_predicted)]  
  F_predicted = [F_predicted]

  if method == "min":
    Y_predicted[0] += Qmin

  return Y_predicted, F_predicted
