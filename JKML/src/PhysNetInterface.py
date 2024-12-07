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
                strs):
#                ,
#                Qenergytradoff,
#                nn_tvv,
#                nn_cutoff,
#                nw,
#                nn_rbf,
#                Qrepresentation,
#                nn_atom_basis,
#                nn_interactions,
#                Qbatch_size,
#                Qlearningrate,
#                parentdir,
#                seed,
#                varsoutfile,
#                Qearlystop,
#                nn_epochs,
#                Qcheckpoint,
#                Qtime):
   
  # PREPARING TRAINING DATABASE 
  np.savez("database.npz", R=np.array([i.get_positions() for i in strs]), Q=Q_charge, D=0.2081943*np.array([np.array(i) for i in D_dipole]), E=27.2107*Y_train, F=27.2107*np.array([np.array(i) for i in F_train]), Z=np.array([i.get_atomic_numbers() for i in strs]), N=np.array([len(i) for i in strs]))
 
  import os,sys
 
  #print('sys.argv[0] =', sys.argv[0])             
  pathname = os.path.dirname(sys.argv[0])        
  #print('path =', pathname)
  #print('full path =', os.path.abspath(pathname)) 
  command  = "python "
  command += os.path.abspath(pathname)+"/src/PhysNet/train.py "
  
  os.system(command)


