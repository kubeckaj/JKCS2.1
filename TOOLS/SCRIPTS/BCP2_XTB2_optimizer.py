#import scipy as sp
import os
import subprocess
import numpy as np
from sys import argv
import shutil
import time

QSCcpu = int(argv[1])
Qparallel = int(argv[2])

def random_folder():
  import random
  testik = 0
  while testik == 0:
    foldername="XTB"+str(random.randint(0,9999))
    if not os.path.exists(foldername):
      testik = 1
  return foldername
  

def cost_function(arr):
    foldername = random_folder()
    os.makedirs(foldername)
    os.chdir(foldername)
    with open('parameter.txt', 'w') as file:
        file.write(' '.join(map(str, arr)))
    process = subprocess.Popen('sh ../XTB3_runXTB.sh '+str(QSCcpu), shell=True)
    process.wait()
    while not os.path.exists('result'):
      time.sleep(1)
    with open('result', 'r') as file:
      number_str = file.readline().strip()
    os.chdir('../')
    shutil.rmtree(foldername)
    return float(number_str)

def parr_cost_function(arr):
    foldername = random_folder()
    os.makedirs(foldername)
    os.chdir(foldername)
    with open('parameter.txt', 'w') as file:
        file.write(' '.join(map(str, arr)))
    process = subprocess.Popen('echo "sbatch -n '+str(Qparallel)+' JKsend bash ../XTB3_runXTB.sh '+str(Qparallel)+' | awk \'{print \$4}\' >> .jobs.txt" > todo.txt; echo echo done >> todo.txt', shell = True)
    process.wait()
    process = subprocess.Popen('manager.sh todo.txt 1 SUB', shell = True)
    process.wait()
    while not os.path.exists('result'):
      time.sleep(1)
    with open('result', 'r') as file:
      number_str = file.readline().strip()
    os.chdir('../')
    shutil.rmtree(foldername)
    return float(number_str)

with open('initial.txt', 'r') as file:
  line = file.readline()
  initial0 = np.array([float(i) for i in line.split()])

#model = sp.optimize.minimize(test, x0=initial, method='L-BFGS-B',options={'ftol':1e-4,'gtol':1e-5,'eps':1e-4})
#print(model)

diff_step = 1e-2
hess_step = diff_step
max_iter = 50
parr = np.array(len(initial0)*[1.0])
if Qparallel > 1:
  from joblib import Parallel, delayed
  import multiprocessing
  num_cores = QSCcpu

#JK optimizer
def GHM(parr,initial0,diff_step,og_cost):
    #l1=[initial0*generate_step(parr,idx,diff_step) for idx in range(len(parr))]
    #l2=[initial0*generate_step(parr,idx,-diff_step) for idx in range(len(parr))]
    #lboth=l1.append(l2)
    lboth=[initial0*generate_step(parr,idx,where) for idx in range(len(parr)) for where in [diff_step,-diff_step]]
    if Qparallel > 1:
      cost_all = Parallel(n_jobs=num_cores)(delayed(parr_cost_function)(sub_i) for sub_i in lboth)
    else:
      cost_all = [cost_function(sub_i) for sub_i in lboth]
    cost_forward = np.array(cost_all)[0:len(parr)]
    cost_backward = np.array(cost_all)[len(parr):2*len(parr)]
    den = (cost_forward-2*og_cost+cost_backward)
    den = den/np.abs(den)*(np.abs(den)+diff_step)
    return ((cost_forward - og_cost) * diff_step) / den

def generate_step(parr,idx,diff_step):
    new_parr = parr.copy()
    new_parr[idx] = new_parr[idx]+diff_step
    return np.array(new_parr)

def lets_optimize(parr,initial0,diff_step,og_cost,hess_step):
    print(f">>> Optimizing for input {initial0}", flush = True)
    print('Original: ',og_cost, " kcal/mol/atom", flush = True)
    rem=og_cost
    for iter in range(max_iter):
        friction = 0.99
        print("####################################")
        print(f'>>> Iteration {iter}', flush = True)
        move = GHM(parr,initial0,diff_step,og_cost)
        if np.linalg.norm(move) > 0.1:
          move = 0.1*move/np.abs(move)
        print('\nMove :', move, flush = True)
        next_test = 0
        frict_red = 0
        while next_test < 2:
          if np.linalg.norm(friction*move) < 1e-10: #diff_step/10:
            print("DONE -- too short step.", flush = True)
            return (parr,og_cost)
          if next_test == 0:
            parr_test = parr-friction*move
            if Qparallel > QSCcpu:
              og_cost_test = parr_cost_function(parr_test*initial0)
            else:
              og_cost_test = cost_function(parr_test*initial0)
            print('\nTest :', og_cost_test, flush = True)
          else:
            parr_test_EXT = parr-friction*move
            if Qparallel > QSCcpu:
              og_cost_test_EXT = parr_cost_function(parr_test_EXT*initial0)
            else:
              og_cost_test_EXT = cost_function(parr_test_EXT*initial0)
            print('\nTest - ext :', og_cost_test_EXT, flush = True)
            if og_cost_test_EXT < og_cost_test:
              parr_test = parr_test_EXT
              og_cost_test = og_cost_test_EXT
              next_test = 1
              friction = 0.66*friction
              continue
            else:
              friction = friction/0.66  
              next_test = 2
              if frict_red == 0:
                friction = friction/0.66
              continue
          if og_cost_test < og_cost:
            next_test = 1
            friction = 0.66*friction
          else:
            friction = 0.66*friction
            frict_red = frict_red + 1
          if frict_red == 100:
            print("DONE -- too much friction used", flush = True)
            return (parr,og_cost)
        parr = parr_test
        if np.abs(og_cost_test-og_cost)*2/(og_cost_test+og_cost) < 0.0001:
          print("DONE -- too small cost function improvements", flush = True)
          return (parr,og_cost)
        og_cost = og_cost_test

        
        
        print('New Parameters:',parr*initial0, flush = True)  
        print('Cost function: ',og_cost, " kcal/mol/atom (original was ",rem,"kcal/mol/atom)", flush = True)
        initial0 = parr*initial0 
        parr= np.array(len(initial0)*[1.0])
        #print(parr*initial0, og_cost, flush = True)
        with open('LAST_parameter.txt', 'w') as file:
          file.write(' '.join(map(str, initial0)))
        process = subprocess.Popen('NUMBA=$(grep "AAA" param_gfn1-xtb.txt -o | wc -l);TEST=($(cat LAST_parameter.txt)); cp param_gfn1-xtb.txt LAST_param_gfn1-xtb.txt; for ((i = 1; i <= NUMBA; i++)); do loop_counter=$((i - 1)); sed -i "s/AAA$i /${TEST[$loop_counter]} /g" LAST_param_gfn1-xtb.txt;done', shell = True)
        process.wait()
  
    return (initial0,og_cost)

if Qparallel > QSCcpu:
  og_cost = parr_cost_function(initial0)
else:
  og_cost = cost_function(initial0)
final, og_cost = lets_optimize(parr,initial0,diff_step,og_cost,hess_step)
print('DONE Cost function: ',og_cost, " kcal/mol/atom ", flush = True)


