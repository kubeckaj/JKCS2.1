import scipy as sp
import os
import subprocess
import numpy as np
from sys import argv
import shutil
import time
import copy
import random
with open('initial.txt', 'r') as file:
  lines = file.readlines()
# Remove any trailing whitespace from each line
lines = [line.strip() for line in lines]
# Convert the first line to a numpy array
#initial0 = np.fromstring(lines[0], sep=' ')
initial0 = np.array([float(i) for i in lines[0].split()])
# Check if there's a second line and if it's not empty
try:
  if len(lines[1]) > 1 and lines[1]:
    #diff_step = np.fromstring(lines[1], sep=' ')
    diff_step = np.array([float(i) for i in lines[1].split()])
  else:
    # If the second line is empty or does not exist, create an array of 0.001s with the same length as array1
    diff_step = np.full_like(initial0, 0.001)
except:
  diff_step = np.full_like(initial0, 0.001)

Qparallel = int(argv[2])
QSCcpu = int(argv[1])
max_iter = 100
parr = np.array(len(initial0)*[1.0])
if Qparallel > 1:
  from joblib import Parallel, delayed
  import multiprocessing
  num_cores = QSCcpu

#folder_num = 0
def random_folder():
  #global folder_num
  testik = 0
  while testik == 0:
    #folder_num += 1
    foldername="XTB"+str(random.randint(1, 9999))+str(random.randint(1, 9999))
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
    try:
      toreturn = float(number_str)
    except:
      toreturn = np.nan
    return toreturn

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
    try:
      toreturn = float(number_str)
    except:
      toreturn = np.nan
    return toreturn



#JK optimizer
def GHM(parr,initial0,diff_step,og_cost):
    #l1=[initial0*generate_step(parr,idx,diff_step) for idx in range(len(parr))]
    #l2=[initial0*generate_step(parr,idx,-diff_step) for idx in range(len(parr))]
    #lboth=l1.append(l2)
    #print("WTF: parr", parr,  flush = True)
    #print("OMG:", [where*diff_step[idx] for where in [1,-1] for idx in range(len(parr))], flush = True)
    lboth=[initial0*generate_step(parr,idx,where*diff_step[idx]) for where in [1,-1] for idx in range(len(parr))]
    if Qparallel > 1:
      cost_all = Parallel(n_jobs=num_cores)(delayed(parr_cost_function)(sub_i) for sub_i in lboth)
    else:
      cost_all = [cost_function(sub_i) for sub_i in lboth]
    print("", flush = True)
    #print(cost_all)
    #print(diff_step)
    cost_forward = np.array(cost_all)[0:len(parr)]
    cost_backward = np.array(cost_all)[len(parr):2*len(parr)]
    den = (cost_forward-2*og_cost+cost_backward)
    res = []
    best_found = og_cost
    best_parr = parr
    for i in range(len(parr)):
      if best_found > cost_backward[i]:
        best_found = np.array(cost_all)[i+len(parr)]
        best_parr = copy.deepcopy(generate_step(parr,i,-diff_step[i]))
      if best_found > cost_forward[i]:
        best_found = np.array(cost_all)[i]
        best_parr = copy.deepcopy(generate_step(parr,i,diff_step[i]))
      if cost_forward[i] == og_cost or cost_backward[i] == og_cost or np.abs(den[i]) < 1e-99 or np.isnan(den[i]):
        if cost_forward[i] < og_cost:
          res.append(-diff_step[i])
          #print([cost_backward[i],og_cost,cost_forward[i],diff_step[i],den[i],-diff_step[i]])
        #elif cost_forward[i] > og_cost:
        #  res.append(diff_step[i])
        #elif cost_backward[i] > og_cost:
        #  res.append(-diff_step[i])
        elif cost_backward[i] < og_cost:
          res.append(diff_step[i])
          #print([cost_backward[i],og_cost,cost_forward[i],diff_step[i],den[i],diff_step[i]])
        elif np.isnan(den[i]):
          res.append(0)
          diff_step[i] = diff_step[i]/4
        else:
          res.append(0)
          #print([cost_backward[i],og_cost,cost_forward[i],diff_step[i],den[i],0])
        diff_step[i] = diff_step[i]*2
      else:
        premove = ((cost_forward[i] - og_cost) * diff_step[i]) / den[i]
        if cost_backward[i] < og_cost and cost_forward[i] < og_cost:
          if cost_backward[i] < cost_forward[i]:
            premove = diff_step[i]
          else: 
            premove = -diff_step[i]
          diff_step[i] = diff_step[i]/2
        elif cost_backward[i] > og_cost and cost_forward[i] > og_cost:
          diff_step[i] = diff_step[i]/2
        elif cost_forward[i] < og_cost:
          if np.abs(premove) > np.abs(diff_step[i]):
            premove = -diff_step[i]
            diff_step[i] = diff_step[i]*2
        elif cost_backward[i] < og_cost:
          if np.abs(premove) > np.abs(diff_step[i]):
            premove = diff_step[i]
            diff_step[i] = diff_step[i]*2
        else:
          print("Oh Jakub probably forgot some condition, tell him you got this message.")

        #THIS SHOULD NOT BE NECESSARY
        #if np.abs(premove) > diff_step[i]:
        #  premove = np.sign(premove)*diff_step[i]
        #  diff_step[i] = diff_step[i]*2
        res.append(premove)
      print("PARR "+str(i)+":",[cost_backward[i],og_cost,cost_forward[i],parr[i],diff_step[i],den[i],res[i]],flush = True)
    #den = den/np.abs(den)*(np.abs(den)+diff_step)
    return np.array(res), diff_step, best_found, best_parr

def generate_step(parr,idx,diff_step):
    new_parr = copy.deepcopy(parr)
    new_parr[idx] = new_parr[idx]+diff_step
    return np.array(new_parr)

#THIS IS NOT USED BUT COULD USEFUL IN FUTURE FOR SOMETHING
def f_jac_hes(x0, delta=1e-4):
    if all(x0==f_jac_hes.lastx):
      return f_jac_hes.lastf

    n = len(x0)
    H = np.zeros((n, n))
    gradient = np.zeros(n)

    # Create list of unique points
    points = [x0]
    for i in range(n):
        x_i_plus = x0.copy()
        x_i_minus = x0.copy()
        x_i_plus[i] += delta
        x_i_minus[i] -= delta
        points.extend([x_i_plus, x_i_minus])
        
        for j in range(i+1, n):
            x_ij_plus = x0.copy()
            x_ij_minus_i = x0.copy()
            x_ij_minus_j = x0.copy()

            x_ij_plus[i] += delta
            x_ij_plus[j] += delta

            x_ij_minus_i[i] += delta
            x_ij_minus_i[j] -= delta

            x_ij_minus_j[i] -= delta
            x_ij_minus_j[j] += delta

            points.extend([x_ij_plus, x_ij_minus_i, x_ij_minus_j])

    # Evaluate function in parallel for all unique points
    #values = parallel_eval(points, f)
    if Qparallel > 1:
      values = Parallel(n_jobs=num_cores)(delayed(parr_cost_function)(sub_i) for sub_i in points)
    else:
      values = [cost_function(sub_i) for sub_i in points]
    values = np.array(values)
    print("\n")
    print(values[0])
    
    # Construct the Hessian and Gradient
    for i in range(n):
        gradient[i] = (values[2*i+1] - values[2*i+2]) / (2*delta)
        H[i, i] = (values[2*i+1] - 2*values[0] + values[2*i+2]) / (delta**2)
        
        for j in range(i+1, n):
            idx_offset = 2*n + 3*(n*(n-1)//2 - (n-i)*(n-i-1)//2)
            H[i, j] = (values[idx_offset] - values[idx_offset+1] - values[idx_offset+2] + values[0]) / (4*delta**2)
            H[j, i] = H[i, j]
    f_jac_hesss.lastx = x
    f_jac_hesss.lastf = H, gradient
    return H, gradient

def lets_optimize(parr,initial0,diff_step,og_cost):
    print(f">>> Optimizing for input {initial0}", flush = True)
    print('Cost function: ',og_cost, " kcal/mol/atom (original)", flush = True)
    rem=og_cost
    friction = 0.5
    for iter in range(max_iter):
        print("####################################")
        print(f'>>> Iteration {iter}', flush = True)
        move, diff_step, best_found, best_parr = GHM(parr,initial0,diff_step,og_cost)
        #print("WTF: best parr", best_parr,  flush = True)
        #if np.linalg.norm(move) > 0.1:
        #  move = 0.1*move/np.linalg.norm(move)
        print('\nMove :', move, flush = True)
        print("Test : starting from", og_cost,  flush = True)
        print("Test : best found", best_found,  flush = True)
        #print('\nMove :', np.linalg.norm(move), flush = True)
        next_test = 0
        frict_red = 0
        prev_test = np.nan
        while next_test < 2:
          if np.linalg.norm(friction*move) < 1e-6: #diff_step/10:
            print("DONE -- too short step.", flush = True)
            parr_test = parr
            og_cost_test = og_cost
            break
            #return (parr,og_cost)
          if next_test == 0:
            parr_test = parr-friction*move
            if Qparallel > QSCcpu:
              og_cost_test = parr_cost_function(parr_test*initial0)
            else:
              og_cost_test = cost_function(parr_test*initial0)
            print('\nTest :', og_cost_test, flush = True)
            if np.isnan(prev_test):
              prev_test = og_cost_test 
            else:
              if np.abs(og_cost_test - prev_test) <= 1e-7:
                og_cost_test = og_cost
                parr_test = parr
                print("DONE -- too small imrovement.", flush = True)
                break
              else:
                prev_test = og_cost_test
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
              friction = friction/2
              continue
            else:
              friction = friction/2  
              next_test = 2
              if frict_red == 0:
                friction = friction/2
              continue
          if og_cost_test < og_cost:
            next_test = 1
            friction = friction/2
          else:
            friction = friction/2
            frict_red = frict_red + 1
          if frict_red == 100:
            print("DONE -- too much friction used", flush = True)
            return (parr,og_cost)
              
        if best_found < og_cost_test and best_found < og_cost:
          og_cost_test = best_found
          #print("MEGA WTF:", best_parr, flush = True)
          parr_test = best_parr.copy()
        parr = parr_test.copy()
        if np.abs(og_cost_test-og_cost)*2/(og_cost_test+og_cost) < 0.0001:
          print("WARNING -- too small cost function improvements", flush = True)
          #return (parr,og_cost)
        og_cost = og_cost_test
        #print("MEGA WTF:", parr, flush = True)
        #print("MEGA WTF2:", initial0, flush = True)

        
        
        print('New Parameters:',parr*initial0, flush = True)  
        print('Cost function: ',og_cost, " kcal/mol/atom (original was ",rem,"kcal/mol/atom)", flush = True)
        #TODO
        #initial0 = parr*initial0 
        #parr= np.array(len(initial0)*[1.0])
        
        #print(parr*initial0, og_cost, flush = True)
        with open('LAST_parameter.txt', 'w') as file:
          file.write(' '.join(map(str, parr*initial0))+"\n")
        process = subprocess.Popen('NUMBA=$(grep "AAA" param_gfn1-xtb.txt -o | wc -l);TEST=($(cat LAST_parameter.txt)); cp param_gfn1-xtb.txt LAST_param_gfn1-xtb.txt; for ((i = 1; i <= NUMBA; i++)); do loop_counter=$((i - 1)); sed -i "s/AAA$i /${TEST[$loop_counter]} /g" LAST_param_gfn1-xtb.txt;done', shell = True)
        process.wait()
        with open('LAST_parameter.txt', 'w') as file:
          file.write(' '.join(map(str, parr*initial0)))
          file.write("\n")
          file.write(' '.join(map(str, diff_step)))
        #print("WTF: " , parr_cost_function(parr*initial0),flush = True)
        friction = friction*4
  
    return (initial0,og_cost)

### HERE WE INITIATE ALL ####
if Qparallel > QSCcpu:
  og_cost = parr_cost_function(initial0)
else:
  og_cost = cost_function(initial0)
final, og_cost = lets_optimize(parr,initial0,diff_step,og_cost)
print('DONE Cost function: ',og_cost, " kcal/mol/atom ", flush = True)

#def test(x):
#  if Qparallel > QSCcpu:
#    return parr_cost_function(x)
#  else:
#    return cost_function(x)

#f_jac_hes.lastx = np.nan*parr
#f_jac = lambda x : f_jac_hes(x)[0]
#f_hes = lambda x : f_jac_hes(x)[1]
#model = sp.optimize.minimize(test, x0=parr, jac=f_jac, hess=f_hes, method='Newton-CG',options={'ftol':1e-4,'gtol':1e-5,'eps':1e-4})
#print(model)
