import argparse as ap
import numpy as np
import os
import subprocess
import time
import shutil
import copy

################## CALLING COST FUNCTION  ####################
def cost_function(guesses,cost_function_file,QSCcpu,QSCparallel):
  #PREPARE FOLDERS
  for guess_i in range(len(guesses)):
    os.makedirs("TRIAL_"+str(guess_i))
    os.chdir("TRIAL_"+str(guess_i))
    with open('parameter.txt', 'w') as file:
      print(' '.join(map(str, guesses[guess_i])), flush=True)
      file.write(' '.join(map(str, guesses[guess_i])))
    os.chdir('../')

  #RUN/SUBMIT JOBS
  if QSCparallel == 0:
    for guess_i in range(len(guesses)):
      os.chdir("TRIAL_"+str(guess_i))
      process = subprocess.Popen('sh '+cost_function_file+' '+str(QSCcpu), shell=True) 
      process.wait()
      os.chdir('../')
  else:
    process = subprocess.Popen('sbatch --time 2:00:00 -n '+str(QSCparallel)+' --array=0-'+str(len(guesses)-1)+' JKsend "cd TRIAL_\$SLURM_ARRAY_TASK_ID; sh '+cost_function_file+' '+str(QSCparallel)+'"', shell = True)
    process.wait()
 
  #WAITING/CHECKING JOBS TO BE FINISHED
  for guess_i in range(len(guesses)-1,-1,-1):
    while not os.path.exists('TRIAL_'+str(guess_i)+'/result'):  
      time.sleep(5)
  time.sleep(1)

  #READING RESULTS
  results = []
  for guess_i in range(len(guesses)):
    with open('TRIAL_'+str(guess_i)+'/result', 'r') as file:
      try:
        result = float(file.readline().strip())
      except:
        result = np.nan
    results.append(result)
    shutil.rmtree('TRIAL_'+str(guess_i))
  
  return np.array(results)

################ APPLY OPERATIONS #######################
def operate(num1, num2, operation, how_many_times = 1, check = 0):
  if operation == "+":
    if check == 1:
      if np.sign(num1) != np.sign(num2) and  np.abs(num1) < np.abs(num2*how_many_times):
        #TODO: if modification does not fulfil criteria, I will return nan!!!
        return np.nan
    if num1 == 0:
      return num2*how_many_times
    else:
      return np.sign(num1)*(np.abs(num1)+num2*how_many_times)
  elif operation == "*":
    if num2 > 0:
      return num1*num2**how_many_times
    else:
      return num1/np.abs(num2)**how_many_times
  else:
    print("Unknown type of operations")
    exit()

############### OPERATE ON ARRAYS ###############
def array_operate(array1,array2,operations,how_many_times = 1,check = 0):
  return np.array([operate(array1[i],array2[i],operations[i],how_many_times,check) for i in range(len(array1))])

############ GENERATE SPECIFIC HESS STEPS ###############
def generate_step(var, var_mod, var_mod_operation, idx, this_ig_var_operation):
  new_var = copy.deepcopy(var)
  if this_ig_var_operation == "+":
    new_var[idx] = operate(new_var[idx], var_mod, var_mod_operation, how_many_times = 1, check = 0)
  else:
    new_var[idx] = operate(new_var[idx], var_mod, var_mod_operation, how_many_times = 1, check = 1)
  return np.array(new_var)

########### OPTIMIZATION DIRECTIONAL HESS MOVE ##########
def mover(initial_guess, var, ig_var_operation, var_mod, var_mod_operation, mod_mod, mod_mod_operation, og_cost, cost_function_file, cpu, parallel):
  #positions of new variables to be tested
  all_positions = [array_operate(initial_guess,generate_step(var, where*var_mod[idx], var_mod_operation[idx], idx, ig_var_operation[idx]),ig_var_operation) for where in [1,-1] for idx in range(len(var))]
  cost_all = cost_function(all_positions, cost_function_file, cpu, parallel)
  print("", flush = True)
  
  #SEPARATE AND PRECALCULATE SOMETHING
  cost_forward = np.array(cost_all)[0:len(var)]
  cost_backward = np.array(cost_all)[len(var):2*len(var)]
  #TODO var_mod_operation will change this for sure FUCKING TOUGH TO DO PROPERLY
  den = (cost_forward-2*og_cost+cost_backward) #denominator
  res = []
  #INTIATE BESTS
  best_found = og_cost
  best_var = copy.deepcopy(var)

  #LOOP OVER ALL VARIABLES AND MODIFY
  for i in range(len(var)):
    #CHECK FOR BEST
    if best_found > cost_backward[i]:
      best_found = np.array(cost_all)[i+len(var)]
      best_var = copy.deepcopy(generate_step(var, -var_mod[i], var_mod_operation[i], i, ig_var_operation[i]))
    if best_found > cost_forward[i]:
      best_found = np.array(cost_all)[i]
      best_var = copy.deepcopy(generate_step(var, var_mod[i], var_mod_operation[i], i, ig_var_operation[i]))
    #TEST ALL VARIANTS
    if cost_forward[i] == og_cost or cost_backward[i] == og_cost or np.abs(den[i]) < 1e-99 or np.isnan(den[i]):
      if cost_forward[i] < og_cost:
        res.append(-var_mod[i])
      elif cost_backward[i] < og_cost:
        res.append(var_mod[i])
      elif np.isnan(den[i]):
        if var_mod_operation[i] == "+":
          res.append(0)
        else:
          res.append(1)
        var_mod[i] = operate(var_mod[i],-mod_mod[i],mod_mod_operation[i],1,1)
        var_mod[i] = operate(var_mod[i],-mod_mod[i],mod_mod_operation[i],1,1)
      else:
        if var_mod_operation[i] == "+":
          res.append(0)
        else:
          res.append(1)
      var_mod[i] = operate(var_mod[i],mod_mod[i],mod_mod_operation[i],1,1)
    else:
      #TODO var_mod_operation will change this for sure FUCKING TOUGH TO DO PROPERLY
      if var_mod_operation[i] == "+" and ((mod_mod_operation[i] == "+" and mod_mod[i] != 0) or (mod_mod_operation[i]) == "*" and mod_mod[i] != 1):
        premove = ((cost_forward[i] - cost_backward[i]) / 2 * var_mod[i]) / den[i]
      else:
        if cost_forward[i] < cost_backward[i]:
          premove = -var_mod[i]
        else:
          premove = var_mod[i]
      if cost_backward[i] < og_cost and cost_forward[i] < og_cost:
        if cost_backward[i] < cost_forward[i]:
          premove = var_mod[i]
        else:
          premove = -var_mod[i]
        var_mod[i] = operate(var_mod[i],-mod_mod[i],mod_mod_operation[i],1,1)
      elif cost_backward[i] > og_cost and cost_forward[i] > og_cost:
        var_mod[i] = operate(var_mod[i],-mod_mod[i],mod_mod_operation[i],1,1)
      elif cost_forward[i] < og_cost:
        if np.abs(premove) > np.abs(var_mod[i]):
          premove = -var_mod[i]
          var_mod[i] = operate(var_mod[i],mod_mod[i],mod_mod_operation[i],1,1)
      elif cost_backward[i] < og_cost:
        if np.abs(premove) > np.abs(var_mod[i]):
          premove = var_mod[i]
          var_mod[i] = operate(var_mod[i],mod_mod[i],mod_mod_operation[i],1,1)
      else:
        print("Oh Jakub probably forgot some condition, tell him you got this message.")

      #THIS SHOULD NOT BE NECESSARY
      #if np.abs(premove) > diff_step[i]:
      #  premove = np.sign(premove)*diff_step[i]
      #  diff_step[i] = diff_step[i]*2
      res.append(premove)
    print("VAR "+str(i)+":",[cost_backward[i],og_cost,cost_forward[i], initial_guess[i], ig_var_operation[i],var[i],var_mod_operation[i],var_mod[i],mod_mod_operation[i],mod_mod[i],den[i],-res[i]],flush = True)
  #den = den/np.abs(den)*(np.abs(den)+diff_step)
  return np.array(res), var_mod, best_found, best_var

########### THE ACTUAL MOVING TESTS ##########
def moving(move, initial_guess, var, ig_var_operation, var_mod, var_mod_operation, mod_mod, mod_mod_operation, og_cost, cost_function_file, cpu, parallel):
  
  testing = "FIRST TESTING"
  friction_applied = 0
  if len(var) > 10:
    in_parallel = 10
  else: 
    in_parallel = len(var)
  prev_test = np.nan
  prev_var = np.nan
  best_test_found = og_cost
  best_test_var = copy.deepcopy(var)
  while_count = 0
  while testing != "END":
    while_count += 1
    if while_count*(in_parallel+1) > 50:
      print("Too much friction")
      testing = "END"
    all_positions = [array_operate(initial_guess,np.array([operate(var[j], array_operate(-move,-mod_mod,mod_mod_operation,friction_applied+i,1)[j],var_mod_operation[j],1, 1 if ig_var_operation[j]=="*" else 0) for j in range(len(var))]),ig_var_operation) for i in range(in_parallel)]
    cost_all = cost_function(all_positions, cost_function_file, cpu, parallel)
    ## TOO SHORT STEP SHOULD BE IRRELEVANT
    #if np.linalg.norm(friction*move) < 1e-6: #diff_step/10:
    #        print("DONE -- too short step.", flush = True)
    #        parr_test = parr
    #        og_cost_test = og_cost
    #        break
    print("")
    for test_i in range(len(cost_all)):
      print("Test - ", cost_all[test_i], flush = True)
      if np.isnan(cost_all[test_i]) and not np.isnan(prev_test) and testing != "END":
        testing = "END"
        print("WTF THIS SHOULD NOT HAPPEN. CONTACT JAKUB")
        continue
      if cost_all[test_i] < best_test_found:
        best_test_found = cost_all[test_i]
        best_test_var = copy.deepcopy(np.array([operate(var[j], array_operate(-move,-mod_mod,mod_mod_operation,friction_applied+test_i,1)[j],var_mod_operation[j],1, 1 if ig_var_operation[j]=="*" else 0) for j in range(len(var))]))
      if not np.isnan(prev_test):
        if np.abs(cost_all[test_i] - prev_test) <= 1e-7:
          print("WARNING -- too small imrovement.", flush = True)
          testing = "END"
        if prev_test < cost_all[test_i]:
          testing = "END"
      prev_test = cost_all[test_i]
      prev_var = copy.deepcopy(np.array([operate(var[j], array_operate(-move,-mod_mod,mod_mod_operation,friction_applied+test_i,1)[j],var_mod_operation[j],1, 1 if ig_var_operation[j]=="*" else 0) for j in range(len(var))]))
    friction_applied=friction_applied+in_parallel
  print("Test : best_test_found", best_test_found,  flush = True)
  return best_test_found, best_test_var

######################## SETUP ########################
def setup(input_file):
  """Read or eventually prepare initial setup"""

  # Read lines
  try:
    with open(input_file, 'r') as file:
      lines = file.readlines()
  except:
    lines = np.array([])
 
  # Remove any trailing whitespace from each line
  lines = [line.strip() for line in lines]
   
  # Convert the first line to a numpy array
  line = 1
  if len(lines) > 0:
    initial_guess = np.array([float(i) for i in lines[line-1].split()])
  else:
    initial_guess = np.array([1.0])

  line = 2
  if len(lines) > line-1 and lines[line-1]:
    var = np.array([float(i) for i in lines[line-1].split()])  
    if len(var) != len(initial_guess):
      print("Wrong input")
      exit()
  else:
    var = np.full_like(initial_guess, 1.0, dtype=np.double)
  
  line = 3
  if len(lines) > line-1 and lines[line-1]:
    ig_var_operation = np.array([str(i) for i in lines[line-1].split()])
    if len(ig_var_operation) != len(initial_guess):
      print("Wrong input")
      exit()
  else:
    ig_var_operation = np.full_like(initial_guess, "*", dtype=str)

  line = 4
  if len(lines) > line-1 and lines[line-1]:
    var_mod = np.array([float(i) for i in lines[line-1].split()])
    if len(var_mod) != len(initial_guess):
      print("Wrong input")
      exit()
  else:
    var_mod = np.full_like(initial_guess, 0.001, dtype=np.double)  

  line = 5
  if len(lines) > line-1 and lines[line-1]:
    var_mod_operation = np.array([str(i) for i in lines[line-1].split()])
    if len(var_mod_operation) != len(initial_guess):
      print("Wrong input")
      exit()
  else:
    var_mod_operation = np.full_like(initial_guess, "+", dtype=str) 

  line = 6
  if len(lines) > line-1 and lines[line-1]:
    mod_mod = np.array([float(i) for i in lines[line-1].split()])
    if len(mod_mod) != len(initial_guess):
      print("Wrong input")
      exit()
  else:
    mod_mod = np.full_like(initial_guess, 2.0, dtype=np.double)

  line = 7
  if len(lines) > line-1 and lines[line-1]:
    mod_mod_operation = np.array([str(i) for i in lines[line-1].split()])
    if len(mod_mod_operation) != len(initial_guess):
      print("Wrong input")
      exit()
  else:
    mod_mod_operation = np.full_like(initial_guess, "*", dtype=str)
   
  return initial_guess, var, ig_var_operation, var_mod, var_mod_operation, mod_mod, mod_mod_operation 

######################### MAIN #########################
def main():
  """Main part of the program"""

  ### READING ARGUMENTS ###
  
  # Define an argument parser object
  parser = ap.ArgumentParser()
  
  # Define the arguments to add
  parser.add_argument("-i", "--input", type=str, required=False, default="initial.txt", help="Input file w max X lines: initial guess, jac/hess step")
  parser.add_argument("-cff", "--cost_function_file", type=str, required=False, default="../XTB3_runXTB.sh", help="Cost function file full path [bash file assumed]")
  parser.add_argument("-cpu", type=int, required=False, default=1, help="Number of CPUs this script can use.")
  parser.add_argument("-parallel", type=int, required=False, default=0, help="Number of CPUs the submitted script can use (def = 0).")
  parser.add_argument("-maxiter", type=int, required=False, default=100, help="Maximum number of iterations (def = 100).")
  args = parser.parse_args()

  ### CHECK INPUT ###
  #TODO

  ### CALLING SETUP ###
  # variable.operation.initial  X  modifier -> variable
  # initial_guess - initial input
  # var - this is what we use to adjust the initial guess with
  # ig_var_operation - how does var adjusts intial guess, e.g.: +,*
  # var_mod = with what we modify var
  # var_mod_operation = how we modify var with var_mod
  # mod_mod = with what we modify mod
  # mod_mod_operation = how we modify var with mod_mod
  initial_guess, var, ig_var_operation, var_mod, var_mod_operation, mod_mod, mod_mod_operation = setup(args.input) 
  #print(initial_guess, ig_var_operation, var, var_mod_operation, var_mod, mod_mod_operation, mod_mod)
  if np.any(ig_var_operation != "*"):
    print("!ig*var: I am not ready for this :-( [EASY TO IMPLEMENT]")
    #exit()
  if np.any(var_mod_operation != "+"):
    print("!var+mod: I am not ready for this :-( [HARD TO IMPLEMENT]")
    #exit()
  if np.any(mod_mod_operation != "*"):
    print("!mod*mod: I am not ready for this :-( [EASY TO IMPLEMENT]")
    #exit()
 
  ### OPTIMIZATION ###
  ## STARTING POINT ## 
  print(f">>> Optimizing for input {initial_guess}", flush = True)
  #print([[initial_guess], args.cost_function_file, args.cpu, args.parallel])
  #original cost function:
  og_cost = cost_function([initial_guess], args.cost_function_file, args.cpu, args.parallel)[0]
  print('Cost function: ',og_cost, " kcal/mol/atom (original)", flush = True)
  if np.isnan(og_cost):
    print("For the initial parameters, XTB returns NaN. Buuuuuuu [EXITING]")
    exit()
  rem_og_cost=og_cost

  ## ITERATIONS ##
  max_iter = args.maxiter
  for iter in range(max_iter):
    print("####################################")
    print(f'>>> Iteration {iter}', flush = True)
    move, var_mod, best_found, best_var = mover(initial_guess, var, ig_var_operation, var_mod, var_mod_operation, mod_mod, mod_mod_operation, og_cost, args.cost_function_file, args.cpu, args.parallel)
    #print('\nMove :', move, flush = True)
    print("Test : starting from", og_cost,  flush = True)
    print("Test : best found", best_found,  flush = True)
    best_test_found, best_test_var = moving(move, initial_guess, var, ig_var_operation, var_mod, var_mod_operation, mod_mod, mod_mod_operation, og_cost, args.cost_function_file, args.cpu, args.parallel)
    if best_test_found < best_found:
      best_found = best_test_found
      best_var =  copy.deepcopy(best_test_var)
    if best_found < og_cost:
      if np.abs(best_found-og_cost) < 0.0001:
        print("WARNING -- too small cost function improvements", flush = True)
      og_cost = best_found
      var =  copy.deepcopy(best_var)

      print('New Parameters:',array_operate(initial_guess,var,ig_var_operation), flush = True)
      print('Cost function: ',og_cost, " kcal/mol/atom (original was ",rem_og_cost,"kcal/mol/atom)", flush = True)
      #TODO not sure if necessary
      #initial_guess = var*initial_guess 
      #var= np.array(len(initial_guess)*[1.0])

      with open('LAST_parameter.txt', 'w') as file:
        file.write(' '.join(map(str, var*initial_guess))+"\n")
      process = subprocess.Popen('NUMBA=$(grep "AAA" param_gfn1-xtb.txt -o | wc -l);TEST=($(cat LAST_parameter.txt)); cp param_gfn1-xtb.txt LAST_param_gfn1-xtb.txt; for ((i = 1; i <= NUMBA; i++)); do loop_counter=$((i - 1)); sed -i "s/AAA$i /${TEST[$loop_counter]} /g" LAST_param_gfn1-xtb.txt;done', shell = True)
      process.wait()
      with open('LAST_parameter.txt', 'w') as file:
        file.write(' '.join(map(str, var*initial_guess)))
        file.write("\n")
        file.write(' '.join(map(str, var_mod)))
      #print("WTF: " , parr_cost_function(parr*initial0),flush = True)
      #TODO: do I want to keep memory on frictions
      #friction = friction*4 
    else:
      #TODO this will be modified by ig_var_operation
      print('Old Parameters:',array_operate(initial_guess,var,ig_var_operation), flush = True)
      print('Cost function: ',og_cost, " kcal/mol/atom (original was ",rem_og_cost,"kcal/mol/atom)", flush = True)
 
 
  return 0

######################### INIT #########################
if __name__ == "__main__":
  main()
