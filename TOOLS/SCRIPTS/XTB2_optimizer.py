#import scipy as sp
import os
import subprocess
import numpy as np

def cost_function(arr):
    with open('parameter.txt', 'w') as file:
        file.write(' '.join(map(str, arr))) 
    process = subprocess.Popen(['bash', 'XTB3_runXTB.sh'])
    process.wait()
    with open('result', 'r') as file:
            number_str = file.readline().strip()
            return float(number_str)

with open('initial.txt', 'r') as file:
  line = file.readline()
  initial0 = [float(i) for i in line.split()]

#model = sp.optimize.minimize(test, x0=initial, method='L-BFGS-B',options={'ftol':1e-4,'gtol':1e-5,'eps':1e-4})
#print(model)

diff_step = 1e-2
hess_step = diff_step
max_iter = 10
initial = len(initial0)*[1]


#JK optimizer
def numeric_gradient(parr,diff_step,og_cost):
    return np.array([cost_function(generate_step(initial0*parr,idx,diff_step))-og_cost for idx in range(len(parr))])

def generate_step(parr,idx,diff_step):
    new_parr = parr.copy()
    new_parr[idx] = new_parr[idx]+diff_step
    return np.array(new_parr)

def OLD_numeric_hessian(gradient_norm,hess_step,parr,diff_step,og_cost):
    cost_stepped = cost_function(initial0*(parr + diff_step * gradient_norm))
    cost_hess_stepped = cost_function(initial0*(parr + (hess_step + diff_step) * gradient_norm))
    cost_hess = cost_function(initial0*(parr + hess_step * gradient_norm))
    if og_cost > cost_stepped or og_cost > cost_hess or og_cost > cost_hess_stepped or cost_stepped > cost_hess or cost_stepped > cost_hess_stepped or cost_hess > cost_hess_stepped:
      print("\nDONE - too large eps or omega", flush = True)
      print([og_cost,cost_stepped,cost_hess,cost_hess_stepped], flush = True)
      return (parr,og_cost)
    if np.abs(cost_hess_stepped - cost_hess - cost_stepped + og_cost) < 1e-99:
      print("\nDONE - minimum found", flush = True)
      return hess_step 
    return ((cost_stepped - og_cost) * hess_step) / (cost_hess_stepped - cost_hess - cost_stepped + og_cost)

def numeric_hessian(gradient_norm,parr,diff_step,og_cost):
    cost_forward = cost_function(initial0*(parr + diff_step * gradient_norm))
    cost_backward = cost_function(initial0*(parr - diff_step * gradient_norm))
    if cost_forward < og_cost or og_cost < cost_backward or cost_forward < cost_backward:
      print("\nDONE - too large eps or omega", flush = True)
      print([cost_backward,og_cost,cost_forward], flush = True)
      return diff_step
    if np.abs(cost_forward-2*og_cost+cost_backward) < 1e-99:
      print("\nDONE - minimum found", flush = True)
      return diff_step
    return ((cost_forward - og_cost) * diff_step) / (cost_forward-2*og_cost+cost_backward)

def lets_optimize(parr,diff_step,og_cost,hess_step):
    print(f">>> Optimizing for input {initial0}", flush = True)
    parr = np.array(parr)
    print('Original: ',og_cost, " kcal/mol/atom", flush = True)
    rem=og_cost
    for iter in range(max_iter):
        friction = 0.99
        print("####################################")
        print(f'>>> Iteration {iter}', flush = True)
        gradient = numeric_gradient(parr,diff_step,og_cost)
        gradient_norm = gradient/np.linalg.norm(gradient) 
        print('\n-Gradient:',-gradient_norm, flush = True)
        move = numeric_hessian(gradient_norm,parr,diff_step,og_cost)
        if np.abs(move) > 0.1:
          move = 0.1*move/np.abs(move)
        print('\nMove :', move, flush = True)
        next_test = 0
        while next_test < 2:
          #if np.abs(friction*move) < diff_step:
          #  print("DONE -- too short step.", flush = True)
          #  return (parr,og_cost)
          if next_test == 0:
            parr_test = parr-friction*move*gradient_norm
            og_cost_test = cost_function(parr_test*initial0)
            print('\nTest :', og_cost_test, flush = True)
          else:
            parr_test_EXT = parr-friction*move*gradient_norm
            og_cost_test_EXT = cost_function(parr_test*initial0)
            print('\nTest - ext :', og_cost_test_EXT, flush = True)
            if og_cost_test_EXT < og_cost_test:
              par_test = parr_test_EXT
              og_cost_test = og_cost_test_EXT
              next_test = 1
              friction = 0.66*friction
              continue
            else: 
              next_test = 2
              continue
          if og_cost_test < og_cost:
            next_test = 1
            friction = 0.66*friction
          else:
            friction = 0.66*friction
          if friction < 0.1:
            print("DONE -- too much break used", flush = True)
            return (parr,og_cost)
        parr = parr_test
        if np.abs(og_cost_test-og_cost)*2/(og_cost_test+og_cost) < 0.0001:
          print("DONE -- too small cost function improvements", flush = True)
          return (parr,og_cost)
        og_cost = og_cost_test
        
        print('New Parameters:',parr*initial0, flush = True)  
        print('Cost function: ',og_cost, " kcal/mol/atom (original was ",rem,"kcal/mol/atom)", flush = True)
        print(parr*initial0, og_cost, flush = True)
    return (parr,og_cost)


og_cost = cost_function(initial0)
final, og_cost = lets_optimize(initial,diff_step,og_cost,hess_step)
print('DONE Cost function: ',og_cost, " kcal/mol/atom ", flush = True)


