from joblib import Parallel, delayed
from os import environ
import sys

print("")
print("################")
print("### ArbAlign ###")
print("################")
print("")
missing = float("nan")

try:
    num_cores = int(environ['SLURM_JOB_CPUS_PER_NODE'])
except:
    from multiprocessing import cpu_count
    num_cores = cpu_count()

import pandas as pd
clusters_df = pd.read_pickle("db_for_rmsd.pkl")

original_length = len(clusters_df)
myNaN = lambda x : missing if x == "NaN" else x

if len(sys.argv) == 2:
  def compare_pair(arg):
    from ArbAlign import compare
    tobecompared = comparepairs[arg]
    AAci = clusters_df.loc[allindexes[tobecompared[0]],("xyz","structure")]
    AAcj = clusters_df.loc[allindexes[tobecompared[1]],("xyz","structure")]
    return compare(AAci,AAcj)
else:
  def compare_pair(arg):
    from qmllib2.representations import generate_fchl18 as gen_repr
    #from qmllib2.representations.fchl import get_local_symmetric_kernels as JKML_symm_kernel    
    from qmllib2.representations.fchl import get_local_kernels as JKML_kernel    
    from qmllib2.solvers import cho_solve
    from numpy import array, eye, dot, log, exp
    tobecompared = comparepairs[arg]
    AAci = clusters_df.loc[allindexes[tobecompared[0]],("xyz","structure")]
    AAcj = clusters_df.loc[allindexes[tobecompared[1]],("xyz","structure")]
    max_atoms = len(AAci.get_atomic_numbers())
    reprA = gen_repr(AAci.get_atomic_numbers(), AAci.get_positions(), max_size = max_atoms, neighbors = max_atoms, cut_distance=5.0)
    reprB = gen_repr(AAcj.get_atomic_numbers(), AAcj.get_positions(), max_size = max_atoms, neighbors = max_atoms, cut_distance=5.0)
    K = JKML_kernel(array([reprA]), array([reprB]), kernel_args = {"sigma": [2.0]})[0][0][0]
    #K0 = JKML_kernel(array([reprA]), array([reprA]), kernel_args = {"sigma": [1.0]})[0][0][0]
    return -log(K)*0.5*2.0**2*10**10
    K = JKML_symm_kernel(array([reprA]), kernel_args = {"sigma": array([1.0])})[0]
    K = K + 1e-7*eye(len(K)) 
    alpha = cho_solve(K, array([1.0]))  
    K = JKML_kernel(array([reprA]), array([reprB]), kernel_args = {"sigma": [1.0]})[0]
    res = (1e4*(dot(K, alpha)[0] - 1))**2
    return res

allindexes = clusters_df.index
comparepairs = []
for AAi in range(len(allindexes)):
  for AAj in range(AAi+1,len(allindexes)):
    pair = [AAi,AAj]
    comparepairs.append(pair)
print("STRUCTURES")
for AAi in range(len(allindexes)):
  print(str(AAi)+"="+clusters_df.loc[allindexes[AAi],("info","file_basename")])
print("")
print("PAIR\tRMSD[Angstrom]")
comparison = Parallel(n_jobs=num_cores)(delayed(compare_pair)(i) for i in range(len(comparepairs)))
for i in range(len(comparepairs)):
  print(str(comparepairs[i])+"\t"+str(comparison[i]))

