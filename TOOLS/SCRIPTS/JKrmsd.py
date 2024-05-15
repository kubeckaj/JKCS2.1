from joblib import Parallel, delayed
from os import environ

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

def compare_pair(arg):
    from ArbAlign import compare
    tobecompared = comparepairs[arg]
    AAci = clusters_df.loc[allindexes[tobecompared[0]],("xyz","structure")]
    AAcj = clusters_df.loc[allindexes[tobecompared[1]],("xyz","structure")]
    return compare(AAci,AAcj)

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

