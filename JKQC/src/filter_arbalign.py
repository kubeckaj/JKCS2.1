def filter_arbalign(clusters_df,Qclustername,Qarbalign,Qout):
  from joblib import Parallel, delayed
  from os import environ

  try:
    num_cores = int(environ['SLURM_JOB_CPUS_PER_NODE'])
  except:
    from multiprocessing import cpu_count
    num_cores = cpu_count()

  if Qclustername != 0:
    from numpy import unique
    uniqueclusters = unique(clusters_df.loc[:,("info","cluster_type")])
  else:
    uniqueclusters = "1"

  original_length = len(clusters_df)
  myNaN = lambda x : missing if x == "NaN" else x
  global preselected_df,comparepairs,allindexes

  def compare_pair(arg):
    from ArbAlign import compare
    tobecompared = comparepairs[arg]
    AAci = preselected_df.loc[allindexes[tobecompared[0]],("xyz","structure")]
    AAcj = preselected_df.loc[allindexes[tobecompared[1]],("xyz","structure")]
    return compare(AAci,AAcj)

  for i in uniqueclusters:
     if Qclustername != 0:
       preselected_df = clusters_df.loc[clusters_df.loc[:,("info","cluster_type")] == i]
     else:
       preselected_df = clusters_df
     allindexes = preselected_df.index
     removedindexes = []
     for AAi in range(len(allindexes)):
       if allindexes[AAi] in removedindexes:
         continue
       comparepairs = []
       for AAj in range(AAi+1,len(allindexes)):
         if allindexes[AAj] in removedindexes:
           continue
         pair = [AAi,AAj]
         comparepairs.append(pair)
       comparison = Parallel(n_jobs=num_cores)(delayed(compare_pair)(i) for i in range(len(comparepairs)))
       for AAc in range(len(comparison)):
         if comparison[AAc] < Qarbalign:
           removedindexes.append(allindexes[comparepairs[AAc][1]])
     clusters_df = clusters_df.drop(removedindexes)

  if Qout >= 1:
    print("ArbAlign: "+str(original_length)+" --> "+str(len(clusters_df)))

  return clusters_df


