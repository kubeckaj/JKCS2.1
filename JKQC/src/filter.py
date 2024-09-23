def filter(clusters_df, Qsort, Qreverse, Qarbalign, QMWarbalign, Quniq, Qsample, Qclustername, Qthreshold, Qcut, Qshuffle, Qselect, Qreacted, bonddistancethreshold, Qout, seed):

  #CHECK WHETHER CLUSTERNAME MAKE SENSE
  if Qclustername == 1:
    if ("info","cluster_type") in clusters_df.columns:
      from pandas import isna
      if isna(clusters_df.loc[:,("info","cluster_type")].values).any():
        #if Qout >= 1:
        #  print("Cluster type with nan present (weird file names, use -noname) -> Qclustername = 0")
        Qclustername = 0
    else:  
      Qclustername = 0

  #IS SORTING NEEDED?
  if str(Qsort) == "0" and str(Qarbalign) != "0":
    if ("log","gibbs_free_energy") in clusters_df.columns:
      if ("out","electronic_energy") in clusters_df.columns:
        Qsort = "gout"
      else:
        Qsort = "g"
    else:
      Qsort = "el"
  if ( str(Qselect) != "0" or ( str(Quniq) != "0" and str(Quniq) != "dup" )) and str(Qsort) == "0":
    if ("log","gibbs_free_energy") in clusters_df.columns:
      if ("out","electronic_energy") in clusters_df.columns:
        Qsort = "gout"
      else:
        Qsort = "g"
    else:
      Qsort = "el"

  #CUTOFF FILTERING
  if Qthreshold != 0:
    from filter_threshold import filter_threshold
    clusters_df = filter_threshold(clusters_df,Qcut,Qclustername,Qout)

  #UNIQUENESS AND SAMPLING
  if str(Quniq) != "0":
    from filter_uniq import filter_uniq
    clusters_df = filter_uniq(clusters_df,Quniq,Qclustername,Qsample,Qout)

  ### REACTED ###
  if Qreacted > 0:
    from filter_reacted import filter_reacted
    clusters_df = filter_reacted(clusters_df,Qclustername,Qreacted,bonddistancethreshold,Qout)

  #SORTING
  if str(Qsort) != "0" and str(Qsort) != "no":
    from filter_sort import filter_sort
    clusters_df = filter_sort(clusters_df,Qsort,Qreverse)

  #ArbAlign
  if Qarbalign > 0 or QMWarbalign > 0:
    from filter_arbalign import filter_arbalign
    clusters_df = filter_arbalign(clusters_df,Qclustername,Qarbalign,QMWarbalign,Qout) 
  
  #SORTING (AGAIN)
  if str(Qsort) != "0" and str(Qsort) != "no" and Qarbalign > 0:
    clusters_df = filter_sort(clusters_df,Qsort,Qreverse)

  #SELECT
  if str(Qselect) != "0":
    from filter_select import filter_select
    clusters_df = filter_select(clusters_df,Qselect,Qclustername,Qout)

  ### SHUFFLE
  if Qshuffle == 1:
    clusters_df = clusters_df.sample(frac=1, random_state=seed)
  
  return clusters_df
