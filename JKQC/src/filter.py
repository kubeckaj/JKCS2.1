def filter(clusters_df, Qsort, Qarbalign, Quniq, Qsample, Qclustername, Qthreshold, Qcut, Qshuffle, Qselect, Qreacted, bonddistancethreshold, Qout):

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
    Qsort = "g"

  #CUTOFF FILTERING
  if Qthreshold != 0:
    from filter_threshold import filter_threshold
    clusters_df = filter_threshold(clusters_df,Qcut,Qclustername,Qout)

  #UNIQUENESS AND SAMPLING
  if str(Quniq) != "0":
    from filter_uniq import filter_uniq
    cluster_df = filter_uniq(clusters_df,Quniq,Qclustername,Qsample,Qout)

  ### REACTED ###
  if Qreacted > 0:
    from filter_reacted import filter_reacted
    clusters_df = filter_reacted(clusters_df,Qclustername,Qreacted,bonddistancethreshold,Qout)
    if Qout == 2:
      print("DONE] Reacted done: "+str(time.time() - start));

  #SORTING
  if str(Qsort) != "0" and str(Qsort) != "no":
    from filter_sort import filter_sort
    clusters_df = filter_sort(clusters_df,Qsort)

  #ArbAlign
  if Qarbalign > 0:
    from filter_arbalign import filter_arbalign
    clusters_df = filter_arbalign(clusters_df,Qclustername,Qarbalign,Qout) 
  
  #TODO: not sure if this sorting is necessary but maybe after uniqueness filtering yes
  #SORTING (AGAIN)
  if str(Qsort) != "0" and str(Qsort) != "no":
    clusters_df = filter_sort(clusters_df,Qsort)

  #SELECT
  if str(Qselect) != "0":
    from filter_select import filter_select
    clusters_df = filter_select(clusters_df,Qselect,Qclustername,Qout)

  ### SHUFFLE
  if Qshuffle == 1:
    clusters_df = clusters_df.sample(frac=1)
  
  return clusters_df
