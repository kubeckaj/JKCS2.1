# IMPORT NECESSARY LIBRARIES
from sys import argv
from os import system

#print command to the output file
cmd="".join(["( echo COMMAND: JKQC "," ".join(argv[1:])," >> output ) 2>/dev/null"])
system(cmd)

#READ ARGUMENTS
from arguments import arguments
locals().update(arguments(argv[1:]))
if Qout == 2:
  from time import time
  start = time()
  print("DONE] Time started: "+str(time() - start));

# Loading input pickles
from load_pickles import load_pickles
clusters_df = load_pickles(input_pkl,Qout,Qid)
if Qout == 2:
  print("DONE] Pickles loading done: "+str(time() - start));

# Loading complement
if Qcomplement != 0:
  from load_complement import load_complement
  clusters_df = load_complement(clusters_df, Qcomplement)
  if Qout == 2:
    print("DONE] Complement done: "+str(time() - start));

# Loading addSP = Adding single-point correction
if len(input_pkl_sp) > 0:
  from load_addsp import load_addsp
  clusters_df = load_addsp(clusters_df,input_pkl_sp,Qout)
  if Qout == 2:
    print("DONE] AddSP done: "+str(time() - start));

#READ QC FILES
if len(files) > 0:
  from read_files import read_files
  clusters_df = read_files(clusters_df, files, orcaextname, orcaext, turbomoleext, Qclustername, Qforces, Qanharm, Qdisp_electronic_energy, Qdisp_forces)
  if Qout == 2:
    print("DONE] Files loaded: "+str(time() - start));

#ADD EXTRA COLUMN
if len(addcolumn) > 0:
  from add_column import add_column
  clusters_df = add_column(clusters_df,addcolumn)
  if Qout == 2:
    print("DONE] Column(s) added: "+str(time() - start));
    
## MODIFY
if Qmodify > 0:
  from data_modification import data_modification
  clusters_df = data_modification(clusters_df,Qunderscore,Qrename,Qclustername,QrenameWHAT,Qiamifo,Qrebasename,Qdrop,Qout2log,Qpresplit,Qindex,seed,Qatomize)
  if Qout == 2:
    print("DONE] Database modified: "+str(time() - start));

## EXTRACT
if Qextract > 0:
  from extract_clusters import extract_clusters
  clusters_df = extract_clusters(clusters_df,Qextract,Pextract,Qclustername,Qout)
  if Qout == 2:
    print("DONE] Extraction done: "+str(time() - start));

## IN ORDER TO SAVE OUTPUT.pkl ##
clusters_df = clusters_df.reset_index(drop=True)
if Qoutpkl > 0:
  original_clusters_df = clusters_df.copy()
  if Qout == 2:
    print("DONE] Original copy done: "+str(time() - start));

## THERMODYNAMICS ##
if Qqha == 1:
  from thermodynamics import thermodynamics
  clusters_df = thermodynamics(clusters_df, Qanh, Qafc, Qfc, Qt, Qdropimg)
  if Qout == 2:
    print("DONE] Data modification done: "+str(time() - start));

## FILTERING ##
from filter import filter
clusters_df = filter(clusters_df, Qsort, Qreverse, Qarbalign, QMWarbalign, Quniq, Qsample, Qclustername, Qthreshold, Qcut, Qshuffle, Qselect, Qreacted,bonddistancethreshold, Qout, seed)

## SAVE OUTPUT.pkl ##
if Qoutpkl > 0:
  from save_pickle import save_pickle
  save_pickle(original_clusters_df.loc[clusters_df.index],output_pkl,Qsplit,Qout) 

## PREPARE DATA PRINT ##
from numpy import array
if Qsolvation != "0" or Qformation != 0:
  if len(Pout) > 0:
    if Pout[0] != "-ct":
      Pout.insert(0,"-ct")
if len(Pout) > 0:
  from print_output import print_output
  output = array(print_output(clusters_df,Qoutpkl,input_pkl,output_pkl,Qsplit,Qclustername,Qt,Qcolumn,Qbonded,Qdistances,Pout,QUenergy,QUentropy))
  if Qout == 2:
    print("DONE] Printing prepared: "+str(time() - start));
else:
  output = array([])
 
## PRINT DATA ##
if not len(output) == 0:
  #TAKING GLOBAL MINIMA ONLY: not needed if sort and select used
  if (Qglob == 1 or Qglob == 2) and (len(clusters_df)>1):
    from take_glob import take_glob
    output = take_glob(output, clusters_df, Qglob)
  
  #TAKING BOLTZMANN AVERAGE OVER ALL MINIMA 
  if Qbavg == 1 or Qbavg == 2:
    from math import isnan
    if isnan(Qt):
      Qt = 298.15
    from take_bavg import take_bavg
    output = take_bavg(output, clusters_df, Qbavg, Qt, QUenergy, Pout, Qclustername)

  #PRINTING (if required)
  if not(Qout < 1 and (Qformation == 1 or Qsolvation != "0")):
    toprint = list(zip(*output)) #[row for row in list(zip(*output))]
    if len(toprint) > 0:
      column_widths = [max(len(str(row[i])) for row in toprint) for i in range(len(toprint[0]))]
      for row in toprint:
        formatted_row = [str(row[i]).ljust(column_widths[i]) for i in range(len(row))]
        print(" ".join(formatted_row),flush = True)
    if Qout == 2:
      print("DONE] Printing done: "+str(time() - start));

#TXT FILE + SOLVATION + FORMATION     
if Qsolvation != "0" or Qformation == 1:
  #LOAD INPUT FILE
  if len(formation_input_file) > 0:
    from load_txt import load_txt
    output = load_txt(output,formation_input_file,QUentropy,QUenergy,Pout)

  #SOLVATION PROPERTIES
  if Qsolvation != "0":
    from print_solvation import print_solvation
    print_solvation(output,Qsolvation,Qt,Qrh,QPsolvent,Pout,QUenergy,Qout)
    if Qout == 2:
      print("DONE] Solvation done: "+str(time() - start));
      
  #FORMATION PROPERTIES
  if Qformation == 1:
    from math import isnan
    from print_formation import print_formation
    if isnan(Qt):
      Qt = 298.15
    if isnan(Qp):
      Qp = 101325
    print_formation(output,Qout,Qt,Qp,Qconc,conc,CNTfactor,QUenergy)
    if Qout == 2:
      print("DONE] Formation done: "+str(time() - start));

# Run AIMNet preparation if requested
if Qaimnet_prep == 1:
  from aimnet_prep import convert_pickle_to_aimnet_h5
  convert_pickle_to_aimnet_h5(clusters_df, output_file="dataset.h5")
  if Qout == 2:
    print("DONE] AIMNet prep done: " + str(time() - start))

if Qout == 2:
  print("DONE] JKQC done: "+str(time() - start));
