####################################################################################################
# Reading input

from os import system, remove, environ, path
from sys import argv
import numpy as np
#import imp
#/home/kubeckaj/Applications/JKCS2.1/JKQC/JKCS/lib64/python3.9/site-packages/numpy
#np = imp.load_module("numpy",None,"/home/kubeckaj/Applications/JKCS2.1/JKQC/JKCS/lib64/python3.9/site-packages/numpy",None)
#import importlib
#np = importlib.import_module('numpy')
#from numpy import nan, array, mod, isnan, linalg, matrix, all, dtype, loadtxt, unique, sqrt, sum, asarray, sort, log, exp, pi, mean, tile, floor, transpose, min, errstate, delete, round, dot, prod, random, apply_along_axis

#Attach command to the output
cmd="".join(["echo COMMAND: JKQC "," ".join(argv[1:])," >> output 2>/dev/null"])
system(cmd)
  

def listToString(s,spaces): 
    # initialize an empty string
    str1 = ""
    # traverse in the string  
    for ele in s: 
        str1 += ele  
        str1 += spaces 
    # return string  
    return str1

# Checking if correct environmnent is loaded
if str(environ.get("VIRTUAL_ENV").split("/")[-1]) != "JKCS":
  print("Trying to load JKCS environment by myself since you were lazy ass and did not do it!.")
  from subprocess import call
  if not call(['/bin/bash', '-i', '-c', "JKpython; python "+listToString(argv," ")]):
    print("This time, JKCS environment was loaded. Please, load it next time by yourself (JKpython).")
  else:
    print("You did not load JKCS environment and I was not able to do it for you.")
  exit()

def print_help():
  print("##################################################")
  print("### python JKQC.py [FILES] [ARGUMENTS] [PRINT] ###")
  print("\nFILES:")
  print(" any .xyz .log .out files")
  print(" any .pkl (you can also use -in/-out)")
  print("\nARGUMENTS:")
  print(" -help       print this help")
  print(" -info       information about a database")
  print(" -in X.pkl   read data from a database")
  print(" -out X.pkl  save data to XX.pkl (-noex = do not print example)")
  print(" -orcaext X  if ORCA has different extension than out")
  print(" -folder X   takes in all X/*.log files (use -R for resursive)")
  print(" -noname     the file names are not analysed (e.g. 1000-10_1.xyz)")
  print(" -rename X Y renames e.g. 1X2sa to 1Y2sa")
  print(" -extract X  prints only selected clusters (e.g. 1sa1w,1sa3-4w or 1sa1w-1_0)")
  print(" -except X   prints only non-selected clusters (e.g. 1sa1w,1sa3-4w or 1sa1w-1_0)")
  print(" -reacted    removes (the minority of) reacted structures")
  print(" -noexample,-noex does not print an example")
  print("\nPRINT:")
  print(" -ct                cluster type (e.g. 1sa,3sa2am)")
  print(" -b                 basename (e.g. 1sa-3, 3sa2am-4_508)")
  print(" -nXYZ,-nLOG,-nOUT  name with its extension")
  print(" -pXYZ,-pLOG,-pOUT  full_path/name with its extension")
  print(" -ePKL              pkl_full_filepath/:EXTRACT:/basename")
  print(" -el       electronic energy from .log [Eh]   -elSP        single-point el.en. (1. in .log) [Eh]")
  print(" -elout    el. energy from .out [Eh]          -rot         RMSD of rotational constants [GHz]") 
  print(" -elc      elc=elout-el [Eh]                  -rots        rotational constants [GHz]")
  print(" -zpec     ZPE correction [Eh]                -mult        muliplicity [-]")
  print(" -zpe      ZPE=el+zpec [Eh]                   -char,-chrg  charge [-electron_charge]")
  print(" -zpeout   ZPEout=elout+zpec [Eh]             -dip         dipole moment [Debye]")
  print(" -uc       energy thermal correction [Eh]     -dips        dipole moments [Debye]")
  print(" -u        internal energy U=el+uc [Eh]       -pol         polarizability [Angstrom^3]")
  print(" -uout     Uout=elout+uc [Eh]                 -templog     temperature used in .log [K]")
  print(" -hc       enthalpy th. corr. hc=uc+kT [Eh]   -preslog     pressure used in .log [atm]")
  print(" -h        enthalpy energy H=el+hc [Eh]       -lf          the lowest vib. freq. [1/cm]")
  print(" -hout     Hout=elout+hc [Eh]                 -f           array of vibration freq. [1/cm]")
  print(" -s        entropy [Eh/K]                     -rsn         rotational symmetry number [-]")
  print(" -gc       Gibbs free energy th. corr. [Eh]   -t,-time     total computational (elapsed) time [mins]")
  print(" -g        Gibbs free energy [Eh]             -termination count lines with \"Normal termination\" status")
  print(" -gout     G with el.en.corr.:Gout=G+elc [Eh] -mi          moments of inertia")
  print(" -mull     Mulliken charges [-el.charge]      -ami         average moment of inertia")
  print(" -xyz      save all XYZ files                 -rg          radius of gyration [Angstrom]")
  print(" -movie    save all XYZs to movie.xyz         -radius      approx. radius of cluster size [Angstrom]")
  print(" -charges  save all Mulliken .charges files   -radius0.5   radius with +0.5 Angstrom correction")
  print(" -mass     prints molecular mass [g/mol]      -level       program version and method used")
  print("\nPOST-CALCULATIONS:")
  print(" -fc [value in cm^-1] frequency cut-off for low-vibrational frequencies CITE: Grimme")
  print(" -temp [value in K]   recalculate for different temperature")
  print(" -v,-as [value]       anharmonicity scaling factor CITE: Grimme")
  print(" -unit                converts units [Eh] -> [kcal/mol] (for entropy: [Eh/K] -> [cal/mol/K])")
  print("\nFILTERING:")
  print(" -sort <str>          sort by: g,gout,el,elout,<full name in database under log>")
  print(" -select <int>        selects <int> best structures from each cluster")
  print(" -uniq,-unique <str>  selects only unique based on, e.g.: rg,el or rg,g or rg,el,dip, or dup=duplicates")
  print(" -sample <int> <str>  selects <INT> distinct based on, e.g.: rg,el or rg,g or rg,el,dip")
  print("                      use e.g. rg2,el0.5 to define threshold as 10**-x [def: rg=2, el/g=3, dip=1]")
  print(" -arbalign <float>    use (modified) ArbAlign program to compare RMSD (by def sort -el). CITE ArbAlign!!")
  print(" -cut/-pass X Y       removes higher/lower values of X=rg,el,g... with cutoff Y (e.g. -cut el -103.45)")
  print(" -cutr/-passr X Y     removes higher/lower rel. values compared to lowest of X=rg,el,g... with cutoff Y (e.g. -cutr g 5)")
  print("\nFORMATION PROPERTIES:")
  print(" -pop                    prints column of population probability")
  print(" -glob OR -globout       prints only values for clusters with the lowest -g OR -gout")
  print(" -bavg OR -bavgout       prints a value that is Boltzmann average over each cluster using -g OR -gout")
  print("                         NOTE: -g/-gout is treated correctly + -s not treated; use (G - H)/T")
  print(" -formation              print values as formation ")
  print(" <input_file> -formation print formations for the input file (no averaging though)")
  print(" -conc sa 0.00001        dG at given conc. [conc in Pa] (or use ppt,ppb,cmmc. e.g. 100ppt")
  print("                         use -cnt for self-consistent dG")
  print("\nOTHERS:")
  print(" -add <column> <file>, -extra <column>, -rebasename, -presplit, -i/-index <int:int>, -imos, -imos_xlsx,")
  print(" -forces [Eh/Ang], -shuffle, -split <int>, -underscore, -addSP <pickle>, -complement <pickle>")
  print(" -column <COL1> <COL2>, -drop <COL>, -out2log, -levels, -atoms, -natoms, -hydration/-solvation <str>")
  print(" -rh <0.0-1.0>, -psolvent <float in Pa>, -anharm, -bonded <float thr.> <element> <element>")

#OTHERS: -imos,-imos_xlsx,-esp,-chargesESP

folder = "./"	
Qcollect = "log" #what extension am I collecting?
Qrecursive = False
files = []  
input_pkl = []
input_pkl_sp = []
output_pkl = "mydatabase.pkl"
Qforces = 0 #should I collect forces? 0/1

addcolumn = []
Qmodify = 0 #Do I want to somehow modify the input dataframe or names?
Qrename = 0 #Do I want to rename some monomer? 
QrenameWHAT = [] #Do I want to rename some monomer? 
Qiamifo = 0 #This is calling renaming only for Ivo Neefjes
Qrebasename = 0 #change names if they are the same 
Qunderscore = 0 #underscore to dash
Qchangeall = 0 #do I want to some column to all
Qcomplement = 0 #subtract these from the list base on basename of the new pickle file
QcolumnDO = 0
Qcolumn = []

Qclustername = 1 #Analyse file names for cluster definition?
Qextract = 0 #Do I want to extarct only some cluster_type(s)?
Pextract = []
Qreacted = 0 #Remove reacted structures?, 1 = yes; 2 = print the reverse 
bonddistancethreshold = 1.75
Qbonded = []

Qout = 0 #Do I want to save output.pkl? 0=NO,1=YES,2=YES but do not print example
Pout = []

missing = np.nan

QUenergy = 1 #if you want to energy 
QUentropy = 1 #if you want to energy 

Qqha = 0 #Run the QHA
Qt = missing
Qp = missing
Qfc = 0 #Run QHA with vib. frequency cutoff
Qanh = 1
Qanharm = 0 #To collect anharmonicities from QC output

Qpresplit = 0 #Do I want to take only part of the data?
Qsplit = 1 #should I split the base on several parts
Qindex = "-1"#
Qdrop = "0" #drop some column
Qout2log = 0 #change column name

Qsort = 0 # 0=no sorting, otherwise string
Qselect = 0 # 0=nothing, otherwise the number of selected structures
Qsample = 0 # The same as uniqueness but adjustable automatically
Quniq = 0 # uniqie based on given arguments
Qarbalign = 0 #use ArbAlign with float parameter
formation_input_file = ""
Qthreshold = 0 #cut/pass something
Qcut = [] #what will be cutted
Qshuffle = 0 #shuffle rows

Qglob = 0 # 1=values for lowest -g, 2=values for lowest -gout
Qbavg = 0 # 1=Boltzmann avg over -g, 2=Boltzmann avg over -gout
Qrh = missing
QPsolvent = missing
Qformation = 0
Qconc = 0
conc = []
CNTfactor = 0 #see Wyslouzil
Qsolvation = "0"

folderMAIN = ""
folderSP = ""

orcaext = "out"
orcaextname = "out"
mrccextname = "out"
turbomoleext = "out"
turbomoleextname = "out"

last = ""
for i in argv[1:]:
  #HELP
  if i == "-help" or i == "--help":
    print_help()
    exit()
  #FOLDER
  if i == "-folder":
    last = "-folder"
    continue
  if last == "-folder":
    last = ""
    if path.exists(i):
      folder = i
      continue
    else:
      print("Folder "+i+" does not exist. [EXITING]")
      exit()
  #COLLECT
  if i == "-collect" or i == "-col" or i == "-coll":
    last = "-collect"
    continue
  if last == "-collect":
    last = ""
    Qcollect = str(i)
    continue
  if i == "-R":
    Qrecursive = True
    continue
  #INPKL
  if i == "-in":
    last = "-in"
    continue
  if last == "-in":
    last = ""
    if path.exists(i):    
      input_pkl.append(i)
      continue
    else:
      print("File "+i+" does not exist. Sorry [EXITING]")
      exit()
  #INPKL SP
  if i == "-addSP":
    last = "-addSP"
    continue
  if last == "-addSP":
    last = ""
    if path.exists(i):
      input_pkl_sp.append(i)
      continue
    else:
      print("File "+i+" does not exist. Sorry [EXITING]")
      exit()
  #OUTPKL
  if i == "-out" or i == "-o":
    last = "-out"
    continue
  if last == "-out":
    last = ""
    if output_pkl != "mydatabase.pkl":
      print("Hey are you trying to use two different outputs? "+output_pkl+"/"+i+" [EXITING]")
      exit()
    else:
      output_pkl = i
      if Qout == 0:
        Qout = 1
      continue
  #Hydration
  if i == "-hydration":
    Qsolvation = "w"
    continue
  #Anharm
  if i == "-anharm":
    Qanharm = 1
    continue
  #Qrh
  if i == "-rh":
    last = "-rh"
    continue
  if last == "-rh":
    last = ""
    Qrh = float(i)
    continue
  #psolvent
  if i == "-psolvent":
    last = "-psolvent"
    continue
  if last == "-psolvent":
    last = ""
    QPsolvent = float(i)
    continue
  #Solvation
  if i == "-solvation":
    last = "-solvation"
    continue
  if last == "-solvation":
    last = ""
    Qsolvation = str(i)
    continue
  #bonded
  if i == "-bonded":
    last = "-bonded"
    rem = []
    continue
  if last == "-bonded":
    last = "-bonded1"
    rem.append(str(i))
    continue
  if last == "-bonded1":
    last = "-bonded2"
    rem.append(str(i))
    continue
  if last == "-bonded2":
    last = "" 
    rem.append(str(i))
    Qbonded.append(rem)
    rem = ""
    Pout.append("-bonded")
    Qbonded_index = -1
    continue
  #IamIfo
  if i == "-IamIfo":
    Qiamifo = 1
    continue
  #FORCES
  if i == "-forces":
    Qforces = 1
    continue
  #ORCA EXTENSION
  if i == "-orcaext" or i == "-orca":
    last = "-orcaext"
    continue
  if last == "-orcaext":
    last = ""
    orcaext = str(i)
    orcaextname = "log"
    continue
  #TURBOMOLE EXTENSION
  if i == "-turbomoleext" or i == "-turbomole":
    last = "-turbomoleext"
    continue
  if last == "-turbomoleext":
    last = ""
    turbomoleext = str(i)
    turbomoleextname = "log"
    continue
  #COLUMN 
  if i == "-column":
    last = "-column"
    QcolumnDO += 1
    Pout.append("-column")
    continue
  if last == "-column":
    last = "-column1"
    QcolumnADD = str(i)
    continue
  if last == "-column1":
    last = ""
    Qcolumn.append([QcolumnADD, str(i)])
    continue
  #LEVELS
  if i == "-levels":
    Pout.append("-levels")
    continue
  #ADD EXTRA COLUMN
  if i == "-add" or i == "-addcolumn":
    last = "-add"
    continue
  if last == "-add":
    last = "-add2"
    columnname = str(i)
    continue
  if last == "-add2":
    last = "" 
    if path.exists(str(i)):
      addcolumn.append([columnname, str(i)])
      continue
    else:
      print("File "+i+" does not exist. Sorry [EXITING]")
      exit()       
  #NONAME
  if i == "-noname":
    Qclustername = 0
    continue
  #PRESPLIT
  if i == "-presplit":
    last = "-presplit"
    continue
  if last == "-presplit":
    last = ""
    Qpresplit = int(i)
    continue
  # DROP
  if i == "-drop":
    last = "-drop"
    continue
  if last == "-drop":
    last = ""
    Qdrop = i
    continue
  #COMPLEMENT
  if i == "-complement":
    last = "-complement"
    continue
  if last == "-complement":
    last = ""
    Qcomplement = str(i)
    continue
  #INDEX
  if i == "-i" or i == "-index":
    last = "-index"
    continue
  if last == "-index":
    last = ""
    Qindex = str(i)
    continue
  #UNDERSCORE
  if i == "-undersocre":
    Qunderscore = 1
    continue
  # OUT 2 LOG
  if i == "-out2log":
    Qout2log = 1
    continue
  #RENAME
  if i == "-rename":
    Qmodify = 1
    Qrename = 1
    last = "-rename"
    continue
  if last == "-rename":
    remember = i
    last = "-rename2"
    continue
  if last == "-rename2":
    last = ""
    QrenameWHAT.append([remember,i])
    continue
  #rebasename
  if i == "-rebasename":
    Qrebasename = 1
    continue
  #REACTED
  if i == "-reacted":
    Qreacted = 1
  #  last = "-reacted"
    continue
  if i == "-reacted2":
    Qreacted = 2
  #  last = "-reacted"
    continue
  #if last == "-reacted":
  #  last = ""
  #  bonddistancethreshold = float(i)
  #  continue  
  #NOEXAMPLE
  if i == "-noexample" or i == "-noex":
    Qout = 2
    continue
  #SPLIT
  if i == "-split":
    last = "-split"
    continue
  if last == "-split":
    Qsplit = int(i)
    last = ""
    continue
  #EXTRACT
  if i == "-extract":
    last = "-extract"
    continue
  if last == "-extract":
    last = ""
    Qextract = 1
    Pextract.append(i)
    continue
  #EXCEPT
  if i == "-except":
    last = "-except"
    continue
  if last == "-except":
    last = ""
    Qextract = 2
    Pextract.append(i)
    continue
  #INPUT FILES
  if len(i) > 3:
    ext = i[-4:]
    if ext == ".xyz" or ext == ".log" or ext == ".out":
      if path.exists(i):
        files.append(i)
        continue
      else:
        print("File "+i+" does not exist. Sorry [EXITING]")
        exit()
    if ext == ".pkl":
      if path.exists(i):
        input_pkl.append(i)
        continue
      else:
        if output_pkl != "mydatabase.pkl":
          print("Hey either one input database does not exist or you are trying to use two different outputs? "+output_pkl+"/"+i+" [EXITING]")
          exit()
        else:
          output_pkl = i
          continue
  ########
  # SORT 
  if i == "-sort" or i == "--sort":
    last = "-sort"
    continue
  if last == "-sort":
    last = ""
    Qsort = str(i)
    continue
  # SELECT
  if i == "-select" or i == "--select":
    last = "-select"
    continue
  if last == "-select":
    last = ""
    Qselect = int(i)
    continue  
  # SAMPLE
  if i == "-sample":
    last = "-sample"
    continue
  if last == "-sample":
    last = "-uniq"
    Qsample = int(i)
    continue
  # UNIQUE
  if i == "-unique" or i == "--unique" or i == "-uniq" or i == "--uniq":
    last = "-uniq"
    continue
  if last == "-uniq":
    last = ""
    Quniq = str(i)
    continue
  if i == "-cut":
    Qthreshold = 1
    last = "-threshold"
    attach = ["<=","absolute"]
    continue
  if i == "-cutr":
    Qthreshold = 1
    last = "-threshold"
    attach = ["<=","relative"]
    continue
  if i == "-pass":
    Qthreshold = 1
    last = "-threshold"
    attach = [">","absolute"]
    continue
  if i == "-passr":
    Qthreshold = 1
    last = "-threshold"
    attach = [">","relative"]
    continue
  if last == "-threshold":
    last = "-threshold2"
    attach.append(i)
    continue
  if last == "-threshold2":
    last = ""
    attach.append(str(i))
    Qcut.append(attach)
    continue 
  #ArbAlign
  if i == "-arbalign" or i == "-ArbAlign":
    last = "-arbalign"
    continue
  if last == "-arbalign":
    last = ""
    Qarbalign = float(i)
    continue
  #SHUFFLE
  if i == "-shuffle":
    Qshuffle = 1
    continue
  ########
  # INFO
  if i == "-info" or i == "--info":
    Pout.append("-info")
    continue
  # CITE
  if i == "-cite" or i == "--cite":
    Pout.append("-cite")
    continue
  # ATOMS
  if i == "-atoms":
    Pout.append("-atoms")
    continue
  # XYZ
  if i == "-xyz" or i == "--xyz" or i == "-XYZ" or i == "--XYZ":
    Pout.append("-xyz")
    continue
  if i == "-imos" or i == "--imos" or i == "-IMOS" or i == "--IMOS" or i == "-IMoS" or i == "--IMoS":
    Pout.append("-imos")
    continue
  if i == "-imos_xlsx" or i == "--imos_xlsx":
    Pout.append("-imos_xlsx")
    continue
  # MOVIE
  if i == "-movie" or i == "--movie":
    Pout.append("-movie")
    continue
  # RG
  if i == "-rg" or i == "--rg" or i == "-Rg" or i == "--Rg":
    Pout.append("-rg")
    continue
  if i == "-radius" or i == "--radius":
    Pout.append("-radius")
    continue
  if i == "-radius0.5" or i == "--radius0.5":
    Pout.append("-radius0.5")
    continue
  # INFO
  if i == "-ct" or i == "-clustertype" or i == "--ct" or i == "--clustertype":
    Pout.append("-ct")
    continue
  if i == "-b" or i == "-basename" or i == "--b" or i == "--basename":
    Pout.append("-b")
    continue
  if i == "-nOUT" or i == "--nOUT":
    Pout.append("-nOUT")
    continue
  if i == "-nLOG" or i == "--nLOG":
    Pout.append("-nLOG")
    continue
  if i == "-nXYZ" or i == "--nXYZ":
    Pout.append("-nXYZ")
    continue
  if i == "-pOUT" or i == "--pOUT":
    Pout.append("-pOUT")
    continue
  if i == "-pLOG" or i == "--pLOG":
    Pout.append("-pLOG")
    continue
  if i == "-pXYZ" or i == "--pXYZ":
    Pout.append("-pXYZ")
    continue
  if i == "-ePKL" or i == "--ePKL":
    Pout.append("-ePKL")
    continue
  # LOG & OUT
  if i == "-elSP" or i == "-elsp" or i == "--elSP" or i == "--elsp":
    Pout.append("-elsp")
    continue
  if i == "-el" or i == "-elen" or i == "--el" or i == "--elen":
    Pout.append("-el")
    continue
  if i == "-elout" or i == "--elout":
    Pout.append("-elout")
    continue
  if i == "-elc" or i == "--elc":
    Pout.append("-elc")
    continue
  if i == "-g" or i == "-gibbs" or i == "--g" or i == "--gibbs":
    Pout.append("-g")
    continue
  if i == "-pop" or i == "--pop":
    Pout.append("-pop")
    continue
  if i == "-natoms":
    Pout.append("-natoms")
    continue
  if i == "-h" or i == "-enthalpy" or i == "--h" or i == "--enthalpy":
    Pout.append("-h")
    continue
  if i == "-s" or i == "-entropy" or i == "--s" or i == "--entropy":
    Pout.append("-s")
    continue
  if i == "-mass":
    Pout.append("-mass")
    continue
  if i == "-u" or i == "--u":
    Pout.append("-u")
    continue
  if i == "-uc" or i == "--uc":
    Pout.append("-uc")
    continue
  if i == "-uout" or i == "--uout":
    Pout.append("-uout")
    continue
  if i == "-hc" or i == "--hc":
    Pout.append("-hc")
    continue
  if i == "-hout" or i == "--hout":
    Pout.append("-hout")
    continue
  if i == "-gc" or i == "--gc":
    Pout.append("-gc")
    continue
  if i == "-g" or i == "--g" or i == "-gibbs" or i == "--gibbs" or i == "-G" or i == "--G":
    Pout.append("-g")
    continue
  if i == "-gout" or i == "--gout" or i == "-GDh" or i == "--GDh" or i == "-GD" or i == "--GD":
    Pout.append("-gout")
    continue
  if i == "-zpec" or i == "-zpecorr" or i == "--zpec" or i == "--zpecorr":
    Pout.append("-zpec")
    continue
  if i == "-zpe" or i == "-ZPE" or i == "--zpe" or i == "--ZPE":
    Pout.append("-zpe")
    continue
  if i == "-zpeout" or i == "-ZPEout" or i == "--zpeout" or i == "--ZPEout":
    Pout.append("-zpeout")
    continue
  if i == "-lf" or i == "-lfreq" or i == "--lf" or i == "--lfreq":
    Pout.append("-lf")
    continue
  if i == "-f" or i == "-freq" or i == "--f" or i == "--freq":
    Pout.append("-f")
    continue
  if i == "-rot" or i == "--rot":
    Pout.append("-rot")
    continue
  if i == "-rots" or i == "--rots":
    Pout.append("-rots")
    continue
  if i == "-mult" or i == "--mult":
    Pout.append("-mult")
    continue
  if i == "-char" or i == "-chrg" or i == "--char" or i == "--chrg":
    Pout.append("-char")
    continue
  if i == "-mull" or i == "--mull":
    Pout.append("-mull")
    continue
  if i == "-esp" or i == "--esp":
    Pout.append("-esp")
    continue
  if i == "-charges" or i == "--charges":
    Pout.append("-charges")
    continue
  if i == "-chargesESP" or i == "--chargesESP":
    Pout.append("-chargesESP")
    continue
  if i == "-dip" or i == "--dip":
    Pout.append("-dip")
    continue
  if i == "-dips" or i == "--dips":
    Pout.append("-dips")
    continue
  if i == "-pol" or i == "--pol":
    Pout.append("-pol")
    continue
  if i == "-templog" or i == "--templog":
    Pout.append("-templog")
    continue
  if i == "-preslog" or i == "--preslog":
    Pout.append("-preslog")
    continue
  if i == "-mi" or i == "--mi":
    Pout.append("-mi")
    continue
  if i == "-ami" or i == "--ami":
    Pout.append("-ami")
    continue
  if i == "-rsn" or i == "--rsn":
    Pout.append("-rsn")
    continue
  if i == "-level" or i == "-lvl":
    Pout.append("-level")
    continue
  if i == "-t" or i == "--t" or i == "-time" or i == "--time":
    Pout.append("-t")
    continue
  if i == "-termination" or i == "--termination":
    Pout.append("-termination")
    continue
  if i == "-extra":
    Pout.append("-extra")
    last = "-extra"
    continue
  if last == "-extra":
    last = ""
    Pout.append(str(i))
    continue
  #PRE_EXTRACT_DATA MANIPULATION
  if i == "-fc" or i == "--fc":
    Qqha = 1
    last = "-fc"
    continue
  if last == "-fc":
    last = ""
    Qfc = float(i)
    continue
  if i == "-temp" or i == "--temp":
    Qqha = 1
    last = "-temp"
    continue
  if last == "-temp":
    last = ""
    Qt = float(i)
    continue
  if i == "-as" or i == "--as" or i == "-v" or i == "--v":
    Qqha = 1
    last = "-as"
    continue
  if last == "-as":
    last = ""
    Qanh = str(i)
    continue
  if i == "-unit" or i == "--unit":
    QUenergy = 627.503
    QUentropy = 627503.0
    continue
  #FORMATION
  if i == "-glob" or i == "--glob":
    Qglob = 1
    continue
  if i == "-globout" or i == "--globout":
    Qglob = 2
    continue
  if i == "-bavg" or i == "--bavg":
    Qbavg = 1
    continue
  if i == "-bavgout" or i == "--bavgout":
    Qbavg = 2
    continue
  if i == "-formation" or i == "--formation":
    Qformation = 1
    continue
  if i == "-conc" or i == "--conc":
    Qconc = 1
    last = "-conc"
    continue
  if last == "-conc":
    last = "-conc2"
    remember = str(i)
    continue
  if last == "-conc2":
    last = ""
    conc.append(np.array([remember, i]))
    continue
  if path.exists(i):
    formation_input_file = i
    continue
  if i == "-cnt" or i == "--cnt":
    CNTfactor = 1
    continue
  #UNKNOWN ARGUMENT   
  print("I am sorry but I do not understand the argument: "+i+" [EXITING]")
  exit()  

####################################################################################################
# Checking for arguments

if len(files) == 0:
  if len(input_pkl) == 0 and len(formation_input_file) == 0:
    from glob import glob
    if Qrecursive:
      files = glob(folder+"/**/*."+Qcollect, recursive = True)
    else:
      files = glob(folder+"/*."+Qcollect)
    if len(files) == 0:
      print("No inputs. No *.log files. [EXITING]")
      exit()
  if len(input_pkl) == 1 and Qsplit != 1 and output_pkl == "mydatabase.pkl":
    output_pkl = input_pkl[0]
    if Qout == 0:
      Qout = 1
  if len(input_pkl) > 1:
    if Qout == 0:
      Qout = 1
else:
  if len(Pout) == 0:
    if Qout==0:
      Qout = 1
    print("Number of files: "+str(len(files)))
  else:
    Qout = 2

import pandas as pd
import re
from ase.io import read,write

####################################################################################################
####################################################################################################
####################################################################################################
## MANIPULATING DATAFRAMES
# THIS IS TAKEN FROM tool.py IN JKCS.py 
### Append to dataframe 
def df_add_append(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    newdataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
    dataframe = dataframe.append(newdataframe)
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )

  return dataframe

### Add to dataframe via iteration
def df_add_iter(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    df_l = dataframe.shape[0]
    var_l = len(variables)
    for i in range(var_l):
      dataframe[label,name][df_l-var_l+i] = variables[i]
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
    for i in range(len(variables)):
      dataframe[label,name].iloc[-1-i] = variables[-1-i]

  return dataframe
####################################################################################################
####################################################################################################
####################################################################################################
## HERE ARE SOME FUNCTIONS THAT I JUST NEEDED SO FAST THAT I EVEN DID NOT COMMENT ON THEM

def seperate_string_number(string):
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i == "_" or previous_character == "_":
            newword += i
        elif i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i
        previous_character = i
        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups

def seperate_string_number2(string):
    ### rg3,el2.4,g -> [rg,3,el,2.4,g]
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i == "_" or previous_character == "_":
            newword += i
        elif i.isalpha() and previous_character.isalpha():
            newword += i
        elif (i.isnumeric() or i == ".") and (previous_character.isnumeric() or previous_character == "."):
            newword += i
        else:
            if previous_character != ",":
              if previous_character.isnumeric() or previous_character == ".":
                groups[-1]=[groups[-1],newword]
              else:
                groups.append(newword)
            if i != ",":  
              newword = i
        previous_character = i
        if x == len(string) - 2:
            if previous_character.isnumeric() or previous_character == ".":
              groups[-1]=[groups[-1],newword]
            else:
              groups.append(newword)
            newword = ''
    return groups

def zeros(input_array):
  output_string = ""
  skip = 0
  for i in range(len(input_array)):
    if skip == 1:
      skip = 0
      continue
    if input_array[i] == "0":
      skip = 1
      continue
    output_string += input_array[i]
  return output_string

def is_nameable(input_array):
  #input_array = seperate_string_number(file_i_BASE.split("-")[0])
  nameable_test = True
  if np.mod(len(input_array),2) == 0:
    for input_array_i in input_array[0::2]:
      if not input_array_i.isnumeric():
        nameable_test = False
        break
    for input_array_i in input_array[1::2]:  
      if input_array_i.isnumeric():
        nameable_test = False
        break
  else:
    nameable_test = False
  return nameable_test

####################################################################################################
####################################################################################################
####################################################################################################

# Loading input pickles
if len(input_pkl) == 0:
  clusters_df = pd.DataFrame()
else:
  for i in range(len(input_pkl)):
    newclusters_df = pd.read_pickle(input_pkl[i])
    newclusters_df.index = [str(j) for j in range(len(newclusters_df))]
    if Qout == 1:
      print("Number of files in "+input_pkl[i]+": "+str(len(newclusters_df)))
    if i == 0: 
      clusters_df = newclusters_df
    else:
      len_clusters_df = len(clusters_df)
      newclusters_df.index = [str(j+len_clusters_df) for j in range(len(newclusters_df))]
      #clusters_df = clusters_df.reset_index().append(newclusters_df.reset_index())
      clusters_df = pd.concat([clusters_df,newclusters_df],axis=0, ignore_index=True)
      #clusters_df = pd.concat([clusters_df,newclusters_df],axis=0,ignore_index=True)
      #clusters_df = clusters_df.append(newclusters_df, ignore_index=True)
  #clusters_df.reindex([str(j) for j in range(len(clusters_df))])
  #print(clusters_df.index)
  #exit()

#Complement
if Qcomplement != 0:
  newclusters_df = pd.read_pickle(Qcomplement)
  clusters_df = clusters_df.merge(newclusters_df, on = [('info', 'file_basename')], how = "outer", indicator = True).loc[lambda x: x['_merge'] == 'left_only'].drop('_merge', axis=1)

#addSP
if len(input_pkl_sp) > 0:
  for i in range(len(input_pkl_sp)):
    newclusters_df_sp = pd.read_pickle(input_pkl_sp[i])
    if Qout == 1:
      print("Number of files in "+input_pkl_sp[i]+": "+str(len(newclusters_df_sp)))
    newclusters_df_sp = newclusters_df_sp.rename(columns={"log": "out"}) 
    #newclusters_df_sp = pd.concat([newclusters_df_sp[("info","file_basename")],newclusters_df_sp["out"]],axis = 1)
    newclusters_df_sp.set_index(('info', 'file_basename'), inplace=True)
    #print(newclusters_df_sp.index)
    #newclusters_df_sp = newclusters_df_sp
    clusters_df.set_index(('info', 'file_basename'), inplace=True)
    #print(clusters_df.index)
    #print(clusters_df.join(newclusters_df_sp,how='left'))
    #print(newclusters_df_sp.loc[:,newclusters_df_sp.columns.get_level_values(0) == 'out'])
    newclusters_df_sp = newclusters_df_sp.loc[:,newclusters_df_sp.columns.get_level_values(0) == 'out']
    clusters_df = pd.concat([clusters_df, newclusters_df_sp], axis=1)
    clusters_df = clusters_df.sort_index()
    newclusters_df_sp = newclusters_df_sp.sort_index()
    #print(clusters_df.merge(newclusters_df_sp, on=('info', 'file_basename')))
    #clusters_df.index = [str(j) for j in range(len(clusters_df))]
    clusters_df = clusters_df.sort_index()
    #print(clusters_df)
    clusters_df.reset_index(inplace=True)

####################################################################################################

# Reading the new files in
rem_orcaextname = orcaextname
for file_i in files:
  folder = path.abspath(file_i)[::-1].split("/",1)[1][::-1]+"/"
  file_i_BASE = file_i[:-4][::-1].split("/",1)[0][::-1]
  file_i_ABC  = folder+file_i_BASE+".log"
  file_i_XTB  = folder+file_i_BASE+".log"
  file_i_XTBengrad  = folder+file_i_BASE+".engrad"
  file_i_G16  = folder+file_i_BASE+".log"
  file_i_XYZ  = folder+file_i_BASE+".xyz"
  #ORCA && MRCC
  file_i_MRCC  = folder+file_i_BASE+".out"
  orcaextname = rem_orcaextname
  if not path.exists(folder+file_i_BASE+".log"):
    file_i_ORCA = folder+file_i_BASE+"."+orcaext
    orcaextname = "log"
    mrccextname = "log"
  else:
    if orcaext == "out":
      if not path.exists(folder+file_i_BASE+"."+orcaext):
        file_i_ORCA = folder+file_i_BASE+".log"
        orcaextname = "log"
      else:
        file_i_ORCA = folder+file_i_BASE+"."+orcaext
    else:
      file_i_ORCA = folder+file_i_BASE+"."+orcaext
      orcaextname = "log"
  ##
  file_i_TURBOMOLE = folder+file_i_BASE+"."+turbomoleext
  file_i_INFO = folder+"info.txt" 
 
  ###############
  ### INFO ######
  ###############
  cluster_id = len(clusters_df)
  clusters_df = df_add_append(clusters_df, "info", "folder_path", [str(cluster_id)], folder)
  clusters_df = df_add_iter(clusters_df, "info", "file_basename", [str(cluster_id)], [file_i_BASE])
  # sorted cluster type
  if Qclustername == 1:
    cluster_type_array = seperate_string_number(file_i_BASE.split("-")[0].split("_")[0]) 
    #if np.mod(len(cluster_type_array),2) == 0:
    if is_nameable(cluster_type_array):
      cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
      cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
      cluster_type_string = zeros(cluster_type_array_sorted)
      clusters_df = df_add_iter(clusters_df, "info", "cluster_type", [str(cluster_id)], [cluster_type_string])
      #
      components = re.split('(\d+)', file_i_BASE.split("-")[0].split("_")[0])[1:][1::2]
      clusters_df = df_add_iter(clusters_df, "info", "components", [str(cluster_id)], [components])
      component_ratio = [int(i) for i in re.split('(\d+)', file_i_BASE.split("-")[0].split("_")[0])[1:][0::2]]
      clusters_df = df_add_iter(clusters_df, "info", "component_ratio", [str(cluster_id)], [component_ratio])
  if path.exists(file_i_INFO):
    file = open(file_i_INFO, "r")  
    for line in file:
      clusters_df = df_add_iter(clusters_df, "info", str(line.split(" ",1)[0]), [str(cluster_id)], [line.split(" ",1)[-1].strip()])    
    file.close()

  #####################
  ### ABC/CREST #######
  #####################
  if path.exists(file_i_ABC):
    file = open(file_i_ABC, "r")
    test = 0
    for i in range(1):
      line=file.readline()
      if re.search("ABC", line):
        test = 1
        break
      if re.search("JXYZ", line):
        test = 1
        break
    if test == 1:
      out_electronic_energy = missing  #1
      for line in file:
        if re.search("ABC energy:", line) or re.search("structure energy:", line): #1
          try:
            out_electronic_energy = float(line.split()[2]) 
          except:
            out_electronic_energy = missing
      clusters_df = df_add_iter(clusters_df, "log", "electronic_energy", [str(cluster_id)], [out_electronic_energy])
    file.close()
      

  ###############
  ### XTB #######
  ###############
  if path.exists(file_i_XTB):
    file = open(file_i_XTB, "r")
    testXTB = 0
    for i in range(5):
      if re.search("\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-", file.readline()):
        testXTB = 1
        break
    if testXTB == 1:
      out_vibrational_frequencies = missing
      out_program = missing                              #PROGRAM
      out_method = missing                               #METHOD
      out_time = missing                                 #TIME
      out_NAtoms = missing
      out_dipole_moment = missing      #0
      out_electronic_energy = missing  #1
      out_enthalpy_energy = missing                      #T8
      out_gibbs_free_energy = missing  #2
      out_entropy = missing                              #E1
      out_mulliken_charges = missing           
      save_something = ""
      for line in file:
        ## PROGRAM
        if re.search("xtb version", line):
          try: 
            out_program = "XTB_" + str(line.split()[3])
          except:
            out_program = missing
        ## METHOD
        if re.search(":  Hamiltonian", line):
          try:
            out_method = str(line.split()[2]).lower()
          except:
            out_method = missing
        ## NAtoms
        if re.search("number of atoms", line):
          try:
            out_NAtoms = int(line.split()[4])
          except:
            out_NAtoms = missing
          continue
        ## MULLIKEN
        if re.search("     #   Z          covCN         q      C6AA      ", line) and not np.isnan(out_NAtoms): #O4
          save_mulliken_charges = out_NAtoms
          save_something = "mulliken_charges"
          try:
            out_mulliken_charges = ["0"]*out_NAtoms
          except:
            out_mulliken_charges = missing
          continue
        #print(line)
        #print(save_something)
        #print(out_NAtoms)
        if save_something == "mulliken_charges":
          try:
            if save_mulliken_charges<=out_NAtoms:
              out_mulliken_charges[out_NAtoms-save_mulliken_charges] = str(line.split()[4])
          except:
            out_mulliken_charges = missing
          save_mulliken_charges-=1
          if save_mulliken_charges == 0:
            save_something = ""
          continue
        ## OTHERS
        if re.search("Debye", line): #0
          try:
            out_dipole_moment = float(line.split()[5])
          except:
            out_dipole_moment = missing
          continue
        if re.search("molecular dipole:", line):
          #print("TEST 1")
          #print(line)
          save_something = "molecular_dipole"
          save_molecular_dipole = 2
          continue
        if save_something == "molecular_dipole": 
          #print("TEST 2")
          #print(line)
          save_molecular_dipole-=1
          if save_molecular_dipole == 0:
            try: 
              out_dipole_moment = float(line.split()[4])
            except:
              out_dipole_moment = missing
            save_something = ""
          continue
        if re.search("TOTAL ENERGY", line): #1
          try:
            out_electronic_energy = float(line.split()[3])
          except:
            out_electronic_energy = missing
          continue
        if re.search("total E", line): #1
          try:
            out_electronic_energy = float(line.split()[3])
          except:
            out_electronic_energy = missing
          continue
        if re.search("^H\(T\)", line): #T8
          try:
            out_enthalpy_energy = out_electronic_energy + float(line.split()[1])
          except:
            out_enthalpy_energy = missing
        if re.search("TOTAL ENTHALPY", line): #T8
          try:
            out_enthalpy_energy = float(line.split()[3])
          except:
            out_enthalpy_energy = missing
        if re.search("^T\*S ", line): #E1
          try:
            out_entropy = float(line.split()[1])*1000*627.503/298.15
          except:
            out_entropy = missing
        if re.search("T\*S/Eh ", line): #E1
          save_something = "entropy"
          save_entropy = 1
          continue
        if save_something == "entropy":
          if save_entropy == 1:
            save_entropy -= 1
            continue
          try:
            out_entropy = float(line.split()[3])*1000*627.503/298.15
          except:
            out_entropy = missing
          save_something = ""
          continue
        if re.search("S\(FR\)/kcal", line):
          save_something = "vibs"
          save_vibs = 0
          out_vibrational_frequencies = []
          continue
        if save_something == "vibs":
          if re.search("----", line):
            save_vibs += 1
            continue
          if save_vibs == 1:
            try:
              out_vibrational_frequencies += [float(line.split()[1])]
            except:
              out_vibrational_frequencies = missing
          if save_vibs == 2:
            save_something = ""
          continue
        #if re.search("^G\(T\)", line): 
        #  try:
        #    out_gibbs_free_energy = out_electronic_energy + float(line.split()[1])
        #  except:
        #    out_gibbs_free_energy = missing
        if re.search("TOTAL FREE ENERGY", line): #2
          try:
            out_gibbs_free_energy = float(line.split()[4])
          except:
            out_gibbs_free_energy = missing
          continue
        if re.search("wall-time", line): #2
          if np.isnan(out_time):
            out_time = 0
          else:
            continue
          try:
            out_time += float(line.split()[2])*24*60+float(line.split()[4])*60+float(line.split()[6])+float(line.split()[8])/60
          except:
            out_time = missing
          continue
      clusters_df = df_add_iter(clusters_df, "log", "time", [str(cluster_id)], [out_time]) #TIME
      clusters_df = df_add_iter(clusters_df, "log", "program", [str(cluster_id)], [out_program]) #PROGRAM
      clusters_df = df_add_iter(clusters_df, "log", "method", [str(cluster_id)], [out_method]) #METHOD
      clusters_df = df_add_iter(clusters_df, "log", "dipole_moment", [str(cluster_id)], [out_dipole_moment])
      clusters_df = df_add_iter(clusters_df, "log", "electronic_energy", [str(cluster_id)], [out_electronic_energy])
      clusters_df = df_add_iter(clusters_df, "log", "gibbs_free_energy", [str(cluster_id)], [out_gibbs_free_energy])
      clusters_df = df_add_iter(clusters_df, "log", "enthalpy_energy", [str(cluster_id)], [out_enthalpy_energy])
      clusters_df = df_add_iter(clusters_df, "log", "entropy", [str(cluster_id)], [out_entropy])
      clusters_df = df_add_iter(clusters_df, "log", "mulliken_charges", [str(cluster_id)], [out_mulliken_charges]) #O4
      clusters_df = df_add_iter(clusters_df, "log", "NAtoms", [str(cluster_id)], [out_NAtoms]) #I3
      clusters_df = df_add_iter(clusters_df, "log", "vibrational_frequencies", [str(cluster_id)], [out_vibrational_frequencies]) 
    file.close()
    if Qforces == 1:
      if path.exists(file_i_XTBengrad):
        file = open(file_i_XTBengrad, "r")
        for gradi in range(11):
          file.readline()
        try:
          save_forces = []
          for gradi in range(3*out_NAtoms):
            save_forces.append(float(file.readline())/0.529177)
          out_forces = [np.array([save_forces[i],save_forces[i+1],save_forces[i+2]]) for i in range(0,len(save_forces),3)]
        except:
          out_forces = missing
        clusters_df = df_add_iter(clusters_df, "extra", "forces", [str(cluster_id)], [out_forces]) #F1
        file.close()
        

  ###############
  ### G16 #######
  ###############
  if path.exists(file_i_G16):
    file = open(file_i_G16, "r")
    testG16 = 0
    for i in range(5):
      if re.search("Gaussian", file.readline()):
        testG16 = 1
        break
    if testG16 == 1:
      out_program = missing                              #PROGRAM
      out_method = missing                               #METHOD
      out_time = missing                                 #TIME
      out_termination = 0                                #TERMINATION
      out_charge = missing                               #I1
      out_multiplicity = missing                         #I2
      out_NAtoms = missing                               #I3
      out_rotational_constants = missing                 #O1
      out_rotational_constant = missing                  #O2
      out_sp_electronic_energy = missing                 #O3x
      out_electronic_energy = missing                    #O3
      out_mulliken_charges = missing                     #O4
      out_esp_charges = missing                          #O4b
      out_dipole_moment = missing                        #O5
      out_dipole_moments = missing                       #O6
      out_polarizability = missing                       #O7
      out_vibrational_frequencies = missing              #V1
      out_temperature = missing                          #V2
      out_pressure = missing                             #V3
      out_moments_of_inertia = missing                   #V4
      out_rotational_symmetry_number = missing           #T1
      out_zero_point_correction = missing                #T2
      out_energy_thermal_correction = missing            #T3
      out_enthalpy_thermal_correction = missing          #T4
      out_gibbs_free_energy_thermal_correction = missing #T5
      out_zero_point_energy = missing                    #T6
      out_internal_energy = missing                      #T7
      out_enthalpy_energy = missing                      #T8
      out_gibbs_free_energy = missing                    #T9
      out_entropy = missing                              #E1
      if Qforces == 1:
        out_forces = missing                             #F1
      if Qanharm == 1:
        out_anharm = missing
      search=-1
      save_anharm=0
      save_mulliken_charges=0
      save_something=""
      for line in file:
        #TIME
        if re.search("Elapsed time",line): 
          if np.isnan(out_time):
            out_time = 0
          try:
            out_time += float(line.split()[2])*24*60+float(line.split()[4])*60+float(line.split()[6])+float(line.split()[8])/60 
          except:
            out_time = missing
          continue
        #TERMINATION
        if re.search("Normal termination",line):
          out_termination += 1
          continue
        #INITIAL SEARCH
        if search==-1:
          #PROGRAM
          if re.search(" Gaussian.*, Revision.*,", line): #PROGRAM
            try: 
              out_program = "G" + str(line.split(",")[0].split()[1]) + "_" + str(line.split(",")[1].split()[1]) 
            except:
              out_program = missing
          #METHOD
          if re.search("^ \#", line): #PROGRAM
            try:
              presaved = line.lower().split()
              out_method = "_".join(presaved)
              #if presaved[1] in ["am1", "pm7", "pm6", "pm3"]:
              #  out_method = "_".join(presaved[0:2])
              #else: 
              #  out_method = "_".join(presaved[0:3])
            except:
              out_method = missing
          #Charge
          if re.search(" Charge = ", line): #I1/I2
            try:
              out_charge = int(line.split()[2])
              out_multiplicity = int(line.split()[5])
            except:
              out_charge = missing
              out_multiplicity = missing
            search += 1
            continue
        #OPTIMIZATION SEARCH
        ## I will convert the strings into floats/ints later !!!!!!!
        if search==0 or search==-1: 
          if np.isnan(out_NAtoms): 
            if re.search(" NAtoms=", line): #I3
              try:
                out_NAtoms = int(line.split()[1])
              except:
                out_NAtoms = missing
              continue
          if re.search("Rotational constants ", line): #O1/O2
            try:
              out_rotational_constants = [str(line.split()[3]),str(line.split()[4]),str(line.split()[5])]
            except:
              out_rotational_constants = missing
            continue
          if re.search("SCF Done", line): #O3
            try:
              out_electronic_energy = float(line.split()[4])
            except:
              out_electronic_energy = missing
            if np.isnan(out_sp_electronic_energy):
              out_sp_electronic_energy = out_electronic_energy
            continue
          ## MULLIKEN
          if re.search(" Mulliken charges:", line): #O4
            save_mulliken_charges = out_NAtoms+1
            save_something = "mulliken_charges"
            try:
              out_mulliken_charges = ["0"]*out_NAtoms
            except:
              out_mulliken_charges = missing
            continue
          if save_something == "mulliken_charges":
            try:
              if save_mulliken_charges<=out_NAtoms:
                out_mulliken_charges[out_NAtoms-save_mulliken_charges] = str(line.split()[2])
            except:
              out_mulliken_charges = missing
            save_mulliken_charges-=1
            if save_mulliken_charges == 0:
              save_something = ""
            continue
          ## DIPOLE
          if re.search('Dipole moment \(field-independent basis, Debye\):', line): #O5/O6
            save_something = "dipole_moment"
            continue
          if save_something == "dipole_moment":
            try:
              out_dipole_moment = float(line.split()[7])
            except:
              out_dipole_moment = missing
            try:
              out_dipole_moments = [float(line.split()[1]), float(line.split()[3]), float(line.split()[5])]
            except:
              out_dipole_moments = missing
            save_something = ""
            continue
          # POLARIZABILITY 
          if re.search("Exact polarizability", line): #O7
            try:
              pol = [float(i) for i in line.split()[2:]]
              pol_mat = np.matrix([[pol[0],pol[1],pol[3]],[pol[1],pol[2],pol[4]],[pol[3],pol[4],pol[5]]])
              out_polarizability = 0.14818471147*sum(np.linalg.eigh(pol_mat)[0])/3.
            except:
              out_polarizability = missing
            search+=1
            #+unit conversion
            try:
              out_rotational_constants = [float(i) for i in out_rotational_constants]
            except:
              out_rotational_constants = missing
            out_rotational_constant = np.linalg.norm(out_rotational_constants)
            out_electronic_energy = float(out_electronic_energy)
            try:
              out_mulliken_charges = [float(i) for i in out_mulliken_charges]
            except:
              out_mulliken_charges = missing
            out_dipole_moment = float(out_dipole_moment)
            try:
              out_dipole_moments = [float(i) for i in out_dipole_moments]
            except:
              out_dipole_moments = missing
            continue       
        if search==1:
          ## ESP
          if re.search(" ESP charges:", line): #O4
            save_esp_charges = out_NAtoms+1
            save_something = "esp_charges"
            try:
              out_esp_charges = ["0"]*out_NAtoms
            except:
              out_esp_charges = missing
            continue
          if save_something == "esp_charges":
            try:
              if save_esp_charges<=out_NAtoms:
                out_esp_charges[out_NAtoms-save_esp_charges] = float(line.split()[2])
            except:
              out_esp_charges = missing
            save_esp_charges-=1
            if save_esp_charges == 0:
              save_something = ""
            continue
          #VIBRATIONAL FREQUENCIES
          if re.search("Frequencies", line): #V1
            if np.all(np.isnan(out_vibrational_frequencies)):
              out_vibrational_frequencies = []
            try:
              out_vibrational_frequencies += [float(i) for i in line.split()[2:]]
            except:
              out_vibrational_frequencies = missing
            continue
          if re.search("Kelvin.  Pressure", line): #V2/V3
            try:
              out_temperature = float(line.split()[1])
              out_pressure = float(line.split()[4])
            except:
              out_temperature = missing
              out_pressure = missing
            continue
          if re.search("Eigenvalues -- ", line): #V4
            if np.isnan(out_moments_of_inertia):
              out_moments_of_inertia = []
            try:
              out_moments_of_inertia += [float(i) for i in line.split()[2:]]
            except:
              out_moments_of_inertia = missing
            search+=1
            continue
          if re.search("Zero-point correction=", line):
            search+=1
        #THERMOCHEMISTRY
        if search==2:
          if re.search("Rotational symmetry number", line): #T1
            try:
              out_rotational_symmetry_number = float(line.split()[3])
            except:
              out_rotational_symmetry_number = missing
            continue
          if re.search("Zero-point correction=", line): #T2
            try:
              out_zero_point_correction = float(line.split()[2])
            except:
              out_zero_point_correction = missing
            continue
          if re.search("Thermal correction to Energy", line): #T3
            try:
              out_energy_thermal_correction = float(line.split()[4])
            except:
              out_energy_thermal_correction = missing
            continue
          if re.search("Thermal correction to Enthalpy", line): #T4
            try:
              out_enthalpy_thermal_correction = float(line.split()[4])
            except:
              out_enthalpy_thermal_correction = missing
            continue
          if re.search("Thermal correction to Gibbs Free Energy", line): #T5
            try:
              out_gibbs_free_energy_thermal_correction = float(line.split()[6])
            except:
              out_gibbs_free_energy_thermal_correction = missing
            continue
          if re.search("Sum of electronic and zero-point Energies", line): #T6
            try:
              out_zero_point_energy = float(line.split()[6])
            except:
              out_zero_point_energy = missing
            continue
          if re.search("Sum of electronic and thermal Energies", line): #T7
            try:
              out_internal_energy = float(line.split()[6])
            except:
              out_internal_energy = missing
            continue
          if re.search("Sum of electronic and thermal Enthalpies", line): #T8
            try:
              out_enthalpy_energy = float(line.split()[6])
            except:
              out_enthalpy_energy = missing
            continue
          if re.search("Sum of electronic and thermal Free Energies", line): #T9
            try:
              out_gibbs_free_energy = float(line.split()[7])
            except:
              out_gibbs_free_energy = missing
            search += 1
            continue
        #ENTROPY
        if search==3:
          if re.search("Total", line): #E1
            try:
              out_entropy = float(line.split()[3])
            except:
              out_entropy = missing
            search += 1 
            continue
        if Qforces == 1:
          if re.search("Center     Atomic                   Forces", line): #F1
            save_something = "forces"
            save_forces = -2
            out_forces = []
            continue
          if save_something == "forces":
            save_forces = save_forces + 1
            if save_forces > 0:
              try:
                out_forces.append(np.array([float(line.split()[2])/0.529177,float(line.split()[3])/0.529177,float(line.split()[4])/0.529177]))
              except:
                out_forces = missing
                save_something = ""
            if save_forces == out_NAtoms:
              save_something = ""
            continue
        if Qanharm == 1:
          if re.search("Fundamental Bands", line) and not save_anharm > 0:
            save_something = "anharm"
            save_anharm = -2
            out_anharm = []
            out_anharm_help = []
            continue
          if save_something == "anharm":
            save_anharm = save_anharm + 1
            if save_anharm > 0:
              try:
                out_anharm.append(float(line[4:].split()[3]))
                out_anharm_help.append(float(line[4:].split()[2]))
              except:
                out_anharm = missing
                save_something = ""
            if save_anharm == len(out_vibrational_frequencies):
              sorted_lists = sorted(zip(out_anharm_help, out_anharm))
              out_anharm = [item[1] for item in sorted_lists]
              save_something = ""
            continue 
      #SAVE
      clusters_df = df_add_iter(clusters_df, "log", "program", [str(cluster_id)], [out_program]) #PROGRAM
      clusters_df = df_add_iter(clusters_df, "log", "method", [str(cluster_id)], [out_method]) #METHOD
      clusters_df = df_add_iter(clusters_df, "log", "time", [str(cluster_id)], [out_time]) #TIME
      clusters_df = df_add_iter(clusters_df, "log", "termination", [str(cluster_id)], [out_termination]) #TERMINATION
      clusters_df = df_add_iter(clusters_df, "log", "charge", [str(cluster_id)], [out_charge]) #I1
      clusters_df = df_add_iter(clusters_df, "log", "multiplicity", [str(cluster_id)], [out_multiplicity]) #I2
      clusters_df = df_add_iter(clusters_df, "log", "NAtoms", [str(cluster_id)], [out_NAtoms]) #I3
      clusters_df = df_add_iter(clusters_df, "log", "rotational_constants", [str(cluster_id)], [out_rotational_constants]) #O1
      clusters_df = df_add_iter(clusters_df, "log", "rotational_constant", [str(cluster_id)], [out_rotational_constant]) #O2
      clusters_df = df_add_iter(clusters_df, "log", "sp_electronic_energy", [str(cluster_id)], [out_sp_electronic_energy]) #O3
      clusters_df = df_add_iter(clusters_df, "log", "electronic_energy", [str(cluster_id)], [out_electronic_energy]) #O3
      clusters_df = df_add_iter(clusters_df, "log", "mulliken_charges", [str(cluster_id)], [out_mulliken_charges]) #O4
      clusters_df = df_add_iter(clusters_df, "log", "esp_charges", [str(cluster_id)], [out_esp_charges]) #O4b
      clusters_df = df_add_iter(clusters_df, "log", "dipole_moment", [str(cluster_id)], [out_dipole_moment]) #O5
      clusters_df = df_add_iter(clusters_df, "log", "dipole_moments", [str(cluster_id)], [out_dipole_moments]) #O6
      clusters_df = df_add_iter(clusters_df, "log", "polarizability", [str(cluster_id)], [out_polarizability]) #O7
      clusters_df = df_add_iter(clusters_df, "log", "vibrational_frequencies", [str(cluster_id)], [out_vibrational_frequencies]) #V1
      clusters_df = df_add_iter(clusters_df, "log", "temperature", [str(cluster_id)], [out_temperature]) #V2
      clusters_df = df_add_iter(clusters_df, "log", "pressure", [str(cluster_id)], [out_pressure]) #V3
      clusters_df = df_add_iter(clusters_df, "log", "moments_of_inertia", [str(cluster_id)], [out_moments_of_inertia]) #V4
      clusters_df = df_add_iter(clusters_df, "log", "rotational_symmetry_number", [str(cluster_id)], [out_rotational_symmetry_number]) #T1
      clusters_df = df_add_iter(clusters_df, "log", "zero_point_correction", [str(cluster_id)], [out_zero_point_correction]) #T2
      clusters_df = df_add_iter(clusters_df, "log", "energy_thermal_correction", [str(cluster_id)], [out_energy_thermal_correction]) #T3
      clusters_df = df_add_iter(clusters_df, "log", "enthalpy_thermal_correction", [str(cluster_id)], [out_enthalpy_thermal_correction]) #T4
      clusters_df = df_add_iter(clusters_df, "log", "gibbs_free_energy_thermal_correction", [str(cluster_id)], [out_gibbs_free_energy_thermal_correction]) #T5
      clusters_df = df_add_iter(clusters_df, "log", "zero_point_energy", [str(cluster_id)], [out_zero_point_energy]) #T6
      clusters_df = df_add_iter(clusters_df, "log", "internal_energy", [str(cluster_id)], [out_internal_energy]) #T7
      clusters_df = df_add_iter(clusters_df, "log", "enthalpy_energy", [str(cluster_id)], [out_enthalpy_energy]) #T8
      clusters_df = df_add_iter(clusters_df, "log", "gibbs_free_energy", [str(cluster_id)], [out_gibbs_free_energy]) #T9
      clusters_df = df_add_iter(clusters_df, "log", "entropy", [str(cluster_id)], [out_entropy]) #E1
      if Qforces == 1:
        clusters_df = df_add_iter(clusters_df, "extra", "forces", [str(cluster_id)], [out_forces]) #F1
      if Qanharm == 1:
        clusters_df = df_add_iter(clusters_df, "extra", "anharm", [str(cluster_id)], [out_anharm]) 
    file.close()

  ###############
  ### XYZ #######
  ###############
  if path.exists(file_i_XYZ):
    try:
      out = read(file_i_XYZ)
    except:
      out = missing
  else:
    out = missing
  clusters_df = df_add_iter(clusters_df, "xyz", "structure", [str(cluster_id)], [out])

  ###############
  ### MRCC ######
  ###############
  if path.exists(file_i_MRCC):
    file = open(file_i_MRCC, "r")
    testMRCC = 0
    for i in range(3):
      line = file.readline()
      if re.search("MRCC program system", line):
        testMRCC = 1
        break
    if testMRCC == 1:
      from datetime import datetime
      out_program = "MRCC"                              #PROGRAM
      out_method = missing                               #METHOD
      out_time = missing                                 #TIME
      out_electronic_energy = missing                    #O3
      ######
      first_time=""
      for line in file:
        #TIME
        if re.search("\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* ",line):
          if first_time == "":
            first_time = line.split()[1]+" "+line.split()[2]
            print(datetime.strptime(first_time,"%Y-%m-%d %H:%M:%S"))
          else:
            try:
              this_time = line.split()[1]+" "+line.split()[2]
              out_time = datetime.strptime(this_time, "%Y-%m-%d %H:%M:%S")-datetime.strptime(first_time, "%Y-%m-%d %H:%M:%S")
              out_time = out_time.total_seconds() / 60
            except:
              out_time = missing
          continue   
        #EL.ENERGY
        if re.search("Total LNO-CCSD\(T\) energy with MP2 corrections", line): #O3
          try:
            out_electronic_energy = float(line.split()[7])
          except:
            out_electronic_energy = missing
          continue 
      clusters_df = df_add_iter(clusters_df, mrccextname, "program", [str(cluster_id)], [out_program]) #PROGRAM
      clusters_df = df_add_iter(clusters_df, mrccextname, "method", [str(cluster_id)], [out_method]) #METHOD
      clusters_df = df_add_iter(clusters_df, mrccextname, "time", [str(cluster_id)], [out_time]) #TIME
      clusters_df = df_add_iter(clusters_df, mrccextname, "electronic_energy", [str(cluster_id)], [out_electronic_energy]) #O3
  
  ###############
  ### ORCA ######
  ###############
  if path.exists(file_i_ORCA):
    file = open(file_i_ORCA, "r")
    testORCA = 0
    for i in range(130):
      line = file.readline()
      if re.search("O   R   C   A", line):
        testORCA = 1
        break
      if re.search("ORCA", line):
        testORCA = 1
        break
    if testORCA == 1:
      out_program = missing                              #PROGRAM
      out_method = missing                               #METHOD
      out_time = missing                                 #TIME
      out_termination = 0                                #TERMINATION
      out_charge = missing                               #I1
      out_multiplicity = missing                         #I2
      out_NAtoms = missing                               #I3
      out_rotational_constants = missing                 #O1
      out_rotational_constant = missing                  #O2
      out_sp_electronic_energy = missing                 #O3x
      out_electronic_energy = missing                    #O3
      out_mulliken_charges = missing                     #O4
      out_dipole_moment = missing                        #O5
      out_dipole_moments = missing                       #O6
      out_polarizability = missing                       #O7
      out_vibrational_frequencies = missing              #V1
      out_temperature = missing                          #V2
      out_pressure = missing                             #V3
      out_moments_of_inertia = missing                   #V4
      out_rotational_symmetry_number = missing           #T1
      out_zero_point_correction = missing                #T2
      out_energy_thermal_correction = missing            #T3
      out_enthalpy_thermal_correction = missing          #T4
      out_gibbs_free_energy_thermal_correction = missing #T5
      out_zero_point_energy = missing                    #T6
      out_internal_energy = missing                      #T7
      out_enthalpy_energy = missing                      #T8
      out_gibbs_free_energy = missing                    #T9
      out_entropy = missing                              #E1
      if Qforces == 1:
        out_forces = missing                             #F1
      search=-1
      save_mulliken_charges=0
      skip_0_freq=1
      save_something=""
      for line in file:
        #TIME
        if re.search("TOTAL RUN TIME",line):
          if np.isnan(out_time):
            out_time = 0
          try:
            out_time += float(line.split()[3])*24*60+float(line.split()[5])*60+float(line.split()[7])+float(line.split()[9])/60+float(line.split()[11])/6000
          except:
            out_time = missing
          continue
        #TERMINATION
        if re.search("ORCA TERMINATED NORMALLY",line):
          out_termination += 1
          continue
        #INITIAL SEARCH
        if search==-1:
          if re.search("Program Version", line): #PROGRAM
            try:
              out_program = "ORCA_" + str(line.split()[2])
            except:
              out_program = missing
            continue
          if re.search(r"^\|.*>.*!", line): #METHOD
            try:
              if pd.isna(out_method):
                out_method = "_".join(sorted(line.split(">")[1].lower().split()))
              else:
                out_method += "_"+"_".join(sorted(line.split(">")[1].lower().split()))
            except:
              out_method = missing
          if re.search("Number of atoms", line): #I3
            try:
              out_NAtoms = int(line.split()[4])
            except:
              out_NAtoms = missing
            continue
          if re.search(" Total Charge", line): #I1/I2
            try:
              out_charge = int(line.split()[4])
            except:
              out_charge = missing
            continue
          if re.search(" Multiplicity", line): #I1/I2
            try:
              out_multiplicity = int(line.split()[3])
            except:
              out_multiplicity = missing
            continue
        #OPTIMIZATION SEARCH
        ## I will convert the strings into floats/ints later !!!!!!!
        if search==0 or search==-1:
          if re.search("FINAL SINGLE POINT ENERGY", line): #O3
            try:
              out_electronic_energy = float(line.split()[4])
            except:
              out_electronic_energy = missing
            if np.isnan(out_sp_electronic_energy):
              out_sp_electronic_energy = out_electronic_energy
              search+=1
            continue
          ## MULLIKEN
          if re.search("MULLIKEN ATOMIC CHARGES", line): #O4
            save_mulliken_charges = out_NAtoms+1
            save_something = "mulliken_charges"
            try:
              out_mulliken_charges = ["0"]*out_NAtoms
            except:
              out_mulliken_charges = missing
            continue
          if save_something == "mulliken_charges":
            try:
              if save_mulliken_charges<=out_NAtoms:
                out_mulliken_charges[out_NAtoms-save_mulliken_charges] = float(line.split()[3])
            except:
              out_mulliken_charges = missing
            save_mulliken_charges-=1
            if save_mulliken_charges == 0:
              save_something = ""
            continue
          ## DIPOLE
          if re.search('Total Dipole Moment', line): #O5/O6
            try:
              out_dipole_moments = [float(line.split()[4])/0.393456, float(line.split()[5])/0.393456, float(line.split()[6])/0.393456]
            except:
              out_dipole_moments = missing
            search+=1
            continue
        if search==1:
          if re.search('Magnitude \(Debye\)', line): #O5/O6
            try:            
              out_dipole_moment = float(line.split()[3])
            except:
              out_dipole_moment = missing
            continue
          if re.search("Rotational constants in MHz", line): #O1/O2
            try:
              out_rotational_constants = [float(line.split()[5])/1000,float(line.split()[6])/1000,float(line.split()[7])/1000]
            except:
              out_rotational_constants = missing
            out_rotational_constant = np.linalg.norm(out_rotational_constants)
            search+=1
            continue
        if search==2:
          #VIBRATIONAL FREQUENCIES
          if re.search("cm\*\*-1", line): #V1
            if re.search("0.00 cm\*\*-1", line) and skip_0_freq == 1: 
              continue
            else:
              skip_0_freq = 0
            if np.all(np.isnan(out_vibrational_frequencies)):
              out_vibrational_frequencies = []
            try:
              out_vibrational_frequencies += [float(line.split()[1])]
            except:
              out_vibrational_frequencies = missing
            continue
          if re.search("NORMAL MODES",line):
            search+=1
            continue
        if search==3:
          if re.search("Temperature         ...", line): #V2
            try:
              out_temperature = float(line.split()[2])
            except:
              out_temperature = missing
            continue
          if re.search("Pressure            ...", line): #V3
            try:
              out_pressure = float(line.split()[2])
            except:
              out_pressure = missing
            continue
          if re.search("Zero point energy                ...", line): #T2
            try:
              out_zero_point_correction = float(line.split()[4])
            except:
              out_zero_point_correction = missing
            continue
          if re.search("Total thermal energy", line): #T7
            try:
              out_internal_energy = float(line.split()[3])
            except:
              out_internal_energy = missing
            continue
          if re.search("Symmetry Number", line): #T1
            try:
              out_rotational_symmetry_number = float(line.split()[5])
            except:
              out_rotational_symmetry_number = missing
            continue
          if re.search("Final entropy term", line): #E1
            try:
              out_entropy = float(line.split()[4])*1000*627.503/out_temperature
            except:
              out_entropy = missing
            continue
          if re.search("Final Gibbs free energy", line): #T9
            try:
              out_gibbs_free_energy = float(line.split()[5])
            except:
              out_gibbs_free_energy = missing
            continue
          if re.search("Total Enthalpy", line): #T8
            try:
              out_enthalpy_energy = float(line.split()[3])
            except:
              out_enthalpy_energy = missing
            continue
        if Qforces == 1: 
          if re.search("CARTESIAN GRADIENT", line): #F1
            save_something = "forces"
            save_forces = -2
            out_forces = []
            continue
          if save_something == "forces":
            save_forces = save_forces + 1
            if save_forces > 0:
              try:
                out_forces.append(np.array([float(line.split()[3])/0.529177,float(line.split()[4])/0.529177,float(line.split()[5])/0.529177]))
              except:
                out_forces = missing
                save_something = ""
            if save_forces == out_NAtoms:
              save_something = ""
            continue
      #FINISH MISSING
      try:
        out_energy_thermal_correction = out_internal_energy - out_electronic_energy
      except:
        out_energy_thermal_correction = missing
      try: 
        out_enthalpy_thermal_correction = out_enthalpy_energy - out_electronic_energy
      except:
        out_enthalpy_thermal_correction = missing
      try:
        out_gibbs_free_energy_thermal_correction = out_gibbs_free_energy - out_electronic_energy
      except:
        out_gibbs_free_energy_thermal_correction = missing
      try:
        out_zero_point_energy = out_zero_point_correction + out_electronic_energy
      except:
        out_zero_point_energy = missing
      out_moments_of_inertia = missing
      out_polarizability = missing
      #SAVE
      clusters_df = df_add_iter(clusters_df, orcaextname, "program", [str(cluster_id)], [out_program]) #PROGRAM
      clusters_df = df_add_iter(clusters_df, orcaextname, "method", [str(cluster_id)], [out_method]) #METHOD
      clusters_df = df_add_iter(clusters_df, orcaextname, "time", [str(cluster_id)], [out_time]) #TIME
      clusters_df = df_add_iter(clusters_df, orcaextname, "termination", [str(cluster_id)], [out_termination]) #TERMINATION
      clusters_df = df_add_iter(clusters_df, orcaextname, "charge", [str(cluster_id)], [out_charge]) #I1
      clusters_df = df_add_iter(clusters_df, orcaextname, "multiplicity", [str(cluster_id)], [out_multiplicity]) #I2
      clusters_df = df_add_iter(clusters_df, orcaextname, "NAtoms", [str(cluster_id)], [out_NAtoms]) #I3
      clusters_df = df_add_iter(clusters_df, orcaextname, "rotational_constants", [str(cluster_id)], [out_rotational_constants]) #O1
      clusters_df = df_add_iter(clusters_df, orcaextname, "rotational_constant", [str(cluster_id)], [out_rotational_constant]) #O2
      clusters_df = df_add_iter(clusters_df, orcaextname, "sp_electronic_energy", [str(cluster_id)], [out_sp_electronic_energy]) #O3
      clusters_df = df_add_iter(clusters_df, orcaextname, "electronic_energy", [str(cluster_id)], [out_electronic_energy]) #O3
      clusters_df = df_add_iter(clusters_df, orcaextname, "mulliken_charges", [str(cluster_id)], [out_mulliken_charges]) #O4
      clusters_df = df_add_iter(clusters_df, orcaextname, "dipole_moment", [str(cluster_id)], [out_dipole_moment]) #O5
      clusters_df = df_add_iter(clusters_df, orcaextname, "dipole_moments", [str(cluster_id)], [out_dipole_moments]) #O6
      clusters_df = df_add_iter(clusters_df, orcaextname, "polarizability", [str(cluster_id)], [out_polarizability]) #O7
      clusters_df = df_add_iter(clusters_df, orcaextname, "vibrational_frequencies", [str(cluster_id)], [out_vibrational_frequencies]) #V1
      clusters_df = df_add_iter(clusters_df, orcaextname, "temperature", [str(cluster_id)], [out_temperature]) #V2
      clusters_df = df_add_iter(clusters_df, orcaextname, "pressure", [str(cluster_id)], [out_pressure]) #V3
      clusters_df = df_add_iter(clusters_df, orcaextname, "moments_of_inertia", [str(cluster_id)], [out_moments_of_inertia]) #V4
      clusters_df = df_add_iter(clusters_df, orcaextname, "rotational_symmetry_number", [str(cluster_id)], [out_rotational_symmetry_number]) #T1
      clusters_df = df_add_iter(clusters_df, orcaextname, "zero_point_correction", [str(cluster_id)], [out_zero_point_correction]) #T2
      clusters_df = df_add_iter(clusters_df, orcaextname, "energy_thermal_correction", [str(cluster_id)], [out_energy_thermal_correction]) #T3
      clusters_df = df_add_iter(clusters_df, orcaextname, "enthalpy_thermal_correction", [str(cluster_id)], [out_enthalpy_thermal_correction]) #T4
      clusters_df = df_add_iter(clusters_df, orcaextname, "gibbs_free_energy_thermal_correction", [str(cluster_id)], [out_gibbs_free_energy_thermal_correction]) #T5
      clusters_df = df_add_iter(clusters_df, orcaextname, "zero_point_energy", [str(cluster_id)], [out_zero_point_energy]) #T6
      clusters_df = df_add_iter(clusters_df, orcaextname, "internal_energy", [str(cluster_id)], [out_internal_energy]) #T7
      clusters_df = df_add_iter(clusters_df, orcaextname, "enthalpy_energy", [str(cluster_id)], [out_enthalpy_energy]) #T8
      clusters_df = df_add_iter(clusters_df, orcaextname, "gibbs_free_energy", [str(cluster_id)], [out_gibbs_free_energy]) #T9
      clusters_df = df_add_iter(clusters_df, orcaextname, "entropy", [str(cluster_id)], [out_entropy]) #E1
      if Qforces == 1:
        clusters_df = df_add_iter(clusters_df, "extra", "forces", [str(cluster_id)], [out_forces]) #F1
    file.close()

  ###############
  ### TURBOMOLE #
  ###############
  if path.exists(file_i_TURBOMOLE):
    file = open(file_i_TURBOMOLE, "r")
    testTURBOMOLE = 0
    for i in range(5):
      if re.search("TURBOMOLE", file.readline()):
        testTURBOMOLE = 1
        break
    if testTURBOMOLE == 1:
      out = missing
      for line in file:
        if re.search("Final CCSD\(F12\*\)\(T\) energy", line):
          out = float(line.split()[5])
        if re.search("Final MP2 energy", line):
          out = float(line.split()[5])
      file.close()
      clusters_df = df_add_iter(clusters_df, turbomoleextname, "electronic_energy", [str(cluster_id)], [out])
####################################################################################################
####################################################################################################

if len(addcolumn) > 0:
  dtype = np.dtype([("label", "U100"), ("value", float)])
  for i in range(len(addcolumn)):
    if addcolumn[i][0] == "SP" and addcolumn[i][1][-3:] == "pkl":
      loadedclustersdf = pd.read_pickle(addcolumn[i][1])
      loadedmatrix = np.array([loadedclustersdf["info"]["file_basename"],loadedclustersdf["log"]["electronic_energy"]])
      loadeddictionary = {loadedmatrix[0][idx] : loadedmatrix[1][idx] for idx in range(len(loadedmatrix[1]))}
      tobeadded = clusters_df["info"]['file_basename'].map(loadeddictionary)
      clusters_df = df_add_iter(clusters_df, "out", "electronic_energy", tobeadded.index, tobeadded.values) 
    else:
      loadedmatrix = np.loadtxt(addcolumn[i][1], usecols=(0, 1), unpack=True, dtype=dtype)
      loadeddictionary = {loadedmatrix[0][idx] : loadedmatrix[1][idx] for idx in range(len(loadedmatrix[1]))}
      tobeadded = clusters_df["info"]['file_basename'].map(loadeddictionary)
      clusters_df = df_add_iter(clusters_df, "extra", addcolumn[i][0], tobeadded.index, tobeadded.values) 
    
####################################################################################################
####################################################################################################

if Qpresplit > 0:
  clusters_df = clusters_df.sample(n=Qpresplit, random_state=42)

if Qindex != "-1":
  try:
    clusters_df = eval(f"clusters_df[{Qindex}]")
  except:
    try: 
      Qindex = Qindex+":"+str(int(Qindex)+1)
      clusters_df = eval(f"clusters_df[{Qindex}]")
    except:
      print("Error selecting given index")
      exit()

####################################################################################################
####################################################################################################
## HERE I MODIFY THE INPUT
if Qunderscore == 1:
  def replace_first_occurrence(string, old_char, new_char):
    index = string.find(old_char)  # Find the index of the first occurrence
    if index != -1:  # If the character is found
     string = string[:index] + new_char + string[index+1:]  # Replace the character
    return string
  
  for cluster_id in clusters_df.index:
    clusters_df.loc[cluster_id]['info']['file_basename'] = replace_first_occurrence(clusters_df.loc[cluster_id]['info']['file_basename'],"_","-")

if Qmodify > 0:
  if Qrename == 1:
    # sorted cluster type
    if Qclustername == 1:
      #print(clusters_df)
      l1 = []
      l2 = []
      l3 = []
      l4 = []
      for cluster_id in clusters_df.index:
        file_i_BASE = clusters_df.loc[cluster_id]['info']['file_basename']
        cluster_type_array = seperate_string_number(file_i_BASE.split("-")[0])
        #if np.mod(len(cluster_type_array),2) == 0:
        if is_nameable(cluster_type_array):
          #print(cluster_id)
          for QrenameONE in QrenameWHAT: 
            cluster_type_array = [w.replace(QrenameONE[0],QrenameONE[1]) if QrenameONE[0] == w else w for w in cluster_type_array]
          l4.append("".join(cluster_type_array)+"-"+file_i_BASE.split("-")[1]) 
          cluster_type_2array_sorted = sorted([cluster_type_array[i:i + 2] for i in range(0, len(cluster_type_array), 2)],key=lambda x: x[1])
          cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
          cluster_type_string = zeros(cluster_type_array_sorted)
          #clusters_df = df_add_iter(clusters_df, "info", "cluster_type", [str(cluster_id)], [cluster_type_string])
          #
          #components = re.split('(\d+)', file_i_BASE.split("-")[0])[1:][1::2]
          components = cluster_type_array[1::2]
          #clusters_df = df_add_iter(clusters_df, "info", "components", [str(cluster_id)], [components])
          component_ratio = [int(i) for i in re.split('(\d+)', file_i_BASE.split("-")[0])[1:][0::2]]
          #clusters_df = df_add_iter(clusters_df, "info", "component_ratio", [str(cluster_id)], [component_ratio])
          l1.append(cluster_type_string)
          l2.append(components)
          l3.append(component_ratio)
      clusters_df =  df_add_iter(clusters_df, "info", "cluster_type", [str(missing)], l1)
      clusters_df =  df_add_iter(clusters_df, "info", "components", [str(missing)], l2)
      clusters_df =  df_add_iter(clusters_df, "info", "component_ratio", [str(missing)], l3)
      clusters_df =  df_add_iter(clusters_df, "info", "file_basename", [str(missing)], l4)
      #print(clusters_df)   

if Qiamifo > 0:
  predashstrings = [i.split("-")[0] for i in clusters_df["info"]["file_basename"].values]
  #uniquestrings = list(set(predashstrings))
  #import string,random
  #randomstrings = ["1"+''.join(random.choices(string.ascii_lowercase, k=4)) for _ in range(len(uniquestrings))]
  #newlist = [missing for i in range(len(predashstrings))]
  #for i in uniquestrings:
  #  positions = [index for index, value in enumerate(predashstrings) if value == i]
  #  for position in positions:
  #    newlist[position] = randomstrings[uniquestrings.index(i)]
  clusters_df = df_add_append(clusters_df, "info", "cluster_type", clusters_df.index, predashstrings)

if Qdrop != "0":
  if Qdrop in clusters_df:
    clusters_df = clusters_df.drop([Qdrop], axis=1, level=0) 

if Qout2log == 1:
  clusters_df = clusters_df.rename({'out': 'log'}, axis='columns')


####################################################################################################
####################################################################################################
## HERE I MANUPILATE WHICH CLUSTERS PASSES THROUGH!! NOT THE SAME AS FILTERING
# EXTRACT
def dash(input_array):
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "-":
      partbefore = input_array[0:element-1]
      partafter = input_array[element+2:]
      output_array_1 = []
      try:
        num1=int(input_array[element-1])
        num2=int(input_array[element+1])+1
        output_array = []
      except:
        break
      for thenumber in range(int(input_array[element-1]),int(input_array[element+1])+1):
        output_array_1.append(partbefore+[str(thenumber)]+partafter)
      for i in output_array_1:
        for j in dash(i): 
          output_array.append(j)
      break
  return output_array

def dash_comment(input_array):
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "-" or input_array[element] == "_":
      output_array = []
      partbefore = input_array[0:element]
      output_array.append(partbefore)
      break
  return output_array

def comma(input_array):
  output_array = []
  output_string = ""
  bracket_val = 0
  for element in range(len(input_array)):
    if input_array[element] == "(":
      bracket_val += 1
    if input_array[element] == ")":
      bracket_val -= 1
    if input_array[element] == ",":
      if bracket_val == 0:
        output_array.append(output_string)
        output_string = ""
        continue 
    output_string += input_array[element]
  output_array.append(output_string)
  return output_array

def comma2(input_array):
  output_array = [input_array]
  for element in range(len(input_array)):
    if input_array[element] == "(":
      elementL = element
    if input_array[element] == ")":
      output_array = []
      elementR = element
      partbefore = input_array[0:elementL]
      partafter = input_array[elementR+1:]
      thenumbers = listToString(input_array[elementL+1:elementR],"").split(",")
      for thenumber in thenumbers:
        for j in comma2(partbefore+[str(thenumber)]+partafter):
          output_array.append(j)
      break
  return output_array

def unique(list1):
    # initialize a null list
    unique_list = []
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # return list
    return unique_list

if Qextract > 0:
  #print(Pextract)
  Pextract_comma = []
  for extract_i in Pextract:
    comma_separated = comma(extract_i)
    for separated_i in range(len(comma_separated)):
      Pextract_comma.append(comma_separated[separated_i])
  #print(Pextract_comma)    
  if Qclustername == 0:
    Pextract_ultimate = unique(Pextract_comma)
  else:
    #print(Pextract_comma)    
    Pextract_dash = []
    for extract_i in Pextract_comma:
      separated = seperate_string_number(extract_i)
      dash_separated = dash(separated)
      for separated_i in range(len(dash_separated)):
        Pextract_dash.append(dash_separated[separated_i])
    #print([listToString(i,"") for i in Pextract_dash])
    Pextract_dash = [ dash_comment(i)[0] for i in Pextract_dash ]
    #print(Pextract_dash)
    #print([listToString(i,"") for i in Pextract_dash])
    Pextract_comma2 = []
    for extract_i in Pextract_dash:
      comma2_separated = comma2(extract_i)
      for separated_i in range(len(comma2_separated)):
        Pextract_comma2.append(comma2_separated[separated_i])
    #print(Pextract_comma2) #[['1', 'sa', '2', '*']]
    ##remove asterix
    #Qasterix=-1
    #
    #for extract_i in Pextract_comma2:
    #    
    #
    ##
    Pextract_sorted = []
    for extract_i in Pextract_comma2:
      if len(extract_i) > 1:
        array_sorted = sorted([extract_i[i:i + 2] for i in range(0, len(extract_i), 2)],key=lambda x: x[1])
        Pextract_sorted.append([item for sublist in array_sorted for item in sublist])
      else:
        Pextract_sorted.append(extract_i)
        Qclustername = 0
    #print(Pextract_sorted)
    Pextract_final = []
    for extract_i in Pextract_sorted:
      corrected = zeros(extract_i)
      if len(corrected) > 0:
        Pextract_final.append(corrected)
    Pextract_ultimate = unique(Pextract_final)
    #print(Pextract_ultimate)
  newclusters_df = pd.DataFrame()
  #print(Pextract_ultimate)
  #exit()
  for extract_i in Pextract_ultimate: 
    #print(extract_i)
    if Qclustername == 0:
      extracted_df = clusters_df[clusters_df["info"]["file_basename"].values == extract_i].copy()
      #print(len(extracted_df))
    else:
      if "*" not in extract_i:
        extracted_df = clusters_df[clusters_df["info"]["cluster_type"].values == extract_i].copy()
      else:
        try:
          asterix_position=seperate_string_number(extract_i)[1::2].index('*')
          new_extract_i=seperate_string_number(extract_i)[:2*asterix_position]+seperate_string_number(extract_i)[asterix_position*2+1+1:]
          try:
            howmany=int(extract_i[asterix_position*2])
          except:
            howmany=-1
          what=["whatever"]
        except:
          what=[]
          while 1==1:
            try:
              asterix_position=seperate_string_number(extract_i)[::2].index('*')
              print(asterix_position)
              try:
                if len(asterix_position)>1:
                  asterix_position=asterix_position[0]
              except:
                howmany=-1
              new_extract_i=seperate_string_number(extract_i)[:2*asterix_position]+seperate_string_number(extract_i)[asterix_position*2+1+1:]
              howmany=-1
              what.append(seperate_string_number(extract_i)[asterix_position*2+1])
              extract_i=new_extract_i
            except:
              break
        def my_special_compare_with_asterix(string1,string2_array,howmany,what):
          string1_array = seperate_string_number(string1)
          #string2_array = seperate_string_number(string2)
          tothowmany=howmany
          for indx in range(int(len(string1_array)/2)):
            #print(string2_array)
            if string1_array[1::2][indx] in string2_array[1::2]:
              found_position = string2_array[1::2].index(string1_array[1::2][indx])
              if string1_array[::2][indx] != string2_array[::2][found_position]:
                return False
              else:
                string2_array=string2_array[:2*found_position]+string2_array[2*found_position+2:] 
                continue
            else:
              if what[0] == "whatever":
                tothowmany-=int(string1_array[::2][indx])
                if howmany>0 and tothowmany<0:
                  return False
                else:
                  continue
              elif string1_array[1::2][indx] in what:
                tothowmany-=int(string1_array[::2][indx])
                continue
              else:
                return False
          if tothowmany == 0 or howmany<0:
            if len(string2_array) > 0:
              return False
            else:
              return True
          else:
            return False
        extracted_df = clusters_df[[my_special_compare_with_asterix(value_i,new_extract_i,howmany,what) for value_i in clusters_df["info"]["cluster_type"].values]].copy()
              
    if len(extracted_df) > 0:
      if len(newclusters_df) == 0:
        newclusters_df = extracted_df.copy()
      else:
        #print(newclusters_df.index)
        #print(newclusters_df["info","file_basename"])
        #newclusters_df = newclusters_df.append(extracted_df.copy())
        newclusters_df = pd.concat([newclusters_df, extracted_df.copy()], axis=0, ignore_index=True)
        # pd.concat([clusters_df,newclusters_df],axis=0, ignore_index=True)
        #try:
        #print(newclusters_df.loc['19268'])
        #print(extracted_df)
        #except: 
        #  print("0")
        #print(newclusters_df)
  prelen=len(clusters_df)
  if Qextract == 2:
    a1 = pd.Index(clusters_df.index)
    a2 = pd.Index(newclusters_df.index)
    clusters_df = clusters_df.iloc[a1.difference(a2, sort=False)]
  else:
    #print("here")
    #print(clusters_df.loc['19268'])
    #print(newclusters_df.loc['19268'])
    clusters_df = newclusters_df.copy()
  if Qout == 1:
    print("Extracting: "+str(prelen)+" --> "+str(len(clusters_df)))
  #clusters_df.reindex([str(j) for j in range(len(clusters_df))])


### REACTED ###
if Qreacted > 0:
  dt = np.dtype(object)
  #Are there some cluster types which I should distinguish?
  if Qclustername != 0:
    unique_cluster_types = np.unique(clusters_df["info"]["cluster_type"].values)
    cluster_subsets = []
    for unique_cluster_type in unique_cluster_types:
      indexes = clusters_df[clusters_df["info"]["cluster_type"]==unique_cluster_type].index.values
      cluster_subsets.append(indexes)
  else:
    cluster_subsets = []
    indexes = clusters_df.index.values
    cluster_subsets.append(indexes)   
  for k0 in range(len(cluster_subsets)):
    #print(k0)
    all_molecules = []
    for k in range(len(cluster_subsets[k0])):
      #print(k)
      b = clusters_df["xyz"]["structure"][cluster_subsets[k0][k]]
      if pd.isna(b):
        all_molecules.append(str(missing))
      else:
        p = b.positions
        symb = np.array(b.get_chemical_symbols())
        ind = [i != 'test' for i in symb]
    
        dist = lambda p1, p2: np.sqrt(np.sum(((p1-p2)**2)))
        dm = np.asarray([[dist(p1, p2) for p2 in p[ind]] for p1 in p[ind]])
    
        def bonded(x,xA,xB):
          if xA == "C" and xB == "N" and x < 1.75:
            return 1
          elif xA == "N" and xB == "C" and x < 1.75:
            return 1
          elif xA == "N" and xB == "N" and x < 1.5:
            return 1
          elif xA == "S" and xB == "O" and x < 1.9:
            return 1
          elif xA == "O" and xB == "S" and x < 1.9:
            return 1
          elif xA == "O" and xB == "N" and x < 1.9:
            return 1
          elif xA == "N" and xB == "O" and x < 1.9:
            return 1
          elif xA == "H" or xB == "H":
            if xA == "H" and xB == "H" and x < 0.8:
              return 1
            else:
              return 0
          elif x < bonddistancethreshold:
            return 1
          else:
            return 0

        bm = np.array([[bonded(dm[i][j],symb[ind][i],symb[ind][j]) for j in range(len(dm[i]))] for i in range(len(dm))])
        #print(bm) 

        test = 0
        choosing_list = range(len(bm))
        #print(choosing_list)
        molecules=[]
        if len(choosing_list) == 0:
          test = 1
        while test == 0:
          selected = [choosing_list[0]]
          #print(selected)
          test_chosen = 0
          j = -1
          while test_chosen == 0:
            j += 1
            #print(str(j)+"/"+str(len(bm)))
            #print(selected)
            chosen = selected[j]
            #print([chosen,j,len(selected)-1,"checking",choosing_list])
            for i in choosing_list:
              
              if bm[chosen][i] == 1 and chosen != i and (not (i in selected)):
                selected.append(i)
                #print([i,j,selected,len(selected),bm[choosing_list[j]][i],choosing_list[j],choosing_list,(i in selected)])
            #print("I am here")
            if len(selected)-1 == j:
              test_chosen = 1
              #print("DONE")
            #print("----------------")
          #print(molecules)
          molecules.append("".join(np.sort([symb[ind][i] for i in selected])))
          #print(np.sort(molecules))
          choosing_list = [i for i in choosing_list if i not in selected]
          if len(choosing_list) == 0:
            test = 1
        #print(str(np.sort(np.array(molecules,dtype = dt))))
#print(molecules)
        all_molecules.append(str(np.sort(np.array(np.sort(molecules),dtype = dt))))
  
    def most_frequent(List):
        return max(set(List), key = List.count)
    
    #print(all_molecules)
    #print(most_frequent(all_molecules))
    mf = most_frequent(all_molecules)
    
    #print(mf)
    #print(all_molecules[0])
    if Qreacted != 2:
      nind = [ i == mf for i in all_molecules]
    else:
      nind = [ i != mf for i in all_molecules]
    #print(all_molecules)
    #print(nind)
    if k0 == 0:
      selclusters_df = clusters_df.loc[cluster_subsets[0]][nind].copy()
    else:
      selclusters_df = selclusters_df.append(clusters_df.loc[cluster_subsets[k0]][nind].copy())
  clusters_df = selclusters_df.copy()
  ### END OF REACTED ###

## SAVE OUTPUT.pkl ##
#print(clusters_df.index)
#exit()
#new_clusters_df = clusters_df.reindex([str(j) for j in range(len(clusters_df))], copy=False)
#print(new_clusters_df)
if Qout > 0:
  original_clusters_df = clusters_df.copy()

## PRE_EXTRACT_DATA MANIPULATION ##
if Qqha == 1:
  h = 6.626176*10**-34 #m^2 kg s^-1
  R = 1.987 #cal/mol/K #=8.31441
  k = 1.380662*10**-23 #m^2 kg s^-2 K^-1
  if Qanh != "1":
    # VIBRATIONAL FREQ MODIFICATION e.g. anharmonicity (vib.freq.,ZPE,ZPEc,U,Uc,H,Hc,S // G,Gc)
    for i in range(len(clusters_df)):
      try:
        test = len(clusters_df["xyz"]["structure"][i].get_atomic_numbers())
      except:
        test = 0
        lf = 0
      if test != 1:
        if test != 0:
          try:
            QtOLD = clusters_df["log","temperature"].values[i]
            if pd.isna(clusters_df["log"]["vibrational_frequencies"].values[i]).any():
              continue
          except:
            QtOLD = clusters_df["log","temperature"].values[i]
          try:
            lf = float(clusters_df["log"]["vibrational_frequencies"].values[i][0])
          except:
            lf = 0
        if lf <= 0:
          clusters_df["log","entropy"][i] = missing
          clusters_df["log","enthalpy_energy"][i] = missing
          clusters_df["log","enthalpy_thermal_correction"][i] = missing
          clusters_df["log","internal_energy"][i] = missing
          clusters_df["log","energy_thermal_correction"][i] = missing
          clusters_df["log","zero_point_correction"][i] = missing
          clusters_df["log","zero_point_energy"][i] = missing
          continue
        try:
          Sv_OLD = np.sum([R*h*vib*2.99793*10**10/k/QtOLD/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K
          Ev_OLD = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
        except:
          Sv_OLD = missing
          Ev_OLD = missing
        #
        if Qanh != "anh" and Qanh != "anh2":
          try:
            clusters_df["log","vibrational_frequencies"][i] = [float(Qanh) * j for j in clusters_df["log","vibrational_frequencies"].values[i]]
          except:
            clusters_df["log","vibrational_frequencies"][i] = [missing]
        else:
          def replace_by_nonnegative(new, orig, q):
            if q == 0:
              mask = np.array(new) > 0
            else:
              mask = [ new[ii] > 0 and new[ii] < orig[ii] for ii in range(len(new)) ]
            orig = np.array(orig)
            new = np.array(new)
            orig[mask] = new[mask]
            return list(orig)
          #print(len(clusters_df["extra","anharm"].values[i]) == len(clusters_df["log","vibrational_frequencies"].values[i]))
          try:
            if Qanh == "anh":
              clusters_df["log","vibrational_frequencies"][i] = replace_by_nonnegative(clusters_df["extra","anharm"].values[i],clusters_df["log","vibrational_frequencies"].values[i],0)
            else:
              clusters_df["log","vibrational_frequencies"][i] = replace_by_nonnegative(clusters_df["extra","anharm"].values[i],clusters_df["log","vibrational_frequencies"].values[i],1)
          except:
            clusters_df["log","vibrational_frequencies"][i] = [missing]
        #
        try:
          Sv = np.sum([R*h*vib*2.99793*10**10/k/QtOLD/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K  
          Ev = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
        except:
          Sv = missing
          Ev = missing
        ###
        clusters_df["log","zero_point_correction"][i] = np.sum([0.5*h*vib*2.99793*10**10 for vib in clusters_df["log","vibrational_frequencies"][i]])*0.00038088*6.022*10**23/1000
        clusters_df["log","zero_point_energy"][i] = clusters_df["log","electronic_energy"][i] + clusters_df["log","zero_point_correction"][i]  
        clusters_df["log","internal_energy"][i] += (Ev - Ev_OLD)/1000/627.503    
        clusters_df["log","energy_thermal_correction"][i] += (Ev - Ev_OLD)/1000/627.503    
        clusters_df["log","enthalpy_energy"][i] += (Ev - Ev_OLD)/1000/627.503
        clusters_df["log","enthalpy_thermal_correction"][i] += (Ev - Ev_OLD)/1000/627.503
        clusters_df["log","entropy"][i] += Sv - Sv_OLD      
        ###
    
  
  # NEW TEMPERATURE (T,S,H,Hc,U,Uc // G,Gc)
  if ~np.isnan(Qt):
    for i in range(len(clusters_df)):
      try:
        if pd.isna(clusters_df["log"]["vibrational_frequencies"].values[i]):
          clusters_df["log","temperature"][i] = Qt
          continue
      except:
        QtOLD = clusters_df["log","temperature"].values[i]
      if Qt != QtOLD:
        try:
          lf = float(clusters_df["log"]["vibrational_frequencies"].values[i][0])
        except:
          lf = 0
        if lf <= 0:
          clusters_df["log","temperature"][i] = Qt
          clusters_df["log","entropy"][i] = missing
          clusters_df["log","enthalpy_energy"][i] = missing
          clusters_df["log","enthalpy_thermal_correction"][i] = missing
          clusters_df["log","internal_energy"][i] = missing
          clusters_df["log","energy_thermal_correction"][i] = missing
          continue
        Sv_OLD = np.sum([R*h*vib*2.99793*10**10/k/QtOLD/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K
        Sv = np.sum([R*h*vib*2.99793*10**10/k/Qt/(np.exp(h*vib*2.99793*10**10/k/Qt)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/Qt)) for vib in clusters_df["log"]["vibrational_frequencies"].values[i]]) #cal/mol/K
        Ev_OLD = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
        Ev = np.sum([R*h*vib*2.99793*10**10/k/(np.exp(h*vib*2.99793*10**10/k/Qt)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df["log"]["vibrational_frequencies"].values[i]])
        ###
        clusters_df["log","temperature"][i] = Qt
        clusters_df["log","entropy"][i] += Sv - Sv_OLD + 4*R*np.log(Qt/QtOLD)
        clusters_df["log","enthalpy_energy"][i] += (Ev - Ev_OLD + 4*R*(Qt-QtOLD))/1000/627.503
        clusters_df["log","enthalpy_thermal_correction"][i] += (Ev - Ev_OLD + 4*R*(Qt-QtOLD))/1000/627.503
        clusters_df["log","internal_energy"][i] += (Ev - Ev_OLD + 3*R*(Qt-QtOLD))/1000/627.503
        clusters_df["log","energy_thermal_correction"][i] += (Ev - Ev_OLD + 3*R*(Qt-QtOLD))/1000/627.503
        ###
  
  # LOW VIBRATIONAL FREQUNECY TREATMENT (S // G,Gc)
  if Qfc > 0:
    for i in range(len(clusters_df)):
      try:
        if pd.isna(clusters_df["log"]["vibrational_frequencies"].values[i]):
          continue
      except:
        lf = 0
      try:
        lf = float(clusters_df["log"]["vibrational_frequencies"].values[i][0])
      except:
        lf = 0
      if lf <= 0:
        clusters_df["log","entropy"][i] = missing
        continue
      Qt = clusters_df["log"]["temperature"].values[i]
      vibs = clusters_df["log"]["vibrational_frequencies"].values[i]
      structure = clusters_df["xyz"]["structure"].values[i]
      
      mu = [float(h/(8*np.pi**2*2.99793*10**10*vib)) for vib in vibs]
      try:
        mi = np.mean(structure.get_moments_of_inertia())
        Sr = [R*(0.5+np.log((8*np.pi**2.99793*(mu[j]*mi/(mu[j]+mi))*k*Qt/h**2)**0.5)) for j in range(len(mu))]  #cal/mol/K
        Sv = [R*h*vib*2.99793*10**10/k/Qt/(np.exp(h*vib*2.99793*10**10/k/Qt)-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/Qt)) for vib in vibs] #cal/mol/K
        w = [1/(1+(Qfc/vib)**4) for vib in vibs]
        Sv_corr = np.sum([w[j]*Sv[j]+(1-w[j])*Sr[j] for j in range(len(w))])
        Sv_each = np.sum(Sv)  #cal/mol/K
        clusters_df["log","entropy"][i] = clusters_df["log","entropy"][i]+(Sv_corr-Sv_each)
      except:
        mi = missing
        clusters_df["log","entropy"][i] = missing
    #print(clusters_df.index)
    #print(clusters_df["log","entropy"])  
    #exit()
      #print(mu)
      #print(mi)
      #print(Sr,flush=True)
      #print(vibs)
      #print(Sv,flush=True)
      #print(w,flush=True)
      #St_each = [R*np.log((2*np.pi*0.001*np.sum(structure.get_masses())*k**2/8.31441*Qt/h**2)**(3/2)*k*Qt/101325)+5/2 for structure in clusters_df["xyz"]["structure"].values]
      ###
      #print(Sv_each)
    ###

#  if Qfc > 0:
#    Qt = clusters_df["log"]["temperature"].values
#    vibs = clusters_df["log"]["vibrational_frequencies"].values
#    structures = clusters_df["xyz"]["structure"].values
#    mu = [[h/(8*np.pi**2*2.99793*10**10*vib) for vib in vibs[i]] for i in range(len(vibs))]
#    mi = [np.mean(structure.get_moments_of_inertia()) for structure in structures]
#    Sr = [[R*(0.5+np.log((8*np.pi**2.99793*(mu[i][j]*mi[i]/(mu[i][j]+mi[i]))*k*Qt[i]/h**2)**0.5)) for j in range(len(mu[i]))] for i in range(len(mu))] #cal/mol/K
#    Sv = [[R*h*vib*2.99793*10**10/k/Qt[i]/(np.exp(h*vib*2.99793*10**10/k/Qt[i])-1)-R*np.log(1-np.exp(-h*vib*2.99793*10**10/k/Qt[i])) for vib in vibs[i]] for i in range(len(vibs))] #cal/mol/K
#    w = [[1/(1+(Qfc/vib)**4) for vib in vibs[i]] for i in range(len(vibs))]
#    Sv_corr = [np.sum([w[i][j]*Sv[i][j]+(1-w[i][j])*Sr[i][j] for j in range(len(w[i]))]) for i in range(len(w))]
#    Sv_each = [np.sum(i) for i in Sv] #cal/mol/K
#    #St_each = [R*np.log((2*np.pi*0.001*np.sum(structure.get_masses())*k**2/8.31441*Qt/h**2)**(3/2)*k*Qt/101325)+5/2 for structure in clusters_df["xyz"]["structure"].values]
#    ###
#    clusters_df["log","entropy"] = [clusters_df["log","entropy"][i]+(Sv_corr[i]-Sv_each[i]) for i in range(len(clusters_df))]
#    ###

  ## CORRECTIONS FOR GIBBS FREE ENERGY
  for i in range(len(clusters_df)):
    try:
      clusters_df["log","gibbs_free_energy"][i] = clusters_df["log","enthalpy_energy"][i] - clusters_df["log","entropy"][i]/1000/627.503 * clusters_df["log","temperature"][i]
    except:
      clusters_df["log","gibbs_free_energy"][i] = missing
    try:
      clusters_df["log","gibbs_free_energy_thermal_correction"][i] = clusters_df["log","gibbs_free_energy"][i] - clusters_df["log","electronic_energy"][i]
    except:
      clusters_df["log","gibbs_free_energy_thermal_correction"][i] = missing

  #  clusters_df["log","gibbs_free_energy"] = [clusters_df["log","enthalpy_energy"][i] - clusters_df["log","entropy"][i]/1000/627.503 * clusters_df["log","temperature"][i] for i in range(len(clusters_df))] 
   #  clusters_df["log","gibbs_free_energy_thermal_correction"] = [clusters_df["log","gibbs_free_energy"][i] - clusters_df["log","electronic_energy"][i] for i in range(len(clusters_df))] 


###### FILTERING ######
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
if str(Qsort) != "0":
  if Qsort == "g":
    Qsort = "gibbs_free_energy"
  if Qsort == "el":
    Qsort = "electronic_energy"
  if Qsort == "gout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = (clusters_df["log"]["gibbs_free_energy"]-clusters_df["log"]["electronic_energy"]+clusters_df["out"]["electronic_energy"]).sort_values().index
      clusters_df = clusters_df.loc[sorted_indices]
  if Qsort == "elout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = clusters_df["out"]["electronic_energy"].sort_values().index
      clusters_df = clusters_df.loc[sorted_indices]
  if str(Qsort) != "no" and str(Qsort) != "gout" and str(Qsort) != "elout":
    clusters_df = clusters_df.sort_values([("log",Qsort)])
if str(Quniq) != "0":
  if Quniq == "dup":
    newclusters_df = clusters_df.copy()
    newclusters_df = newclusters_df.drop_duplicates(subset=[("info","file_basename")])
  else:
    if Qclustername != 0:
      uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
    else:
      uniqueclusters = "1"
    newclusters_df = []
    myNaN = lambda x : missing if x == "NaN" else x
    for i in uniqueclusters:
       if Qclustername != 0:
         preselected_df = clusters_df[clusters_df["info"]["cluster_type"] == i]
       else:
         preselected_df = clusters_df
       separated_inputs = seperate_string_number2(str(Quniq))
       compare_list = []
       compare_list_num = []
       for separated_input in separated_inputs:
         if isinstance(separated_input,list):
           if separated_input[0] == "rg":
             compare_list.append("rg")
           elif separated_input[0] == "el":
             compare_list.append("electronic_energy")         
           elif separated_input[0] == "g":
             compare_list.append("gibbs_free_energy")
           elif separated_input[0] == "d" or separated_input[0] == "dip":
             compare_list.append("dipole_moment")
           else:
             compare_list.append(separated_input[0])
           compare_list_num.append(float(separated_input[1]))
         else:
           #compare_list.append(separated_input)
           if separated_input == "rg":
             compare_list.append("rg")
             compare_list_num.append(2)
           elif separated_input == "el":
             compare_list.append("electronic_energy")
             compare_list_num.append(3)
           elif separated_input == "g":
             compare_list.append("gibbs_free_energy")
             compare_list_num.append(3)
           elif separated_input == "d" or separated_input == "dip":
             compare_list.append("dipole_moment")
             compare_list_num.append(1)
           else:
             compare_list.append(separated_input)
             compare_list_num.append(1)
       #if str(Quniq) == "rg,g":
       #  compare_list = ["rg","gibbs_free_energy"]
       #elif str(Quniq) == "rg":
       #  compare_list = ["rg"]
       #else:
       #  compare_list = ["rg","electronic_energy"]
       #print(compare_list)
       #print(compare_list_num)
       scale_test=1
       scale=1.0
       scale_scale = 0.1
       scale_iter = 0
       while scale_test == 1:
         scale_iter = scale_iter + 1 
         tocompare = []
         for compare_element_i in range(len(compare_list)):
           j = compare_list[compare_element_i]
           jj = scale*compare_list_num[compare_element_i]
           if j == "rg":
             rg = []
             for aseCL in preselected_df["xyz"]["structure"]:
               try:
                 rg.append((np.sum(np.sum((aseCL.positions-np.tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/np.sum(aseCL.get_masses()))**0.5)
               except:
                 rg.append(missing)
             values = [np.floor(myNaN(o)*10**jj) for o in rg]
           elif j == "gout":
             gout = []
             for Gouti in range(len(preselected_df)):
               try:
                 gout.append(preselected_df["log"]["gibbs_free_energy"].values[Gouti]-preselected_df["log"]["electronic_energy"].values[Gouti]+preselected_df["out"]["electronic_energy"].values[Gouti])
               except:
                 gout.append(missing)
             values = [np.floor(myNaN(o)*10**jj) for o in gout]
           else:  
             values = [np.floor(float(myNaN(o))*10**jj) for o in preselected_df["log"][j].values]
           tocompare.append(values)
         tocompare = np.transpose(tocompare)
         uniqueindexes = np.unique(tocompare,axis = 0,return_index=True)[1]
         if Qsample > 0 and len(uniqueindexes) > Qsample and scale_iter < 100:
           if scale_scale > 0:
             scale_scale = -0.9*scale_scale
           scale = scale + scale_scale
         elif Qsample > 0 and len(uniqueindexes) < Qsample and len(compare_list) >= Qsample and scale_iter < 100:
           if scale_scale < 0:
             scale_scale = -0.9*scale_scale
           scale = scale + scale_scale  
         else:
           scale_test = 0
       selected_df = preselected_df.iloc[uniqueindexes]
       if len(newclusters_df) == 0:
         newclusters_df = selected_df
       else:
         #print(newclusters_df)
         newclusters_df = newclusters_df.append(selected_df)
         #print(newclusters_df)
  if Qout == 1:
    if Qsample > 0:
      print("Sampled: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
    else:
      print("Uniqueness: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
  clusters_df = newclusters_df
if Qarbalign > 0:
  import ArbAlign
  from joblib import Parallel, delayed
  import multiprocessing
  num_cores = multiprocessing.cpu_count()
  def compare_pair(arg):
    tobecompared = comparepairs[arg]
    AAci = preselected_df.loc[allindexes[tobecompared[0]]]
    AAcj = preselected_df.loc[allindexes[tobecompared[1]]]
    return ArbAlign.compare(AAci["xyz"]["structure"],AAcj["xyz"]["structure"])
  if Qclustername != 0:
    uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
  else:
    uniqueclusters = "1"
  newclusters_df = []
  myNaN = lambda x : missing if x == "NaN" else x
  for i in uniqueclusters:
     if Qclustername != 0:
       preselected_df = clusters_df[clusters_df["info"]["cluster_type"] == i]
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
           #TODO try:
           #print(comparepairs[AAc][1]) 
           removedindexes.append(allindexes[comparepairs[AAc][1]])
     clusters_df = clusters_df.drop(removedindexes)
#not sure if this sorting is necessary but maybe after uniqueness filtering yes
if str(Qsort) != "0" and str(Qsort) != "no":
  if Qsort == "gout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = (clusters_df["log"]["gibbs_free_energy"]-clusters_df["log"]["electronic_energy"]+clusters_df["out"]["electronic_energy"]).sort_values().index
      clusters_df = clusters_df.loc[sorted_indices]
  elif Qsort == "elout":
    if ("out","electronic_energy") in clusters_df.columns:
      sorted_indices = clusters_df["out"]["electronic_energy"].sort_values().index
      clusters_df = clusters_df.loc[sorted_indices]
  else:
    clusters_df = clusters_df.sort_values([("log",Qsort)])
if Qthreshold != 0:
  for i in range(len(Qcut)):
    if Qcut[i][2] == "el":
      what = 627.503*clusters_df["log"]["electronic_energy"].values
    elif Qcut[i][2] == "g":
      what = 627.503*clusters_df["log"]["gibbs_free_energy"].values
    elif Qcut[i][2] == "elout":
      what = 627.503*clusters_df["out"]["electronic_energy"].values
    elif Qcut[i][2] == "lf":
      #print([len(i) for i in clusters_df["log"]["vibrational_frequencies"].values])
      #what = np.array([np.array([ii]) if pd.isna([ii]).any() else np.array(ii) for ii in clusters_df["log"]["vibrational_frequencies"].values])[:,0]
      what = np.array([np.array(ii) if pd.isna([ii]).any() else np.array(ii[0]) for ii in clusters_df["log"]["vibrational_frequencies"].values])
    elif Qcut[i][2] == "rg":
      rg = []
      for aseCL in clusters_df["xyz"]["structure"]:
        try:
          rg.append((np.sum(np.sum((aseCL.positions-np.tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/np.sum(aseCL.get_masses()))**0.5)
        except:
          rg.append(missing)
      what = np.array(rg)
    else:
      what = clusters_df["log"][Qcut[i][2]].values
    if Qcut[i][1] == "relative":
      min = np.min(what)
    else:
      min = 0
    if Qcut[i][3] == "nan" or Qcut[i][3] == "NA" or Qcut[i][3] == "na" or Qcut[i][3] == "NaN":   
      if Qcut[i][0] == ">":
        clusters_df = clusters_df[pd.isna(what-min)]       
      else:
        clusters_df = clusters_df.drop(index=clusters_df[pd.isna(what-min)].index)
    else:
      with np.errstate(invalid='ignore'):
        if Qcut[i][0] == ">":
          clusters_df = clusters_df[what-min > float(Qcut[i][3])]       
        else:
          clusters_df = clusters_df[what-min <= float(Qcut[i][3])]
if str(Qselect) != "0":
  if Qclustername != 0:
    uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
  else:
    uniqueclusters = "1"
  newclusters_df = []
  for i in uniqueclusters:
     if Qclustername != 0:
       selected_df = clusters_df[clusters_df["info"]["cluster_type"] == i][0:Qselect] 
     else:
       selected_df = clusters_df[0:Qselect]
     if len(newclusters_df) == 0:
       newclusters_df = selected_df
     else:
       #print(newclusters_df)
       newclusters_df = newclusters_df.append(selected_df)
       #print(newclusters_df)
  if Qout == 1:
    print("Selecting/Sampling: "+str(len(clusters_df))+" --> "+str(len(newclusters_df)))
  clusters_df = newclusters_df
### SHUFFLE
if Qshuffle == 1:
  clusters_df = clusters_df.sample(frac=1)

#x = clusters_df["info"]["file_basename"].astype("category")
#print(x.astype("category").categories)
if Qrebasename == 1:
  values=clusters_df["info"]["file_basename"].values
  for i in range(len(clusters_df)):
    if values[i] in np.delete(values, i, axis=0):
      version=1
      values=clusters_df["info"]["file_basename"].values
      while values[i]+"-v"+str(version) in values:
        version+=1
      clusters_df["info","file_basename"][i] = values[i]+"-v"+str(version)
      original_clusters_df.loc[clusters_df.index[i]]["info","file_basename"] = values[i]+"-v"+str(version)

## SAVE OUTPUT.pkl ##
if Qout > 0:
  #print(original_clusters_df.iloc[0])
  tosave = original_clusters_df.loc[clusters_df.index]
  #print(clusters_df.index[0])
  #print(original_clusters_df.loc['19268'])  
  #print(len(clusters_df.index))
  #print(len(clusters_df.index))
  #print(len(tosave))
  if Qsplit == 1:
    try:
      tosave.to_pickle(output_pkl)
    except:
      print("Pickle was not written down due to an error.")
    if Qout == 1:
      if len(tosave) == 0:
        print(tosave)
        print("No files in the input!")
      else:
        print("Example output:")
        print(tosave.iloc[0])
        print("Number of files in "+output_pkl+": "+str(len(tosave)))
  else:
    if len(tosave) == 0:
      if Qout ==1:
        print(tosave)
        print("No files in the input!")
    else:
      lengths = -(-len(tosave)//Qsplit)
      for split in range(Qsplit):
        output_pkl_split = output_pkl[:-4]+"_s"+str(split+1)+".pkl"
        start=split*lengths
        end=(split+1)*lengths
        if end > len(tosave):
          end = len(tosave)
        tosave[start:end].to_pickle(output_pkl_split)
        if Qout == 1:
          print("Number of files in "+output_pkl_split+": "+str(len(tosave[start:end])))
    

## EXTRACT DATA ##
if Qsolvation != "0" or Qformation != 0:
  if len(Pout) > 0:
    if Pout[0] != "-b" and Pout[0] != "-ct":
      Pout.insert(0,"-ct")
output = []
last = ''
for i in Pout:
  if i == "-info":
    print(clusters_df.info())
    continue
  if i == "-levels":
    pd.set_option('display.max_colwidth', None)
    if not pd.isna(clusters_df["log"]["program"]).all():
      print("# LOG #")
      print(clusters_df["log"][["program","method"]].drop_duplicates())
      print("#######")   
    if "out" in clusters_df:
      if not pd.isna(clusters_df["out"]["program"]).all():
        print("# OUT #")
        print(clusters_df["out"][["program","method"]].drop_duplicates())
        print("#######")   
    continue
  if i == "-extra":
    last = "-extra"
    continue
  if last == "-extra":
    last = ""
    output.append(clusters_df["extra"][i].values)
    continue
  if i == "-cite":
    try:
      output.append(clusters_df["info"]["citation"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  #XYZ
  if i == "-xyz":
    for ind in clusters_df.index:
      try:
        write(clusters_df["info"]["file_basename"][ind]+".xyz",clusters_df["xyz"]["structure"][ind])
      except:
        print("Corrupted structure saved for "+clusters_df["info"]["file_basename"][ind])
    continue
  #PDB IMOS
  if i == "-imos":
    from ase.io.proteindatabank import write_proteindatabank
    for ind in clusters_df.index:
      write_proteindatabank(".JKQChelp.pdb",clusters_df["xyz"]["structure"][ind])
      f1=open(".JKQChelp.pdb", "r")
      f2=open(clusters_df["info"]["file_basename"][ind]+".pdb", "w")
      f2.write(f1.readline())
      for j in range(len(clusters_df["xyz"]["structure"][ind].get_chemical_symbols())):
        line=f1.readline()
        ch4="{:+.2f}".format(clusters_df["log"]["esp_charges"][ind][j])
        #print(line)
        #print(line[:56]+ch4+line[61:])
        #print("ATOM      2  H1  GLY     1      52.826  13.287  -6.429  0.1642   1.300       H")
        #line[52:62]="0.012164968"
        f2.write(line[:56]+ch4+line[61:])
      f2.write(f1.readline())
      f1.close()
      f2.close()
      remove(".JKQChelp.pdb")
    continue
  #XLSX IMOS
  if i == "-imos_xlsx":
    from ase.data.vdw import vdw_radii
    import xlsxwriter
    workbook = xlsxwriter.Workbook('imos.xlsx')
    bold = workbook.add_format({'bold': True,'fg_color': '#FFFF00', 'border':1})
    for ind in clusters_df.index:
      clustername = clusters_df["info"]["file_basename"][ind]
      if len(clustername) > 30:
        clustername = clustername[:20]+'-TBC'+str(ind)
      worksheet = workbook.add_worksheet(clustername)
      
      pos=clusters_df["xyz"]["structure"][ind].get_positions()
      mass=clusters_df["xyz"]["structure"][ind].get_masses()
      at=clusters_df["xyz"]["structure"][ind].get_atomic_numbers()
      ch=clusters_df["log"]["esp_charges"][ind]
      row = 0
      for j in range(len(pos)):
        worksheet.write(row, 0, pos[j][0])
        worksheet.write(row, 1, pos[j][1])
        worksheet.write(row, 2, pos[j][2])
        worksheet.write(row, 3, vdw_radii[at[j]])
        worksheet.write(row, 4, ch[j])
        if at[j] == 18:
          worksheet.write(row, 6, 41) #argon
        elif at[j] == 84:
          worksheet.write(row, 6, 212) #polonium
        else:
          worksheet.write(row, 6, np.round(mass[j]))
        row += 1
      
      worksheet.write(0, 5, 'TOTAL z',bold)
      worksheet.write(1, 5, clusters_df["log"]["charge"][ind])
      worksheet.write(2, 5, 'Totalmass',bold)
      worksheet.write(3, 5, np.sum(clusters_df["xyz"]["structure"][ind].get_masses()))
      worksheet.write(4, 5, 'Crossection',bold)
      
    workbook.close()
    continue
  #CHARGES
  if i == "-chargesESP":
    for ind in clusters_df.index:
      f = open(clusters_df["info"]["file_basename"][ind]+".charges","w")
      f.write("\n".join([str(tt) for tt in clusters_df["log"]["esp_charges"][ind]])+"\n")
      f.close()
    continue
  if i == "-charges":
    for ind in clusters_df.index:
      f = open(clusters_df["info"]["file_basename"][ind]+".charges","w")
      f.write("\n".join([str(tt) for tt in clusters_df["log"]["mulliken_charges"][ind]])+"\n")
      f.close()
    continue
  #MOVIE
  if i == "-movie":
    f = open("movie.xyz","w")
    for ind in clusters_df.index:
      try:
        write(".movie.xyz",clusters_df["xyz"]["structure"][ind])
        with open(".movie.xyz", 'r') as f2:
          lines = f2.readlines()
          lines[1] = clusters_df["info"]["file_basename"][ind]+"\n"
        f.writelines(lines)
        f2.close()
      except:
        continue
    f.close()
    remove(".movie.xyz")
    continue
  #Atoms
  if i == "-atoms":
    atoms = []
    for ind in clusters_df.index:
      try:
        aseCL=clusters_df["xyz"]["structure"][ind]
        #print(aseCL.get_atomic_numbers())
        atoms.extend(aseCL.get_atomic_numbers())
      except:
        continue
    print(" ".join([str(i) for i in np.unique(atoms)]))
    continue
  #bonded
  if i == "-bonded":
    bonded = []
    Qbonded_index += 1
    for ind in clusters_df.index:
      try:
        aseCL=clusters_df["xyz"]["structure"][ind] 
        positions = aseCL.positions
        symb = np.array(aseCL.get_chemical_symbols())
        dist = lambda p1, p2: np.sqrt(np.sum(((p1-p2)**2)))
        symb_ind = np.array(aseCL.get_chemical_symbols())
        mask1 = symb_ind == Qbonded[Qbonded_index][1]
        mask2 = symb_ind == Qbonded[Qbonded_index][2]
        dm = [dist(p1, p2) for p2 in positions[mask1] for p1 in positions[mask2]]
        bonds = sum(test_i <= float(Qbonded[Qbonded_index][0]) for test_i in dm)
        if Qbonded[Qbonded_index][1] == Qbonded[Qbonded_index][2]:
          bonds = (bonds-sum(mask1))/2
        bonded.append(str(int(bonds)))
      except:
        bonded.append(missing)
    output.append(bonded)
    continue
  #Rg
  if i == "-rg":
    rg = []
    for ind in clusters_df.index:
      try:
        aseCL=clusters_df["xyz"]["structure"][ind]
        rg.append((np.sum(np.sum((aseCL.positions-np.tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/np.sum(aseCL.get_masses()))**0.5)
      except:
        rg.append(missing)
    output.append(rg)
    continue
  if i == "-radius" or i == "-radius0.5":
    radius = []
    for aseCL in clusters_df["xyz"]["structure"]:
      try:
        dist = lambda p1, p2: np.sqrt(np.sum((p1-p2)**2))
        centered = aseCL.positions-aseCL.positions.mean(axis = 0)
        ratios = np.sqrt(np.linalg.eigvalsh(np.dot(centered.transpose(),centered)/len(centered)))
        maxdist = np.max(np.asarray([[dist(p1, p2) for p2 in aseCL.positions] for p1 in aseCL.positions]))
        if i == "-radius0.5":
          if ratios[0] < 0.5:
            ratios[0] = 0.5
          if ratios[1] < 0.5:
            ratios[1] = 0.5
          if ratios[2] < 0.5:
            ratios[2] = 0.5
          if maxdist < 0.5:
            maxdist = 0.5
        if np.max(ratios) > 0:
          ratios = ratios / np.max(ratios)
        else:
          ratios = missing
        radius.append((np.prod(maxdist*ratios))**(1/3))	
      except:
        radius.append(missing)
    output.append(radius)
    continue
  #INFO
  if i == "-ct":
    try: 
      output.append(clusters_df["info"]["cluster_type"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-b": 
    try:
      output.append(clusters_df["info"]["file_basename"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-nOUT":
    try:
      output.append(clusters_df["info"]["file_basename"].values+".out")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-nLOG":
    try:
      output.append(clusters_df["info"]["file_basename"].values+".log")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-nXYZ":
    try:
      output.append(clusters_df["info"]["file_basename"].values+".xyz")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pOUT":
    try:
      output.append(clusters_df["info"]["folder_path"].values+clusters_df["info"]["file_basename"].values+".out")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pLOG":
    try:
      output.append(clusters_df["info"]["folder_path"].values+clusters_df["info"]["file_basename"].values+".log")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pXYZ":
    try:
      output.append(clusters_df["info"]["folder_path"].values+clusters_df["info"]["file_basename"].values+".xyz")
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-ePKL":
    if Qout > 0 or len(input_pkl) == 1:
      if Qout > 0:
        try:
          output.append(path.abspath(output_pkl)+"/:EXTRACT:/"+clusters_df["info"]["file_basename"].values)
        except:
          output.append([missing]*len(clusters_df))
      else:
        try:
          output.append(path.abspath(input_pkl[0])+"/:EXTRACT:/"+clusters_df["info"]["file_basename"].values)
        except:
          output.append([missing]*len(clusters_df))
    else:
      print("Sorry but it seems to me that you are taking this from more file or no file is formed (yet)")
      exit()
    continue
  if i == "-mass":
    masses = []  
    for ind in clusters_df.index:
      try:
        masses.append(str(np.sum(clusters_df["xyz"]["structure"][ind].get_masses())))
      except:
        masses.append(missing)
    output.append(masses)
    continue
  if i == "-natoms":
    natoms = []
    for ind in clusters_df.index:
      try:
        natoms.append(str(len(clusters_df["xyz"]["structure"][ind].get_atomic_numbers())))
      except:
        natoms.append(missing)
    output.append(natoms)
    continue
  if i == "-elsp":
    try:
      output.append(QUenergy*clusters_df["log"]["sp_electronic_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-el":
    try:
      output.append(QUenergy*clusters_df["log"]["electronic_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-elout":
    try:
      output.append(QUenergy*clusters_df["out"]["electronic_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-elc":
    try:
      output.append(QUenergy*(clusters_df["out"]["electronic_energy"].values-clusters_df["log"]["electronic_energy"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-uc":
    try:
      output.append(QUenergy*clusters_df["log"]["energy_thermal_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-u":
    try:
      output.append(QUenergy*(clusters_df["log"]["electronic_energy"].values+clusters_df["log"]["energy_thermal_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-uout":
    try:
      output.append(QUenergy*(clusters_df["out"]["electronic_energy"].values+clusters_df["log"]["energy_thermal_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-zpec": 
    try:
      output.append(QUenergy*clusters_df["log"]["zero_point_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-zpe":
    try:
      output.append(QUenergy*(clusters_df["log"]["electronic_energy"].values+clusters_df["log"]["zero_point_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-zpeout":
    try:
      output.append(QUenergy*(clusters_df["out"]["electronic_energy"].values+clusters_df["log"]["zero_point_correction"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-g":
    try:
      output.append(QUenergy*clusters_df["log"]["gibbs_free_energy"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pop":
    population = []
    R = 1.987 #cal/mol/K #=8.31441
    for pop in range(len(clusters_df)):
      try:
        if Qclustername != 0:
          allsameclusters = clusters_df[clusters_df["info"]["cluster_type"]==clusters_df["info"]["cluster_type"].values[pop]]
        else: 
          allsameclusters = clusters_df
        try:
          Qt = float(clusters_df["log","temperature"].values[pop])
        except:
          Qt = 298.15
        me=clusters_df["log"]["gibbs_free_energy"].values[pop] #*627.503*1000/Qt/R
        partition_function = np.sum([np.exp(iii) for iii in -(allsameclusters["log"]["gibbs_free_energy"].values-me)*627.503*1000/Qt/R])
        ratio=np.exp(-0*627.503*1000/Qt/R)/partition_function
        population.append(f"{ratio:8f}")
      except:
        population.append(missing)
    output.append(population)
    continue
  if i == "-gc":
    try:
      output.append(QUenergy*clusters_df["log"]["gibbs_free_energy_thermal_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-gout":
    try:
      output.append(QUenergy*(clusters_df["log"]["gibbs_free_energy"].values+clusters_df["out"]["electronic_energy"].values-clusters_df["log"]["electronic_energy"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-h":
    try:
      entalpies = clusters_df["log"]["enthalpy_energy"].values
      entalpies[entalpies == "NaN"] = missing 
      output.append(QUenergy*entalpies)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-hc":
    try:
      output.append(QUenergy*clusters_df["log"]["enthalpy_thermal_correction"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-hout":
    try:
      output.append(QUenergy*(clusters_df["log"]["enthalpy_energy"].values+clusters_df["out"]["electronic_energy"].values-clusters_df["log"]["electronic_energy"].values))
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-s":
    try:
      entropies = clusters_df["log"]["entropy"].values
      #entropies[entropies == "NaN"] = missing
      entropies = np.array([missing if x == "NaN" else x for x in entropies])
      output.append(QUentropy*entropies/1000/627.503)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-lf": 
    lowestfreq = []
    for aseCL in clusters_df["log"]["vibrational_frequencies"]:
      try:
        lowestfreq.append(aseCL[0])
      except:
        lowestfreq.append(missing)
    output.append(lowestfreq)
    continue
  if i == "-level":
    levels = []
    for ind in clusters_df.index:
      try:
        if ("log","method") in clusters_df and ("log","program") in clusters_df:
          if ("out","method") in clusters_df and ("out","program") in clusters_df:
            if not pd.isna(clusters_df["out"]["program"][ind]):
              levels.append(clusters_df["log"]["program"][ind]+"_"+clusters_df["log"]["method"][ind]+"__"+clusters_df["out"]["program"][ind]+"_"+clusters_df["out"]["method"][ind])
            else:
              levels.append(clusters_df["log"]["program"][ind]+"_"+clusters_df["log"]["method"][ind])
          else:
            levels.append(clusters_df["log"]["program"][ind]+"_"+clusters_df["log"]["method"][ind])  
        else:
          if ("out","method") in clusters_df and ("out","program") in clusters_df:
            levels.append(clusters_df["out"]["program"][ind]+"_"+clusters_df["out"]["method"][ind])
          else:
            levels.append(missing)
      except:
        levels.append(missing)
    output.append(levels)
    continue
  if i == "-f":
    try:
      output.append(clusters_df["log"]["vibrational_frequencies"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-rot":
    try:
      output.append(clusters_df["log"]["rotational_constant"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-rots":
    try:
      output.append(clusters_df["log"]["rotational_constants"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-mult":
    try:
      output.append(clusters_df["log"]["multiplicity"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-char":
    try:
      output.append(clusters_df["log"]["charge"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-esp":
    try:
      output.append(clusters_df["log"]["esp_charges"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-mull":
    try:
      output.append(clusters_df["log"]["mulliken_charges"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-dip":
    try:
      output.append(clusters_df["log"]["dipole_moment"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-dips":
    try:
      output.append(clusters_df["log"]["dipole_moments"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-pol":
    try:
      output.append(clusters_df["log"]["polarizability"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-templog":
    try:
      output.append(clusters_df["log"]["temperature"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-preslog":
    try:
      output.append(clusters_df["log"]["pressure"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-mi":
    try:
      output.append([structure.get_moments_of_inertia() for structure in clusters_df["xyz"]["structure"].values])
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-ami":
    try:
      output.append([np.mean(structure.get_moments_of_inertia()) for structure in clusters_df["xyz"]["structure"].values])
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-rsn":
    try:
      output.append(clusters_df["log"]["rotational_symmetry_number"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-t":
    try:
      output.append(clusters_df["log"]["time"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-termination":
    try: 
      output.append(clusters_df["log"]["termination"].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  if i == "-column":
    try:
      output.append(clusters_df[Qcolumn[0][0]][Qcolumn[0][1]].values)
    except:
      output.append([missing]*len(clusters_df))
    continue
  output.append(["UNKNOWN_ARGUMENT"]*len(clusters_df))

## PRINT DATA ##
if not len(output) == 0:
  #print(output)
  output = np.array(output)

  #TAKING GLOBAL MINIMA ONLY
  if Qglob == 1 or Qglob == 2:
    uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
    indexes = []
    for i in uniqueclusters:
      if Qglob == 1:
        GFE = clusters_df["log","gibbs_free_energy"].values[clusters_df["info"]["cluster_type"] == i]
      elif Qglob == 2:
        GFE = clusters_df["log","gibbs_free_energy"].values + clusters_df["out"]["electronic_energy"].values - clusters_df["log"]["electronic_energy"].values 
        GFE = np.array(GFE)[clusters_df["info"]["cluster_type"] == i]
      else:
        print("Qglob error. [EXITING]")
        exit()
      GFE[GFE == "NaN"] = missing
      if len(GFE[~pd.isna(GFE)]) != 0:
        globindex = clusters_df.index.values[clusters_df["info"]["cluster_type"] == i][GFE == np.min(GFE[~pd.isna(GFE)])][0]
        indexes.append(int(np.array(range(len(clusters_df)))[clusters_df.index.values == globindex][0]))
    
    newoutput = []
    for j in range(output.shape[0]):
      toappend=[]
      for i in indexes:
        toappend.append(output[j][i])
      newoutput.append(np.array(toappend,dtype=object)) 
    #output = [[output[j][i] for i in indexes] for j in range(output.shape[0])]
    output = np.array(newoutput)
    #output = np.array(output)
  
  #TAKING BOLTZMANN AVERAGE OVER ALL MINIMA 
  if Qbavg == 1 or Qbavg == 2:
    k = 1.380662*10**-23 # [J/K]
    uniqueclusters = np.unique(clusters_df["info"]["cluster_type"].values)
    portions = []
    freeenergies = [] 
    entropies = []
    indexes = []
    for i in uniqueclusters:
      if Qbavg == 1:
        GFE = clusters_df["log","gibbs_free_energy"].values[clusters_df["info"]["cluster_type"] == i]
        temps = clusters_df["log","temperature"].values[clusters_df["info"]["cluster_type"] == i]
      elif Qbavg == 2:
        GFE = clusters_df["log","gibbs_free_energy"].values + clusters_df["out"]["electronic_energy"].values - clusters_df["log"]["electronic_energy"].values
        GFE = np.array(GFE)[clusters_df["info"]["cluster_type"] == i]
        temps = clusters_df["log","temperature"].values[clusters_df["info"]["cluster_type"] == i]
      else:
        print("Qglob error. [EXITING]")
        exit()
      nonans = ~pd.isna(GFE)
      GFE = GFE[nonans]
      try:
        minimum = np.min(GFE)
      except:
        minimum = missing
      GFE = GFE-minimum
      #print(GFE)     
 
      preportions = [np.exp(-GFE[i]*43.60*10**-19/k/temps[nonans][i]) for i in range(GFE.shape[0])]
      #print(preportions)
      toaddfreeenergies = []
      try:
        addfreeenergy = QUenergy*(minimum - 1/43.60/10**-19*k*temps[nonans][0]*np.log(np.sum([np.exp(-GFE[i]*43.60*10**-19/k/temps[nonans][i]) for i in range(GFE.shape[0])])))
      except:
        addfreeenergy = missing
      freeenergies.append(addfreeenergy)
      #freeenergies.append(QUenergy*(minimum - 1/43.60/10**-19*k*temps[nonans][0]*np.log(np.sum([np.exp(-GFE[i]*43.60*10**-19/k/temps[nonans][i]) for i in range(GFE.shape[0])]))))
      #print(freeenergies)
      sumpreportions = np.sum(preportions)
      portions.append(preportions/sumpreportions)
      #print(portions)
      indexes.append(np.array(range(len(clusters_df)))[clusters_df["info"]["cluster_type"] == i][nonans])
    
    def is_averagable(input_array): # :-D
      try:
        [float(i) for i in input_array]
        test = 0
        for i in range(len(input_array)):
          if input_array[i] != input_array[0]:
            test = 1
            break
        if test == 0:
          return False
        else:
          return True
      except ValueError:
        return False
    def is_the_same(input_array):
      test = 0
      for i in range(len(input_array)):
        if input_array[i] != input_array[0]:
          test = 1
          break
      if test == 0:
        return input_array[0]
      else:
        return missing

    def myif(cond,opt1,opt2,opt3):
      if Pout[cond] == "-g" or Pout[cond] == "-gout":
        return opt2
      elif Pout[cond] == "-s":
        return opt3
      else:
        return opt1
   
    #print(output)
    newoutput = []
    skip = []
    #print(output)
    for l in range(output.shape[0]):
      toappend = []
      for i in range(len(portions)):
        if i in skip:
          continue
        try:
          appendable = myif(l,np.sum([float(portions[i][j])*float(output[l,indexes[i][j]]) for j in range(len(portions[i]))]),freeenergies[i],missing) if is_averagable(output[l][indexes[i]]) else is_the_same(output[l,indexes[i]])
        except: 
          appendable = missing
        if l == 0:
          if pd.isna(appendable):
            skip.append(i)
            continue
        #print(appendable)
        toappend.append(appendable)
      #print(toappend)
      newoutput.append(np.array(toappend, dtype=object))
    output = np.array(newoutput)
    #output = [[myif(l,np.sum([float(portions[i][j])*float(output[l,indexes[i][j]]) for j in range(len(portions[i]))]),freeenergies[i],missing) if is_averagable(output[l][indexes[i]]) else is_the_same(output[l,indexes[i]]) for i in range(len(portions))] for l in range(output.shape[0])]   
    #print(output)

  #fn = ".help"+str(np.random.randint(100000,size=1)[0])+".txt" 
  #f = open(fn, "w")
  #[f.write(" ".join(map(str,row))+"\n") for row in list(zip(*output))]
  #f.close()
  ##TODO can you make this working using only python?
  #if Qout != 2 or Qformation == 0:
  #  system("cat "+fn+" | column -t")
  #  remove(fn)
  if not(Qout == 2 and Qformation==1):
    toprint = list(zip(*output)) #[row for row in list(zip(*output))]
    if len(toprint) > 0:
      column_widths = [max(len(str(row[i])) for row in toprint) for i in range(len(toprint[0]))]
      for row in toprint:
        formatted_row = [str(row[i]).ljust(column_widths[i]) for i in range(len(row))]
        print(" ".join(formatted_row),flush = True)

#print(output)
#print(type(output[3][3]))
def myFunc(e):
  try:
    numero = np.sum([int(i) for i in seperate_string_number(dash_comment(e[0])[0])[0::2]])
  except:
    numero = 0
  return numero+len(e[0])/100.
dt = np.dtype(object)

#LOAD INPUT FILE
if Qsolvation != "0" or Qformation == 1:
  if len(formation_input_file) > 0:
    f = open(formation_input_file, "r")
    output = []
    def mytofloat(string):
      try:
        return QUenergy*float(string)
      except:
        return string
    for line in f.readlines():
      output.append([mytofloat(i) for i in line.split()])
    output = np.array(output,dtype=dt).transpose()
    f.close()

if Qsolvation != "0":
  if Qout != 2:
    print("#####################################",flush = True)
    print("##########  SOLVATION  ##############",flush = True)
    print("#####################################",flush = True)
  if len(Pout)>1:
    if Pout[1] != "-g" and Pout[1] != "-gout":
      print("Please use -ct -g or -ct -gout as the first two arguments.")
      exit()
  #SORT OUTPUT1
  #print(output) 
  output = np.transpose(np.transpose(output)[np.apply_along_axis(myFunc, axis=1, arr=np.transpose(output)).argsort()])
  #print(output) 
  cluster_types = [dash_comment(seperate_string_number(i))[0] for i in output[0]]
  no_solvent_cluster_types = []
  solvent_content = []
  for i in cluster_types:
    solvent_content_i = 0
    no_solvent_cluster_types_i = []
    for j in range(int(len(i)/2)):
      if i[2*j+1] == Qsolvation:
        solvent_content_i = i[2*j]
      else:
        no_solvent_cluster_types_i.append(i[2*j]) 
        no_solvent_cluster_types_i.append(i[2*j+1]) 
    no_solvent_cluster_types.append(no_solvent_cluster_types_i)
    solvent_content.append(solvent_content_i)
  cluster_types = ["".join(i) for i in cluster_types]
  no_solvent_cluster_types = ["".join(i) for i in no_solvent_cluster_types]
  #unique_clusters = [x for x in np.array(np.unique(no_solvent_cluster_types), dtype=object) if x != []]
  unique_clusters = [x for x in np.array(np.unique(no_solvent_cluster_types), dtype=object) if x != ""]
  #unique_clusters = np.unique(no_solvent_cluster_types)
  #unique_clusters = filter(None, unique_clusters)
  free_energies = output[1]
  
  #print(cluster_types)
  #print(no_solvent_cluster_types)
  #print(solvent_content)
  #print(unique_clusters)
  index_of_solvent = -1
  for i in range(len(cluster_types)):
    if cluster_types[i] == "1"+Qsolvation:
      index_of_solvent = i
  if index_of_solvent == -1:
    print("Missing solvent")
    exit()
  #print(index_of_solvent)
  #solvent_free_energy = 

  #TODO
  if not pd.isna(Qt):
    Temp = Qt
  else:
    Temp = 298.15
  print(f"Temp = %.2f K; "%(Temp), end = "")
  if not pd.isna(QPsolvent):
    p_solvent = QPsolvent
  else:
    if not pd.isna(Qrh):
      rh = Qrh
    else:
      rh = 1.0
    print(f"RH = %.2f %%; "%(rh*100), end = "") 
    #p_solvent = rh*10**(8.14019-1810.9/(244.485+Temp-273.15))*133.322
    p_solvent = rh*100*6.1094*np.exp(17.625*(Temp-273.15)/(Temp-273.15+234.04))
  print(f"p_solvent = %.2f Pa "%(p_solvent))
  #p_solvent = 100  
  p_ref = 101325
  R = 1.987 #cal/mol/K #=8.31441
  #Antoine equation
  #https://www.omnicalculator.com/chemistry/vapour-pressure-of-water
  #print("P_solvent")
  #print(p_solvent) 

  #print(output)
  new_output = []  
  for i in unique_clusters:
    indexes = []
    for j in range(len(no_solvent_cluster_types)):
      if i == no_solvent_cluster_types[j]:
        indexes.append(j)
    #print("Indexes:")
    #print(indexes)
    tot_conc = 0
    nominators = []
    free_energies_i = []
    for j in indexes:
      free_energies_i.append(output[1][j]-float(solvent_content[j])*output[1][index_of_solvent])
    #print("Free energies:")
    #print(free_energies_i)
    minimum = np.min(free_energies_i)
    free_energies_i = free_energies_i - minimum
    #print("Minimum:")
    #print(minimum)
    #print("Free energies:")
    #print(627.503/QUenergy*free_energies_i/R*1000/Temp)
    #print(indexes)
    #print(1/R*1000/Temp)
    for j in range(len(indexes)):
      #print(-QUenergy*free_energies_i[j]/R*1000/Temp)
      nominator = (p_solvent/p_ref)**float(solvent_content[indexes[j]])*np.exp(-627.503/QUenergy*free_energies_i[j]/R*1000/Temp)
      nominators.append(nominator)
      #print(cluster_types[indexes[j]])
    denominator = np.sum(nominators)
    #print(nominators/denominator)
    
    #print(output[0])
    for j in range(len(indexes)):
      print(f"%8s "%(cluster_types[indexes[j]]), end="")
    print("")
    for j in range(len(indexes)):
      print(f"%8.1f "%(nominators[j]/denominator*100),end="")
    print("")
    print("----------------------------------------")
    #new_output
    

if Qformation == 1:
  if Qout != 2:
    print("#####################################",flush = True)
    print("##########  FORMATION  ##############",flush = True)
    print("#####################################",flush = True)
  #SORT OUTPUT
  output = np.transpose(np.transpose(output)[np.apply_along_axis(myFunc, axis=1, arr=np.transpose(output)).argsort()])
  #
  #print(seperate_string_number(output[0][0]))
  #print(dash_comment(seperate_string_number(output[0][0])))
  cluster_types = [dash_comment(seperate_string_number(i))[0] for i in output[0]]
  #print(cluster_types)
  ######## SOLVING PROTONATION
  for i in range(len(cluster_types)):
    chosen_cluster_type = cluster_types[i]
    if "p" in chosen_cluster_type:
      protonated_base = ""
      for base in ["gd","eda","tma","dma","am","buta","dbma","dea","dhexa","dmea","dpenta","dpropa","ibuta","IIebuta","ipropa","ipropea","mea","nona","propa","sbuta","tbuta","tea","tibuta","tpropa","w"]:
        if base in chosen_cluster_type:
          protonated_base = base
          break
      new_cluster_type = []
      for j in range(1,len(chosen_cluster_type),2):
        if chosen_cluster_type[j] == protonated_base:
          subsctracted = int(chosen_cluster_type[j-1]) - 1
          if subsctracted > 0:
            new_cluster_type.append(str(subsctracted))
            new_cluster_type.append(chosen_cluster_type[j])
        elif chosen_cluster_type[j] == "p":
          if int(chosen_cluster_type[j-1]) > 1:
            new_cluster_type.append(str(int(chosen_cluster_type[j-1]) - 1))
            new_cluster_type.append(chosen_cluster_type[j])
          else:
            continue
        else:
          new_cluster_type.append(chosen_cluster_type[j-1])
          new_cluster_type.append(chosen_cluster_type[j])
      new_cluster_type.append("1")
      new_cluster_type.append("protonated_"+protonated_base)       
      cluster_types[i] = new_cluster_type 
  ########################        
  #print(cluster_types)
  cluster_types_sorted = [sorted([extract_i[i:i + 2] for i in range(0, len(extract_i), 2)],key=lambda x: x[1]) for extract_i in cluster_types]
  #print(cluster_types_sorted)
  monomers = np.array([np.sum([int(j) for j in np.array(i)[:,0]]) for i in cluster_types_sorted]) == 1
  dt = np.dtype(object)
  #print(np.array(cluster_types,dtype=dt)[monomers])
  monomer_nonunique_types = [i[1] for i in np.array(cluster_types,dtype=dt)[monomers]]
  monomer_types = list(np.unique(monomer_nonunique_types))
  #monomers = []
  for i in monomer_types:
    found = 0
    for j in range(len(monomer_nonunique_types)):
      if i == monomer_nonunique_types[j]:
        found += 1
        if found > 1:
          monomers[j] = False
  #    if i == monomer_nonunique_types[j]:
  #      monomers.append(PRE_monomers[j])
  #      break
  #monomers = np.array(monomers)    
  #print(monomers)
  if Qout != 2:
    print("MONOMERS: " + " ".join(monomer_types),flush = True)
  f = open(".help.txt", "w")
  #TODO
  
  for i in range(len(output[0])):
    #print(i)
      

    #continue
    #TODO
    line = np.array(output)[:,i]
    cluster_total_number = np.sum([int(sel[0]) for sel in cluster_types_sorted[i]])
    #print(cluster_types_sorted[i])
    for j in range(len(cluster_types_sorted[i])):
      cluster_molecule = cluster_types_sorted[i][j][1]
      cluster_molecule_number = cluster_types_sorted[i][j][0]
      #print(np.array(output)[:,monomers])
      test_monomer = 0
      for k in range(len(np.array(output)[:,monomers][0])):
        #selected_monomer = seperate_string_number(np.array(output)[:,monomers][0,k])[1]
        selected_monomer = np.array(cluster_types_sorted,dtype=dt)[monomers][k][0][1]
        #print("--------------------------")
        #print(selected_monomer)
        #print(cluster_molecule)
        if cluster_molecule == selected_monomer:
          for line_i in range(1,len(line)):
            if type(line[line_i]) != type("str"):
              try:
                line[line_i] = float(line[line_i]) - float(cluster_molecule_number) * float(np.array(output)[:,monomers][line_i,k])
                if Qconc > 0:
                  for conc_j in range(len(conc)):
                    if conc[conc_j][0] == selected_monomer:         
                      R = 1.987 #cal/mol/K #=8.31441
                      if np.isnan(Qt):
                        try:
                          Qt = clusters_df["log","temperature"].values[i]
                        except:
                          Qt = 298.15
                      if np.isnan(Qp):
                        try:
                          Qp = 101325.0*float(clusters_df["log","pressure"].values[i])
                        except:
                          Qp = 101325.0
                      conc_mon=float(eval(conc[conc_j][1].replace("ppt","*10**-9*"+str(Qp)).replace("ppb","*10**-6*"+str(Qp)).replace("^","**").replace("cmmc","*10**6*1.380662*10**-23*"+str(Qt)) ))
                      line[line_i] = float(line[line_i]) - QUenergy*(float(cluster_molecule_number) - CNTfactor*float(cluster_molecule_number)/cluster_total_number) * R/1000/627.503 * Qt * np.log( conc_mon / Qp)
              except:
                line[line_i] = missing
          test_monomer = 1
      if test_monomer == 0:
        line[1:] = missing 
    f.write(" ".join(map(str,line))+"\n")
  f.close()
  system("cat .help.txt | column -t")
  remove(".help.txt")

