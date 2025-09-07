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
  nameable_test = True
  if len(input_array) % 2 == 0:
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

def adjustnames(file_basename, QINFOfile_basename, QINFOcluster_type, QINFOcomponents, QINFOcomponent_ratio):
  from re import split
  file_basename_split = file_basename.split("-")[0].split("_")[0]
  split_numbers_letters = split('(\d+)',file_basename_split)[1:]
  #cluster_type_array = seperate_string_number(file_basename_split)
  if is_nameable(split_numbers_letters):
    components = split_numbers_letters[1::2]
    component_ratio = [int(i) for i in split_numbers_letters[0::2]]
    for i in range(len(components)):
      wasthere = False
      for j in range(len(QINFOcomponents)):
        if components[i] == QINFOcomponents[j]:
          QINFOcomponent_ratio[j] += component_ratio[i]
          wasthere = True
          break
      if wasthere == False:
        QINFOcomponents.append(components[i])
        QINFOcomponent_ratio.append(component_ratio[i])
    split_numbers_letters = [str(item) for pair in zip(QINFOcomponent_ratio, QINFOcomponents) for item in pair]
    cluster_type_2array_sorted = sorted([split_numbers_letters[i:i + 2] for i in range(0, len(split_numbers_letters), 2)],key=lambda x: x[1])
    cluster_type_array_sorted = [item for sublist in cluster_type_2array_sorted for item in sublist]
    cluster_type = zeros(cluster_type_array_sorted)
    QINFOcomponents = split_numbers_letters[1::2]
    QINFOcomponent_ratio = [int(i) for i in split_numbers_letters[0::2]]
    QINFOcluster_type = cluster_type
    return QINFOcluster_type, QINFOcluster_type, QINFOcomponents, QINFOcomponent_ratio
  else:
    return "str", None, [], []


def print_help():
  print("""###################################
    JKMD [ SIMULATION_ARGUMENTS ] [-follow [OTHER_SIMULATION_ARGUMENTS]]

  SIMULATION_ARGUMENTS:
    [SPECIE] [[SPECIE(S) CONSTRAINTS]] [[OTHER SPECIE]] [CALCULATOR] [THERMOSTAT] [SIMULATION_SETUP]

  SPECIE:
    -index <int>           index of structure to be taken from pickle file [default = -1 (i.e., last)]
    -indexrange <int>      random index from the range from <int> to -1 (e.g., -10)
    <file>                 xyz or pickle file (last structure taken)
    -char <int>            charge [default = 0]
    -mult <int>            multiplicity [default = 1]
    -recenter              move to [0,0,0]
    -move,-moveby <array>  move center of mass by [x,y,z] vector (e.g., [5,0,0])
    -moveto <array>        move center of mass to [x,y,z] coordinates (e.g., [5,0,0])
    -moveto2 <float> <pos> move center of mass to [x,y,z] coordinates while true x=max(float,x) (e.g., 6 [5,0,0] moves to [6,0,0])
    -mb <int>              initiate vel. from Maxwell-Boltzmann distribution ex. at <int> K
    -setvel <0/1>          0) removes COM velocites, 1) removes all velocities
    -vel <array>           adds velocities as a vector [x,y,z] in Angstrom/fs
    -box <float>           set cell of size LxLxL Angstrom with PBC (must be set for all species)

  CONSTRAINTS:
    -fix                          fixed position of atoms
    -EF_h_A <float>               ext. force field in the form of harmonic potential 
                                  0.5*k*|COM-[x,y,z]|^2 
    -EF_h_A_xyz <pos> <float>     -EF_h_A with centrum of harminc in pos = [x0,y0,z0]
    -EF_fbh_A <float> <float>     ext. force field in the form of flat bott harmonic potential 
                                  0.5*k*(|COM-[x,y,z]|-r0)^2*Heaviside(|COM-[x,y,z]|,r0) 
                                  (e.g., <10> [Ang] <1> kcal/mol/A^2)
    -EF_fbh_A_xyz <pos> <float>   -EF_fbh_A with centrum of flat bott harmonic in pos = [x0,y0,z0]
    -EF_c_COM <array>             constant ext. force on COM (e.g., [1,0,0]) 
    -EF_h_COM_COM <float> <float> harmonic potential between COMs of last two molecules [harm k_bias]
    -EF_h_RMSD <float> <float>    <testing reaction RMSD constrain> [harmonic k_bias]
    UMBRELLA SAMPLING:
    -harm <float>      add harmonic potential COM <float> distance constrain [2 species]
    -k_bias <float>    strength of the biasing harmonic potential in kcal/mol/A^2 [e.g., 100]
    -slow <int>        linearly increases the US potential in <int> steps [default = 0]
    -ha,-heavyatoms    apply only to non-hydrogen atoms       
    -QMMM              the last structure is treated as QM(GFN1-xTB) assuming you use MM(GFNFF)
    -rmsd <int> <file> RMSD to be reached between simulaiton and file structure
 
  CALCULATOR
    -xtb1             GFN1-xTB {TBlite} [set as default]
    -xtb2             GFN2-xTB {TBlite}
    -xtb "<str>"      xtb method (e.g., GFN1-xTB, GFNFF) {XTB}
    -nn_model <path>  Neural Network model defined by path {SchNetPack}
     -EFD               use EFD module {SchNetPack}
     -cutoff <float>    cutoff radius [Angstrom] {SchNetPack} [default = 10.0]
     -dl,-deltalearning use delta learning +GFN1-xTB {SchNetPack}
    -pn_model <path>  Neural Network model defined by path {PhysNet_DER} [requires input.inp too]    -aiment <path>    AIMNet model defined by path {AIMNet}
    -orca "<str>"     QC method (e.g., "XTB1" or "B97-3c") {ORCA}
                      -additional setup might be required!!!
    -max_iter <int>   maximum number of SCF iterations (e.g., 50) [default = 250] {TBlite}
    -mix_damp <float> mixer damping (e.g., 0.1) [default = 0.4] {TBlite}

  THERMOSTAT:
    -vv                  Velocity Verlet [set as default]
    -langevin <float>    Langevin thermostat with friction <float> 1/fs (e.g., 0.01)
    -nose_hoover <float> Nose-Hoover NVT thermostat with <float> time constant (e.g., 25)
    -csvr,bussi <float>  CSVR, stochastic velocity rescaling alg., def. by time const (e.g., 25)
    -andersen <float>    Andersen thermostat with  random collision probability [typically, 0.0001-0.1]
    -fix_COM             the center of mass is not fixed by default
    -temp <float>        temperature (for langevin and nose_hoover) [default = 300]

  SIMULATION SETUP:
    -dt <float>       time step [in fs, default = 1] 
    -ns,-steps <int>  number of steps [default = 1000]
    -dump <int>       dumping properties every <int> step [0 means no dump, default = 1]

  OTHER:
    -nf <str>         folder where the simulation will be performed
    -distout          save distance between two molecules
    -follow           takes the last structure and performes subsequent simulation
    -test(-test2)     see (very) detailed output
    -noex             minimize print on screen
    -nopickle         do not store and save structures
    -seed             set random seed (use at the begginning before specie setup) [TESTING]
    -repeat <int>     repeat x times
    --                use for looping [m=minus] (e.g. 0--10, m1.1--0.1--m0.1, or 0--0.2--5)
    
  EXAMPLES:
      ### Water equilibration followed by longer simulation. 
       JKMD w.xyz -mb 300 -langevin 0.01 -dump 0 -follow -ns 10000 -dump 10  
      ### same as this:
       JKMD w.xyz -mb 300 -char 0                                 #SYSTEM\\\\
                    -dt 1 -ns 1000  -langevin 0.01 -dump 0  -xtb1 #SIMULATION 1\\\\
            -follow -dt 1 -ns 10000 -langevin 0.01 -dump 10 -xtb1 #SIMULATION 2

      ### Construction and Simulation of 1sa1w with NN
       JKMD w.xyz -recenter -add sa.xyz -recenter -move [4,0,0] -nn_model model.pkl -langevin 0.01

      ### Umbrella sampling step
       JKMD w.xyz -mb 300 -langevin 0.01 -nf EQ1
       JKMD simEQ1.pkl -langevin 0.01 -nf EQ2 -ns 100
       JKMD simEQ1.pkl -recenter simEQ2.pkl -moveto [0--0.2--10,0,0] -nose_hoover 25 -ns 10000 -harm 0--0.2--10 -nf US_w_w_NPT

""")


def arguments(argument_list = [], species_from_previous_run = [], charge_from_previous_run = 0, multiplicity_from_previous_run = 1, QINFOcluster_type = "", QINFOcomponents = [], QINFOcomponent_ratio = []):
  from os import path
  missing = float("nan")

  #SPECIES
  Qindex = -1
  Qcharge = charge_from_previous_run
  Qmultiplicity = multiplicity_from_previous_run
  Qtry = 0
  if len(species_from_previous_run) == 0:
    Qindex_of_specie = -1
    species = []
    Qconstraints = 0 #0 = none, 1 = umbrella sampling, 2 = external forces
    Qdistance = 0
    
    Qsavepickle = 1
    Qseed = 42
    from numpy import random
    Qrng = random
    Qout = 1 #output level. 0=only neccessary,1=yes,2=rich print
    Qfolder = ""
  
    #CONSTRAINTS
    Qharm = 10
    Qk_bias = 100
    Qslow = 0
    Qheavyatoms = 0
    QEF = []         #h_COM_COM, h_A, c_COM
    QEF_par = []
    QEF_systems = []
    QMMM = []
  
    #CALCULATOR
    Qcalculator = "XTB1"
    Qcalculator_input = ""
    Qcalculator_max_iterations = 250
    Qmixer_damping = 0.4
    Qcutoff = 10.0
    Qnn_EFD = 0
   
    #THERMOSTAT AND SIMULATION
    Qdt = 1    #timestep
    Qns = 1000 #number of steps
    Qdump = 1  #dump every
    Qsave = -1
    Qtemp = 300
    
    Qthermostat = "VV" #VV = Velocity Verlet, L = Langevin, NH = Nose-Hoover, B = Bussi
    Qthermostat_L = 0.01
    Qthermostat_NH = 25
    Qthermostat_B = 25
    Qfixcm = 0
 
    QINFOfile_basename = "" #IF SOMETHING GETS FUCKED UP, I WILL USE "str" as file_basename
    QINFOcluster_type = ""
    QINFOcomponents = []
    QINFOcomponent_ratio = []

  else:
    Qindex_of_specie = 0
    species = [species_from_previous_run]

  Qfollow_activated = 0
  Qfollow = []

  last = ""
  for i in argument_list:
    ## I have no idea what is this for
    if "--" in i[1:]:
      print("WARNING: Range replaced!")    
      i = "123456789"
    
    #HELP
    if i == "-help" or i == "--help":
      print_help()
      exit()
  
    #PRINT 
    if i == "-print":
      last = "-print"
      continue
    if last == "-print":
      last = ""
      Qout = i
      continue
    #TEST
    if i == "-test":
      Qout = 2
      continue
    if i == "-test2":
      Qout = 3
      continue

    #FOLDER
    if i == "-nf":
      last = "-nf"
      continue
    if last == "-nf":
      last = ""
      Qfolder = i
      continue 

    #NOPICKLE
    if i == "-nopickle":
      Qsavepickle = 0
      continue   

    #TRY
    if i == "-try":
      Qtry = 1
      continue

    ##########################
    ### FOLLOW ###############
    ##########################

    if Qfollow_activated == 1:
      Qfollow.append(i)
      continue

    if i == "-follow":
      Qfollow_activated = 1
      continue

    #NOEXAMPLE
    if i == "-noexample" or i == "-noex":
      Qout = 0
      continue
    #SLOW for umbrella sampling
    if i == "-slow":
      last = "-slow"
      continue
    if last == "-slow":
      Qslow = int(i)
      last = ""
      continue
    #-heavyatoms for umbrella sampling
    if i == "-heavyatoms" or i == "-ha":
      Qheavyatoms = 1
      continue

    #INDEX
    if i == "-index":
      last = "-index"
      continue
    if last == "-index":
      last = ""
      Qindex = int(i)
      continue
    if i == "-indexrange":
      last = "-indexrange"
      continue
    if last == "-indexrange":
      from random import randint
      last = ""
      range_end = abs(int(i))
      print("STRUCTURE USED: -"+str(range_end)+":-1", flush=True)
      Qindex = randint(-range_end, -1)
      continue

    #FIX COM
    if i == "-fix_COM":
      Qfixcm = 1
      continue

    #CALCULATOR
    if i == "-xtb1":
      Qcalculator = "XTB1"  
      continue
    if i == "-xtb2":
      Qcalculator = "XTB2"
      continue
    if i == "-xtb":
      Qcalculator = "XTB"
      last = "-xtb"
      continue
    if last == "-xtb":
      last = ""
      Qcalculator_input = i
      continue
    if i == "-orca":
      Qcalculator = "ORCA"
      last = "-orca"
      continue
    if last == "-orca":
      last = ""
      Qcalculator_input = i
      continue
    if i == "-nn_model":
      import sys, os, glob
      try:
        thepath=glob.glob(os.path.dirname(os.path.abspath(__file__))+'/../../JKQC/JKCS/SCHNETPACK/lib/py*/site-packages/')[0]
      except:
        print("SCHNETPACK was not set properly during setup (run: sh setup.sh -nn -up grendel). [EXITING]")
        exit()
      sys.path.append(thepath)
      Qcalculator = "NN"
      last = "-nn_model"
      continue
    if last == "-nn_model":
      last = ""
      Qcalculator_input = i
      continue
    if i == "-EFD":
      Qnn_EFD = 1
      continue
    if i == "-pn_model":
      Qcalculator = "PhysNet"
      last = "-nn_model"
      continue
    if i == "-aimnet2_model" or i == "-aimnet":
      import sys, os, glob
      try:
        thepath=glob.glob(os.path.dirname(os.path.abspath(__file__))+'/../../JKQC/JKCS/AIMNET/lib/py*/site-packages/')[0]
      except:
        print("AIMNET was not set properly during setup (run: sh setup.sh -aimnet -up grendel). [EXITING]")
        exit()
      sys.path.append(thepath)
      Qcalculator = "AIMNET2"
      last = "-aimnet2_model"
      continue
    if last == "-aimnet2_model":
      last = ""
      Qcalculator_input = i
      continue
    if i == "-max_iter":
      last = "-max_iter"
      continue
    if last == "-max_iter":
      last = ""
      Qcalculator_max_iterations = int(i)
      continue
    if i == "-mix_damp":
      last = "-mix_damp"
      continue
    if last == "-mix_damp":
      last = ""
      Qmixer_damping = float(i)
      continue
    if i == "-cutoff":
      last = "-cutoff"
      continue
    if last == "-cutoff":
      last = ""
      Qcutoff = float(i)
      continue

    #SPECIES
    if i[-4:] == ".xyz" and last == "":
      from ase.io import read
      if not QINFOfile_basename == "str": 
        file_basename = i[:-4][::-1].split("/",1)[0][::-1]
        QINFOfile_basename, QINFOcluster_type, QINFOcomponents, QINFOcomponent_ratio = adjustnames(file_basename, QINFOfile_basename, QINFOcluster_type, QINFOcomponents, QINFOcomponent_ratio)
      species.append(read(i,"-2"))
      
      if Qindex_of_specie == -1:
        Qlenfirst = len(species[0])
      Qindex_of_specie = len(species) - 1
      Qindex = -1
      continue
    if i[-4:] == ".pkl" and last == "":
      from pandas import read_pickle
      if not QINFOfile_basename == "str":
        file_basename = read_pickle(i).iloc[Qindex][("info","file_basename")]
        QINFOfile_basename, QINFOcluster_type, QINFOcomponents, QINFOcomponent_ratio = adjustnames(file_basename, QINFOfile_basename, QINFOcluster_type, QINFOcomponents, QINFOcomponent_ratio)
      species.append(read_pickle(i).iloc[Qindex][("xyz","structure")])
      if Qindex_of_specie == -1:
        Qlenfirst = len(species[0])
      Qindex_of_specie = len(species) - 1
      Qindex = -1
      continue

    # SPLIT
    if i == "-split":
      last = "-split"
      continue
    if last == "-split":
      last = ""
      split = int(i)
      if split > len(species[Qindex_of_specie]):
        print("You are trying to split somehing that has only length "+len(species[Qindex_of_specie])+".")
        exit()
      species.append(species[Qindex_of_specie][split:])
      species[Qindex_of_specie] = species[Qindex_of_specie][0:split]
      if Qindex_of_specie == 0:
        Qlenfirst = len(species[0])
      Qindex_of_specie = len(species) - 1
      Qindex = -1
      continue
 
    # CHARGE
    if i == "-chrg" or i == "-char":
      last = "-char"
      continue
    if last == "-char":
      Qcharge = Qcharge + int(i)
      last = ""
      continue
    # MULITPLICITY
    if i == "-mult":
      last = "-mult"
      continue
    if last == "-mult":
      Qmultiplicity = int(i)
      last = ""
      continue
 
    # SEED
    if i == "-seed":
      last = "-seed"
      continue
    if last == "-seed":
      Qseed = int(i)
      from numpy import random
      from random import seed
      seed(Qseed)
      Qrng = random.default_rng(Qseed)
      last = ""
      continue

    #RECENTER
    if i == "-recenter":
      species[Qindex_of_specie].translate(-species[Qindex_of_specie].get_center_of_mass())
      continue
  
    #MOVE
    if i == "-move" or i == "-moveby":
      last = "-move"
      continue
    if last == "-move":
      last = ""
      from ast import literal_eval
      species[Qindex_of_specie].translate(literal_eval(i))
      continue
    if i == "-moveto":
      last = "-moveto"
      continue
    if last == "-moveto":
      last = ""
      from ast import literal_eval
      themove=literal_eval(i)
      #if themove[0]<=6:
      #  themove[0]=6
      #ix = species[Qindex_of_specie]
      #x[1:3].translate(-x[1:3].get_center_of_mass()+themove)
      #species[Qindex_of_specie] = x
      species[Qindex_of_specie].translate(-species[Qindex_of_specie].get_center_of_mass()+themove)
      continue
    if i == "-moveto2":
      last = "-moveto2"
      continue
    if last == "-moveto2":
      last = "-moveto2b"
      moveTHR=float(i)
      continue
    if last == "-moveto2b":
      last = ""
      from ast import literal_eval
      themove=literal_eval(i)
      if themove[0]<=moveTHR:
        themove[0]=moveTHR
      species[Qindex_of_specie].translate(-species[Qindex_of_specie].get_center_of_mass()+themove)
      continue
  
    #SELECT
    if i == "-select":
      last = "-select"
      continue
    if last == "-select":
      last = ""
      if int(i) >= len(species[Qindex_of_specie]) or int(i) < 0:
        print("Very weird try of selecting a specie that does not exist. Check your -select!")
        exit()
      Qindex_of_specie = int(i)
      continue

    #INITIATE VELOCITIES
    if i == "-mb":
      last = "-mb"
      continue
    if last == "-mb":
      last = ""
      from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
      MaxwellBoltzmannDistribution(species[Qindex_of_specie], temperature_K = int(i), force_temp = True, rng=Qrng)
      continue
    #SETVEL 0/1
    if i == "-setvel":
      last = "-setvel"
      continue
    if last == "-setvel":
      last = ""
      if int(i) == 0:
        from ase.md.velocitydistribution import Stationary
        Stationary(species[Qindex_of_specie])
        continue
      elif int(i) == 1:
        species[Qindex_of_specie].set_velocities(0*species[Qindex_of_specie].get_velocities())
        continue
      else:
        print("WTF")
        exit()
    #vel 
    if i == "-vel":
      last = "-vel"
      continue
    if last == "-vel":
      last = ""
      from ast import literal_eval
      from ase.units import Ang, fs
      from numpy import array
      species[Qindex_of_specie].set_velocities(species[Qindex_of_specie].get_velocities()+array(literal_eval(i))*Ang/fs)
      continue

    #BOX
    if i == "-box":
      last = "-box"
      continue
    if last == "-box":
      last = ""
      species[Qindex_of_specie].set_cell([float(i),float(i),float(i)])
      species[Qindex_of_specie].set_pbc([1,1,1])
      continue
 
    #TIMESTEP
    if i == "-dt":
      last = "-dt"
      continue
    if last == "-dt":
      last = ""
      Qdt = float(i)
      continue
  
    #NUMBER OF STEPS 
    if i == "-ns" or i == "-steps":
      last = "-ns"
      continue
    if last == "-ns":
      last = ""
      Qns = int(i)
      continue
  
    # DUMPING
    if i == "-dump":
      last = "-dump"
      continue
    if last == "-dump":
      last = ""
      Qdump = int(i)
      continue
    # SAVE
    if i == "-save":
      last = "-save"
      continue
    if last == "-save":
      last = ""
      Qsave = int(i)
      continue

    #TEMPERATURE
    if i == "-temp":
      last = "-temp"
      continue
    if last == "-temp":
      last = ""
      Qtemp = float(i)
      continue

    #THERMOSTATS
    if i == "-vv":
      Qthermostat = "VV"
      continue
    if i == "-langevin":
      last = "-langevin"
      continue
    if last == "-langevin":
      last = ""
      Qthermostat = "L"
      Qthermostat_L = float(i)
      continue
    if i == "-nose_hoover":
      last = "-nose_hoover"
      continue
    if last == "-nose_hoover":
      last = ""
      Qthermostat = "NH"
      Qthermostat_NH = float(i)
      continue
    if i == "-csvr" or i == "-bussi":
      last = "-csvr"
      continue
    if last == "-csvr":
      last = ""
      Qthermostat = "B"
      Qthermostat_B = float(i)
      continue
    if i == "-andersen":
      last = "-andersen"
      continue
    if last == "-andersen":
      last = ""
      Qthermostat = "A"
      Qthermostat_A = float(i)
      continue
    if i == "-opt":
      Qthermostat = "OPT"
      continue
 
    #CONSTRAINTS
    if i == "-rmsd":
      last = "-rmsd"
      continue
    if last == "-rmsd":
      last = "-rmsd2"
      Qconstraints = 2
      Qrmsddiff = float(i)
      continue
    if last == "-rmsd2":
      last = ""
      Qrmsdreffile = str(i)
      continue
    if i == "-harm":
      last = "-harm" 
      Qconstraints = 1
      continue
    if last == "-harm":
      last = ""
      Qharm = float(i)
      continue
    if i == "-k_bias":
      last = "-k_bias"
      continue
    if last == "-k_bias":
      last = ""
      Qk_bias = float(i)
      continue
    if i == "-distout":
      Qdistance = 1
      continue
    if i == "-qmmm":
      mfrom = 0
      for j in range(Qindex_of_specie):
        mfrom=mfrom+len(species[j])
      mto=mfrom+len(species[Qindex_of_specie])
      QMMM.append([mfrom,mto])
      continue
    if i == "-dl" or i == "-deltalearning":
      QEF_systems.append("all")
      QEF.append("deltalearning")
      QEF_par.append("GFN1-xTB")
      continue

    ### EXTERNAL FORCES ###
    ## HARMONIC POTENTIAL ##
    if i == "-EF_h_A":
      last = "-EF_h_A"
      QEF.append("h_A")
      rem = None
      continue
    if i == "-EF_h_A_xyz":
      last = "-EF_h_A_xyz"
      QEF.append("h_A_xyz")
      continue
    if last == "-EF_h_A_xyz":
      last = "-EF_h_A"
      rem = literal_eval(i)
      continue
    if last == "-EF_h_A":
      last = ""
      mfrom = 0
      for j in range(Qindex_of_specie):
        mfrom=mfrom+len(species[j])
      mto=mfrom+len(species[Qindex_of_specie])
      QEF_systems.append([mfrom,mto])
      QEF_par.append([rem,float(i)])
      continue
    ## FLAT BOTTOM HARMONIC POTENTIAL ##
    if i == "-EF_fbh_A":
      last = "-EF_fbh_A"
      QEF.append("fbh_A")
      rem = None
      continue
    if i == "-EF_fbh_A_xyz":
      last = "-EF_fbh_A_xyz"
      QEF.append("fbh_A_xyz")
      continue
    if last == "-EF_fbh_A_xyz":
      last = "-EF_fbh_A"
      rem = literal_eval(i)
      continue
    if last == "-EF_fbh_A":
      last = "-EF_fbh_A_2"
      mfrom = 0
      for j in range(Qindex_of_specie):
        mfrom=mfrom+len(species[j])
      mto=mfrom+len(species[Qindex_of_specie])
      QEF_systems.append([mfrom,mto])
      rem2 = float(i)
      continue
    if last == "-EF_fbh_A_2":
      last = ""
      QEF_par.append([rem,rem2,float(i)])
      continue
    ## CONSTANT EXTERNAL FORCE ##
    if i == "-EF_c_COM":
      last = "-EF_c_COM"
      continue
    if last == "-EF_c_COM":
      last = ""
      QEF.append("c_COM")
      mfrom = 0
      for j in range(Qindex_of_specie):
        mfrom=mfrom+len(species[j])
      mto=mfrom+len(species[Qindex_of_specie])
      QEF_systems.append([mfrom,mto])
      QEF_par.append(literal_eval(i))
      #TODO
      continue
    ## COM COM DISTANCE ##
    if i == "-EF_h_COM_COM":
      last = "-EF_h_COM_COM"
      continue
    #TODO THIS IS A BIT WEIRD
    if last == "-EF_h_COM_COM":
      last = "-EF_h_COM_COM_2"
      QEF.append("h_COM_COM")
      continue
    if last == "-EF_h_COM_COM_2":
      last = "-EF_h_COM_COM_3"
      x=float(i)
      continue
    if last == "-EF_h_COM_COM_3":
      last = ""
      QEF_par.append([x,float(i)])
      QEF_systems.append([Qspecies-1,Qspecies])
      continue
    ## RMSD HARMONIC ##
    if i == "-EF_h_RMSD":
      last = "-EF_h_RMSD"
      QEF.append("h_RMSD")
      continue
    if last == "-EF_h_RMSD":
      last = "-EF_h_RMSD_2"
      x=float(i)
      continue
    if last == "-EF_h_RMSD_2":
      last = ""
      QEF_par.append([x,float(i)])
      QEF_systems.append([])
      continue
    if i == "-fix":
      from ase.constraints import FixAtoms
      species[Qindex_of_specie].set_constraint(FixAtoms(indices=range(len(species[Qindex_of_specie]))))
      continue
    ##########################
 
    #UNKNOWN ARGUMENT
    print("I am sorry but I do not understand the argument: "+i+" [EXITING]")
    exit()

  if last != "":
    print(last)
    print("Hey looser, the last argument is incomplete")
    exit()

  print("==============")
  print("Total charge: " + str(Qcharge))
  print("Total multiplicity: " + str(Qmultiplicity))
  print("Number of species: " + str(len(species)))
  print("==============")
  for i in range(len(species)):
    print("SPECIE "+str(i)+": "+str(species[i]))
    print("COM: "+str(species[i].get_center_of_mass()))
    print("COORDINATES:")
    print(species[i].get_positions())
    print("VELOCITIES:")
    print(species[i].get_velocities())
    print("==============")
  print("NAME: " + str(QINFOcluster_type))
  print("FOLLOW: ")
  print(Qfollow)
  print("==============", flush=True)
 
  if Qtry == 1: 
    exit()

  return locals()
