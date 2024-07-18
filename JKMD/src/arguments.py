def print_help():
  print("""###################################
    JKMD [ SIMULATION_ARGUMENTS ] [-follow [OTHER_SIMULATION_ARGUMENTS]]

  SIMULATION_ARGUMENTS:
    [ SPECIE ] [[OTHER SPECIE]] [ CALCULATOR ] [ THERMOSTAT ] [ SIMULATION_SETUP ] [[ CONSTRAINS ]]

  SPECIE:
    -index <int>          index of structure to be taken from pickle file [default = -1 (i.e., last)]
    <file>                xyz or pickle file (last structure taken)
    -char <int>           charge [default = 0]
    -recenter             move to [0,0,0]
    -move,-moveby <array> move center of mass by [x,y,z] vector (e.g., [5,0,0])
    -moveto <array>       move center of mass to [x,y,z] coordinates (e.g., [5,0,0])
    -mb <int>             initiate vel. from Maxwell-Boltzmann distribution ex. at <int> K
    -setvel <0/1>         0) removes COM velocites, 1) removes all velocities
    -vel <array>          adds velocities as a vector [x,y,z] in Angstrom/fs

  CALCULATOR
    -xtb1            GFN1-xTB {XTBlite} [set as default]
    -xtb2            GFN2-xTB {XTBlite}
    -xtb "<str>"     xtb method (e.g., GFN1-xTB) {XTB}
    -nn_model <path> Neural Network model defined by path {SchNetPack}
    -orca "<str>"    QC method (e.g., "XTB1" or "B97-3c") {ORCA}
                     -additional setup might be required!!!

  THERMOSTAT:
    -vv                  Velocity Verlet [set as default]
    -langevin <float>    Langevin thermostat with friction <float> 1/fs (e.g., 0.01)
    -nose_hoover <float> Nose-Hoover NVT thermostat with <float> time constant (e.g., 25)
    -csvr,bussi <float>  CSVR, stochastic velocity rescaling alg., def. by time const (e.g., 25)
    -fix_COM             the center of mass is not fixed by default
    -temp <float>        temperature (for langevin and nose_hoover) [default = 300]

  SIMULATION SETUP:
    -dt <float>       time step [in fs, default = 1] 
    -ns,-steps <int>  number of steps [default = 1000]
    -dump <int>       dumping properties every <int> step [0 means no dump, default = 1]

  CONSTRAIN:
    UMBRELLA SAMPLING:
    -harm <float>     add harmonic potential COM <float> distance constrain [2 species]
    -k_bias <float>   strength of the biasing harmonic potential in kcal/mol/A^2 [e.g., 100]
 
  OTHER:
    -nf <str>         folder where the simulation will be performed
    -follow           takes the last structure and performes subsequent simulation
    -test             see detailed output
    -noex             minimize print on screen
    -repeat <int>     repeat x times
    --                use for looping (e.g. 0--10, or 0--0.2--5)
    
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


def arguments(argument_list = [], species_from_previous_run = [], charge_from_previous_run = 0):
  from os import path
  missing = float("nan")

  #SPECIES
  Qindex = -1
  Qcharge = charge_from_previous_run
  if len(species_from_previous_run) == 0:
    Qindex_of_specie = -1
    species = []
    Qconstraints = 0

    Qseed = 42
    Qout = 1 #output level. 0=only neccessary,1=yes,2=rich print
    Qfolder = ""
  
    #CONSTRAINTS
    Qharm = 10
    Qk_bias = 100
  
    #CALCULATOR
    Qcalculator = "XTB1"
    Qcalculator_input = ""
   
    #THERMOSTAT AND SIMULATION
    Qdt = 1    #timestep
    Qns = 1000 #number of steps
    Qdump = 1  #dump every
    Qtemp = 300
    
    Qthermostat = "VV" #VV = Velocity Verlet, L = Langevin, NH = Nose-Hoover, B = Bussi
    Qthermostat_L = 0.01
    Qthermostat_NH = 25
    Qthermostat_B = 25
    Qfixcm = 0

  else:
    Qindex_of_specie = 0
    species = [species_from_previous_run]

  Qfollow_activated = 0
  Qfollow = []

  last = ""
  for i in argument_list:
    #HELP
    if i == "-help" or i == "--help":
      print_help()
      exit()
  
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
    #TEST
    if i == "-test":
      Qout = 2
      continue

    #INDEX
    if i == "-index":
      last = "-index"
      continue
    if last == "-index":
      last = ""
      Qindex = int(i)
      continue

    #FOLDER
    if i == "-nf":
      last = "-nf"
      continue
    if last == "-nf":
      last = ""
      Qfolder = i
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
      Qcalculator = "NN"
      last = "-nn_model"
      continue
    if last == "-nn_model":
      last = ""
      Qcalculator_input = i
      continue

    #SPECIES
    if i[-4:] == ".xyz":
      from ase.io import read
      species.append(read(i,"-2"))
      if Qindex_of_specie == -1:
        Qlenfirst = len(species[0])
      Qindex_of_specie = Qindex_of_specie + 1
      Qindex = -1
      continue
    if i[-4:] == ".pkl":
      from pandas import read_pickle
      species.append(read_pickle(i).iloc[Qindex][("xyz","structure")])
      if Qindex_of_specie == -1:
        Qlenfirst = len(species[0])
      Qindex_of_specie = Qindex_of_specie + 1
      Qindex = -1
      continue
 
    # CHARGE
    if i == "-chrg" or i == "char":
      last = "-char"
      continue
    if last == "-char":
      Qcharge = Qcharge + int(i)
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
      species[Qindex_of_specie].translate(-species[Qindex_of_specie].get_center_of_mass()+literal_eval(i))
      continue

    #INITIATE VELOCITIES
    if i == "-mb":
      last = "-mb"
      continue
    if last == "-mb":
      last = ""
      from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
      MaxwellBoltzmannDistribution(species[Qindex_of_specie], temperature_K = int(i), force_temp = True)
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
      species[Qindex_of_specie].set_velocities(species[Qindex_of_specie].get_velocities()+literal_eval(i))
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
 
    #CONSTRAINTS
    if i == "-harm":
      last = "-harm" 
      continue
    if last == "-harm":
      last = ""
      Qconstraints = 1
      Qharm = float(i)
      continue
    if i == "-k_bias":
      last = "-k_bias"
      continue
    if last == "-k_bias":
      last = ""
      Qconstraints = 1
      Qk_bias = float(i)
      continue

    #UNKNOWN ARGUMENT
    print("I am sorry but I do not understand the argument: "+i+" [EXITING]")
    exit()

  if last != "":
    print(last)
    print("Hey looser, the last argument is incomplete")
    exit()

  return locals()
