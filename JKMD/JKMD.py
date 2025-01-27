# IMPORT NECESSARY LIBRARIES
from sys import argv
from os import system

#print command to the output file
#cmd="".join(["( echo COMMAND: JKQC "," ".join(argv[1:])," >> output ) 2>/dev/null"])
#system(cmd)

import os
import psutil
#os.system("if ! command -v module &> /dev/null; then source /com/bin/modules.sh; fi; module load intel; module load openmpi;")
#os.system("if ! command -v module &> /dev/null; then source /com/bin/modules.sh; fi; module load gcc openmpi mkl")
#os.environ['OMP_STACKSIZE'] = '4G'
try:
  os.environ['OMP_NUM_THREADS'] = f'{len(psutil.Process().cpu_affinity())}'
except:
  os.environ['OMP_NUM_THREADS'] = str(1)
                               #f'{len(psutil.Process().cpu_affinity())},1'
#os.environ['OMP_MAX_ACTIVE_LEVELS'] = '1'
#import resource
#resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


print("""
   JJJJJJ KK   KK MMM    MMM DDDDD
   JJ  JJ KK  KK  MMMM  MMMM DDDDDD
       JJ KK KK   MM MMMM MM DD   DD
   JJ  JJ KKKK    MM  MM  MM DD   DD 
   JJ  JJ KKKKK   MM      MM DD   DD 
   JJJJJJ KK KKK  MM      MM DDDDDD
    JJJJ  KK  KKK MM      MM DDDDD
""")



current_time = 0
current_step = 0

#READ ARGUMENTS
from arguments import arguments
Qfollow_activated = -1
QEF_applied = 0
while not Qfollow_activated == 0:
  if Qfollow_activated == -1:
    locals().update(arguments(argv[1:]))
    if Qfollow_activated == 1 and Qout == 2:
      print([Qfollow,Qfollow_activated])
  else:
    if Qout == 2:
      print("Next -follow")
    locals().update(arguments(Qfollow,all_species,charge_from_previous_run=Qcharge))
  if Qout == 2:
    from time import time
    start = time()
    print("DONE] Time started: "+str(time() - start));
  
  #CREATE THE SYSTEM
  if Qout == 2:
    print("Combining species.")
  all_species = species[0]
  for i in range(1,len(species)):
    all_species = all_species + species[i]
  if Qout == 2:
    print(all_species)
  
  #CONSTRAINTS
  constraints = []
  if Qconstraints == 1:
    if len(species) != 2:
      print("Nice try. The umbrella sampling part of JKMD is yet not ready for your jokes.")
      exit()
    from umbrellaconstraint import UmbrellaConstraint
    constraints.append(UmbrellaConstraint(all_species,Qk_bias,len(species[0]),Qharm,Qslow))
    Qconstraints = 2
  if Qdistance == 1:
    if len(species) != 2:
      print("Nice try. The -distout part of JKMD is yet not ready for your jokes.")
      exit()
  #EXTERNAL FORCE
  if len(QEF) > 0:
    for i in range(QEF_applied,len(QEF)):
      QEF_applied += 1
      if QEF[i] == "h_A" or QEF[i] == "fbh_A" or QEF[i] == "c_COM":
        from externalforce import ExternalForce
        constraints.append(ExternalForce(QEF[i],QEF_par[i],QEF_systems[i]))
        if Qout == 2:
           print("External FF applied")
           print("External Force: "+QEF[i]+" on "+str(QEF_systems[i][0])+" to "+str(QEF_systems[i][1])+" with parameters "+str(QEF_par[i]))
      if QEF[i] == "h_COM_COM":
        if len(species) != 2:
          print("Nice try. The umbrella sampling part of JKMD is yet not ready for your jokes.")
          exit()
        from umbrellaconstraint import UmbrellaConstraint
        constraints.append(UmbrellaConstraint(all_species,QEF_par[i][0],len(species[0]),QEF_par[i][1],Qslow))
  #SET CONSTRAINTS
  if len(constraints) > 0:
    all_species.set_constraint(constraints)
  
  #SET CALCULATOR
  if Qout == 2:
    print("Setting calculator.")
  def call_calculator():
    from calculator import calculator
    all_species.calc = calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations, Qcharge, Qout, all_species)
  call_calculator()
  if Qout == 2:
    print(all_species)
    #print(all_species.get_positions())
    #print(all_species.get_potential_energy())
    #exit()
  
  #THERMOSTAT
  if Qthermostat == "VV":
    from ase import units
    from ase.md.verlet import VelocityVerlet
    dyn = VelocityVerlet(all_species, Qdt * units.fs)
    if Qfixcm == 1:
      print("Sorry did not check yet how to fix COM for VV")
      exit()
    #TODO
    #dyn.zero_center_of_mass_momentum(verbose = 1)
  elif Qthermostat == "L":
    from ase import units
    from ase.md.langevin import Langevin
    dyn = Langevin(all_species, Qdt * units.fs, temperature_K = Qtemp, friction = Qthermostat_L / units.fs, fixcm = Qfixcm, rng = Qrng)
  elif Qthermostat == "NH":
    from ase import units
    from ase.md.npt import NPT
    from numpy import identity
    all_species.set_cell(2 * identity(3))
    #TODO the npt library is outdated and I adjusted it manually to make this working
      #line 148
      #if externalstress is not None:
      #  self.set_stress(externalstress)
    dyn = NPT(atoms = all_species, timestep = Qdt * units.fs, temperature_K = Qtemp, ttime = Qthermostat_NH * units.fs, externalstress = None)
    #from ase.md.nptberendsen import NPTBerendsen
    #dyn = NPTBerendsen(all_species, timestep=0.1 * units.fs, temperature_K=300,
    #               taut=100 * units.fs, pressure_au=1.01325 * units.bar,
    #               taup=1000 * units.fs, compressibility_au=4.57e-5 / units.bar)
    if Qfixcm == 1:
      dyn.zero_center_of_mass_momentum(verbose = 1)
  elif Qthermostat == "B":
    from ase import units
    #from ase.md.bussi import Bussi
    from ase_bussi import Bussi
    dyn = Bussi(all_species, Qdt * units.fs, temperature_K = Qtemp, taut = Qthermostat_NH * units.fs, rng = Qrng)
  else:
    print("Some weird thermostat.")
    exit()
 
  #DUMPING
  if Qdump != 0:
    from print_properties import print_properties, init
    def mergeDictionary(dict_1, dict_2):
      if len(dict_1) == 0:
        return dict_2
      dict_3 = {**dict_1, **dict_2}
      for key, value in dict_3.items():
        if key in dict_1 and key in dict_2:
          dict_3[key] = value + dict_1[key]
        elif key in dict_1:
          dict_3[key] = [float("nan")] + value
        else:
          ml = max([len(x) for x in dict_1.values()])
          dict_3[key] = value + ml*[float("nan")]
      return dict_3
      #dict1 = dist3

    init(current_time,current_step)
    #global cluster_dic
    if current_step == 0:
      cluster_dic = {}
    def save(fail = False):
      global cluster_dic,current_time,current_step
      toupdate,current_time,current_step = print_properties(species = all_species, timestep = Qdt, interval = Qdump, Qconstraints = Qconstraints, Qdistance = Qdistance, split = Qlenfirst, fail = fail)
      if Qsavepickle == 1:
        toupdate.update({("log","method"):[" ".join(argv[1:])],("log","program"):["Python"]})
        if Qconstraints != 0:
          toupdate.update({("log","k_bias"):[min(current_step/max(Qslow,0.0000001),1)*Qk_bias],("log","harm_distance"):[Qharm]})
        cluster_dic = mergeDictionary(cluster_dic, toupdate)
    dyn.attach(save, interval = Qdump)
    if Qcalculator == "PhysNet": 
      from calculator import calculator
      def updatephysnet():
        all_species.calc = calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations, Qcharge, Qout, all_species)
      dyn.attach(updatephysnet, interval = 1)
    #dyn.attach(mergeDictionary(cluster_dic, print_properties(species = all_species, timestep = Qdt, interval = Qdump)), interval = Qdump) 

  #SIMULATION
  if 'current_step' in globals():
    steps_before_sim = current_step
  else:
    steps_before_sim = 0
  sim_errors = 0
  sim_last_error = -1
  while current_step - steps_before_sim <= Qns:
    try:
      dyn.run(Qns - (current_step - steps_before_sim))
    except Exception as e:
      from numpy import random
      positions = all_species.get_positions()
      noise = random.normal(scale=0.01, size=positions.shape)
      all_species.set_positions(positions + noise)
      call_calculator()
      print(str(e))
      print("Something got screwed up within the dyn.run(Qns). Small adjustment to the structure!!!")
      if current_step == sim_last_error:
        sim_errors += 1
      else:
        sim_errors = 1
        sim_last_error = current_step
      if sim_errors == 4:
        if "cluster_dic" not in globals():
          print("I have nothing to save as the run failed immediately.")
        else:
          print("I have saved the error structure in the appropriate folder.")
          Qsavepickle = 1
          save(fail = True)
          print(Qfolder)
          from pandas import DataFrame
          for key in cluster_dic:
            cluster_dic[key] = cluster_dic[key][::-1]
          clusters_df = DataFrame(cluster_dic)
          clusters_df.to_pickle(Qfolder+"/error.pkl")
        exit()
    if Qdump == 0:
      current_time = current_time + Qdt*Qns
      current_step = current_step + Qns

  if Qout == 2:
    print("Simulation round done.")

if Qsavepickle == 1:
  if Qout == 2:
    print("Done and now just saving pickle.")
  from pandas import DataFrame
  for key in cluster_dic:
    cluster_dic[key] = cluster_dic[key][::-1]
  clusters_df = DataFrame(cluster_dic) #, index = range(len(cluster_dic)))
  try:
    clusters_df.to_pickle(Qfolder+"/../sim"+Qfolder.split("/")[-1]+".pkl")
    print("The sim"+Qfolder.split("/")[-1]+".pkl has been hopefully created.")
  except:
    print("Something got fucked up.")

print("Done.")
