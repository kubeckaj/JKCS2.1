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
#os.environ['MKL_STACKSIZE'] = '4G'
#try:
#  os.environ['OMP_NUM_THREADS'] = f'{len(psutil.Process().cpu_affinity())}'
#  os.environ['MKL_NUM_THREADS'] = f'{len(psutil.Process().cpu_affinity())}'
#  os.environ['OPENBLAS_NUM_THREADS'] = f'{len(psutil.Process().cpu_affinity())}'
#except:
os.environ['OMP_NUM_THREADS'] = str(1)
os.environ['MKL_NUM_THREADS'] = str(1)
os.environ['OPENBLAS_NUM_THREADS'] = str(1)
os.environ['NUMEXPR_MAX_THREADS'] = str(1)

#                               #f'{len(psutil.Process().cpu_affinity())},1'
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

def savepickle():
  global cluster_dic
  if "cluster_dic" in globals() or "cluster_dic1" in globals() or "cluster_dic2" in globals():
    from pandas import DataFrame
    if Qconstraints == 4 and len(species) == 2:
      global cluster_dic1,cluster_dic2
      for key in cluster_dic1:
        cluster_dic2[key] = cluster_dic2[key][::-1]
      cluster_dic = mergeDictionary(cluster_dic1, cluster_dic2)
    else: 
      for key in cluster_dic:
        cluster_dic[key] = cluster_dic[key][::-1]
    clusters_df = DataFrame(cluster_dic) #, index = range(len(cluster_dic)))
    global Qfolder
    try:
      clusters_df.to_pickle(Qfolder+"/../sim"+Qfolder.split("/")[-1]+".pkl")
      print("The sim"+Qfolder.split("/")[-1]+".pkl has been hopefully created.")
    except:
      print("Something got fucked up.")

#print(f"Using {os.environ['OMP_NUM_THREADS']} threads.")
current_time = 0
current_step = 0


#READ ARGUMENTS
from arguments import arguments
Qfollow_activated = -1
QEF_applied = 0
while not Qfollow_activated == 0:
  if Qfollow_activated == -1:
    locals().update(arguments(argv[1:]))
    if Qfollow_activated == 1 and Qout > 1:
      print([Qfollow,Qfollow_activated])
  else:
    if Qout > 1:
      print("Next -follow")
    locals().update(arguments(Qfollow,all_species,charge_from_previous_run=Qcharge,multiplicity_from_previous_run=Qmultiplicity, QINFOcluster_type = QINFOcluster_type, QINFOcomponents = QINFOcomponents, QINFOcomponent_ratio = QINFOcomponent_ratio))
  if Qout > 1:
    from time import time
    start = time()
    print("DONE] Time started: "+str(time() - start));
  
  #CREATE THE SYSTEM
  if Qout > 1:
    print("Combining species.")
  all_species = species[0]
  for i in range(1,len(species)):
    all_species = all_species + species[i]
  if Qout > 1:
    print(all_species)
  
  #CONSTRAINTS
  constraints = []
  if Qconstraints == 1:
    if len(species) != 2:
      print("Nice try. The umbrella sampling part of JKMD is yet not ready for your jokes.")
      exit()
    from umbrellaconstraint import UmbrellaConstraint
    constraints.append(UmbrellaConstraint(all_species,Qk_bias,len(species[0]),Qharm,Qslow,Qheavyatoms))
    Qconstraints = 3
  if Qconstraints == 2:
    if not ( len(species) == 1 or len(species) == 2):
      print("Nice try. The RMSD part of JKMD is yet not ready for your jokes.")
      exit()
    from rmsdconstraint import RMSDConstraint
    from energyconstrain import EnergyConstraint
    if len(species) == 2:
      #constraints.append(RMSDConstraint(species[0],Qk_bias,Qrmsddiff,species[0],Qslow))
      #constraints.append(RMSDConstraint(species[0],Qk_bias,Qrmsddiff,species[1],Qslow))
      #constrain1 = RMSDConstraint(species[0],Qk_bias,Qrmsddiff,species[1],Qslow)
      #constrain2 = RMSDConstraint(species[1],Qk_bias,Qrmsddiff,species[0],Qslow)
      #species[0].set_constraint(constrain1)
      from ArbAlign import compare
      #print(species[0].get_chemical_symbols())
      #print(species[1].get_chemical_symbols())
      species[0], _ = compare(species[1], species[0], Qreturn_geometry = 1)
      species[1], _ = compare(species[0], species[1], Qreturn_geometry = 1)
      #print(species[0].get_chemical_symbols())
      #print(species[1].get_chemical_symbols())
      
      def call_rmsd_constrain1():
        species[0].set_constraint([EnergyConstraint(species[0].get_potential_energy(),species[1].get_potential_energy()),RMSDConstraint(species[0],Qk_bias,Qrmsddiff,species[1],Qslow)])
        #species[0].set_constraint(RMSDConstraint(species[0],Qk_bias,Qrmsddiff,species[1],Qslow))
        #constrain1.ref = species[1]
      def call_rmsd_constrain2():
        #constrain2.ref = species[0]
        species[1].set_constraint([EnergyConstraint(species[1].get_potential_energy(),species[0].get_potential_energy()),RMSDConstraint(species[1],Qk_bias,Qrmsddiff,species[0],Qslow)])  
        #species[1].set_constraint(RMSDConstraint(species[1],Qk_bias,Qrmsddiff,species[0],Qslow))
      #call_rmsd_constrain1()
      #call_rmsd_constrain2()
      print(species)
    else:
      constraints.append(RMSDConstraint(all_species,Qk_bias,Qrmsddiff,Qrmsdreffile,Qslow))
    Qconstraints = 4
  if len(QMMM) > 0:
    from src.QMMM import QMMM as QMMMcalc  
    for i in range(len(QMMM)):
      constraints.append(QMMMcalc(QMMM[i]))
  if Qdistance == 1:
    if len(species) != 2:
      print("Nice try. The -distout part of JKMD is yet not ready for your jokes.")
      exit()
  #EXTERNAL FORCE
  if len(QEF) > 0:
    for i in range(QEF_applied,len(QEF)):
      QEF_applied += 1
      if QEF[i] == "h_A" or QEF[i] == "h_A_xyz" or QEF[i] == "fbh_A" or QEF[i] == "fbh_A_xyz" or QEF[i] == "c_COM":
        from externalforce import ExternalForce
        constraints.append(ExternalForce(QEF[i],QEF_par[i],QEF_systems[i]))
        if Qout > 1:
           print("External FF applied")
           print("External Force: "+QEF[i]+" on "+str(QEF_systems[i][0])+" to "+str(QEF_systems[i][1])+" with parameters "+str(QEF_par[i]))
      if QEF[i] == "h_COM_COM":
        if len(species) != 2:
          print("Nice try. The umbrella sampling part of JKMD is yet not ready for your jokes.")
          exit()
        from umbrellaconstraint import UmbrellaConstraint
        constraints.append(UmbrellaConstraint(all_species,QEF_par[i][0],len(species[0]),QEF_par[i][1],Qslow))
      if QEF[i] == "deltalearning":
        from deltalearning import DeltaLearning
        constraints.append(DeltaLearning(QEF[i],QEF_par[i],QEF_systems[i]))
        if Qout > 1:
          print("Delta learning applied")
      if QEF[i] == "h_RMSD":
        from umbrellaRMSDconstraint import UmbrellaConstraint
        Qdistance = 2
        constraints.append(UmbrellaConstraint(all_species,QEF_par[i][1],QEF_par[i][0],Qslow))
  #SET CONSTRAINTS
  if len(constraints) > 0:
    all_species.set_constraint(constraints)
    if Qout > 1:
      print("CONSTRAINTS: "+str(constraints), flush = True)
  
  #SET CALCULATOR
  if Qout > 1:
    print("Setting calculator.", flush = True)
  def call_calculator():
    from calculator import calculator
    if Qconstraints == 4 and len(species) == 2:
      species[0].calc = calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations, Qcharge, Qmultiplicity, Qout, species[0],Qmixer_damping,Qcutoff)
      species[1].calc = calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations, Qcharge, Qmultiplicity, Qout, species[1],Qmixer_damping,Qcutoff)
    else:
      all_species.calc = calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations, Qcharge, Qmultiplicity, Qout, all_species,Qmixer_damping,Qcutoff)
  if Qout > 1:
    print("Calling calculator.", flush = True)
    print("  Qconstraints: "+str(Qconstraints), flush = True)
  call_calculator()
  if Qout > 1:
    print(all_species)
    #print(all_species.get_positions())
    #print(all_species.get_potential_energy())
    #print(all_species.get_forces())
    #exit()
  
  #THERMOSTAT
  if Qout > 1:
    print("Setting thermostat.", flush = True)
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
    if Qconstraints == 4 and len(species) == 2:
      dyn1 = Langevin(species[0], Qdt * units.fs, temperature_K = Qtemp, friction = Qthermostat_L / units.fs, fixcm = Qfixcm, rng = Qrng)
      dyn2 = Langevin(species[1], Qdt * units.fs, temperature_K = Qtemp, friction = Qthermostat_L / units.fs, fixcm = Qfixcm, rng = Qrng)
      dyn1.attach(call_rmsd_constrain1, interval = 1)
      dyn2.attach(call_rmsd_constrain2, interval = 1)
    else:
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
  elif Qthermostat == "A":
    from ase import units
    from ase.md.andersen import Andersen 
    dyn = Andersen(all_species, Qdt * units.fs, temperature_K = Qtemp, andersen_prob = Qthermostat_A * units.fs, rng = Qrng)
  elif Qthermostat == "OPT":
    from ase.optimize import BFGS
    dyn = BFGS(all_species)
  else:
    print("Some weird thermostat.")
    exit()
  if Qout > 1:
    print("Finished setting Qthermostat: "+str(Qthermostat), flush = True)
 
  #DUMPING
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
  if Qconstraints == 4 and len(species) == 2:
    cluster_dic1 = {}
    cluster_dic2 = {}
  else:
    if "cluster_dic" not in globals(): #current_step == 0:
      cluster_dic = {}
  def save(fail = False):
    global current_time,current_step
    if Qconstraints == 4 and len(species) == 2:
      global cluster_dic1,cluster_dic2
      toupdate1,current_time,current_step = print_properties(species = species[0], timestep = 0, interval = 0, Qconstraints = Qconstraints, Qdistance = Qdistance, split = Qlenfirst, fail = fail, QINFOfile_basename = QINFOfile_basename, QINFOcluster_type = QINFOcluster_type, QINFOcomponents = QINFOcomponents, QINFOcomponent_ratio = QINFOcomponent_ratio, heavyatoms = Qheavyatoms)
      toupdate2,current_time,current_step = print_properties(species = species[1], timestep = Qdt, interval = Qdump, Qconstraints = Qconstraints, Qdistance = Qdistance, split = Qlenfirst, fail = fail, QINFOfile_basename = QINFOfile_basename, QINFOcluster_type = QINFOcluster_type, QINFOcomponents = QINFOcomponents, QINFOcomponent_ratio = QINFOcomponent_ratio, heavyatoms = Qheavyatoms)
    else:
      global cluster_dic
      toupdate,current_time,current_step = print_properties(species = all_species, timestep = Qdt, interval = Qdump, Qconstraints = Qconstraints, Qdistance = Qdistance, split = Qlenfirst, fail = fail, QINFOfile_basename = QINFOfile_basename, QINFOcluster_type = QINFOcluster_type, QINFOcomponents = QINFOcomponents, QINFOcomponent_ratio = QINFOcomponent_ratio, heavyatoms =Qheavyatoms)
    if Qsavepickle == 1:
      if Qconstraints == 4 and len(species) == 2:
        toupdate1.update({("log","method"):[" ".join(argv[1:])],("log","program"):["Python"]})
        toupdate2.update({("log","method"):[" ".join(argv[1:])],("log","program"):["Python"]})
      else:
        toupdate.update({("log","method"):[" ".join(argv[1:])],("log","program"):["Python"]})
        if Qconstraints == 3:
          toupdate.update({("log","k_bias"):[min(current_step/max(Qslow,0.0000001),1)*Qk_bias],("log","harm_distance"):[Qharm]})
      if Qconstraints == 4 and len(species) == 2:
        cluster_dic1 = mergeDictionary(cluster_dic1, toupdate1)
        cluster_dic2 = mergeDictionary(cluster_dic2, toupdate2)
      else:
        cluster_dic = mergeDictionary(cluster_dic, toupdate)
  if Qdump == 0:
    save()
  else:
    if Qconstraints == 4 and len(species) == 2:
      #dyn1.attach(save, interval = 1)
      dyn2.attach(save, interval = 1)
    else:
      dyn.attach(save, interval = Qdump)
  if Qsave >= 0:
    dyn.attach(savepickle, interval = Qsave)
  if Qcalculator == "PhysNet": 
    #from calculator import calculator
    def updatephysnet():
      call_calculator()
      #all_species.calc = calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations, Qcharge, Qout, all_species)
    dyn.attach(updatephysnet, interval = 1)
  if Qcalculator == "XTB":
    if Qcalculator_input == "GFNFF":
      def removetopologyfiles():
        import os
        os.system("ls")
        print(os.system("rm -rf gfnff_adjacency"))
        os.system("rm -rf gfnff_topo")
        os.system("ls")
        call_calculator()
      dyn.attach(removetopologyfiles, interval = 1)
  if Qnn_EFD == 1:
    from src.JKelectrostatics import compute_energies_forces
    from src.JKdispersions import compute_d4_energy_forces as compute_dispersions
    #from src.JKdispersions import compute_d3bj_energy_forces as compute_dispersions
    def updateElDisp():
      internal_E = species.get_potential_energy()
      internal_F = species.get_forces()
      from tblite.ase import TBLite
      CS = species.constraints
      del species.constraints
      tmpatoms = species.copy()
      species.set_constraint(CS)
      tmpatoms.calc = TBLite(method="GFN1-xTB", cache_api=True, charge=float(0), verbosity = 0, max_iterations = 300, accuracy = 1.0)
      Q_charges = array([tmpatoms.get_charges()]) #.transpose()
      #Q_charges = spk_calc.model_results['partial_charges'].detach().numpy()
      electrostatics_E, electrostatics_F = compute_energies_forces(species.get_positions(), Q_charges)
      dispersions_E, dispersions_F = compute_dispersions(species.get_positions(), symbols = array(species.get_chemical_symbols()), totalcharge = 0)
      print(f"JKML(SchNetPack): {0.0367493 * internal_E} {electrostatics_E} {dispersions_E}")
      species.set_potential_energy(0.0367493 * internal_E + electrostatics_E + dispersions_E)
      species.set_forces(0.0367493 * internal_F + electrostatics_F + dispersions_F)
  stepsmade = 0
  def stepsmadeadd():
    global stepsmade
    stepsmade += 1
  if Qconstraints == 4 and len(species) == 2:
    dyn1.attach(stepsmadeadd, interval = 1)
    dyn2.attach(stepsmadeadd, interval = 1)
  else:
    dyn.attach(stepsmadeadd, interval = 1)

  #SIMULATION
  #print(all_species)
  if Qout > 1:
    print("Starting simulation.", flush = True)
  if 'current_step' in globals():
    steps_before_sim = current_step
  else:
    steps_before_sim = 0
  sim_errors = 0
  sim_tot_errors = 0
  sim_last_error = -1
  while stepsmade < Qns:
    #if Qout > 1:
    #  dyn.run(1)
    #  print("1 done")
    try:
      if Qconstraints == 4 and len(species) == 2:
        dyn1.run(1)
        dyn2.run(1)
      elif Qthermostat == "OPT":
        dyn.run(fmax=0.05)
      else:
        dyn.run(Qns - stepsmade)
      if Qdump == 0:
        current_time = current_time + Qdt*Qns
        current_step = current_step + Qns
      else:
        current_time = current_time - Qdt*Qdump
        current_step = current_step - Qdump
    except Exception as e:
      print("Something got screwed up within the dyn.run(Qns). Small adjustment to the structure!!!", flush = True)
      print(str(e), flush = True)
      from numpy import random
      positions = all_species.get_positions()
      noise = random.normal(scale=0.01, size=positions.shape)
      all_species.set_positions(positions + noise)
      call_calculator()
      print(all_species)
      sim_tot_errors += 1
      if current_step == sim_last_error:
        sim_errors += 1
      else:
        sim_errors = 1
        sim_last_error = current_step
      if sim_errors == 4 or sim_tot_errors > 20:
        if sim_tot_errors > 20:
          print("Too many errors in the simulation. Exiting.")
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

  if Qout > 1:
    print("Simulation round done.")

if Qsavepickle == 1:
  if Qout > 1:
    print("Done and now just saving pickle.")
  savepickle()
  #save(Qfolder+"/../sim"+Qfolder.split("/")[-1]+".pkl")
 # if "cluster_dic" in globals() or "cluster_dic1" in globals() or "cluster_dic2" in globals():
 #   from pandas import DataFrame
 #   if Qconstraints == 4 and len(species) == 2:
 #     for key in cluster_dic1:
 #       cluster_dic2[key] = cluster_dic2[key][::-1]
 #     cluster_dic = mergeDictionary(cluster_dic1, cluster_dic2) 
 #   else:
 #     for key in cluster_dic:
 #       cluster_dic[key] = cluster_dic[key][::-1]
 #   clusters_df = DataFrame(cluster_dic) #, index = range(len(cluster_dic)))
 #   try:
 #     clusters_df.to_pickle(Qfolder+"/../sim"+Qfolder.split("/")[-1]+".pkl")
 #     print("The sim"+Qfolder.split("/")[-1]+".pkl has been hopefully created.")
 #   except:
 #     print("Something got fucked up.")

print("Done.")
