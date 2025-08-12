def print_properties(species , timestep = 1, interval = 1, Qconstraints = 0, Qdistance = 0, split = None, fail = False, QINFOfile_basename = "str", QINFOcluster_type = "", QINFOcomponents = [], QINFOcomponent_ratio = []):
  from ase.md.velocitydistribution import Stationary
  from ase.md.velocitydistribution import ZeroRotation
  from ase import units
  #from ase.io import read, write
  global current_time, current_step
 
  ### """Function to print the potential, kinetic and total energy"""
  if not fail:
    epot = species.get_potential_energy()/0.043364115308770496 #kcal/mol
    ekin = species.get_kinetic_energy()/0.043364115308770496    #kca/lmol
  else:
    epot = float("nan")
    ekin = float("nan")

  ### CONSTRAINTS
  CS = species.constraints
  del species.constraints
  species_copy = species.copy()
  species_copy_for_saving = species.copy()
  species.set_constraint(CS)

  ### DISTANCE
  if Qconstraints == 3 or Qdistance == 1:
    from numpy import sqrt, sum, ones, array
    if 1==0:
      mask1 = ones(len(species_copy[0:split]), dtype=bool)
      mask2 = ones(len(species_copy[split:]), dtype=bool)
    else:
      mask1 = array(species_copy[0:split].symbols) != 'H'
      mask2 = array(species_copy[split:].symbols) != 'H'
    dist_n = sqrt(sum(((species_copy[0:split][mask1].get_center_of_mass()-species_copy[split:][mask2].get_center_of_mass())**2)))
    spread_a = species_copy[0:split][mask1].get_all_distances().max()
    spread_b = species_copy[split:][mask2].get_all_distances().max()
  elif Qdistance == 2:
    from umbrellaRMSDconstraint import RMSD3
    dist_n, spread_a, spread_b = RMSD3(species_copy)
  else:
    dist_n = 0.0
    spread_a = 0.0
    spread_b = 0.0

  ### """ TEMPERATURES """
  T_temp = species_copy.get_temperature()
  Stationary(species_copy, False)
  T_com = species_copy.get_temperature()
  #TODO
  ZeroRotation(species_copy, False)
  T_rotate = species_copy.get_temperature()
  T_tr =  len(species_copy)*(T_temp-T_com)
  T_rot = len(species_copy)*(T_com-T_rotate)
  if len(species_copy) > 2:
    T_vib = len(species_copy)/(len(species_copy)-2)*T_rotate
  elif len(species_copy) == 2:
    T_vib = 3*len(species_copy)/(3*len(species_copy)-5)*T_rotate
  else:
    T_vib = 0
  if Qconstraints > 0 or Qdistance > 0:
    if current_step == 0:
      print('      STEP_[-] TIME_[fs] | Et[kcal/mol] Ep[kcal/mol] Ek[kcal/mol] | T_[K] Tt[K] Tr[K] Tv[K] | COMd_[A] MaxA_[A] MaxB_[A]', flush=True)
    print('JKMD: %-*i %-*.1f | %-*.3f %-*.3f %-*.3f | %-*.0f %-*.0f %-*.0f %-*.0f | %-8.4f %-8.2f %-8.2f' % (8,current_step, 9,current_time, 12,epot + ekin, 12,epot, 12,ekin, 5,T_temp, 5,T_tr, 5,T_rot, 5,T_vib, dist_n, spread_a, spread_b), flush=True)
  else:
    if current_step == 0:
      print('      STEP_[-] TIME_[fs] | Et[kcal/mol] Ep[kcal/mol] Ek[kcal/mol] | T_[K] Tt[K] Tr[K] Tv[K]', flush=True)
    print('JKMD: %-*i %-*.1f | %-*.3f %-*.3f %-*.3f | %-*.0f %-*.0f %-*.0f %-*.0f' % (8,current_step, 9,current_time, 12,epot + ekin, 12,epot, 12,ekin, 5,T_temp, 5,T_tr, 5,T_rot, 5,T_vib), flush=True)


  from os import path
  folder_path = path.abspath("./test")[::-1].split("/",1)[1][::-1]+"/"
  dic = {("info","folder_path"):[folder_path]}
  dic.update({("info","file_basename"):[QINFOfile_basename+"-"+str(current_step)]})
  if QINFOfile_basename != "str":
    dic.update({("info","cluster_type"):[QINFOcluster_type]})
    dic.update({("info","components"):[QINFOcomponents]})
    dic.update({("info","component_ratio"):[QINFOcomponent_ratio]})
  dic.update({("xyz","structure"):[species_copy_for_saving]})
  dic.update({("log","md_time"):[current_time]})
  dic.update({("log","md_step"):[current_step]})
  dic.update({("log","electronic_energy"):[epot/627.503]})
  dic.update({("log","kinetic_energy"):[ekin/627.503]})
  dic.update({("log","total_energy"):[(epot+ekin)/627.503]})
  dic.update({("log","temperature"):[T_temp]})
  dic.update({("log","translational_temperature"):[T_tr]})
  dic.update({("log","rotational_temperature"):[T_rot]})
  dic.update({("log","vibrational_temperature"):[T_vib]})
  if Qconstraints != 0 or Qdistance == 1: 
    dic.update({("log","COM_distance"):[dist_n]})
    dic.update({("log","maxA_distance"):[spread_a]})
    dic.update({("log","maxB_distance"):[spread_b]})
  #write("traj.xyz", a, append = True)
  
  current_time = current_time + interval*timestep
  current_step = current_step + interval
  
  return dic,current_time,current_step

def init(Qtime,Qstep):
  global current_step, current_time
  current_time = Qtime
  current_step = Qstep
