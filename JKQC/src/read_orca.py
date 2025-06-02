###############
### ORCA ######
###############
def read_orca_init(Qforces = 0, Qanharm = 0, Qdisp_forces = 0):
  from re import compile
  global PATTERN_ORCA_out_method,PATTERN_ORCA_out_vibrational_frequencies,PATTERN_ORCA_out_mulliken_charges
  PATTERN_ORCA_out_method = compile(rb"\n\|.*>.*!.*")
  PATTERN_ORCA_out_vibrational_frequencies = compile(rb'VIBRATIONAL FREQUENCIES.*\n.*-{2,}\n.*\nScaling factor for frequencies.*\n.*\n((?:\s*\d+:\s*[-]?\d+\.\d+\s*cm\*\*-1.*\n)+)\n')
  PATTERN_ORCA_out_mulliken_charges = compile(rb'MULLIKEN ATOMIC CHARGES\s*-{2,}\n((?:\s+\d+\s+\w+\s*:\s*[-+]?\d*\.\d+\n)+)Sum') 
  if Qanharm == 1:
    global PATTERN_ORCA_anharm
    PATTERN_ORCA_out_vibrational_frequencies = compile(rb'Fundamental transitions.*\n.*-{2,}\n.*Mode.*\n.*-{2,}.*\n((?:\s*\d+\s+[-]?\d+.\d+\s+[-]?\d+.\d+\s+[-]?\d+\.\d+\s*\n)+).*-{2,}.*\n.*\n')
    PATTERN_ORCA_anharm = compile(rb'Anharmonic constants.*\n.*-{2,}\n.*r.*\n.*-{2,}\n((?:\s*\d+\s*\d+\s*[-+]?\d+\.\d+\s.*\n)+).*-{2,}\n')
  if Qforces == 1:
    global PATTERN_ORCA_out_forces
    PATTERN_ORCA_out_forces = compile(rb'CARTESIAN GRADIENT.*\w*\n*-{2,}\n*((?:\s+\d+\s+\w+\s*:\s*[-+]?\d*\.\d+\s*[-+]?\d*\.\d+\s*[-+]?\d*\.\d+\n)+)\s*\n[D,N]') 
  if Qdisp_forces == 1:
    global PATTERN_ORCA_out_dispersion_forces
    PATTERN_ORCA_out_dispersion_forces = compile(rb'DISPERSION GRADIENT\s*\n-{2,}\n((?:\s+\d+\s+\w+\s*:\s*[-+]?\d*\.\d+(?:[eE][-+]?\d+)?\s+[-+]?\d*\.\d+(?:[eE][-+]?\d+)?\s+[-+]?\d*\.\d+(?:[eE][-+]?\d+)?\n)+)')

def find_line(bytes_string,take_first = 0,idx = 0):
  mm.seek(idx)
  if take_first:
    start_index = mm.find(bytes_string)
  else:
    start_index = mm.rfind(bytes_string)
  if start_index != -1:
    mm.seek(start_index)
    line = mm.readline()
    return line.decode("utf-8").strip(), start_index
  else:
    return None, idx

def read_orca(mmm, orcaextname, Qforces = 0, Qanharm = 0, Qdisp_electronic_energy = 0, Qdisp_forces = 0):
  from numpy import array,sqrt
  missing = float("nan")

  columns = ["program","method","time","termination","charge","multiplicity","NAtoms","rotational_constants","rotational_constant","sp_electronic_energy","electronic_energy","mulliken_charges","dipole_moment","dipole_moments","polarizability","vibrational_frequencies","temperature","pressure","moments_of_inertia","rotational_symmetry_number","zero_point_correction","energy_thermal_correction","enthalpy_thermal_correction","gibbs_free_energy_thermal_correction","zero_point_energy","internal_energy","enthalpy_energy","gibbs_free_energy","entropy","scf_energy"]

  global mm
  mm = mmm

  #TIME
  try:
    line,idx = find_line(rb"TOTAL RUN TIME", 0, 0)
    out_time = float(line.split()[3])*24*60+float(line.split()[5])*60+float(line.split()[7])+float(line.split()[9])/60+float(line.split()[11])/6000
  except:
    out_time = missing

  #ELECTRONIC ENERGY
  try:
    line,idx = find_line(rb'FINAL SINGLE POINT ENERGY', 1, 0)
    out_sp_electronic_energy = float(line.split()[4])
    line,idx = find_line(rb'FINAL SINGLE POINT ENERGY', 0, 0)
    out_electronic_energy = float(line.split()[4])
  except:
    out_sp_electronic_energy = missing 
    out_electronic_energy = missing 

  #VIBRATIONAL FREQUENCIES
  try:
    lines = PATTERN_ORCA_out_vibrational_frequencies.findall(mm)[-1].decode("utf-8").split("\n")
    if Qanharm == 1:
      out_vibrational_frequencies = [float(line.split()[1]) for line in lines[0:-1]]
      try:
        lines = PATTERN_ORCA_anharm.findall(mm)[-1].decode("utf-8").split("\n")
        corrections = array([array([float(number) for number in line.split()]) for line in lines[0:-1]])
        out_anharmonicties = []
        for i in range(len(out_vibrational_frequencies)):
          corr = 0
          for j in range(len(corrections)):
            if corrections[j,0] == i and corrections[j,1] == i:
              corr = corr + 2*corrections[j,2]
            elif corrections[j,0] == i:
              corr = corr + 0.5*corrections[j,2]
            elif corrections[j,1] == i:
              corr = corr + 0.5*corrections[j,2]
          out_anharmonicties.append(corr+out_vibrational_frequencies[i])  
        out_anharm = out_anharmonicties
      except:
        out_anharm = missing
    else:
      out_vibrational_frequencies = [float(line.split()[1]) for line in lines[6:-1]]
  except:
    out_vibrational_frequencies = missing
    if Qanharm == 1:
      out_anharm = missing
 
  #METHOD
  try:
    out_method = "_".join(PATTERN_ORCA_out_method.findall(mm)[-1].decode("utf-8").split("> ")[1].split())
  except:
    out_method = missing

  #MULLIKEN ATOMIC CHARGES 
  try: 
    out_mulliken_charges = [float(line.split()[3]) for line in PATTERN_ORCA_out_mulliken_charges.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]     
  except:
    out_mulliken_charges = missing
     
  #TERMINATION
  try:
    #out_termination = len(PATTERN_ORCA_out_termination.findall(mm))
    if mm.find(rb"ORCA TERMINATED NORMALLY") > 0:
      out_termination = 1
    else:
      out_termination = 0
  except:
    out_termination = missing

  #PROGRAM VERSION
  try:
    line, idx = find_line(rb'Program Version')
    out_program = "ORCA_" + line.split()[2]
  except:
    out_program = missing

  #Number of atoms
  try:
    line, idx = find_line(rb"Number of atoms", 1, idx)
    out_NAtoms = int(line.split()[4])
  except:
    out_NAtoms = missing

  #Total Charge
  try:
    line, idx = find_line(rb"Total Charge", 1, idx)
    out_charge = int(line.split()[4])
  except:
    out_charge = missing

  #Multiplicity
  try:
    line, idx = find_line(rb' Multiplicity', 1, idx)
    out_multiplicity = int(line.split()[3])
  except:
    out_multiplicity = missing

  #DIPOLE MOMENTS
  try:
    line, idx = find_line(rb'Total Dipole Moment', 0)
    out_dipole_moments = [float(line.split()[4])/0.393456, float(line.split()[5])/0.393456, float(line.split()[6])/0.393456]
  except:
    out_dipole_moments = missing

  #DIPOLE MOMENT
  try:
    line, idx = find_line(rb'Magnitude (Debye)', 1, idx)
    out_dipole_moment = float(line.split()[3])
  except:
    out_dipole_moment = missing

  #ROTATIONAL CONSTANT(S)
  try:
    line, idx = find_line(rb'Rotational constants in MHz', 1, idx)
    out_rotational_constants = [float(line.split()[5])/1000,float(line.split()[6])/1000,float(line.split()[7])/1000]
  except:
    out_rotational_constants = [missing]
  out_rotational_constant = sqrt(sum(array(out_rotational_constants)**2))

  #TEMPERATURE
  try:
    line, idx = find_line(rb'Temperature         ...', 1, idx)
    out_temperature = float(line.split()[2])
  except:
    out_temperature = missing

  #PRESSURE
  try:
    line, idx = find_line(rb'Pressure            ...', 1, idx)
    out_pressure = float(line.split()[2])
  except:
    out_pressure = missing    

  #ZERO POINT ENERGY
  try:
    line, idx = find_line(rb'Zero point energy                ...', 1, idx)
    out_zero_point_correction = float(line.split()[4])
  except:
    out_zero_point_correction = missing

  #TOTAL THERMAL ENERGY
  try:
    line, idx = find_line(rb'Total thermal energy', 1, idx)
    out_internal_energy = float(line.split()[3])
  except:
    out_internal_energy = missing

  #SYMMETRY NUMBER
  try:
    line, idx = find_line(rb'Symmetry Number', 1, idx)
    out_rotational_symmetry_number = float(line.split()[2])
  except:
    out_rotational_symmetry_number = missing

  #FINAL ENTROPY
  try:
    line, idx = find_line(rb'Final entropy term', 1, idx)
    out_entropy = float(line.split()[4])*1000*627.503/out_temperature
  except:
    out_entropy = missing

  #GIBBS FREE ENERGY
  try:
    line, idx = find_line(rb'Final Gibbs free energy', 1, idx)
    out_gibbs_free_energy = float(line.split()[5])
  except:
    out_gibbs_free_energy = missing

  #ENTHALPY
  try:
    line, idx = find_line(rb'Total enthalpy', 1, 0)
    out_enthalpy_energy = float(line.split()[3])
  except:
    out_enthalpy_energy = missing

  #SCF ENERGY
  try:
    line, idx = find_line(rb'Total Energy', 1, 0)
    out_scf_energy = float(line.split()[3])
  except:
    out_scf_energy = missing

  #FORCES
  if Qforces == 1:
    try:
      out_forces = [-array([float(line.split()[3]),float(line.split()[4]),float(line.split()[5])])/0.529177 for line in PATTERN_ORCA_out_forces.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]
    except:
      out_forces = missing

  #ELECTRONIC ENERGY DISPERSION CORRECTION
  if Qdisp_electronic_energy == 1:
    try:
      line,idx = find_line(rb'Dispersion correction', 1, 0)
      out_disp_electronic_energy = float(line.split()[-1])
    except:
      out_disp_electronic_energy = missing
  
  #FORCES DISPERSION CORRECTION
  if Qdisp_forces == 1:
    try:
      out_disp_forces = [-array([float(line.split()[3]),float(line.split()[4]),float(line.split()[5])]) / 0.529177 for line in PATTERN_ORCA_out_dispersion_forces.findall(mm)[-1].decode("utf-8").split("\n")[1:-1]]
    except:
      out_disp_forces = missing

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
  mm.close()
  all_locals = locals()
  dic = {(orcaextname,column):[all_locals.get("out_"+column)] for column in columns}
  if Qforces == 1:
    dic.update({("extra","forces"):[out_forces]})
  if Qanharm == 1:
    dic.update({("extra","anharm"):[out_anharm]})
  if Qdisp_electronic_energy == 1:
    dic.update({("extra", "dispersion_electronic_energy"): [out_disp_electronic_energy]})
  if Qdisp_forces == 1:
    dic.update({("extra", "dispersion_forces"): [out_disp_forces]})

  #ANHARMONIC FREQS
  return dic

#else:
#  dic = {(orcaextname,column):[missing] for column in columns}
#  if Qforces == 1:
#    dic.update({("extra","forces"):[missing]})
#  return dic

