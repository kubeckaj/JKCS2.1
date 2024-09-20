###############
###  G16  #####
###############
def read_g16_init(Qforces = 0, Qanharm = 0):
  from re import compile
  global PATTERN_G16_out_program,PATTERN_G16_out_method,PATTERN_G16_out_mulliken_charges,PATTERN_G16_out_dipole_moment,PATTERN_G16_out_esp_charges
  PATTERN_G16_out_program = compile(rb" Gaussian.*, Revision.*,")
  PATTERN_G16_out_method = compile(rb"\n #.*")
  PATTERN_G16_out_mulliken_charges = compile(rb' Mulliken charges:.*\n.*\n((?:\s+\d+\s+\w+\s*[-+]?\d*\.\d+.*\n)+).*Sum') 
  PATTERN_G16_out_esp_charges = compile(rb' ESP charges:.*\n.*\n((?:\s+\d+\s+\w+\s*[-+]?\d*\.\d+.*\n)+).*Sum') 
  PATTERN_G16_out_dipole_moment = compile(rb'Dipole moment \(field-independent basis, Debye\):.*\n.*\n')
  if Qforces == 1:
    global PATTERN_G16_out_forces
    PATTERN_G16_out_forces = compile(rb'Center     Atomic                   Forces.*\n.*Number.*\n\s*[-]*\s*\n((?:\s+\d+\s+.*\n)+)\s*[-]*\s*\n')
  if Qanharm == 1:
    global PATTERN_G16_out_anharm
    PATTERN_G16_out_anharm = compile(rb'Fundamental Bands.*\n.*\n.*Mode.*\n((?:\s*\w*\s*\d+\(.\)\s+\w+\s+[-+]?\d+\.\d+.*\n)+).*\n.*Overtones\n')

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

def find_lines(bytes_string,take_first = 0,idx = 0):
  lines = []
  line = ""
  while line is not None:
    line, idx = find_line(bytes_string,take_first,idx+1)
    if line is not None:
      lines.append(line)
  return lines, idx

def read_g16(mmm, Qforces = 0, Qanharm = 0):
  from numpy import array
  from math import sqrt
  missing = float("nan")
  
  columns = ["program","method","time","termination","charge","multiplicity","NAtoms","rotational_constants","rotational_constant","sp_electronic_energy","electronic_energy","mulliken_charges","dipole_moment","dipole_moments","polarizability","vibrational_frequencies","temperature","pressure","moments_of_inertia","rotational_symmetry_number","zero_point_correction","energy_thermal_correction","enthalpy_thermal_correction","gibbs_free_energy_thermal_correction","zero_point_energy","internal_energy","enthalpy_energy","gibbs_free_energy","entropy"]
  columns.append("esp_charges")
  columns.append("moments_of_inertia")

  global mm
  mm = mmm
      
  #TIME
  try:
    lines,idx = find_lines(rb"Elapsed time", 0, 0)
    out_time = sum([float(line.split()[2])*24*60+float(line.split()[4])*60+float(line.split()[6])+float(line.split()[8])/60 for line in lines])
  except:
    out_time = missing

  #TERMINATION
  try:
    if mm.find(rb"Normal termination") > 0:
      out_termination = 1
    else:
      out_termination = 0
  except:
    out_termination = missing

  #MULLIKEN ATOMIC CHARGES 
  try: 
    out_mulliken_charges = [float(line.split()[2]) for line in PATTERN_G16_out_mulliken_charges.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]     
  except:
    out_mulliken_charges = missing

  #ESP CHARGES
  try:
    out_esp_charges = [float(line.split()[2]) for line in PATTERN_G16_out_esp_charges.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]
  except:
    out_esp_charges = missing
     
  #PROGRAM VERSION
  try:
    line = PATTERN_G16_out_program.findall(mm)[0].decode("utf-8").split("\n")[-1]
    out_program = "G" + str(line.split(",")[0].split()[1]) + "_" + str(line.split(",")[1].split()[1])
  except:
    out_program = missing

  #METHOD
  try:
    line = PATTERN_G16_out_method.findall(mm)[-1].decode("utf-8").split("\n")[-1]
    out_method = "_".join(line.lower().split())
  except:
    out_method = missing
  
  #DIPOLE MOMENTS
  try:
    line = PATTERN_G16_out_dipole_moment.findall(mm)[-1].decode("utf-8").split("\n")[-2]
    out_dipole_moments = [float(line.split()[1]), float(line.split()[3]), float(line.split()[5])]
    out_dipole_moment = float(line.split()[7])
  except:
    out_dipole_moments = [missing]
    out_dipole_moment = missing

  #Total Charge & Multiplicity
  try:
    line, idx = find_line(rb" Charge = ", 1, 0)
    out_charge = int(line.split()[2])
    out_multiplicity = int(line.split()[5])
  except:
    out_charge = missing
    out_multiplicity = missing

  #Number of atoms
  try:
    line, idx = find_line(rb"NAtoms=", 1, idx)
    out_NAtoms = int(line.split()[1])
  except:
    out_NAtoms = missing

  #ROTATIONAL CONSTANT(S)
  try:
    line, idx = find_line(rb'Rotational constants ', 0, idx)
    out_rotational_constants = [float(line.split()[3]),float(line.split()[4]),float(line.split()[5])]
    out_rotational_constant = sqrt(sum(array(out_rotational_constants)**2))
  except:
    out_rotational_constants = [missing]
    out_rotational_constant = missing

  #POLARIZABILITY
  try:
    line, idx = find_line(rb'Exact polarizability', 0, 0)
    pol = [float(i) for i in line.split()[2:]]
    from numpy import matrix
    from numpy.linalg import eigh
    pol_mat = matrix([[pol[0],pol[1],pol[3]],[pol[1],pol[2],pol[4]],[pol[3],pol[4],pol[5]]])
    out_polarizability = 0.14818471147*sum(eigh(pol_mat)[0])/3.
  except:
    out_polarizability = missing

  #ELECTRONIC ENERGY
  try:
    line,idx = find_line(rb'SCF Done', 1, 0)
    out_sp_electronic_energy = float(line.split()[4])
    line,idx = find_line(rb'SCF Done', 0, 0)
    out_electronic_energy = float(line.split()[4])
  except:
    out_sp_electronic_energy = missing 
    out_electronic_energy = missing 

  #VIBRATIONAL FREQUENCIES
  try:
    lines,idx = find_lines(rb" Frequencies -- ", 1, 0)
    if Qanharm == 1:
      lines = lines[0:int(len(lines)/3)]
    out_vibrational_frequencies = [float(element) for line in lines for element in line.split()[2:]]
    if len(out_vibrational_frequencies) == 0:
      out_vibrational_frequencies = [missing]
  except:
    out_vibrational_frequencies = [missing]
 
  #TEMPERATURE & PRESSURE
  try:
    line, idx = find_line(rb' Temperature ', 1, idx)
    out_temperature = float(line.split()[1])
    out_pressure = float(line.split()[4])
  except:
    out_temperature = missing
    out_pressure = missing    

  #MOMENT OF INERTIA
  try:
    line, idx = find_line(rb'Eigenvalues -- ', 1, idx)
    out_moments_of_inertia = [float(i) for i in line.split()[2:]]
  except:
    out_moments_of_inertia = missing

  #SYMMETRY NUMBER
  try:
    line, idx = find_line(rb'Rotational symmetry number', 1, idx)
    out_rotational_symmetry_number = float(line.split()[3])
  except:
    out_rotational_symmetry_number = missing

  #ZERO POINT ENERGY CORRECTION
  try:
    line, idx = find_line(rb'Zero-point correction=', 1, idx)
    out_zero_point_correction = float(line.split()[2])
  except:
    out_zero_point_correction = missing

  #TOTAL THERMAL ENERGY CORRECTION
  try:
    line, idx = find_line(rb'Thermal correction to Energy', 1, idx)
    out_energy_thermal_correction = float(line.split()[4])
  except:
    out_energy_thermal_correction = missing

  #TOTAL THERMAL CORRECTION TO ENHALPY
  try:
    line, idx = find_line(rb'Thermal correction to Enthalpy', 1, idx)
    out_enthalpy_thermal_correction = float(line.split()[4])
  except:
    out_enthalpy_thermal_correction = missing

  #TOTAL THERMAL CORRECTION TO GIBBS FREE ENERGY
  try:
    line, idx = find_line(rb'Thermal correction to Gibbs Free Energy', 1, idx)
    out_gibbs_free_energy_thermal_correction = float(line.split()[6])
  except:
    out_gibbs_free_energy_thermal_correction = missing

  #ZERO POINT ENERGY
  try:
    line, idx = find_line(rb'Sum of electronic and zero-point Energies', 1, idx)
    out_zero_point_energy= float(line.split()[6])
  except:
    out_zero_point_energy = missing

  #TOTAL THERMAL ENERGY
  try:
    line, idx = find_line(rb'Sum of electronic and thermal Energies', 0, 0)
    out_internal_energy = float(line.split()[6])
  except:
    out_internal_energy = missing

  #ENTHALPY
  try:
    line, idx = find_line(rb'Sum of electronic and thermal Enthalpies', 0, idx)
    out_enthalpy_energy = float(line.split()[6])
  except:
    out_enthalpy_energy = missing

  #GIBBS FREE ENERGY
  try:
    line, idx = find_line(rb'Sum of electronic and thermal Free Energies', 1, idx)
    out_gibbs_free_energy = float(line.split()[7])
  except:
    out_gibbs_free_energy = missing

  #FINAL ENTROPY
  try:
    line, idx = find_line(rb'Total', 1, idx)
    out_entropy = float(line.split()[3])
  except:
    out_entropy = missing

  #FORCES
  if Qforces == 1:
    try:
      lines = PATTERN_G16_out_forces.findall(mm)[-1].decode("utf-8").split("\n")[:-2]
      out_forces = [array([float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])/0.529177 for line in PATTERN_G16_out_forces.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]
    except:
      out_forces = missing

  #ANHARMONIC FREQS
  if Qanharm == 1:
    try:
      lines = PATTERN_G16_out_anharm.findall(mm)[-1].decode("utf-8").split("\n")[:-1]
      out_anharm = [float(line.split()[-4]) for line in lines][-1::-1]
    except:
      out_anharm = missing

  #SAVE
  mm.close()
  all_locals = locals()
  dic = {("log",column):[all_locals.get("out_"+column)] for column in columns}
  if Qforces == 1:
    dic.update({("extra","forces"):[out_forces]})
  if Qanharm == 1:
    dic.update({("extra","anharm"):[out_anharm]})
  return dic

#else:
#  dic = {(orcaextname,column):[missing] for column in columns}
#  if Qforces == 1:
#    dic.update({("extra","forces"):[missing]})
#  if Qanharm == 1:
#    dic.update({("extra","anharm"):[missing]})
#  return dic

