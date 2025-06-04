##############
### XTB ######
##############
def read_xtb_init():
  from re import compile
  global PATTERN_XTB_out_dipole_moment,PATTERN_XTB_out_dipole_moment2,PATTERN_XTB_out_vibrational_frequencies,PATTERN_XTB_out_mulliken_charges,PATTERN_XTB_out_mulliken_charges2
  PATTERN_XTB_out_dipole_moment = compile(rb"dipole moment from electron density .*\n.*\n.*total .*Debye.*\n")
  PATTERN_XTB_out_dipole_moment2 = compile(rb"molecular dipole:.*\n.*tot .*Debye.*\n.*q only:.*\n.*full.*\n")
  PATTERN_XTB_out_vibrational_frequencies = compile(rb'projected vibrational frequencies.*\n((?:.*eigval\s*:.*\n)+).*reduced masses.*\n')
  PATTERN_XTB_out_mulliken_charges = compile(rb'Mulliken.*\n((?:\s+\d+\w+\s+.*\n)+)\n') 
  PATTERN_XTB_out_mulliken_charges2 = compile(rb'Z\s+covCN\s+q\s+C6AA.*\n((?:\s*\d*\s+\d+\s+\w\s+.*\n)+)\n')

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

def read_xtb(mmm):
  global mm
  mm = mmm
  
  from numpy import array,sqrt
  from io import open
  missing = float("nan")
  
  columns = ["program","method","time","NAtoms","electronic_energy","mulliken_charges","dipole_moment","dipole_moments","vibrational_frequencies","enthalpy_energy","gibbs_free_energy","entropy","zero_point_correction","zero_point_energy","rotational_symmetry_number"]

  #PROGRAM VERSION
  try:
    line, idx = find_line(rb'* xtb version', 1, 0)
    out_program = "XTB_" + line.split()[3]
  except:
    out_program = missing

  #METHOD
  try:
    line, idx = find_line(rb':  Hamiltonian', 1, 0)
    out_method = str(line.split()[2]).lower()
  except:
    out_method = missing

  #Number of atoms
  try:
    line, idx = find_line(rb"number of atoms", 1, 0)
    out_NAtoms = int(line.split()[4])
  except:
    out_NAtoms = missing

  #MULLIKEN ATOMIC CHARGES 
  try:
    out_mulliken_charges = [float(line.split()[1]) for line in PATTERN_XTB_out_mulliken_charges.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]     
  except:
    try:
      out_mulliken_charges = [float(line.split()[4]) for line in PATTERN_XTB_out_mulliken_charges2.findall(mm)[-1].decode("utf-8").split("\n")[:-1]]     
    except:
      out_mulliken_charges = missing
     
  #DIPOLE MOMENT
  try:
    lines = PATTERN_XTB_out_dipole_moment.findall(mm)[-1].decode("utf-8").split("\n")
    out_dipole_moment = float(lines[-2].split()[-1])
    out_dipole_moments = [float(i)/0.393456 for i in lines[-2].split()[0:3]]
  except:
    try:
      lines = PATTERN_XTB_out_dipole_moment2.findall(mm)[-1].decode("utf-8").split("\n")
      out_dipole_moment = float(lines[-2].split()[-1])
      out_dipole_moments = [float(i)/0.393456 for i in lines[-2].split()[1:4]]
    except:
      out_dipole_moment = missing
      out_dipole_moments = [missing]

  #rotational number
  try:
    line,idx = find_line(rb':  rotational number', 0, 0)
    out_rotational_symmetry_number = float(line.split()[3])
  except:
    out_rotational_symmetry_number = missing

  #ELECTRONIC ENERGY
  try:
    line,idx = find_line(rb'| TOTAL ENERGY', 0, 0)
    out_electronic_energy = float(line.split()[3])
  except:
    try:
      line,idx = find_line(rb'total E', 0, 0)
      out_electronic_energy = float(line.split()[3])
    except:
      out_electronic_energy = missing 

  #ENTHALPY "^H\(T\)"
  try:
    line, idx = find_line(rb'H(T)', 0, 0)
    out_enthalpy_energy = out_electronic_energy + float(line.split()[1])
  except:
    try:
      line, idx = find_line(rb'| TOTAL ENTHALPY', 0, 0)
      out_enthalpy_energy = float(line.split()[3])          
    except:
      out_enthalpy_energy = missing

  #FINAL ENTROPY
  try:
    line, idx = find_line(rb'        TOT                      ', 0, 0)
    out_entropy = float(line.split()[3])#/1000/627.503
  except:
    out_entropy = missing

  #VIBRATIONAL FREQUENCIES
  try:
    lines = PATTERN_XTB_out_vibrational_frequencies.findall(mm)[-1].decode("utf-8").split("\n")
    out_vibrational_frequencies = [float(item) for line in lines[1:-1] for item in line.split()[2:]]
  except:
    out_vibrational_frequencies = missing

  #ZERO POINT ENERGY CORRECTION + ENERGY
  try:
    line, idx = find_line(rb':: zero point energy', 1, idx)
    out_zero_point_correction = float(line.split()[4])
    out_zero_point_energy = out_electronic_energy + out_zero_point_correction
  except:
    out_zero_point_correction = missing  
    out_zero_point_energy = missing

  #GIBBS FREE ENERGY
  try:
    line, idx = find_line(rb'G(T)', 0, 0)
    out_gibbs_free_energy = out_electronic_energy + float(line.split()[1])
  except:
    try:
      line, idx = find_line(rb'| TOTAL FREE ENERGY', 0, 0)
      out_gibbs_free_energy = float(line.split()[4])
    except:
      out_gibbs_free_energy = missing

  #TIME
  try:
    line,idx = find_line(rb"* wall-time", 1, 0)
    out_time = float(line.split()[2])*24*60+float(line.split()[4])*60+float(line.split()[6])+float(line.split()[8])/60
  except:
    out_time = missing

  #SAVE
  mm.close()
  all_locals = locals()
  dic = {("log",column):[all_locals.get("out_"+column)] for column in columns}
  return dic

#else:
#  dic = {("log",column):[missing] for column in columns}
#  return dic
    
