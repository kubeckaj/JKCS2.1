##############
### MRCC #####
##############
def read_mrcc_init():
  #from re import compile
  return

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

def read_mrcc(mmm):
  global mm
  mm = mmm
  
  from datetime import datetime
  from io import open
  missing = float("nan")

  columns = ["program","method","time","electronic_energy","scf_energy","correlation_energy"]

  #PROGRAM VERSION
  #try:
  #  line, idx = find_line(rb'* xtb version', 1, 0)
  #  out_program = "XTB_" + line.split()[3]
  #except:
  #  out_program = missing
  out_program = "MRCC"

  #METHOD
  #try:
  #  line, idx = find_line(rb':  Hamiltonian', 1, 0)
  #  out_method = str(line.split()[2]).lower()
  #except:
  #  out_method = missing
  out_method = missing

  #TIME
  try:
    lines,idx = find_lines(rb"************************ ", 1, 0) 
    out_time = datetime.strptime(lines[-1].split()[1]+" "+lines[-1].split()[2], "%Y-%m-%d %H:%M:%S")-datetime.strptime(lines[0].split()[1]+" "+lines[0].split()[2], "%Y-%m-%d %H:%M:%S")
    out_time =out_time.total_seconds() / 60
  except:
    out_time = missing
  #NOTE: MRCC for some reason has unique text for each type of calculations, these will only work for LNO-CCSD(T) methods...
  #ELECTRONIC ENERGY
  try:
    line,idx = find_line(rb'Total LNO-CCSD(T) energy with MP2 corrections', 0, 0)
    out_electronic_energy = float(line.split()[7])
  except:
    out_electronic_energy = missing
  #SCF ENERGY
  try:
    line,idx = find_line(rb'Reference energy', 0, 0)
    out_scf_energy = float(line.split()[3])
  except:
    out_scf_energy = missing
  #CORRELATION ENERGY
  try:
    line,idx = find_line(rb'CCSD(T) correlation energy + MP2 corrections', 0, 0)
    out_correlation_energy = float(line.split()[7])
  except:
    out_correlation_energy = missing
  #SAVE
  mm.close()
  all_locals = locals()
  dic = {("log",column):[all_locals.get("out_"+column)] for column in columns}
  return dic

