######################
#### ABC/CREST #######
######################
def read_abc_crest_init():
  pass

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

def read_abc_crest(mmm):
  global mm
  mm = mmm
  
  from numpy import array,sqrt
  from io import open
  missing = float("nan")
  
  columns = ["electronic_energy"]

  #ELECTRONIC ENERGY
  try:
    line,idx = find_line(rb'ABC energy:', 0, 0)
    out_electronic_energy = float(line.split()[2])
  except:
    try:
      line,idx = find_line(rb'structure energy:', 0, 0)
      out_electronic_energy = float(line.split()[2])
    except:
      out_electronic_energy = missing 

  #SAVE
  mm.close()
  all_locals = locals()
  dic = {("log",column):[all_locals.get("out_"+column)] for column in columns}
  return dic

#else:
#  dic = {("log",column):[missing] for column in columns}
#  return dic
    
