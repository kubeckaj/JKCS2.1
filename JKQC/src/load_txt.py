def mytofloat(string,element,QUentropy,QUenergy,Pout):
    try:
      try:
        if Pout[element] == "-s":
          multiplier = QUentropy
        elif Pout[element] == "-pXYZ":
          multiplier = 1
        else:
          multiplier = QUenergy
      except:
        multiplier = QUenergy
      return multiplier*float(string)
    except:
      return string

def load_txt(output,formation_input_file,QUentropy,QUenergy,Pout):
  from numpy import array
  f = open(formation_input_file, "r")
  output = list(output.transpose())
  for line in f.readlines():
    if not line.strip(): #allows for blank lines
      continue
    if str(line[0]) != "#":
      line_elements = line.split()
      output.append([mytofloat(line_elements[i],i,QUentropy,QUenergy,Pout) for i in range(len(line_elements))])
  output = array(output,dtype=object).transpose()
  f.close()
  return output
