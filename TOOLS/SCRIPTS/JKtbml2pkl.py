#https://wiki.fysik.dtu.dk/ase/_modules/ase/io/turbomole.html
from ase import Atoms
from os import path
from pandas import DataFrame
import os
import sys

def check_file_exists(file_path):
  if not os.path.isfile(file_path):
    print(f"The file '{file_path}' does not exist.")
    exit()

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

def find_next_occurrence(arr, target):
  result = []
  n = len(arr)
  
  for i in range(n):
    found = False
    for j in range(i + 1, n):
      if arr[j] == target:
        result.append(j)
        found = True
        break
    if not found:
      result.append(-1)
  
  return result

if __name__ == "__main__":
  if len(sys.argv) < 2:
    print("Please provide a file path as an argument.")
  else:
    file_path = sys.argv[1]
    if file_path == "-help":
      print("Just run: JKtmbl2pkl <file>")
      exit()
    check_file_exists(file_path) 

  f = open(file_path, "r")
  file_data = f.read()
  f.close()
  folder_path = path.abspath(file_path)[::-1].split("/",1)[1][::-1]+"/"
  file_data_splitted = file_data.split("\n")
  ii = 0
  clusters_dict = {}
  next_end = find_next_occurrence(file_data_splitted,"$end")  

  i = 0
  while i < len(file_data_splitted)-1:
  #for i in range(0,len(file_data_splitted)-9,9):
    #print([i,next_end[i]])
    i_end = next_end[i]
    if i_end == -1:
      print("Oh, something weird is going on with the file format.")
    n_atoms = int((i_end - i - 2)/2)
    ii = ii + 1
    #EXTRACT NEW DATA
    energy = float(file_data_splitted[i+1].split()[6].replace("D", "E"))
    atoms = "".join([ j.split()[3].capitalize() for j in file_data_splitted[i+2:i+2+n_atoms] ])
    xyz = Atoms(atoms,[[0.529177*float(iii.replace("D", "E")) for iii in j.split()[0:3]] for j in file_data_splitted[i+2:i+2+n_atoms]])
    forces = [[-1/0.529177*float(iii.replace("D", "E")) for iii in j.split()[0:3]] for j in file_data_splitted[i+2+n_atoms:i+2+n_atoms+n_atoms]]
  
    #ADD TO DICTIONARY 
    dic = {("info","folder_path"):[folder_path]}
    dic.update({("info","file_basename"):["str-"+str(ii)]})
    dic.update({("xyz","structure"):[xyz]})
    dic.update({("log","electronic_energy"):[energy]})
    dic.update({("extra","forces"):[forces]})
   
    clusters_dict = mergeDictionary(clusters_dict,dic)
    i = i_end + 1
    
  clusters_df = DataFrame(clusters_dict,index=range(ii))
  
  try:
    clusters_df.to_pickle("mydatabase.pkl")
    print("The mydatabase.pkl has been hopefully created.")
  except:
    print("Something got fucked up.")
  
  
