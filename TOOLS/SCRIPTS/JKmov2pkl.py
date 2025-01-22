from ase.io import read
from sys import argv
from os import path
from pandas import DataFrame, to_pickle
from read_files import mergeDictionary

file = argv[1]
if not path.exists(file):
  print(f"The file '{file}' does not exist.", flush=True)
  exit()

folder_path         = path.abspath(file)[::-1].split("/",1)[1][::-1]+"/"
file_basename_0     = file[:-4][::-1].split("/",1)[0][::-1]
strs = read(file,index=':')
columns = ["folder_path","file_basename"]
columns = columns + ["cluster_type","components","component_ratio"]
cluster_type = float("nan")
components = float("nan")
component_ratio = float("nan")
all_locals = locals()

clusters_dict = {}
for i in range(len(strs)):
  file_basename = file_basename_0+"-"+str(i)
  dic = {("info",column):[all_locals.get(column)] for column in columns}
  dic.update({("xyz","structure"):[strs[i]]})
  dic.update({("xyz","id1"):[float("nan")]})
  try:
    en = float(strs[i].info['Energy'])
    dic.update({("log","electronic_energy"):[en]})
  except:
    dic.update({("log","electronic_energy"):[float("nan")]})
 
  clusters_dict = mergeDictionary(clusters_dict,dic)

clusters_df = DataFrame(clusters_dict,index=range(0,len(strs)))
to_pickle(clusters_df,file_basename_0+".pkl")
