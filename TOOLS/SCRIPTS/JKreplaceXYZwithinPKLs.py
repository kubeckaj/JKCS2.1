import pandas as pd
from pandas import read_pickle
import sys

if len(sys.argv) < 3:
	print("Usage: JKreplaceXYZwithinPKLs.py <input pkl file> <input pkl file with the xyz> <output pkl file>")
	sys.exit()

a = read_pickle(sys.argv[1]).sort_values([('info','file_basename')])
b = read_pickle(sys.argv[2]).sort_values([('info','file_basename')])
if len(a) != len(b):
	print("Files are not the same length")
	sys.exit()

a['xyz'] = b['xyz']
a.to_pickle(sys.argv[3])

