#!/usr/bin/python
import numpy as np
import sys
import glob, os

# use either
#     python JKbonds.py S O 1.5 0.xyz
# or
#     python JKbonds.py S - 1.5 0.xyz

atom1=sys.argv[1]
atom2=sys.argv[2]
atom3=sys.argv[3]
threshold=float(sys.argv[4])

# I actually assume that just one file is given
if len(sys.argv) > 5:
  selectedfiles=sys.argv[5:]
else:
  selectedfiles=glob.glob("*.xyz")


for files in selectedfiles:
  #f1 = open(files,'r')
  #l1 = f1.readline()
  #l2 = f1.readline()
  #l3 = f1.read()
  #f1.close
  #if 'Gyration_radius:' in l2:
  #  continue

  try:
    data_raw=np.genfromtxt(files,dtype=("|S10", float, float, float),skip_header=2)
  except ValueError:
    print("JKrg.py: ERROR to read file "+files)
    continue

  totcount=0
  totcount2=0
  for i in range(len(data_raw)):
    if data_raw[i][0]==atom1:
      count=0
      count2=0
      for j in range(len(data_raw)):
        if data_raw[j][0]==atom2 or '-'==atom2:
          if j!=i:
            distance=((data_raw[i][1]-data_raw[j][1])**2+(data_raw[i][2]-data_raw[j][2])**2+(data_raw[i][3]-data_raw[j][3])**2)**(0.5)
            #print(i,distance,threshold,count)
            if distance<threshold:
              count=count+1
              if '-'!=atom3:
                for k in range(len(data_raw)):
                  if data_raw[k][0]==atom3:
                    if k!=i and k!=j:
                      distance=((data_raw[i][1]-data_raw[k][1])**2+(data_raw[i][2]-data_raw[k][2])**2+(data_raw[i][3]-data_raw[k][3])**2)**(0.5)
		      if distance<threshold:
                        count2=count2+1
      totcount=totcount+count
      totcount2=totcount2+count2
  if '-'==atom3:
    print totcount  
  else:
    print totcount2

        
