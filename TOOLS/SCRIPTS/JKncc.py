#!/usr/bin/python2
import numpy as np
import sys
import glob, os

if len(sys.argv) > 1:
  selectedfiles=sys.argv[1:]
else:
  selectedfiles=glob.glob("*.xyz")

for files in selectedfiles:
  f1 = open(files,'r')
  l1 = f1.readline()
  l2 = f1.readline()
  l3 = f1.read()
  f1.close
  if 'NCC:' in l2:
    continue

  try: 
    data_raw=np.genfromtxt(files,dtype=("|S10", float, float, float),skip_header=2)
  except ValueError:
    print("JKrg.py: ERROR to read file "+files)
    continue
  
  data_with_mass = np.zeros(len(data_raw),dtype=[('Atom','|S10'),
  						 ('x','float'),
  						 ('y','float'),	
  						 ('z','float'),
  						 ('Mass','float')] )
  M=0.0 #total mass
  for i in range(len(data_raw)):
  	data_with_mass[i][0]=data_raw[i][0]
  	data_with_mass[i][1]=data_raw[i][1]
  	data_with_mass[i][2]=data_raw[i][2]
  	data_with_mass[i][3]=data_raw[i][3]
  	if data_raw[i][0]=='AR' or data_raw[i][0]=='Ar':
  		data_with_mass[i][4]=40
  		M+=40
  	if data_raw[i][0]=='C':
  		data_with_mass[i][4]=12
  		M+=12
  	if data_raw[i][0]=='H':
  		data_with_mass[i][4]=1
  		M+=1
  	if data_raw[i][0]=='O':
  		data_with_mass[i][4]=16
  		M+=16
  	if data_raw[i][0]=='N':
  		data_with_mass[i][4]=14
  		M+=14
  	if data_raw[i][0]=='I':
  		data_with_mass[i][4]=127
  		M+=127
  	if data_raw[i][0]=='S':
  		data_with_mass[i][4]=32
  		M+=32
  	if data_raw[i][0]=='B':
  		data_with_mass[i][4]=11
  		M+=11
  
  ### JK 
  e1=0.0
  e2=0.0
  for i in range(len(data_with_mass)):
        if data_raw[i][0]=='H': 
		min1=100
		min2=100
		for j in range(len(data_with_mass)):
            		if data_raw[j][0]=='O':
				dist=((data_with_mass[i][1]-data_with_mass[j][1])**2+(data_with_mass[i][2]-data_with_mass[j][2])**2+(data_with_mass[i][3]-data_with_mass[j][3])**2)**0.5
  				e1 += 1/dist
            		if data_raw[j][0]=='N':
				dist=((data_with_mass[i][1]-data_with_mass[j][1])**2+(data_with_mass[i][2]-data_with_mass[j][2])**2+(data_with_mass[i][3]-data_with_mass[j][3])**2)**0.5
  				e2 += 1/dist

  #print 'NCC = ', e1/e2
  ### np.savetxt('F_conf.dat',np.c_[F[:,0],F[:,1]])
  
  f2 = open(files,'w')
  f2.write(l1)
  f2.write(l2.rstrip('\n')+'  NCC:  '+str(e1/(e2+e1))+'\n ')
  f2.write(l3)
  f2.close
####
