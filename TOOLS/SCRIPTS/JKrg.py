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
  if 'Gyration_radius:' in l2:
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
        if data_raw[i][0]=='AR' or data_raw[i][0]=='Ar' or data_raw[i][0]==b'AR' or data_raw[i][0]==b'Ar':
                data_with_mass[i][4]=40
                M+=40
        if data_raw[i][0]=='C' or data_raw[i][0]==b'C':
                data_with_mass[i][4]=12
                M+=12
        if data_raw[i][0]=='H' or data_raw[i][0]==b'H':
                data_with_mass[i][4]=1
                M+=1
        if data_raw[i][0]=='O' or data_raw[i][0]==b'O':
                data_with_mass[i][4]=16
                M+=16
        if data_raw[i][0]=='N' or data_raw[i][0]==b'N':
                data_with_mass[i][4]=14
                M+=14
        if data_raw[i][0]=='I' or data_raw[i][0]==b'I':
                data_with_mass[i][4]=127
                M+=127
        if data_raw[i][0]=='S' or data_raw[i][0]==b'S':
                data_with_mass[i][4]=32
                M+=32
        if data_raw[i][0]=='B' or data_raw[i][0]==b'B':
                data_with_mass[i][4]=11
                M+=11
        if data_raw[i][0]=='Na' or data_raw[i][0]==b'Na':
                data_with_mass[i][4]=11
                M+=11
        if data_raw[i][0]=='Cl' or data_raw[i][0]==b'Cl':
                data_with_mass[i][4]=11
                M+=11
  
  ### JK move to center of mass
  mixi=0.0 #sum of m_i x_i
  miyi=0.0
  mizi=0.0
  for i in range(len(data_with_mass)):
  	mixi += data_with_mass[i][4]*data_with_mass[i][1]
  	miyi += data_with_mass[i][4]*data_with_mass[i][2]
  	mizi += data_with_mass[i][4]*data_with_mass[i][3]
  
  for i in range(len(data_with_mass)):
        data_with_mass[i][1]=data_with_mass[i][1]-mixi/M
        data_with_mass[i][2]=data_with_mass[i][2]-miyi/M
        data_with_mass[i][3]=data_with_mass[i][3]-mizi/M 

  # gyration radius
  mixi=0.0 #sum of m_i x_i
  miyi=0.0
  mizi=0.0
  for i in range(len(data_with_mass)):
          mixi += data_with_mass[i][4]*data_with_mass[i][1]**2
          miyi += data_with_mass[i][4]*data_with_mass[i][2]**2
          mizi += data_with_mass[i][4]*data_with_mass[i][3]**2
  #print 'Gyration_radius [Angstrom] = ', math.sqrt((mixi+miyi+mizi)/M)
  ### np.savetxt('F_conf.dat',np.c_[F[:,0],F[:,1]])
  
  f2 = open(files,'w')
  f2.write(l1)
  f2.write(l2.rstrip('\n')+'  Gyration_radius:  '+str(((mixi+miyi+mizi)/M)**0.5)+'\n ')
  f2.write(l3)
  f2.close
####
