#!/appl/opt/python/2.7.13-gcc540-shared/bin/python
###
### Loosly based on B.B. Chaudhuri, Pattern Recognition Letters, Volume 15, Issue 9, 1994, Pages 893-99
###

###
### Chooses from a 'resultXTB.dat' file a representative subset of data based on a manually chosen parameter
###

###		 ### 
### Initializing ###
###              ###

import sys

def printhelp():
  print('THIS SCRIPT WAS MAINLY WRITTEN BY VITTUSAATANA :-D')
  print('   Use this file like this:')
  print('           DataFilter.py [file] {OPTIONS}')
  print(' OPTIONS:')
  print('  *** UNIQUENESS ***')
  print('  -u or -uniq  uniqueness')
  print('  -u1   2      uniqueness grid on 2 decimals for Rg (column 1) [default=2]')
  print('  -u2   5      uniqueness grid on 5 decimals for E  (column 2) [default=3]')
  print('  -u3   4      uniqueness grid on 4 decimals for     column 3  [default=3]')
  print('  -sort 4      sort values with respect to column 4 [default=3]')
  print('  *** FILTERING ***')
  print('  -rg  10.0    filter out all files with gyr. rad. greater than 10.0 a.u.')
  print('  -rgm  4.0    filter out all files with gyr. rad. greatet than 4.0*(#ofMolecules) a.u.')
  print('  -d   15      filter out structures with energy higher than 15 kcal/mol THAN GLOBAL MINNIMUM')
  print('  -dm   5      filter out structures with energy higher than 5*(#ofMolecules) kcal/mol THAN GLOBAL MINNIMUM')
  print('  -en -10.5    filter out all files with energy higher than -10.5 a.u.')
  print('  -bonds 4     filter out if amount of bonds for specified atom is greater than 4')
  print('  *** SAMPLING/SELECTION ***')
  print('  -s   20      will select aprox. 20 structures')
  print('  -sm  10      will select aprox. 10*(#ofMolecules) structures')
  print('  -r    0.05   select points just from area of radius 0.05 reduced units')
  print('  ###########')
  print('  -c1   2      1. column in analysed file (Rg)       [default = 2]')
  print('  -c2   4      2. column in analysed file (Energy)   [default = 3]             ! GIBBS FREE ENERGY IS IN 4th column')
  print('  -c3   4      3. column in analysed file (Whatever) [default = NO 4th COLUMN] ! GIBBS FREE ENERGY IS IN 4th column')
  print('  ')
  print('EXAMPLE:')
  #print('   DataFilter.py resultsXTB.dat -r 0.05 -rg 15 -s 50') 
  #print('   DataFilter.py resultsDFT_HIGH_freq.dat -r 0.05 -s 5 -d 3') 
  #print(' or')
  print('   JKCS7_filter resultsXTB.dat -u -u1 2 -u2 4 -u3 1 -sort 3') 
  print('   JKCS7_filter resultsXTB.dat -rgm 4 -dm 4 -sm 20') 
  print('   JKCS7_filter resultsXTB.dat -u3 1 -c3 4 -rgm 4 -dm 4 -s 100') 
  print('   JKCS7_filter resultsDFT_HIGH.dat -s 10 -dm 1.2') 
  print('   JKCS7_filter resultsDFT_HIGH_freq.dat -c2 4 -d 1') 

#################################################################################################################
#################################################################################################################
#################################################################################################################
c1=2            # AHHHHHHHHHHHHH Columns are hardcoded here
c2=3		#Switchable from on to two dimensions??
c3=0
Xrg=100
Xen=0

Qen=0
Qrg=0
Qbonds=0
Qfilter=0
Qdelete=0
Qselect=0

last=''
Narg=len(sys.argv)
for i in range(1,Narg):
  arg=sys.argv[i]
  if str(arg) == "-help":
    printhelp()
    exit()
  # c1
  if str(arg) == "-c1":
    last='-c1'
    continue
  if last == '-c1':
    c1=int(arg)
    last=''
    continue
  # c2
  if str(arg) == "-c2":
    last='-c2'
    continue
  if last == '-c2':
    c2=int(arg)
    last=''
    continue
  # c3
  if str(arg) == "-c3":
    last='-c3'
    continue
  if last == '-c3':
    c3=int(arg)
    last=''
    continue
  # select
  if str(arg) == "-s":
    Qselect=1
    last='-s'
    continue
  if last == '-s':
    select=float(arg)
    last=''
    continue
  # select(Molecules)
  if str(arg) == "-sm":
    Qselect=1
    last='-sm'
    continue
  if last == '-sm':
    molecules=int(open("parameters.txt").readline().rstrip().split()[5])
    select=float(arg)*molecules
    last=''
    continue
  # radius
  if str(arg) == "-r":
    Qfilter=1
    last='-r'
    continue
  if last == '-r':
    radius=float(arg)
    last=''
    continue
  #delete
  if str(arg) == "-d":
    Qdelete=1
    last='-d'
    continue
  if last == '-d':
    delete=float(arg)
    last=''
    continue
  #delete(Molecules)
  if str(arg) == "-dm":
    Qdelete=1
    last='-dm'
    continue
  if last == '-dm':
    molecules=int(open("parameters.txt").readline().rstrip().split()[5])
    delete=float(arg)*molecules
    last=''
    continue
  # remove excesed energies
  if str(arg) == "-en":
    Qen=1
    last='-en'
    continue
  if last == '-en':
    deleteEN=float(arg)
    last=''
    continue
  #remove excesed rg
  if str(arg) == "-rg":
    Qrg=1
    last='-rg'
    continue
  if last == '-rg':
    deleteRG=float(arg)
    last=''
    continue
  #remove excesed bonds
  if str(arg) == "-bonds":
    Qbonds=1
    last='-bonds'
    continue
  if last == '-bonds':
    deleteBONDS=float(arg)
    last=''
    continue
  #remove excesed rg(Molecules)
  if str(arg) == "-rgm":
    Qrg=1
    last='-rgm'
    continue
  if last == '-rgm':
    molecules=int(open("parameters.txt").readline().rstrip().split()[5])
    deleteRG=float(arg)*molecules
    last=''
    continue
  #file
  file_in=arg
#################################################################################################################
#################################################################################################################
#checking some inputs

if ( file_in == "" ):
  print('DataFilter.py: Missing input file. EXITING')
  exit()

# check amount of rows
count = 0
for line in open(file_in): 
  count += 1
  if ( count > 1 ):
    break 
if ( count == 0 ):
  print('DataFilter.py: Empty input file. EXITING')
  exit()
if ( count == 1 ):
  print('DataFilter.py: Just one structure is in the input file (no filtering/sampling is needed). EXITING')
  exit()

#if nothing then at least Qfilter
if ( ( Qdelete == 0 ) and ( Qfilter == 0 ) and ( Qen == 0) and ( Qrg == 0) ):
  Qfilter=1
  radius=0.05

#if just delete
if ( ( Qselect == 0 ) and ( Qfilter == 0 ) ):
  Qfilter=1
  radius=0.0


#if nothing then at least Qfilter
if ( ( Qselect == 1 ) and ( Qfilter == 0) ):
  Qfilter=1
  radius=0.05

c1=c1-1      #More hardcode
c2=c2-1
c3=c3-1
#################################################################################################################
#################################################################################################################
#################################################################################################################
import numpy
from random import randint
import math
print('DataFilter.py: Start.')
#file to be printed  
file_out=file_in[:-4]+'_FILTERED.dat'

############### Introduce number for number of columns to be taken
#colNum=3 # Maybe array with one element per column?
#colNumArr=arange(3)  #Gives me array([0,1,2])
numpy.set_printoptions(threshold=numpy.nan)

#read filei
if Qbonds == 0:
  size=3
  if c3 == -1:
    data=numpy.loadtxt(file_in,usecols=(c1,c2))
    data=numpy.insert(data,2,0,axis=1)
  else:
    data=numpy.loadtxt(file_in,usecols=(c1,c2,c3))  #Need variable or array in usecols=(???) Looping, reading in he columns singly and the appending them?
else:
  size=4
  if c3 == -1:
    data=numpy.loadtxt(file_in,usecols=(c1,c2,-1))
    data=numpy.insert(data,2,0,axis=1)
  else:
    data=numpy.loadtxt(file_in,usecols=(c1,c2,c3,-1))  #Need variable or array in usecols=(???) Looping, reading in he columns singly and the appending them?

with open(file_in) as f:			# IF c3 is empty, there isnone, let c3 run through with only zeros!
    data2 = f.read().splitlines()

if data.size == 2:  ####?????????????? I hope this does not matter
  thefile = open(file_out, 'w')
  print>>thefile, data2[0]
  exit()
  
l=len(data)
#print(data)
xdata0=data.reshape(size*l,1)[0::size]
ydata0=data.reshape(size*l,1)[1::size]
zdata0=data.reshape(size*l,1)[2::size]

######################################
### PREPARING - removing stupid values
######################################
QdeletestupidE=0		#I do not exactly understand what is happening here - I ignore it fornow
QdeletestupidRG=1
if QdeletestupidE == 1:
  data= data[data[:, 1].argsort()]
  kickout=l
  for i in range(0,l):
    if data[i,1] >= Xen:
      kickout=i
      break
  data=data[0:kickout]
  l=len(data)

if QdeletestupidRG == 0:
  data= data[data[:, 0].argsort()]
  kickout=l
  for i in range(0,l):
    if data[i,0] >= Xrg:
      kickout=i
      break
  data=data[0:kickout]
  l=len(data)

if Qrg == 1 or Qen == 1 or Qdelete == 1 or Qfilter == 1:
  print('DataFilter.py: FILTERING')

######################################
### Qen + Qrg / cutting out some energies and Rg
######################################
if Qrg == 1:
  data= data[data[:, 0].argsort()]
  kickout=l
  for i in range(0,l):
    if data[i,0] >= deleteRG:
      kickout=i
      break
  data=data[0:kickout]
  l=len(data)

if Qen == 1:
  data= data[data[:, 1].argsort()]
  kickout=l
  for i in range(0,l):
    if data[i,1] >= deleteEN:
      kickout=i
      break
  data=data[0:kickout]
  l=len(data)

if Qbonds == 1:
  data= data[data[:, -1].argsort()]
  kickout=l
  for i in range(0,l):
    if data[i,-1] > deleteBONDS:
      kickout=i
      break
  data=data[0:kickout]
  l=len(data)
  #change to normal size
  data=data[:,0:size-1]
 
####################
### FILTERING/DELETING
###################

if Qdelete == 1:
  #l=len(data)
  data= data[data[:, 1].argsort()]
  energyminimum=data[0,1]
  kickout=l
  if  file_in == 'resultsABC.dat' :
    energyconvert=1.0/0.239006
  else:
    energyconvert=1.0/627.5
  for i in range(0,l):
    if data[i,1] >= energyminimum + float(delete)*energyconvert:
      kickout=i
      break
  data=data[0:kickout]
  l=len(data)

################################
### SAMPLING/SELECTION
################################
if Qfilter == 1:
  l=len(data)
  xdata=data.reshape(3*l,1)[0::3]
  ydata=data.reshape(3*l,1)[1::3]
  zdata=data.reshape(3*l,1)[2::3]
  lx=len(xdata)
 #For Normalization
  xmin=min(xdata)
  xmax=max(xdata)
  dx=xmax-xmin
  ymin=min(ydata)
  ymax=max(ydata)
  dy=ymax-ymin
  zmin=min(zdata)
  zmax=max(zdata)
  dz=zmax-zmin
  #print('   ')
  #print(zmin)
  #print(zmax)
  #print(dz)
  if ( dx == 0 ):
    dx=numpy.array([1])

  if ( dy == 0 ):
    dy=numpy.array([1])
 
  if ( dz == 0 ):
    dz=numpy.array([1])

  xdataNorm=(xdata-xmin)/dx
  ydataNorm=(ydata-ymin)/dy
  zdataNorm=(zdata-zmin)/dz
  
#  print(ymin)
#  print(ymax)
#  print(dy)
#  print(xmin)
#  print(xmax)
#  print(dx)
  #drChosen=0.01
  #dEChosen=0.001

  ###                                       ###
  ### Calculate Densities for all Elements  ###
  ###                                       ###

 #JKsorry d= numpy.zeros(lx) 
 #JKsorry print('calc densities')
 #JKsorry for i in range(l):
 #JKsorry   for j in range(l):
 #JKsorry     diff= 1.0*math.sqrt((xdataNorm[i]-xdataNorm[j])**2.+(ydataNorm[i]-ydataNorm[j])**2.)
 #JKsorry     if ( diff < 0.1 ): #0.1 is just a manually chosen factor for size of density
 #JKsorry       d[i] = 1.*d[i]+1

  #Assign those densities to the according data points
  #JKsorry d = (1)*d.reshape(lx,1)
  #JKsorry numpy.set_printoptions(threshold='nan')
  #JKsorry all = numpy.column_stack((d, xdataNorm, ydataNorm))
  #JKsorry all = all[all[:, 0].argsort()]
  # JK let sort with respect to energy
  all = numpy.column_stack(( xdataNorm, ydataNorm, zdataNorm))
  all = all[all[:, 1].argsort()]  
   ###	            ###
  ### The algorithm ###
  ###		    ###
  #print(all)
  if Qselect == 0:                #Is this case happening at all?
    print('DataFilter.py: SAMPLING/SELECTING')
    c=[10,10,10] #The 0,0,0 Element is removed in the end automatically by the formatting loop
    while (0 < len(all)):
      #print('new round')
      #print(all[0:10])
      #print(all[0])
      
      c = numpy.vstack((c, all[0])) 
      j = 1
      while (j < len(all)):       #Not sure what is going on here?
        if ( math.sqrt((all[0,0] - all[j,0])**2 +(all[0,1]-all[j,1])**2 + (all[0,2]-all[j,2])**2 ) < radius ):
          all = numpy.delete(all, j, 0)
        else:	
          j = j + 1
      all = numpy.delete(all, 0, 0)	
  else:
    print('DataFilter.py: SAMPLING/SELECTING')
    allBCKP=all
    test=0
    factor=0.05
    state=0 #many=2 or low=1
    while ( test == 0 ): 
      if len(allBCKP) <= select+1:
        radius=0.0
        test=1 
      all=allBCKP
      c=[10000,100000,100000] #The 0,0 Element is removed in the end automatically by the formatting loop
      while (0 < len(all)):
        c = numpy.vstack((c, all[0]))
        j = 1
        while (j < len(all)):
          if ( math.sqrt((all[0,0] - all[j,0])**2 +(all[0,1]-all[j,1])**2 + (all[0,2]-all[j,2])**2) < radius ):
            all = numpy.delete(all, j, 0)
          else:
            j = j + 1
        all = numpy.delete(all, 0, 0)
      print('DataFilter.py:    selected: '+str(len(c)-1)+', radius: '+str(radius))
      
      if len(c) == select+1:
        test=1
      else:
        if len(c) < select+1:
          #too large radius
          if state == 1:
            factor=factor/2
          state=2
          radius=radius-factor
        else:
          #too small radius
          if state == 2:
            factor=factor/2
          state=1
          radius=radius+factor
      if radius < 0.0:
        radius=0.0
      if factor < 0.000001:
        test=1
      
 

  ### 							###
  ### Format the chosen Data points and write result file ###
  ###							###
  #print(c)
  p=len(c)
  xChosen = c.reshape(3*p,1)[0::3]*dx+xmin
  yChosen = c.reshape(3*p,1)[1::3]*dy+ymin
  zChosen = c.reshape(3*p,1)[2::3]*dz+zmin
  #print(zChosen) 
  # Possibly iterating through all lines of XTB data for every line  of xChosen/yChosen to find matches and print the structures after
  k = 0
  tolerance1=0.0000001
  tolerance2=0.0000001
  tolerance3=0.0000001
  cFormat=[]
  while (k < len(c)):
    i = 0	
    while (i < len(xdata0) ): #i should start at 0
      if( abs(xChosen[k]-xdata0[i]) <= tolerance1 and abs(yChosen[k]-ydata0[i]) <= tolerance2 and abs(zChosen[k]-zdata0[i] <= tolerance3)) :
        cFormat.append(data2[i])
        break
      i = i + 1
    k = k + 1
   
  #print(str(data2))
#################################################

#################################################
# print out
#################################################
f = open("FILTER.txt", "a")
f.write(file_in + ' : radius = ' + str(radius) + ' xmin/dx/ymin/dy/zmin/dz ' + str(xmin[0]) + ' ' + str(dx[0]) + ' ' + str(ymin[0]) + ' ' + str(dy[0]) + ' ' + str(zmin[0]) + ' ' + str(dz[0]) + '\n')
f.close()
thefile = open(file_out, 'w')
for item in cFormat:
  print>>thefile, item



print('DataFilter.py: Done.')






