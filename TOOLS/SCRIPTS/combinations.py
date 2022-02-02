import numpy
import sys
#import math

#numpy: asarray, array
def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = list(map(tuple, args)) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def productJK(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = list(map(tuple, args)) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool if sum(x)+y<=max(max(args))]
    for prod in result:
        yield tuple(prod)

#sys.argv(1)
maximum=list(map(int,sys.argv[2].split(',')[1:]))
structures=list(map(int,sys.argv[1].split(',')[1:]))
charges=list(map(int,sys.argv[3].split(',')[1:]))
totcharge=int(sys.argv[4])
npcharges=numpy.array(charges)

#maximum=map(int,maximum)
#print("maximum")
#print(maximum)
#print(structures)
#print(type(structures))
#print("len(maximum)")
#print(len(maximum))
#print(maximum[0])
#print(size(npcharges))

list2=[]
listch=[]
ch1=0
for j in range(len(maximum)):
  #print("LOOP"+str(j))
  comb = productJK(range(maximum[j]+1),repeat=structures[j])
  list1=numpy.empty([0,structures[j]], int)
  #print("NUMPY SIZE "+str(numpy.size(list1,1)))
  #print("comb:")
  #print(comb)
  #print("list comb:")
  #print(list(comb))
  for i in list(comb): 
    #print("LOOP"+str(j)+str(i))
    #print(numpy.size(numpy.asarray(i,int)))
    toappend=numpy.asarray(i,int)
    #print(toappend)
    #print(numpy.sum(toappend))
    if numpy.sum(toappend)==maximum[j]:
      list1=numpy.concatenate((list1,[toappend]),axis=0)
  #print(list1)
  list2.append(list1)
  ch2=ch1+structures[j]
  listch.append(charges[ch1:ch2])
  ch1=ch2
  #del list1

#for j in range(len(maximum)):
#  print("LOOP"+str(j))
#  comb = product(range(maximum[j]+1),repeat=structures[j])
#  list1=numpy.empty([0,structures[j]], int)
#  print("NUMPY SIZE "+str(numpy.size(list1,1)))
#  #print("comb:")
#  #print(comb)
#  print("list comb:")
#  print(list(comb))
#  for i in list(comb): 
#    print("LOOP"+str(j)+str(i))
#    print(numpy.size(numpy.asarray(i,int)))
#    toappend=numpy.asarray(i,int)
#    print(toappend)
#    print(numpy.sum(toappend))
#    if numpy.sum(toappend)==maximum[j]:
#      list1=numpy.concatenate((list1,[toappend]),axis=0)
#  print(list1)
#  list2.append(list1)
#  ch2=ch1+structures[j]
#  listch.append(charges[ch1:ch2])
#  ch1=ch2
#  #del list1

#print(list2)
#print("DONE")
listJK=[]
for i in product(*list2): 
  toappend=numpy.concatenate(i)
  if numpy.sum(toappend*npcharges)==totcharge:
    listJK.append(toappend)

#print(listJK)
for i in listJK:
  text=""
  for j in i:
    add=numpy.array2string(j)
    if text=="":
      text=add
    else:
      text=text+"_"+add
  print(text)
