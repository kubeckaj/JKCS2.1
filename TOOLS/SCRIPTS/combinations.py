import numpy
import sys

#numpy: asarray, array
def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)


#sys.argv(1)
maximum=map(int,sys.argv[2].split(',')[1:])
structures=map(int,sys.argv[1].split(',')[1:])
charges=map(int,sys.argv[3].split(',')[1:])
totcharge=int(sys.argv[4])
npcharges=numpy.array(charges)

#maximum=map(int,maximum)
#print(maximum)
#print(structures)
#print(type(structures))
#print(len(maximum))
#print(maximum[0])
#print(size(npcharges))

list2=[]
listch=[]
ch1=0
for j in range(len(maximum)):
  comb = product(range(maximum[j]+1),repeat=structures[j])
  list1=numpy.empty([0,structures[j]], int)
  #print(numpy.size(list1,1))
  for i in list(comb): 
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

#print(list2)
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
