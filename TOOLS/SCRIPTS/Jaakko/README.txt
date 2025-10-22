Here's the updated cluster filtering code. It needs to be run in a folder with input.txt or parameters.txt .
For example:
ClusterFilter results.dat -h 0H -stat -t -iso isomers.pkl
"-h 0H" select clusters with as many H bonds as there are H-donors. (changing the number shifts the cutoff) 
"-stat" prints out how hydrogen bond counts are distributed
"-t" selects clusters with unique H-bonding topologies (for large OOM clusters this can be more practical than -uniq rg,el,dip)
"-iso isomers.pkl" selects only clusters with containing isomers given in the isomers.pkl (useful if checking if cis-trans isometry has changed)

Also, nn-optimizer.py shows how I'm using the UMA model. ("turbo" mode makes it >2x faster, but it takes time to compile, so there's added code which loads/saves the compilation cache.)

