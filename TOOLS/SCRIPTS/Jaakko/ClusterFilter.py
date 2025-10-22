
import os, sys

def printhelp():
  print('   Use this file like this:')
  print('           ClusterFilter.py [FILE] {OPTIONS}')
  print(' OPTIONS:')
  print('  -i/-itol 0.3     filter out if any intermolecular bond length has changed more than 30% [default=0.2]')
  print('  -d/-dmol 2.5     filter out if smallest distance between molecules greater than 2.5 [default=2.7]')
  print('  -h 2(D/H/A)      filter clusters based on the number of hydrogen bonds')
  print('  -t               filter the lowest unique topologies')
  print('  -r/-cutr 5       filter clusters with electronic energies above a relative cutoff')
  print('  -s/-select 5     select at most this many clusters per cluster type')
  print('  -maxf/minf 0.01  filter (un)converged clusters with max force limit of 0.01 eV/Ang')
  print('  -except sa,am    filter out subclusters consisting only specified molecules (separated with ",")')
  print('  -subcl 2         minimum size of returned subclusters (1 returns monomers) [default=2]')
  print('  -iso iso.txt     monomer isomers to select')
  print('  -sort el/g/gout  sort energies by X [default=el]')
  print('  ')

#######################################################################

internal_tol = 0.2      # largest allowed relative error
max_dist = 2.7          # maximum distance for H-bonding between molecules
min_dist = 1.4          # minimum allowed distance between molecules (otherwise a reaction has occured)
print_stats = False
subcl_size = 2
Hbonding = None
topology = 0
relative_cutoff = -1
select = 0
isomer_file = None
to_extract = None
to_except = None
sort_by=('log','electronic_energy')

maxf = None
converged = None

# Look for inputs
file_in = ""
file_not = None
mol_file = None
last = ""
rejected_types = [] # e.g. sa,a,so4
n_arg = len(sys.argv)
for i in range(1,n_arg):
    arg=sys.argv[i]
    if str(arg) == "--help" or str(arg) == "-help":
        printhelp()
        exit()
    if arg == "-itol" or arg == '-i':
        last = "-itol"
        continue
    if last == "-itol":
        last = ""
        internal_tol = float(arg)
        continue
    if arg == "-dmol" or arg == '-d':
        last = "-dmol"
        continue
    if last == "-dmol":
        last = ""
        max_dist = float(arg)
        continue
    if arg == "-maxf" or arg == '-d':
        last = "-maxf"
        continue
    if last == "-maxf":
        last = ""
        maxf = float(arg)
        converged = True
        continue
    if arg == "-minf" or arg == '-d':
        last = "-minf"
        continue
    if last == "-minf":
        last = ""
        maxf = float(arg)
        converged = False
        continue
    if arg == "-h":
        last = "-h"
        continue
    if last == "-h":
        last = ""
        if arg[-1].isdigit():
            Hbonding = [int(arg), None]
        else:
            Hbonding = [int(arg[:-1]), arg[-1].lower()]
        continue
    # TODO: el,g,etc...
    if arg == "-r" or arg == '-cutr':
        last = "-r"
        continue
    if last == "-r":
        if arg.isalpha():
            continue
        last = ""
        relative_cutoff = float(arg) # in kcal/mol
        continue
    if arg == "-s" or arg == '-select':
        last = "-s"
        continue
    if last == "-s":
        last = ""
        select = int(arg)
        continue
    if arg == "-mol_file" or arg=='--mol_file':
        last = "-mol_file"
        continue
    if last == "-mol_file":
        last = ""
        mol_file = arg
        continue
    if arg == "-iso" or arg=='-isomers':
        last = "-iso"
        continue
    if last == "-iso":
        last = ""
        isomer_file = arg
        if topology==0:
            topology=2
        continue
    if arg == "-extract":
        last = "-extract"
        continue
    if last == "-extract":
        last = ""
        to_extract = arg
        continue
    if arg == "-except":
        last = "-except"
        continue
    if last == "-except":
        last = ""
        to_except = arg
        continue
    if arg == "-sort":
        last = "-sort"
        continue
    if last == "-sort":
        last = ""
        sort_by = arg
        continue
    if arg == "-subcl":
        last = "-subcl"
        continue
    if arg == "-subcl":
        last = ""
        subcl_size = int(arg)
        continue
    if arg == "-stat" or arg == "-stats":
        last = ""
        print_stats = True
        continue
    if last == "-minimums":
        last = ""
        print_minimums = False
        count_minimums = int(arg)
        continue
    if arg == "-t" or arg == '-topo':
        last = "-t"
        topology = 1
        continue
    if last == "-t" and arg.isdigit():
        last = ""
        topology = int(arg)
        continue
    if arg == "-not":
        last = "-not"
        continue
    if last == "-not":
        last = ""
        file_not = arg
        continue
    file_in = arg
    
if ( file_in == "" ):
  print('ClusterFilter.py: Missing input file. EXITING')
  exit()

# check amount of rows
count = 0
for line in open(file_in): 
    count += 1
    if ( count > 1 ):
        break 
if ( count == 0 ):
    print('ClusterFilter.py: Empty input file. EXITING')
    exit()

###############################################################################
from filter import ClusterFilter
from time import time

t = time()
print('ClusterFilter.py: Start.')
print('Loading data...')

cf = ClusterFilter(file_in, mol_file=mol_file)
total_size = cf.get_filtered_length()


# sort by
column_dict = {'el': ('log','electronic_energy'), 'g': ('log','gibbs_free_energy'),
            'elout': ('out','electronic_energy'), 'gout': ('out','gibbs_free_energy'),
            'rg': (), 'dip': ('log', 'dipole_moment')}
el,g,elout,gout = list(column_dict.values())[:4]
if sort_by in column_dict.keys():
    sort_by = column_dict[sort_by]
if (sort_by!=gout) and (sort_by not in cf.clusters_df.columns):
    print('Error! Cannot sort: '+str(sort_by)+' column missing in pickled data.')

if sort_by == gout:
    for col in [el,g,elout]:
        if col not in cf.clusters_df.columns:
            print('Error! Cannot calculate '+str(sort_by)+' column');exit()
    cf.clusters_df[gout] = cf.clusters_df[g] - cf.clusters_df[el] + cf.clusters_df[elout] 


n_mol = len(cf.mol_df)
print("Found", n_mol, "molecules in parameters.txt.")

print('Running filters...')

if isinstance(to_extract, str):
    print('  Extracting: '+to_extract)
    cf.extract_clusters(to_extract.split(','))
if isinstance(to_except, str):
    print('  Excepting: '+to_except)
    cf.except_clusters(to_except.split(','))

if isinstance(maxf, float):
    print('  Convergence: '+('<' if converged else '>')+str(maxf)+' eV/Ang')
    cf.converged(maxf, not converged)

if topology > 0 and Hbonding == None:
    Hbonding = [100, "X"]

if Hbonding == None:
    # simple filters
    cf.reacted(itol=internal_tol)
    cf.distance(minH=min_dist,max=max_dist)
else:
    # H-bond filtering
    print('  H-bonds...')
    n,H = Hbonding
    stats = cf.Hbonded(n,H,internal_tol,return_stats=print_stats)
    if print_stats:
        print(stats)

if topology:
    print('  Topologies...')
    cf.topology(topology, isomer_file, file_not, sort_by)

if relative_cutoff >= 0:
    print('  Energy Threshold: '+str(relative_cutoff)[:6]+' kcal/mol')
    relative_cutoff = relative_cutoff / 627.503 # in Hartree
    cf.cutr(relative_cutoff, sort_by)
if select > 0:
    print('  Selecting up to '+str(select))
    cf.select(select, sort_by)

t = time() - t
print('Filtered '+str(cf.get_filtered_length())+'/'+str(total_size)+' clusters')
print("Time taken:", round(t,4), "seconds")

file_out=file_in[:-4]+'_FILTERED.dat'
cf.save_to(file_out)
