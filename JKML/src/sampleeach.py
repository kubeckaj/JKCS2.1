def compare_mbtr(structa, structb):
  '''
  Expects two xyz pickle structures
  Returns a float equal the sum of the average mean square deviation
  between each bond type for the twombtr inputs.
  Is 0 if equal, lower values equal similiarty.
  '''
  from numpy import sum, square
  return sum(square(structa-structb))


def task_train(arg):
  #mbtr_train[arg] = cr_mbtr_k2(train_database[arg])
  #return cr_mbtr_k2(train_database[arg])
  return cr_mbtr_k2(arg)

def task_test(arg):
  #mbtr_test[arg] = cr_mbtr_k2(test_database[arg])
  #return cr_mbtr_k2(test_database[arg])
  return cr_mbtr_k2(arg)

def flatten(l):
  return [item for sublist in l for item in sublist]

def call_mbtr(chemsyms_uniques):
  from dscribe.descriptors import MBTR
  k2min, k2max, k2n = 0.7, 2.0, 100
  if 1==1:
    import warnings
    warnings.filterwarnings(
        "ignore", ".*Please use atoms.calc.*"
    )
  return MBTR(
        species=chemsyms_uniques,
        #ATOMIC NUMBER
        #k1={
        #    "geometry": {"function": "atomic_number"},
        #    "grid": {"min": 0, "max": 16, "n": 100, "sigma": 0.1},
        #},
        #DISTANCE OR INVERSE DISTANCE
        #k2={
        geometry={"function": "distance"},
        grid={"min": k2min, "max": k2max, "n": k2n, "sigma": 0.000000001},
        weighting={"function": "exp", "scale": 0.5, "threshold": 3e-3},
        #"weighting": {"function": "unity"},#, "scale": 0.5, "threshold": 3e-3},
        #},
        #ANGLE OR COSINE(ANGLE)
        #TODO to be tested as well
        #k3={
        #    "geometry": {"function": "cosine"},
        #    "grid": {"min": -1, "max": 1, "n": 100, "sigma": 0.000000001},
        #    "weighting": {"function": "exp", "scale": 0.5, "threshold": 3e-3},
        #    #"weighting": {"function": "unity"},# "scale": 0.5, "threshold": 3e-3},
        #},
        periodic=False,
        normalization="l2",
        #flatten=False  # without this it's gonna be a dict and then should find out which part corresponds to which stuff
  )

def sampleeach_mbtr(train_high_database, test_high_database):
  from joblib import Parallel, delayed
  import multiprocessing
  from dscribe.descriptors import MBTR


  def cr_mbtr_k2(struct):
    '''
    Expects xyz pickle structure
    Returns the mbtr for bond lenghts. 
    '''
    import warnings
    warnings.filterwarnings("ignore", ".*Please use atoms.calc.*")
    return mbtr.create(struct)#["k2"][:][:]

  def task_train(arg):
    #mbtr_train[arg] = cr_mbtr_k2(train_database[arg])
    #return cr_mbtr_k2(train_database[arg])
    return cr_mbtr_k2(arg)
  
  def task_test(arg):
    #mbtr_test[arg] = cr_mbtr_k2(test_database[arg])
    #return cr_mbtr_k2(test_database[arg])
    return cr_mbtr_k2(arg)

  print("This likely does not work anymore !!!!!!!!", flush = True)
  chemsyms_uniques = list(set(flatten([i.get_chemical_symbols() for i in train_high_database["xyz"]["structure"]])))
  #max_atoms=max([len(i.get_chemical_symbols()) for i in pd.read_pickle(TRAIN_HIGH).sort_values([('info','file_basename')])["xyz"]["structure"]])
  mbtr = call_mbtr(chemsyms_uniques)
  print("MBTR must be calculated (just once) for both train and test datasets", flush = True)

  num_cores = multiprocessing.cpu_count()
  print("Trying to use "+str(num_cores)+" CPUs for MBTR. (If less are available I hope that nothing gets fucked up.)")

  mbtr_train = Parallel(n_jobs=num_cores)(delayed(task_train)(i) for i in train_high_database["xyz"]["structure"])
  mbtr_test = Parallel(n_jobs=num_cores)(delayed(task_train)(i) for i in test_high_database["xyz"]["structure"])
  print("MBTR done", flush = True)
  sampleeach_all = range(len(mbtr_test))
  
  return sampleeach_all, mbtr_train, mbtr_test

def sampleeach_fchl(train_high_database, test_high_database, krr_cutoff):

  #SAMLEEACH = SIMILARITY based on FCHL 
  from joblib import Parallel, delayed
  import multiprocessing
  print("This likely does not work anymore !!!!!!!!", flush = True)
  def task(arg):
    return generate_representation(arg.get_positions(),arg.get_atomic_numbers(),max_size = max_atoms, neighbors = max_atoms, cut_distance=krr_cutoff)

  num_cores = 2 #multiprocessing.cpu_count()
  print("Trying to use "+str(num_cores)+" CPUs for FCHL. (If less are available I hope that nothing gets fucked up.)")
  max_atoms = max([len(i.get_atomic_numbers()) for i in train_high_database["xyz"]["structure"]])
  max_atoms2 = max([len(i.get_atomic_numbers()) for i in test_high_database["xyz"]["structure"]])
  if max_atoms2 > max_atoms:
    max_atoms = max_atoms2
  fchl_train = Parallel(n_jobs=num_cores)(delayed(task)(i) for i in train_high_database["xyz"]["structure"])
  fchl_test = Parallel(n_jobs=num_cores)(delayed(task)(i) for i in test_high_database["xyz"]["structure"])
  sampleeach_all = range(len(fchl_test))
  print("FCHL done", flush = True)
  return sampleeach_all, fchl_train, fchl_test

def sampledist_mbtr(sample_train, sample_test, sampleeach_i, Qsampleeach):
  from numpy import array, shape
  dist = array([compare_mbtr(sample_train[i],sample_test[sampleeach_i]) for i in range(shape(sample_train)[0])])
  sampledist = dist.argsort()[:Qsampleeach]
  return sampledist

def sampledist_fchl(sample_train, sample_test, sampleeach_i, Qsampleeach, Qkernel, kernel_args, train_high_database):
  if Qkernel == "Gaussian":
    from qmllib.representations.fchl import get_local_symmetric_kernels as JKML_sym_kernel
    from qmllib.representations.fchl import get_local_kernels as JKML_kernel
  else:
    from qmllib.representations.fchl import laplacian_kernel_symmetric as JKML_sym_kernel
    from qmllib.representations.fchl import laplacian_kernel as JKML_kernel
  from numpy import array, shape, sqrt

  simil = JKML_kernel(array([m for m in [sample_test[sampleeach_i]]]), array([m for m in sample_train]), [0.001], **kernel_args)[0][0]
  simil = [ simil[i]/sqrt(JKML_kernel(array([m for m in [sample_train[i]]]),array([m for m in [sample_train[i]]]))[0][0][0]) for i in range(len(train_high_database))]
  dist = array([-m for m in simil])
  sampledist = dist.argsort()[:-Qsampleeach]
  return sampledist

