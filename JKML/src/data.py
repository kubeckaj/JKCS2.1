def substract_monomers(the_clusters_df, the_monomers_df, Qmonomers, column_name_1, column_name_2, Qifeldisp):
    if Qmonomers == 2:  # no difference from monomers at all
        the_ens_correction = [0] * len(the_clusters_df)
    else:
        from numpy import sum, transpose
        if Qmonomers == 1:
            clusters_df0 = the_monomers_df
            ins_monomers = [sum(clusters_df0["info"]["component_ratio"].values[i]) == 1 for i in range(len(clusters_df0))]
            monomers_df = clusters_df0[ins_monomers]
        else:
            ins_monomers = [sum(the_clusters_df["info"]["component_ratio"].values[i]) == 1 for i in range(len(the_clusters_df))]
            monomers_df = the_clusters_df[ins_monomers]
        the_ens_correction = [0] * len(the_clusters_df)
        n = 0
        for i in transpose([the_clusters_df["info"]["components"].values, the_clusters_df["info"]["component_ratio"].values]):
            nn = 0
            if type(i[0]) == type([]):
                all_mons = i[0]
                all_mons_ratios = i[1]
            else:
                all_mons = [i[0]]
                all_mons_ratios = [i[1]]
            for j in all_mons:
                test = 0
                for k in range(len(monomers_df)):
                    monomer_k = monomers_df.iloc[k]
                    monomer_k_name = monomer_k["info"]["components"]
                    if type(monomer_k_name) == type([]):
                        monomer_k_name = monomer_k_name[0]
                    if j == monomer_k_name:
                        mon_energy = monomer_k[column_name_1][column_name_2]
                        #print("here", flush=True)
                        if Qifeldisp == 1 and ("extra", "dispersion_electronic_energy") in monomer_k and monomer_k[("extra", "dispersion_electronic_energy")] is not None:
                            #print(monomer_k)
                            mon_energy -= monomer_k[("extra", "dispersion_electronic_energy")]
                        the_ens_correction[n] += mon_energy * all_mons_ratios[nn]
                        test = 1
                        break
                if test == 0:
                    print("JKML(data): OMG; monomer " + j + " was not found in:")
                    print(monomers_df["info"]["components"].values)
                    exit()
                nn += 1
            n += 1
    return the_ens_correction

def prepare_data_for_training(train_high_database, monomers_high_database, train_low_database, monomers_low_database, seed, size, method, column_name_1, column_name_2, Qmin, Qifforces, Qifcharges, Qifdipole, Qsampleeach, Qforcemonomers, sampledist, Qmonomers, Qifeldisp, Qifforcedisp):

  from sklearn.model_selection import train_test_split    

  ## The high level of theory
  clusters_df = train_high_database.sample(frac = 1, random_state = seed)
  #clusters_df = train_high_database
  if Qsampleeach != 0:
    clusters_df = clusters_df.iloc[sampledist]
  if Qforcemonomers == 1:
    clusters_df0 = monomers_high_database
    from pandas import concat
    clusters_df = concat([clusters_df, clusters_df0.copy()], ignore_index=True)
    #clusters_df = clusters_df.append(clusters_df0, ignore_index=True)
  ## The low level of theory
  if method == "delta":
    clusters_df2 = train_low_database.sample(frac = 1, random_state = seed)
    if Qsampleeach != 0:
      clusters_df2 = clusters_df2.iloc[sampledist]
    if Qforcemonomers == 1:
      clusters_df0l = monomers_low_database
      from pandas import concat
      clusters_df2 = concat([clusters_df2, clusters_df0l.copy()], ignore_index=True)
      #clusters_df2 = clusters_df2.append(clusters_df0l, ignore_index=True)

  #Do we take only subset for training?
  if size != "full":
    if len(clusters_df) <= int(size):
      size = "full"
  if size != "full":
    clusters_df, clusters_df_trash, idx, idx_trash = train_test_split(clusters_df, range(len(clusters_df)), test_size=(len(clusters_df)-size)/len(clusters_df), random_state=seed)
    if method == "delta":
      clusters_df2 = clusters_df2.iloc[idx]
    size = "full" #IT IS BECAUSE I DO NOT WANT TO MAKE MY TEST SET SMALLER

  ### ENERGIES = VARIABLES / STRUCTURES
  #TODO print names
  #print(clusters_df["info"]["file_basename"].values)
  ens = (clusters_df[column_name_1][column_name_2]).values.astype("float")
  strs = clusters_df["xyz"]["structure"]
  print("JKML(data): data length = "+str(ens.shape[0]), flush = True)
  if method == "delta":
    ens2 = (clusters_df2[column_name_1][column_name_2]).values.astype("float")
    print("JKML: low data length = "+str(ens2.shape[0]), flush = True)
    #str2 should be the same as str by principle
    if ens.shape != ens2.shape:
      print("JKML(data): The LOW and HIGH method train sizes do not match. [EXITING]")
      exit()
  if method == "min":
    ens -= Qmin

  ### NUMBERS
  from ase import Atoms
  from ase.io import read
  Z_atoms = []
  N_atoms = []
  
  for s in strs:
    atoms = s if isinstance(s, Atoms) else read(s)
    atomic_numbers = atoms.get_atomic_numbers()
    Z_atoms.append(atomic_numbers)
    N_atoms.append(len(atomic_numbers))

  ### BASENAME
  file_basenames = clusters_df["info"]["file_basename"].values

  ### FORCES
  if ("extra","forces") in clusters_df.columns and Qifforces == 1:
    if method == "delta":
      from numpy import array
      F1 = array([array(i) for i in clusters_df["extra"]["forces"].values], dtype=object)
      F2 = array([array(i) for i in clusters_df2["extra"]["forces"].values], dtype=object)
      F_train = F1 - F2
      #F_train = array([i.tolist() for i in F_train])
    else:
      F_train = clusters_df["extra"]["forces"].values
    Qforces = 1
  else:
    F_train = float("nan")
    Qforces = 0

  ### DISPERSION CORRECTIONS
  if ("extra", "dispersion_electronic_energy") in clusters_df.columns and Qifeldisp == 1:
    disp_energies = clusters_df["extra"]["dispersion_electronic_energy"].values.astype(float)
    ens -= disp_energies  
  if ("extra", "dispersion_forces") in clusters_df.columns and Qifforcedisp == 1:
    from numpy import array
    disp_forces = clusters_df["extra"]["dispersion_forces"].values
    F_train = [array(f) - array(df) for f, df in zip(F_train, disp_forces)]

  ### CHARGE
  if ("log","charge") in clusters_df.columns and ("log","mulliken_charges") in clusters_df.columns and Qifcharges == 1:
    Q_charge = clusters_df["log"]["charge"].values
    Q_charges = clusters_df["log"]["mulliken_charges"].values
    Qcharge = 1
  elif ("log","charge") in clusters_df.columns and Qifcharges == 1:
    Q_charge = clusters_df["log"]["charge"].values
    Q_charges = None
    Qcharge = 1
  elif ("log","mulliken_charges") in clusters_df.columns and Qifcharges == 1:
    from numpy import array
    Q_charge = array([0]*len(clusters_df))
    Q_charges = clusters_df["log"]["mulliken_charges"].values
    Qcharge = 1
  elif Qifcharges == 1:
    from numpy import array
    Q_charge = array([0]*len(clusters_df))
    Q_charges = None
    Qcharge = 1
  else:
    Q_charges = None
    Q_charge = None
    Qcharge = 0

  ### DIPOLE MOMENT
  if ("log","dipole_moments") in clusters_df.columns and Qifdipole == 1:
    D_dipole = clusters_df["log"]["dipole_moments"].values
    Qdipole = 1
  elif Qifdipole == 1:
    print("JKML(data): Dipoles are fucked up.")
    exit()
  else:
    D_dipole = None
    Qdipole = 0

  ### BINDING PROPERTIES CALCULATION (i.e. relative to monomers) ###
  #HIGH LEVEL
  ens_correction = substract_monomers(clusters_df,monomers_high_database,Qmonomers,column_name_1,column_name_2,Qifeldisp)
  #LOW LEVEL
  if method == "delta":
    ens2_correction = substract_monomers(clusters_df2,monomers_low_database,Qmonomers,column_name_1,column_name_2,Qifeldisp)

  #The binding (formation) energy calculation (or final property)
  form_ens = ens - ens_correction
  #print(form_ens, flush = True)
  if method == "delta":
    form_ens2 = ens2 - ens2_correction
    #print(form_ens2, flush = True)
  if method == "delta":
    Y_train = form_ens - form_ens2
  else:
    Y_train = form_ens

  return locals()
  #return strs, Y_train, F_train, Qforces, Q_charge, Q_charges, Qcharge, D_dipole, Qdipole, size

###########################################################################################
###########################################################################################
###########################################################################################

def prepare_data_for_testing(test_high_database,test_low_database,monomers_high_database,monomers_low_database,Qsampleeach,sampleeach_i,method,size,seed,column_name_1,column_name_2,Qeval,Qifforces,Qmonomers,Qmin,Qprintforces,Qifcharges,Qifeldisp):

  from sklearn.model_selection import train_test_split

  ### DATABASE LOADING ###
  ## The high level of theory
  clusters_df = test_high_database
  if Qsampleeach != 0:
    clusters_df = clusters_df.iloc[[sampleeach_i]]
  if method == "delta":
    clusters_df2 = test_low_database
    if Qsampleeach != 0:
      clusters_df2 = clusters_df2.iloc[[sampleeach_i]]

  #Do we take only subset for testing?
  if size != "full":
    if len(clusters_df) <= int(size):
      size = "full"
  if size != "full":
    clusters_df_trash, clusters_df, idx_trash, idx = train_test_split(clusters_df, range(len(clusters_df)), test_size=size/len(clusters_df), random_state=seed)
    if method == "delta":
      clusters_df2 = clusters_df2.iloc[idx]
  clustersout_df = clusters_df.copy()

  if Qeval == 2:
    try:
      ens = (clusters_df[column_name_1][column_name_2]).values.astype("float")
      if method == "min":
        ens -= Qmin
      #Qeval = 2 #2=Compare ML prediction with QC
    except:
      Qeval = 1 #1=Only predicts
      ens = None
  else:
    ens = None
  if method == "delta":
    ens2 = (clusters_df2[column_name_1][column_name_2]).values.astype("float")
  else:
    ens2 = None
  strs = clusters_df["xyz"]["structure"]
  #str2 should be the same as str by princip

  ### FORCES
  if ("extra","forces") in clusters_df.columns and Qifforces == 1:
    if method == "delta":
      from numpy import array
      F1 = array([array(i) for i in clusters_df["extra"]["forces"].values])
      F2 = array([array(i) for i in clusters_df2["extra"]["forces"].values])
      F_test = F1 - F2
    else:
      F_test = clusters_df["extra"]["forces"].values
    Qforces = 1
  elif Qifforces == 1:
    F_test = None
    Qforces = 1
  else:
    F_test = None
    Qforces = 0

  ### NUMBERS
  from ase import Atoms
  from ase.io import read
  Z_atoms = []
  N_atoms = []

  for s in strs:
    atoms = s if isinstance(s, Atoms) else read(s)
    atomic_numbers = atoms.get_atomic_numbers()
    Z_atoms.append(atomic_numbers)
    N_atoms.append(len(atomic_numbers))

  ### CHARGE
  if ("log","charge") in clusters_df.columns and ("log","mulliken_charges") in clusters_df.columns and Qifcharges == 1:
    Q_charge = clusters_df["log"]["charge"].values
    Q_charges = clusters_df["log"]["mulliken_charges"].values
    Qcharge = 1
  elif ("log","charge") in clusters_df.columns and Qifcharges == 1:
    Q_charge = clusters_df["log"]["charge"].values
    Q_charges = None
    Qcharge = 0
  elif ("log","mulliken_charges") in clusters_df.columns and Qifcharges == 1:
    from numpy import array
    Q_charge = array([0]*len(clusters_df))
    Q_charges = clusters_df["log"]["mulliken_charges"].values
    Qcharge = 1
  elif Qifcharges == 1:
    from numpy import array
    Q_charge = array([0]*len(clusters_df))
    Q_charges = None
    Qcharge = 0
  else:
    Q_charges = None
    Q_charge = None
    Qcharge = 0

  if Qeval == 2:
    print("JKML(data): data length = "+str(ens.shape), flush = True)
  if method == "delta":
    print("JKML(data): low data length = "+str(ens2.shape), flush = True)
    if Qeval == 2:
      if ens.shape != ens2.shape:
        print("JKML(data): The LOW and HIGH method test sizes do not match. [EXITING]", flush = True)
        exit()

  ### BINDING PROPERTIES CALCULATION (i.e. relative to monomers) ###
  #HIGH LEVEL
  ens_correction = substract_monomers(clusters_df,monomers_high_database,Qmonomers,column_name_1,column_name_2,Qifeldisp)
  #LOW LEVEL
  if method == "delta":
    ens2_correction = substract_monomers(clusters_df2,monomers_low_database,Qmonomers,column_name_1,column_name_2,Qifeldisp)
  else:
    ens2_correction = None
  #TODO clustername given as argument
  #monomers_df = pd.read_pickle("/home/kubeckaj/ML_SA_B/DATABASES/database_XTB_monomers_at_DFT.pkl")
  #import re
  #clustername=np.array(re.split(r'(\d+)', sys.argv[3])[1:])
  #clustername=clustername.reshape(int(len(clustername)/2),2)
  #clustername=[clustername[:,(1,0)]]*len(clusters_df)
  #print(clustername)
  #print(clustername[0][0])
  #for i in clustername: #np.transpose([clusters_df["info"]["components"].values,clusters_df["info"]["component_ratio"].values]):

  #The binding (formation) energy calculation
  if Qeval == 2:
    form_ens = ens - ens_correction
    #print(form_ens, flush = True)
    if method == "direct":
      Y_validation = form_ens
      form_ens2 = None
    elif method == "delta":
      form_ens2 = ens2 - ens2_correction
      Y_validation = form_ens - form_ens2
      #print(form_ens2, flush = True)
    else:
      Y_validation = form_ens - Qmin
  else:
    form_ens = None
    form_ens2 = None
    Y_validation = None
   
  return locals()
