def replace_by_nonnegative(new, orig, q):
  from numpy import array
  if q == 0:
    mask = array(new) > 0
  else:
    mask = [ new[ii] > 0 and new[ii] < orig[ii] for ii in range(len(new)) ]
  orig = array(orig)
  new = array(new)
  orig[mask] = new[mask]
  return list(orig)

def thermodynamics(clusters_df, Qanh, Qfc, Qt):
  from numpy import array, sum, log, exp, isnan, mean, pi
  from pandas import isna
  missing = float("nan")

  h = 6.626176*10**-34 #m^2 kg s^-1
  R = 1.987 #cal/mol/K #=8.31441
  k = 1.380662*10**-23 #m^2 kg s^-2 K^-1

  if Qanh != "1":
    # VIBRATIONAL FREQ MODIFICATION e.g. anharmonicity (vib.freq.,ZPE,ZPEc,U,Uc,H,Hc,S // G,Gc)
    for i in range(len(clusters_df)):
      try:
        structure = clusters_df.at[i,("xyz","structure")]
        if isna(structure): 
          print("Structure is missing.")
        test = len(structure.get_atomic_numbers())
      except:
        test = 0
        lf = 0
      if test != 1:
        if test != 0:
          try:
            QtOLD = clusters_df.at[i,("log","temperature")]
          except:
            QtOLD = 298.15
          try: 
            if isna(clusters_df.at[i,("log","vibrational_frequencies")]).any():
              continue
          except:
            continue
          try:
            lf = float(clusters_df.at[i,("log","vibrational_frequencies")][0])
          except:
            lf = 0
        if lf <= 0:
          clusters_df.at[i,("log","entropy")] = missing
          clusters_df.at[i,("log","enthalpy_energy")] = missing
          clusters_df.at[i,("log","enthalpy_thermal_correction")] = missing
          clusters_df.at[i,("log","internal_energy")] = missing
          clusters_df.at[i,("log","energy_thermal_correction")] = missing
          clusters_df.at[i,("log","zero_point_correction")] = missing
          clusters_df.at[i,("log","zero_point_energy")] = missing
          continue
        try:
          Sv_OLD = sum([R*h*vib*2.99793*10**10/k/QtOLD/(exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*log(1-exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df.at[i,("log","vibrational_frequencies")]]) #cal/mol/K
          Ev_OLD = sum([R*h*vib*2.99793*10**10/k/(exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df.at[i,("log","vibrational_frequencies")]])
        except:
          Sv_OLD = missing
          Ev_OLD = missing
        #
        if Qanh != "anh" and Qanh != "anh2":
          try:
            clusters_df.at[i,("log","vibrational_frequencies")] = [float(Qanh) * j for j in clusters_df.at[i,("log","vibrational_frequencies")]]
          except:
            clusters_df.at[i,("log","vibrational_frequencies")] = [missing]
        else:
          try:
            if Qanh == "anh":
              clusters_df.at[i,("log","vibrational_frequencies")] = replace_by_nonnegative(clusters_df.at[i,("extra","anharm")],clusters_df.at[i,("log","vibrational_frequencies")],0)
            else:
              clusters_df.at[i,("log","vibrational_frequencies")] = replace_by_nonnegative(clusters_df.at[i,("extra","anharm")],clusters_df.at[i,("log","vibrational_frequencies")],1)
          except:
            clusters_df.at[i,("log","vibrational_frequencies")] = [missing]
        #
        try:
          Sv = sum([R*h*vib*2.99793*10**10/k/QtOLD/(exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*log(1-exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df.at[i,("log","vibrational_frequencies")]]) #cal/mol/K  
          Ev = sum([R*h*vib*2.99793*10**10/k/(exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df.at[i,("log","vibrational_frequencies")]])
        except:
          Sv = missing
          Ev = missing
        ###
        clusters_df.at[i,("log","zero_point_correction")] = sum([0.5*h*vib*2.99793*10**10 for vib in clusters_df.at[i,("log","vibrational_frequencies")]])*0.00038088*6.022*10**23/1000
        clusters_df.at[i,("log","zero_point_energy")] = clusters_df.at[i,("log","electronic_energy")] + clusters_df.at[i,("log","zero_point_correction")]  
        clusters_df.at[i,("log","internal_energy")] += (Ev - Ev_OLD)/1000/627.503    
        clusters_df.at[i,("log","energy_thermal_correction")] += (Ev - Ev_OLD)/1000/627.503    
        clusters_df.at[i,("log","enthalpy_energy")] += (Ev - Ev_OLD)/1000/627.503
        clusters_df.at[i,("log","enthalpy_thermal_correction")] += (Ev - Ev_OLD)/1000/627.503
        clusters_df.at[i,("log","entropy")] += Sv - Sv_OLD      
        ###
  
  # NEW TEMPERATURE (T,S,H,Hc,U,Uc // G,Gc)
  if ~isnan(Qt):
    for i in range(len(clusters_df)):
      try:
        if isna(clusters_df.at[i,("log","vibrational_frequencies")]).any():
          continue
      except:
        continue
      try:
        QtOLD = clusters_df.at[i,("log","temperature")]
        clusters_df.at[i,("log","temperature")] = Qt
      except: 
        QtOLD = 298.15
      if Qt != QtOLD:
        try:
          lf = float(clusters_df.at[i,("log","vibrational_frequencies")][0])
        except:
          lf = 0
        if lf <= 0:
          clusters_df.at[i,("log","temperature")] = Qt
          clusters_df.at[i,("log","entropy")] = missing
          clusters_df.at[i,("log","enthalpy_energy")] = missing
          clusters_df.at[i,("log","enthalpy_thermal_correction")] = missing
          clusters_df.at[i,("log","internal_energy")] = missing
          clusters_df.at[i,("log","energy_thermal_correction")] = missing
          continue
        Sv_OLD = sum([R*h*vib*2.99793*10**10/k/QtOLD/(exp(h*vib*2.99793*10**10/k/QtOLD)-1)-R*log(1-exp(-h*vib*2.99793*10**10/k/QtOLD)) for vib in clusters_df.at[i,("log","vibrational_frequencies")]]) #cal/mol/K
        Sv = sum([R*h*vib*2.99793*10**10/k/Qt/(exp(h*vib*2.99793*10**10/k/Qt)-1)-R*log(1-exp(-h*vib*2.99793*10**10/k/Qt)) for vib in clusters_df.at[i,("log","vibrational_frequencies")]]) #cal/mol/K
        Ev_OLD = sum([R*h*vib*2.99793*10**10/k/(exp(h*vib*2.99793*10**10/k/QtOLD)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df.at[i,("log","vibrational_frequencies")]])
        Ev = sum([R*h*vib*2.99793*10**10/k/(exp(h*vib*2.99793*10**10/k/Qt)-1)+R*h*vib*2.99793*10**10/k*0.5 for vib in clusters_df.at[i,("log","vibrational_frequencies")]])
        ###
        clusters_df.at[i,("log","temperature")] = Qt
        clusters_df.at[i,("log","entropy")] += Sv - Sv_OLD + 4*R*log(Qt/QtOLD)
        clusters_df.at[i,("log","enthalpy_energy")] += (Ev - Ev_OLD + 4*R*(Qt-QtOLD))/1000/627.503
        clusters_df.at[i,("log","enthalpy_thermal_correction")] += (Ev - Ev_OLD + 4*R*(Qt-QtOLD))/1000/627.503
        clusters_df.at[i,("log","internal_energy")] += (Ev - Ev_OLD + 3*R*(Qt-QtOLD))/1000/627.503
        clusters_df.at[i,("log","energy_thermal_correction")] += (Ev - Ev_OLD + 3*R*(Qt-QtOLD))/1000/627.503
        ###
  
  # LOW VIBRATIONAL FREQUNECY TREATMENT (S // G,Gc)
  if Qfc > 0:
    for i in range(len(clusters_df)):
      try:
        if isna(clusters_df.at[i,("log","vibrational_frequencies")]):
          continue
      except:
        lf = 0
      try:
        lf = float(clusters_df.at[i,("log","vibrational_frequencies")][0])
      except:
        lf = 0
      if lf <= 0:
        clusters_df.at[i,("log","entropy")] = missing
        continue
      if isnan(Qt):
        Qt = clusters_df.at[i,("log","temperature")]
      vibs = clusters_df.at[i,("log","vibrational_frequencies")]
      structure = clusters_df.at[i,("xyz","structure")]
      
      mu = [float(h/(8*pi**2*2.99793*10**10*vib)) for vib in vibs]
      try:
        mi = mean(structure.get_moments_of_inertia())
        Sr = [R*(0.5+log((8*pi**2.99793*(mu[j]*mi/(mu[j]+mi))*k*Qt/h**2)**0.5)) for j in range(len(mu))]  #cal/mol/K
        Sv = [R*h*vib*2.99793*10**10/k/Qt/(exp(h*vib*2.99793*10**10/k/Qt)-1)-R*log(1-exp(-h*vib*2.99793*10**10/k/Qt)) for vib in vibs] #cal/mol/K
        w = [1/(1+(Qfc/vib)**4) for vib in vibs]
        Sv_corr = sum([w[j]*Sv[j]+(1-w[j])*Sr[j] for j in range(len(w))])
        Sv_each = sum(Sv)  #cal/mol/K
        clusters_df.at[i,("log","entropy")] = clusters_df.at[i,("log","entropy")]+(Sv_corr-Sv_each)
      except:
        mi = missing
        clusters_df.at[i,("log","entropy")] = missing
    ###

  ## CORRECTIONS FOR GIBBS FREE ENERGY
  for i in range(len(clusters_df)):
    try:
      clusters_df.at[i,("log","gibbs_free_energy")] = clusters_df.at[i,("log","enthalpy_energy")] - clusters_df.at[i,("log","entropy")]/1000/627.503 * clusters_df.at[i,("log","temperature")]
    except:
      clusters_df.at[i,("log","gibbs_free_energy")] = missing
    try:
      clusters_df.at[i,("log","gibbs_free_energy_thermal_correction")] = clusters_df.at[i,("log","gibbs_free_energy")] - clusters_df.at[i,("log","electronic_energy")]
    except:
      clusters_df.at[i,("log","gibbs_free_energy_thermal_correction")] = missing

  return clusters_df
