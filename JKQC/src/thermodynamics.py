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

def thermodynamics(clusters_df, Qanh, Qafc, Qfc, Qt, Qdropimg):
  from numpy import array, sum, log, exp, isnan, mean, pi
  from pandas import isna
  missing = float("nan")

  h = 6.626176*10**-34 #m^2 kg s^-1
  R = 1.987 #cal/mol/K #=8.31441
  k = 1.380662*10**-23 #m^2 kg s^-2 K^-1

  if Qdropimg != 0:
    for i in range(len(clusters_df)):
      vib = clusters_df.at[i,("log","vibrational_frequencies")]
      if type(vib) == type(missing):
        continue
      clusters_df.at[i,("log","vibrational_frequencies")] = [item for item in vib if item >= 0]

  ########################################################
  # LOW VIBRATIONAL FREQUNECY ANTITREATMENT (S // G,Gc) ##
  ########################################################
  if Qafc > 0:
    for i in range(len(clusters_df)):
      try:
        lf = float(clusters_df.at[i,("log","vibrational_frequencies")][0])
      except:
        lf = 0
      if lf <= 0:
        clusters_df.at[i,("log","entropy")] = missing
        continue
      if isnan(Qt):
        try:
          Qt = clusters_df.at[i,("log","temperature")]
        except:
          Qt = 298.15
      vibs = clusters_df.at[i,("log","vibrational_frequencies")]
      structure = clusters_df.at[i,("xyz","structure")]

      mu = [float(h/(8*pi**2*2.99793*10**10*vib)) for vib in vibs]
      try:
        mi = mean(structure.get_moments_of_inertia())
        Sr = [R*(0.5+log((8*pi**2.99793*(mu[j]*mi/(mu[j]+mi))*k*Qt/h**2)**0.5)) for j in range(len(mu))]  #cal/mol/K
        Sv = [R*h*vib*2.99793*10**10/k/Qt/(exp(h*vib*2.99793*10**10/k/Qt)-1)-R*log(1-exp(-h*vib*2.99793*10**10/k/Qt)) for vib in vibs] #cal/mol/K
        w = [1/(1+(Qafc/vib)**4) for vib in vibs]
        Sv_corr = sum([w[j]*Sv[j]+(1-w[j])*Sr[j] for j in range(len(w))])
        Sv_each = sum(Sv)  #cal/mol/K
        clusters_df.at[i,("log","entropy")] = clusters_df.at[i,("log","entropy")]-(Sv_corr-Sv_each)
      except:
        mi = missing
        clusters_df.at[i,("log","entropy")] = missing
    ###
 
  #########################
  ## VIBRATIONAL SCALING ##
  #########################
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
            lf = float(clusters_df.at[i,("log","vibrational_frequencies")][0])
          except:
            lf = 0
          try: 
            if isna(clusters_df.at[i,("log","vibrational_frequencies")]).any():
              lf = 0
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
        if Qanh != "anh" and Qanh != "anh2" and Qanh != "wB97X-3c" and Qanh != "wB97X-D":
          try:
            clusters_df.at[i,("log","vibrational_frequencies")] = [float(Qanh) * j for j in clusters_df.at[i,("log","vibrational_frequencies")]]
          except:
            clusters_df.at[i,("log","vibrational_frequencies")] = [missing]
        else:
          try:
            if Qanh == "anh":
              clusters_df.at[i,("log","vibrational_frequencies")] = replace_by_nonnegative(clusters_df.at[i,("extra","anharm")],clusters_df.at[i,("log","vibrational_frequencies")],0)
            elif Qanh == "wB97X-3c":
              from numpy import piecewise
              def anh_corr(x):
                if 1==1:
                  #return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), x>=3900],[0.861409, 0.885911, 0.912831, 0.933906, 0.936296, 0.950458, 0.935704, 0.933336, 0.957378, 0.96616, 0.977959, 0.975587, 0.973703, 0.974143, 0.975156, 0.972287, 0.970533, 0.966922, 0.972591, 0.972591, 0.236802, 0.34666, 0.388443, 0.642876, 0.606038, 0.719738, 0.753388, 0.804467, 0.850579, 0.870367, 0.948011, 0.961714, 0.958506, 0.911795, 0.943901, 0.956063, 0.954999, 0.951046, 0.951589, 0.955226])*x
                  return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), (x>=3900)&(x<4000),x>4000],[0.582728, 0.6407, 0.759964, 0.825341, 0.861125, 0.900091, 0.89788, 0.925491, 0.928046, 0.951405, 0.961447, 0.966595, 0.966203, 0.966312, 0.967825, 0.962537, 0.959766, 0.957327, 0.965982, 0.48524, 0.194003, 0.323196, 0.529488, 0.508194, 0.452896, 0.716601, 0.777274, 0.806458, 0.850253, 0.938257, 0.957628, 0.95651, 0.955162, 0.918533, 0.944945, 0.946553, 0.949621, 0.946313, 0.9493, 0.954287,0.954287])*x
                if 1==0:
                  return (1.0204762161801728 - 0.000016154970938850175*x - 50.4350857894416/(239.21983191091954 + x))*x
                if x < 2000:
                  return 0.968*x
                else:
                  return 0.948*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-D":
              from numpy import piecewise
              def anh_corr(x):
                xx=3900
                return 1.0204762161801728 - 0.000016154970938850175*x - 50.4350857894416/(239.21983191091954 + x)
                print("ERROR")
                if x < xx+100 and x > xx-100:
                   return 1.0204762161801728 - 0.000016154970938850175*x - 50.4350857894416/(239.21983191091954 + x)
                   #return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), (x>=3900)&(x<4000),x>4000],[0.582728, 0.6407, 0.759964, 0.825341, 0.861125, 0.900091, 0.89788, 0.925491, 0.928046, 0.951405, 0.961447, 0.966595, 0.966203, 0.966312, 0.967825, 0.962537, 0.959766, 0.957327, 0.965982, 0.48524, 0.194003, 0.323196, 0.529488, 0.508194, 0.452896, 0.716601, 0.777274, 0.806458, 0.850253, 0.938257, 0.957628, 0.95651, 0.955162, 0.918533, 0.944945, 0.946553, 0.949621, 0.946313, 0.9493, 0.954287,0.954287])
                   #if x<2000:
                   #  return 0.950
                   #else:
                   #  return 0.943
                   #return 0.944
                else:
                   return 1 
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(j)*j for j in clusters_df.at[i,("log","vibrational_frequencies")]]
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
  
  #################################################
  #### NEW TEMPERATURE (T,S,H,Hc,U,Uc // G,Gc) ####
  #################################################
  if ~isnan(Qt):
    for i in range(len(clusters_df)):
      #try:
      #  if isna(clusters_df.at[i,("log","vibrational_frequencies")]).any():
      #    lf = 0
      #  try:
      #    if isna(clusters_df.at[i,("log","vibrational_frequencies")]):
      #      lf = 0
      #  except:
      #    lf = 0 #
      #except:
      #  lf = 0
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
  
  ####################################################
  # LOW VIBRATIONAL FREQUNECY TREATMENT (S // G,Gc) ##
  ####################################################
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
        try:
          Qt = clusters_df.at[i,("log","temperature")]
        except:
          Qt = 298.15
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
      if isnan(Qt):
        try:
          Qt = clusters_df.at[i,("log","temperature")]
        except:
          Qt = 298.15
      clusters_df.at[i,("log","gibbs_free_energy")] = clusters_df.at[i,("log","enthalpy_energy")] - clusters_df.at[i,("log","entropy")]/1000/627.503 * Qt
    except:
      clusters_df.at[i,("log","gibbs_free_energy")] = missing
    try:
      clusters_df.at[i,("log","gibbs_free_energy_thermal_correction")] = clusters_df.at[i,("log","gibbs_free_energy")] - clusters_df.at[i,("log","electronic_energy")]
    except:
      clusters_df.at[i,("log","gibbs_free_energy_thermal_correction")] = missing

  return clusters_df
