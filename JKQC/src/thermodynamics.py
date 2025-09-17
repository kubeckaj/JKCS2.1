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
        if str(Qanh) not in {"anh","anh2","B97-3c_1","B97-3c_2","B97-3c_sf","B97-3c_mult","r2SCAN-3c_1","r2SCAN-3c_2","r2SCAN-3c_sf","r2SCAN-3c_mult","wb97-3c_1","wb97-3c_2","wb97-3c_sf","wb97-3c_mult","wB97X-D_1", "wB97X-D_2","wB97X-D_sf","wB97X-D_mult"}:
          try:
            clusters_df.at[i,("log","vibrational_frequencies")] = [float(Qanh) * j for j in clusters_df.at[i,("log","vibrational_frequencies")]]
          except:
            clusters_df.at[i,("log","vibrational_frequencies")] = [missing]
        else:
          try:
            if Qanh == "anh":
              clusters_df.at[i,("log","vibrational_frequencies")] = replace_by_nonnegative(clusters_df.at[i,("extra","anharm")],clusters_df.at[i,("log","vibrational_frequencies")],0)
            ########################################################
            elif Qanh == "B97-3c_1":
              def anh_corr(x):
                return 0.944*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "B97-3c_2":
              def anh_corr(x):
                if x < 2000:
                  return 0.967*x
                else:
                  return 0.937*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "B97-3c_sf":
              def anh_corr(x):
                return (0.969507 - 8.55527*10**-6*x + 1.99602/(-0.0447777 + x))*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "B97-3c_mult":
              from numpy import piecewise
              def anh_corr(x):
                scaling = [1.00058, 0.999846, 0.935579, 0.967753, 0.964498, 0.964641, 0.948904, 0.96143, 0.962068, 0.965688, 0.972093, 0.971545, 0.971685, 0.973273, 0.971551, 0.969923, 0.968905, 0.962642, 0.974881, 0.326225, 0.325892, 0.41384, 0.441502, 0.526397, 0.687563, 0.758214, 0.803279, 0.806293, 0.844447, 0.941828, 0.961305, 0.95736, 0.955739, 0.919227, 0.947026, 0.951936, 0.950952, 0.945497, 0.951473, 0.954571, 0.954571]
                return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), (x>=3900)&(x<4000),x>4000],scaling)*x
               
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            ########################################################
            elif Qanh == "r2SCAN-3c_1":
              def anh_corr(x):
                return 0.950*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "r2SCAN_2":
              def anh_corr(x):
                if x < 2000:
                  return 0.969*x
                else:
                  return 0.944*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "r2SCAN_sf":
              def anh_corr(x):
                return (1.0211 - 1.59745*10**-5*x - 102.482/(1517.15 + x))*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "r2SCAN_mult":
              from numpy import piecewise
              def anh_corr(x):
                scaling = [0.993869, 0.929885, 0.940778, 0.962749, 0.970478, 0.967294, 0.953732, 0.964591, 0.966077, 0.973015, 0.976342, 0.972922, 0.972205, 0.974236, 0.973522, 0.97198, 0.972194, 0.970429, 0.946279, 0.450775, 0.429606, 0.533808, 0.593157, 0.657998, 0.747984, 0.78632, 0.802993, 0.852336, 0.879202, 0.947577, 0.961449, 0.954227, 0.91328, 0.925493, 0.956013, 0.953386, 0.949111, 0.951456, 0.954239, 0.956302, 0.956302]
                return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), (x>=3900)&(x<4000),x>4000],scaling)*x

              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            ########################################################
            elif Qanh == "wB97X-3c_1":
              def anh_corr(x):
                return 0.954*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-3c_2":
              def anh_corr(x):
                if x < 2000:
                  return 0.971*x
                else:
                  return 0.949*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-3c_sf":
              def anh_corr(x):
                return (1.07555 - 2.42816*10**-5*x - 188.46/(1179.03 + x))*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-3c_mult":
              from numpy import piecewise
              def anh_corr(x):
                scaling = [0.938768, 0.905829, 0.921018, 0.942186, 0.959706, 0.964414, 0.950646, 0.955143, 0.964944, 0.968226, 0.979111, 0.976256, 0.974098, 0.97433, 0.975238, 0.973679, 0.970807, 0.968801, 0.974567, 0.978848, 0.30976, 0.342098, 0.420479, 0.553954, 0.693097, 0.733632, 0.76931, 0.804974, 0.83931, 0.876033, 0.950491, 0.962892, 0.958473, 0.915792, 0.946232, 0.956013, 0.955388, 0.9524, 0.952433, 0.955203, 0.955203]
                return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), (x>=3900)&(x<4000),x>4000],scaling)*x

              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            ########################################################
            elif Qanh == "wB97X-D_1":
              def anh_corr(x):
                return 0.950*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-D_2":
              def anh_corr(x):
                if x < 2000:
                  return 0.967*x
                else:
                  return 0.945*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-D_sf":
              def anh_corr(x):
                return (0.961752 - 3.89697*10**-6*x + 2.814/(5.5432 + x))*x
              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            elif Qanh == "wB97X-D_mult":
              from numpy import piecewise
              def anh_corr(x):
                scaling = [1.00058, 0.999846, 0.935579, 0.967753, 0.964498, 0.964641, 0.948904, 0.96143, 0.962068, 0.965688, 0.972093, 0.971545, 0.971685, 0.973273, 0.971551, 0.969923, 0.968905, 0.962642, 0.974881, 0.326225, 0.325892, 0.41384, 0.441502, 0.526397, 0.687563, 0.758214, 0.803279, 0.806293, 0.844447, 0.941828, 0.961305, 0.95736, 0.955739, 0.919227, 0.947026, 0.951936, 0.950952, 0.945497, 0.951473, 0.954571, 0.954571]
                return piecewise(x,[(x>=0)&(x<100), (x>=100)&(x<200), (x>=200)&(x<300), (x>=300)&(x<400), (x>=400)&(x<500), (x>=500)&(x<600), (x>=600)&(x<700), (x>=700)&(x<800), (x>=800)&(x<900), (x>=900)&(x<1000), (x>=1000)&(x<1100), (x>=1100)&(x<1200), (x>=1200)&(x<1300), (x>=1300)&(x<1400), (x>=1400)&(x<1500), (x>=1500)&(x<1600), (x>=1600)&(x<1700), (x>=1700)&(x<1800), (x>=1800)&(x<1900), (x>=1900)&(x<2000), (x>=2000)&(x<2100), (x>=2100)&(x<2200), (x>=2200)&(x<2300), (x>=2300)&(x<2400), (x>=2400)&(x<2500), (x>=2500)&(x<2600), (x>=2600)&(x<2700), (x>=2700)&(x<2800), (x>=2800)&(x<2900), (x>=2900)&(x<3000), (x>=3000)&(x<3100), (x>=3100)&(x<3200), (x>=3200)&(x<3300), (x>=3300)&(x<3400), (x>=3400)&(x<3500), (x>=3500)&(x<3600), (x>=3600)&(x<3700), (x>=3700)&(x<3800), (x>=3800)&(x<3900), (x>=3900)&(x<4000),x>4000],scaling)*x

              clusters_df.at[i,("log","vibrational_frequencies")] = [anh_corr(float(j)) for j in clusters_df.at[i,("log","vibrational_frequencies")]]
            ########################################################
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
