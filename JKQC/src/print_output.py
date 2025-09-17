def print_output(clusters_df, Qoutpkl, input_pkl, output_pkl, Qsplit, Qclustername, Qt, Qcolumn, Qbonded, Pout = [], QUenergy = 1, QUentropy = 1):
  """Print output from JKQC
  clusters_df = Pandas Dataframe
  Pout = list of outputs
  """
  missing = float("nan")
 
  Qbonded_index = -1 
  Qcolumn_i = 0
  output = []
  last = ''
  for i in Pout:
    #INFO
    if i == "-info":
      print(clusters_df.info())
      continue
    #LEVELS
    if i == "-levels":
      from pandas import set_option, isna
      try:
        set_option('display.max_colwidth', None)
        if not isna(clusters_df["log"]["program"]).all():
          print("# LOG #")
          print(clusters_df["log"][["program","method"]].drop_duplicates())
          print("#######")   
        if "out" in clusters_df:
          if not isna(clusters_df["out"]["program"]).all():
            print("# OUT #")
            print(clusters_df["out"][["program","method"]].drop_duplicates())
            print("#######")  
      except:
        print("For some reason, I cannot print programs/methods.") 
      continue
    #EXTRA
    if i == "-extra":
      last = "-extra"
      continue
    if last == "-extra":
      last = ""
      try:
        output.append(clusters_df["extra"][i].values)
      except:
        output.append([missing]*len(clusters_df))
      continue
    #CITE
    if i == "-cite":
      try:
        output.append(clusters_df["info"]["citation"].values)
      except:
        output.append([missing]*len(clusters_df))
      continue
    #ID 
    if i == "-id1":
      try:
        output.append(clusters_df["xyz"]["id1"].values)
      except:
        output.append([missing]*len(clusters_df))
      continue
    #XYZ
    if i == "-xyz":
      from ase.io import write
      for ind in clusters_df.index:
        try:
          write(clusters_df.loc[ind,("info","file_basename")]+".xyz",clusters_df.loc[ind,("xyz","structure")], format='xyz')
        except:
          print("Corrupted structure saved for "+clusters_df["info"]["file_basename"][ind])
      continue
    #PDB IMOS
    if i == "-imos":
      from ase.io.proteindatabank import write_proteindatabank
      from os import remove
      from pandas import isna
      for ind in clusters_df.index:
        if isna(clusters_df.loc[ind,("log","esp_charges")]):
          print("Missing esp charges for "+clusters_df.loc[ind,("info","file_basename")])
          continue
        write_proteindatabank(".JKQChelp.pdb",clusters_df.loc[ind,("xyz","structure")])
        f1=open(".JKQChelp.pdb", "r")
        f2=open(clusters_df.loc[ind,("info","file_basename")]+".pdb", "w")
        f2.write(f1.readline())
        for j in range(len(clusters_df.loc[ind,("xyz","structure")].get_chemical_symbols())):
          line=f1.readline()
          ch4="{:+.2f}".format(clusters_df.loc[ind,("log","esp_charges")][j])
          f2.write(line[:56]+ch4+line[61:])
        f2.write(f1.readline())
        f1.close()
        f2.close()
        remove(".JKQChelp.pdb")
      continue
    #XLSX IMOS
    if i == "-imos_xlsx":
      from ase.data.vdw import vdw_radii
      import xlsxwriter
      from pandas import isna
      from numpy import sum
      workbook = xlsxwriter.Workbook('imos.xlsx')
      bold = workbook.add_format({'bold': True,'fg_color': '#FFFF00', 'border':1})
      for ind in clusters_df.index:
        if ("extra","esp_charges") in clusters_df.columns:
          Qcol1 = "extra"
          Qcol2 = "esp_charges"
        elif ("extra","chelpg") in clusters_df.columns:
          Qcol1 = "extra"
          Qcol2 = "chelpg"
        elif ("log","esp_charges") in clusters_df.columns: 
          if not isna(clusters_df.loc[ind,("log","esp_charges")]):
            Qcol1 = "log"
            Qcol2 = "esp_charges"
          else:
            print("Missing esp charges for "+clusters_df.loc[ind,("info","file_basename")])
            continue
        else:
          print("Missing eps_charges")
          continue
        clustername = clusters_df["info"]["file_basename"][ind]
        if len(clustername) > 30:
          clustername = clustername[:20]+'-TBC'+str(ind)
        worksheet = workbook.add_worksheet(clustername)
        
        pos=clusters_df["xyz"]["structure"][ind].get_positions()
        mass=clusters_df["xyz"]["structure"][ind].get_masses()
        at=clusters_df["xyz"]["structure"][ind].get_atomic_numbers()
        ch=clusters_df[Qcol1][Qcol2][ind]
        if type(ch) == str:
          import ast
          ch = list(map(float, ast.literal_eval(ch)))
        row = 0
        for j in range(len(pos)):
          worksheet.write(row, 0, pos[j][0])
          worksheet.write(row, 1, pos[j][1])
          worksheet.write(row, 2, pos[j][2])
          worksheet.write(row, 3, vdw_radii[at[j]])
          worksheet.write(row, 4, ch[j])
          if at[j] == 18:
            worksheet.write(row, 6, 41) #argon
          elif at[j] == 84:
            worksheet.write(row, 6, 212) #polonium
          else:
            worksheet.write(row, 6, round(mass[j]))
          row += 1
        
        worksheet.write(0, 5, 'TOTAL z',bold)
        worksheet.write(1, 5, clusters_df["log"]["charge"][ind])
        worksheet.write(2, 5, 'Totalmass',bold)
        worksheet.write(3, 5, sum(clusters_df["xyz"]["structure"][ind].get_masses()))
        worksheet.write(4, 5, 'Crossection',bold)
        
      workbook.close()
      continue
    #CHARGES
    if i == "-chargesESP":
      for ind in clusters_df.index:
        f = open(clusters_df["info"]["file_basename"][ind]+".charges","w")
        f.write("\n".join([str(tt) for tt in clusters_df["log"]["esp_charges"][ind]])+"\n")
        f.close()
      continue
    if i == "-charges":
      for ind in clusters_df.index:
        f = open(clusters_df["info"]["file_basename"][ind]+".charges","w")
        f.write("\n".join([str(tt) for tt in clusters_df["log"]["mulliken_charges"][ind]])+"\n")
        f.close()
      continue
    #MOVIE
    if i == "-movie":
      #TODO: this could be done faster just with ASE command 
      from ase.io import write 
      from os import remove
      f = open("movie.xyz","w")
      for ind in clusters_df.index:
        try:
          atoms = clusters_df.loc[ind,("xyz","structure")]
          atoms.wrap()
          write(".movie.xyz",atoms, format='xyz')
          with open(".movie.xyz", 'r') as f2:
            lines = f2.readlines()
            try:
              lines[1] = clusters_df.loc[ind,("info","file_basename")]+f" Total Energy: {clusters_df.loc[ind,('log','electronic_energy')]}\n"
            except:
              lines[1] = clusters_df.loc[ind,("info","file_basename")]+"\n"
          f.writelines(lines)
          f2.close()
        except:
          continue
      f.close()
      remove(".movie.xyz")
      continue
    if i == "-gif":
      import mogli
      import imageio
      #from PIL import Image
      #molecules = mogli.read('movie.xyz')
      #images = []
      #for i in range(len(molecules)):
      #  mogli.export(molecules[i], 'movie.png', width=1920, height=1080,
      #               bonds_param=1.4, camera=((20, 0, 0), (0, 0, 0),(0, 1, 0)))
      #  img = Image.open("movie.png").convert("RGBA")
      #  images.append(img)
      #imageio.mimsave("movie.gif", images, duration=100, loop=0)
      #continue
      from ase.visualize.plot import plot_atoms
      import matplotlib.pyplot as plt
      from numpy import array
      images = []
      i=0;
      for ind in clusters_df.index:
        i=i+1
        print(str(i)+"/"+str(len(clusters_df)))
        fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
        ax.set_facecolor("white")
        ax.set_axis_off()
        plot_atoms(clusters_df.loc[ind,("xyz","structure")], ax=ax, radii=0.5, show_unit_cell=False)
        fig.canvas.draw()
        image = array(fig.canvas.renderer.buffer_rgba())
        images.append(image)
        plt.close(fig)
        imageio.mimsave("movie.gif", images, duration=0.1) 
      continue
    #Atoms
    if i == "-atoms":
      atoms = []
      for ind in clusters_df.index:
        try:
          aseCL=clusters_df.loc[ind,("xyz","structure")]
          atoms.extend(aseCL.get_atomic_numbers())
        except:
          continue
      print(" ".join([str(i) for i in list(set(atoms))]))
      continue
    #bonded
    if i == "-bonded":
      from numpy import array,sum
      bonded = []
      Qbonded_index += 1
      for ind in clusters_df.index:
        try:
          aseCL = clusters_df.loc[ind,("xyz","structure")]
          positions = aseCL.positions
          symb = aseCL.get_chemical_symbols()
          dist = lambda p1, p2: sum((p1-p2)**2)**0.5
          symb_ind = array(aseCL.get_chemical_symbols())
          mask1 = symb_ind == Qbonded[Qbonded_index][1]
          mask2 = symb_ind == Qbonded[Qbonded_index][2]
          dm = [dist(p1, p2) for p2 in positions[mask1] for p1 in positions[mask2]]
          bonds = sum(test_i <= float(Qbonded[Qbonded_index][0]) for test_i in dm)
          if Qbonded[Qbonded_index][1] == Qbonded[Qbonded_index][2]:
            bonds = (bonds-sum(mask1))/2
          bonded.append(str(int(bonds)))
        except:
          bonded.append(missing)
      output.append(bonded)
      continue
    #Rg
    if i == "-rg":
      from numpy import tile,sum
      import warnings
      warnings.filterwarnings("ignore", category=RuntimeWarning)
      rg = []
      for ind in clusters_df.index:
        try:
          aseCL=clusters_df.loc[ind,("xyz","structure")]
          rg.append((sum(sum((aseCL.positions-tile(aseCL.get_center_of_mass().transpose(),(len(aseCL.positions),1)))**2,axis=-1)*aseCL.get_masses())/sum(aseCL.get_masses()))**0.5)
        except:
          rg.append(missing)
      output.append(rg)
      continue
    #mean force
    if i == "-meanforce":
      from numpy import mean,sum,array
      meanforce = []
      for ind in clusters_df.index:
        try:
          aseCL=clusters_df.loc[ind,("extra","forces")]
          #print(aseCL)
          force = lambda p: sum(array(p)**2)**0.5
          meanforce.append(mean([force(array(p)) for p in aseCL]))
        except:
          meanforce.append(missing)
      output.append(meanforce)
      continue
    #max dist
    if i == "-maxdist":
      from numpy import max,sum
      maxdist = []
      for ind in clusters_df.index:
        try:
          aseCL=clusters_df.loc[ind,("xyz","structure")]
          dist = lambda p1, p2: sum((p1-p2)**2)**0.5
          maxdist.append(max([dist(p1, p2) for p1 in aseCL.get_positions() for p2 in aseCL.get_positions()]))
        except:
          maxdist.append(missing)
      output.append(maxdist)
      continue
    #ERRPA
    if i == "-errpa":
      err = []
      for ind in clusters_df.index:
        try:
          aseCL = clusters_df.loc[ind,("xyz","structure")]
          symb = aseCL.get_chemical_symbols()
          er = clusters_df.loc[ind,("extra","error")]
          atoms = len(symb)
          err.append(float(er/atoms))
        except:
          err.append(float("nan"))
      output.append(err)
      continue
    #Radius
    if i == "-radius" or i == "-radius0.5":
      from numpy import sqrt,linalg,asarray,dot,prod,max
      radius = []
      for aseCL in clusters_df["xyz"]["structure"]:
        try:
          dist = lambda p1, p2: sqrt(((p1-p2)**2).sum())
          centered = aseCL.positions-aseCL.positions.mean(axis = 0)
          ratios = sqrt(linalg.eigvalsh(dot(centered.transpose(),centered)/len(centered)))
          maxdist = max(asarray([[dist(p1, p2) for p2 in aseCL.positions] for p1 in aseCL.positions]))
          if i == "-radius0.5":
            if ratios[0] < 0.5:
              ratios[0] = 0.5
            if ratios[1] < 0.5:
              ratios[1] = 0.5
            if ratios[2] < 0.5:
              ratios[2] = 0.5
            if maxdist < 0.5:
              maxdist = 0.5
          if max(ratios) > 0:
            ratios = ratios / max(ratios)
          else:
            ratios = missing
          radius.append((prod(maxdist*ratios))**(1/3)/2)
        except:
          radius.append(missing)
      output.append(radius)
      continue
    ### INFO
    # cluster type
    if i == "-ct":
      try: 
        output.append(clusters_df.loc[:,("info","cluster_type")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # basename
    if i == "-b": 
      try:
        output.append(clusters_df.loc[:,("info","file_basename")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # name.out
    if i == "-nOUT":
      try:
        output.append(clusters_df.loc[:,("info","file_basename")]+".out")
      except:
        output.append([missing]*len(clusters_df))
      continue
    # name.log
    if i == "-nLOG":
      try:
        output.append(clusters_df.loc[:,("info","file_basename")]+".log")
      except:
        output.append([missing]*len(clusters_df))
      continue
    # name.xyz
    if i == "-nXYZ":
      try:
        output.append(clusters_df.loc[:,("info","file_basename")]+".xyz")
      except:
        output.append([missing]*len(clusters_df))
      continue
    # path/name.out
    if i == "-pOUT":
      try:
        output.append(clusters_df.loc[:,("info","folder_path")]+clusters_df.loc[:,("info","file_basename")]+".out")
      except:
        output.append([missing]*len(clusters_df))
      continue
    # path/name.log
    if i == "-pLOG":
      try:
        output.append(clusters_df.loc[:,("info","folder_path")]+clusters_df.loc[:,("info","file_basename")]+".log")
      except:
        output.append([missing]*len(clusters_df))
      continue
    # path/name.xyz
    if i == "-pXYZ":
      try:
        output.append(clusters_df.loc[:,("info","folder_path")]+clusters_df.loc[:,("info","file_basename")]+".xyz")
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ePKL
    if i == "-ePKL":
      from os import path
      filebasenames = clusters_df.loc[:,("info","file_basename")]
      if len(filebasenames) != len(list(set(filebasenames))):
        print("Sorry, there are non-unique clusters. I do not allow you to live.[HINT: -rebasename]")
        exit()
      if Qoutpkl == 1 and Qsplit == 1:
        try:
          output.append(path.abspath(output_pkl)+"/:EXTRACT:/"+filebasenames)
        except:
          output.append([missing]*len(clusters_df))
      elif Qoutpkl == 0 and len(input_pkl) == 1:
        try:
          output.append(path.abspath(input_pkl[0])+"/:EXTRACT:/"+filebasenames)
        except:
          output.append([missing]*len(clusters_df))
      else:
        print("Sorry but it seems you are taking this from more pickles or no file is formed")
        exit()
      continue
    #MASS
    if i == "-mass":
      masses = []  
      for ind in clusters_df.index:
        try:
          masses.append(str((clusters_df.loc[ind,("xyz","structure")].get_masses()).sum()))
        except:
          masses.append(missing)
      output.append(masses)
      continue
    #Natoms
    if i == "-natoms":
      natoms = []
      for ind in clusters_df.index:
        try:
          natoms.append(str(len(clusters_df.loc[ind,("xyz","structure")].get_atomic_numbers())))
        except:
          natoms.append(missing)
      output.append(natoms)
      continue
    if i == "-nel":
      nels = []
      for ind in clusters_df.index:
        try:
          from numpy import sum
          nels.append(str(sum(clusters_df.loc[ind,("xyz","structure")].get_atomic_numbers())))
        except:
          nels.append(missing)
      output.append(nels)
      continue
    # FORCES
    if i == "-maxf":
      fs = []
      for ind in clusters_df.index:
        try:
          from numpy import max,abs
          fs.append(str(QUenergy*max(abs(clusters_df.loc[ind,("extra","forces")]))))
        except:
          fs.append(missing)
      output.append(fs)
      continue
    # First = SP el. energy
    if i == "-elsp":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","sp_electronic_energy")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # el
    if i == "-el":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","electronic_energy")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # elout
    if i == "-elout":
      try:
        output.append(QUenergy*clusters_df.loc[:,("out","electronic_energy")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # el corr
    if i == "-elc":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("out","electronic_energy")]-clusters_df.loc[:,("log","electronic_energy")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # el SCF 
    if i == "-elscf":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("log","scf_energy")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # actual el corr??
    if i == "-elcorr":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("log","correlation_energy")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # energy th corr.
    if i == "-uc":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","energy_thermal_correction")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # en
    if i == "-u":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("log","electronic_energy")]+clusters_df.loc[:,("log","energy_thermal_correction")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # en out
    if i == "-uout":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("out","electronic_energy")]+clusters_df.loc[:,("log","energy_thermal_correction")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ZPE corr.
    if i == "-zpec": 
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","zero_point_correction")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ZPE
    if i == "-zpe":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("log","electronic_energy")]+clusters_df.loc[:,("log","zero_point_correction")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ZPE out
    if i == "-zpeout":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("out","electronic_energy")]+clusters_df.loc[:,("log","zero_point_correction")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # G
    if i == "-g":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","gibbs_free_energy")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # pop
    if i == "-pop":
      from math import isnan, exp
      from numpy import sum
      population = []
      R = 1.987 #cal/mol/K #=8.31441
      for pop in range(len(clusters_df)):
        try:
          if Qclustername != 0:
            allsameclusters = clusters_df.loc[clusters_df.loc[:,("info","cluster_type")]==clusters_df.iloc[pop].loc[("info","cluster_type")]]
          else: 
            allsameclusters = clusters_df
          if isnan(Qt):
            try:
              Qt = float(clusters_df.iloc[pop].loc[("log","temperature")])
            except:
              Qt = 298.15
          iii = clusters_df.iloc[pop]
          me=clusters_df.iloc[pop].loc[("log","gibbs_free_energy")] #*627.503*1000/Qt/R
          partition_function = sum([exp(iii) for iii in -(allsameclusters.loc[:,("log","gibbs_free_energy")]-me)*627.503*1000/Qt/R])
          ratio=exp(-0*627.503*1000/Qt/R)/partition_function
          population.append(f"{ratio:8f}")
        except:
          population.append(missing)
      output.append(population)
      continue
    # popEL
    if i == "-popEL":
      from math import isnan, exp
      from numpy import sum
      population = []
      R = 1.987 #cal/mol/K #=8.31441
      for pop in range(len(clusters_df)):
        try:
          if Qclustername != 0:
            allsameclusters = clusters_df.loc[clusters_df.loc[:,("info","cluster_type")]==clusters_df.loc[pop,("info","cluster_type")]]
          else:
            allsameclusters = clusters_df
          if isnan(Qt):
            try:
              Qt = float(clusters_df.loc[pop,("log","temperature")])
            except:
              Qt = 298.15
          me=clusters_df.iloc[pop].loc[("log","zero_point_energy")] #*627.503*1000/Qt/R
          everything=clusters_df.loc[:,("log","zero_point_energy")].values #*627.503*1000/Qt/R
          minimum=min(everything)
          partition_function = sum([exp(iii) for iii in -(everything-minimum)*627.503*1000/Qt/R])
          ratio=exp(-(me-minimum)*627.503*1000/Qt/R)
          population.append(f"{ratio:8f}")
        except:
          population.append(missing)
      output.append(population)
      continue
    # G corr
    if i == "-gc":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","gibbs_free_energy_thermal_correction")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # G out
    if i == "-gout":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("log","gibbs_free_energy")]+clusters_df.loc[:,("out","electronic_energy")]-clusters_df.loc[:,("log","electronic_energy")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # H
    if i == "-h":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","enthalpy_energy")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # H corr
    if i == "-hc":
      try:
        output.append(QUenergy*clusters_df.loc[:,("log","enthalpy_thermal_correction")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # H out
    if i == "-hout":
      try:
        output.append(QUenergy*(clusters_df.loc[:,("log","enthalpy_energy")]+clusters_df.loc[:,("out","electronic_energy")]-clusters_df.loc[:,("log","electronic_energy")]))
      except:
        output.append([missing]*len(clusters_df))
      continue
    # Entropy
    if i == "-s":
      try:
        output.append(QUentropy*clusters_df.loc[:,("log","entropy")]/1000/627.503)
      except:
        output.append([missing]*len(clusters_df))
      continue
    # Lowest freq.
    if i == "-lf": 
      lowestfreq = []
      for aseCL in clusters_df.loc[:,("log","vibrational_frequencies")]:
        try:
          lowestfreq.append(aseCL[0])
        except:
          lowestfreq.append(missing)
      output.append(lowestfreq)
      continue
    #LEVEL
    if i == "-level":
      from pandas import isna
      levels = []
      for ind in clusters_df.index:
        try:
          level = ""
          if ("log","method") in clusters_df and ("log","program") in clusters_df:
            if not isna(clusters_df.loc[ind,("log","program")]):
              level += clusters_df.at[ind,("log","program")]+"_"+clusters_df.at[ind,("log","method")]
            else:
              level += "nan"
            if ("out","method") in clusters_df and ("out","program") in clusters_df:
                level += "__"
          if ("out","method") in clusters_df and ("out","program") in clusters_df:
            if not isna(clusters_df.loc[ind,("out","program")]):
              level += clusters_df.at[ind,("out","program")]+"_"+clusters_df.at[ind,("out","method")]
            else:
              level += "nan"
          if level == "":
            level = "nan"
        except: 
          level = "nan"
        levels.append(level.replace(" ","_"))
      output.append(levels)
      continue
    # FREQs
    if i == "-f":
      try:
        output.append(clusters_df.loc[:,("log","vibrational_frequencies")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ROT CONST.
    if i == "-rot":
      try:
        output.append(clusters_df.loc[:,("log","rotational_constant")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ROT DIAG
    if i == "-rots":
      try:
        output.append(clusters_df.loc[:,("log","rotational_constants")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # Multiplicity
    if i == "-mult":
      try:
        output.append(clusters_df.loc[:,("log","multiplicity")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    #CHARGE
    if i == "-char":
      try:
        output.append(clusters_df.loc[:,("log","charge")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ESP CHARGE
    if i == "-esp":
      try:
        output.append(clusters_df.loc[:,("log","esp_charges")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # MILLIKEN CHARGE
    if i == "-mull":
      try:
        output.append(clusters_df.loc[:,("log","mulliken_charges")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # DIPOLE MOMENT
    if i == "-dip":
      try:
        output.append(clusters_df.loc[:,("log","dipole_moment")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # DIPS
    if i == "-dips":
      try:
        output.append(clusters_df.loc[:,("log","dipole_moments")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # POLARIZABILITY
    if i == "-pol":
      try:
        output.append(clusters_df.loc[:,("log","polarizability")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # TEMP in LOG
    if i == "-templog":
      try:
        output.append(clusters_df.loc[:,("log","temperature")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # PRESSURE in LOG
    if i == "-preslog":
      try:
        output.append(clusters_df.loc[:,("log","pressure")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # MOMENT OF INERTIA
    if i == "-mi":
      try:
        output.append([structure.get_moments_of_inertia() for structure in clusters_df.loc[:,("xyz","structure")]])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ????
    if i == "-ami":
      from numpy import mean
      try:
        output.append([mean(structure.get_moments_of_inertia()) for structure in clusters_df.loc[:,("xyz","structure")]])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # ROT SYMM NUMBER
    if i == "-rsn":
      try:
        output.append(clusters_df.loc[:,("log","rotational_symmetry_number")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # TIME
    if i == "-t":
      try:
        output.append(clusters_df.loc[:,("log","time")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # TERMINATION
    if i == "-termination":
      try: 
        output.append(clusters_df.loc[:,("log","termination")])
      except:
        output.append([missing]*len(clusters_df))
      continue
    # SPECIFIC COLUMN
    if i == "-column":
      try:
        output.append(clusters_df.loc[:,(Qcolumn[Qcolumn_i][0],Qcolumn[Qcolumn_i][1])])
        Qcolumn_i = Qcolumn_i + 1
      except:
        output.append([missing]*len(clusters_df))
      continue
    output.append(["UNKNOWN_ARGUMENT"]*len(clusters_df))
  return output
  
