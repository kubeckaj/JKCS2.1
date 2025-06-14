def calculator(Qcalculator, Qcalculator_input, Qcalculator_max_iterations = 300, Qcharge = 0, Qmultiplicity = 1, Qout = 1, 
               atoms=None, Qmixer_damping=0.4, Qcutoff=10.0):

# TODO
### 792 /home/kubeckaj/Applications/JKCS2.1/JKQC/JKCS/lib/python3.9/site-packages/ase/calculators/calculator.py
#        if atoms is not None:
#            CS = atoms.constraints
#            del atoms.constraints
#            self.atoms = atoms.copy()
#            atoms.set_constraint(CS)
  if Qout > 0:
    Qprint=Qout-1
  else:
    Qprint=Qout
   
  #### XTB-Lite ###
  if Qcalculator == "XTB1":
    #if Qcharge != 0:
    #  print("Oh sorry, the charge has not been implemented yet, ask Jakub.")
    #  exit()
    from tblite.ase import TBLite
    return TBLite(method="GFN1-xTB", cache_api=True, charge=float(Qcharge), multiplicity = Qmultiplicity, verbosity = Qprint, max_iterations = Qcalculator_max_iterations, accuracy = 1.0, mixer_damping = Qmixer_damping) #, initial_guess = "eeq")

  elif Qcalculator == "XTB2":
    #if Qcharge != 0:
    #  print("Oh sorry, the charge has not been implemented yet, ask Jakub.")
    #  exit()
    from tblite.ase import TBLite
    return TBLite(method="GFN2-xTB", cache_api=True, charge=float(Qcharge), multiplicity = Qmultiplicity, verbosity = Qprint, max_iterations = Qcalculator_max_iterations, accuracy = 1.0, mixer_damping = Qmixer_damping)

  ### XTB ###
  elif Qcalculator == "XTB":
    #if Qcharge != 0:
    #  print("Oh sorry, the charge has not been implemented yet, ask Jakub.")
    #  exit()
    if Qmultiplicity != 1:
      print("Oh sorry, the multiplicity has not been implemented yet, ask Jakub.")
      exit()
    from xtb.ase.calculator import XTB
    return XTB(method=Qcalculator_input, charge = float(Qcharge)) ;#+" --chrg "+str(Qcharge))#, charge=Qcharge)

  ### ORCA ###
  elif Qcalculator == "ORCA":
    from ase.calculators.orca import ORCA
    import psutil
    import os
    #os.system("if ! command -v module &> /dev/null; then source /com/bin/modules.sh; fi; module load intel; module load openmpi;")
    #os.system("if ! command -v module &> /dev/null; then source /com/bin/modules.sh; fi; module load gcc openmpi mkl")
    #os.environ['OMP_STACKSIZE'] = '4G'
    #os.environ['OMP_NUM_THREADS'] = f'{len(psutil.Process().cpu_affinity())}'
                               #f'{len(psutil.Process().cpu_affinity())},1'
    #os.environ['OMP_MAX_ACTIVE_LEVELS'] = '1'
    #https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html#calculator-configuration
    #os.environ['ASE_CONFIG_PATH'] =  "/home/kubeckaj/Applications/ORCA6.0/orca_6_0_1_linux_x86-64_shared_openmpi416/"
    #profile="/home/kubeckaj/Applications/ORCA6.0/orca_6_0_1_linux_x86-64_shared_openmpi416/orca",
    from ase.calculators.orca import OrcaProfile
    profile = OrcaProfile(command = "/home/kubeckaj/Applications/ORCA6.0/orca_6_0_1_linux_x86-64_shared_openmpi416/orca")
    cpus = f'{len(psutil.Process().cpu_affinity())}'
    return ORCA(profile = profile, 
                charge=Qcharge,
                mult=Qmultiplicity,
                orcasimpleinput=Qcalculator_input+" engrad",
                orcablocks='%pal nprocs '+cpus+' end')

  ### SchNetPack ###
  elif Qcalculator == "NN":
    from schnetpack.interfaces import SpkCalculator
    import schnetpack as spk
    if 1==1:
      import warnings
      warnings.filterwarnings(
          "ignore",
          ".*which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling.*"
      )
    return SpkCalculator(
        model_file = Qcalculator_input,
        device="cpu",
        neighbor_list=spk.transform.ASENeighborList(cutoff=Qcutoff),
        #neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0), #,cutoff_shell=2.0),
        #transforms=spk.transform.atomistic.SubtractCenterOfMass(),
        energy_key='energy',
        force_key='forces',
        energy_unit="eV",#not sure about units --- maybe Hartree
        force_unit="eV/Ang",#not sure about units
        position_unit="Ang",
        )

  ### PhysNet ###
  elif Qcalculator == "PhysNet":
    import os,sys
    pathname = os.path.dirname(sys.argv[0])
    sys.path.append(pathname + "/../JKML/src/PhysNet_DER/")
    from PNcalculator import PhysNetCalculator  
    if 1==1:
      import warnings
      warnings.filterwarnings(
          "ignore",
          ".*which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling.*"
      ) 
    from torch import cuda
    if cuda.is_available():
      os.system("sed -i 's/--device=cpu/--device=cuda/g' input.inp")
    else:
      os.system("sed -i 's/--device=cuda/--device=cpu/g' input.inp")
  
    return PhysNetCalculator(
      checkpoint=Qcalculator_input,
      atoms=atoms,
      charge=Qcharge,
      config='input.inp')

  ### AIMNET2 ###
  elif Qcalculator == "AIMNET2":
    #from aimnet2calc import AIMNet2ASE
    from aimnet.calculators.aimnet2ase import AIMNet2ASE
    return AIMNet2ASE(Qcalculator_input, charge=Qcharge, mult=Qmultiplicity)

  else:
    print("Unknown calculator.")
    exit()
