#####################################################################################################
#####################################################################################################

def training(Qforces,Y_train,F_train,Qenergytradoff,strs,nn_tvv,nn_cutoff,nw,nn_rbf,Qrepresentation,nn_atom_basis,nn_interactions,Qbatch_size,Qlearningrate,parentdir,seed,varsoutfile,Qearlystop,nn_epochs,Qcheckpoint,Qtime,Qifcharges,Q_charges,Qifdipole,D_dipole):

  print("JKML(SchNetPack): Load quite heavy SchNetPack libraries.")
  from schnetpack.data import ASEAtomsData
  from schnetpack.data import AtomsDataModule
  import logging
  import schnetpack.transform as trn
  import torchmetrics
  import pytorch_lightning as pl
  from torch import cuda, nn, optim, manual_seed, backends
  import schnetpack as spk
  import os
  import numpy as np
  import time

  #os.environ['PYTHONHASHSEED'] = str(seed)
  #np.random.seed(seed)
  ## TORCH
  manual_seed(seed)
  cuda.manual_seed(seed)
  cuda.manual_seed_all(seed)
  backends.cudnn.deterministic = True
  backends.cudnn.benchmark = False
  pl.seed_everything(seed)


  print("JKML(SchNetPack): Setting up SchNetPack.", flush=True)

  if 1==1:
    import warnings
    warnings.filterwarnings(
        "ignore", ".*Trying to infer the `batch_size` from an ambiguous collection.*"
    )
    warnings.filterwarnings(
        "ignore",
        ".*Attribute 'model' is an instance of `nn.Module` and is already saved during checkpointing. It is recommended to ignore them using.*"
    )
    warnings.filterwarnings(
        "ignore",
        ".*The dataloader, .*_dataloader, does not have many workers which may be a bottleneck. Consider increasing the.*"
    )
    warnings.filterwarnings(
        "ignore",
        ".*which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling.*"
    )
    warnings.filterwarnings(
        "ignore",
        ".*The verbose parameter is deprecated. Please use get.*"
    )

  # PREPARING TRAINING DATABASE
  print("JKML(SchNetPack): Adjusting training database for SchNetPack.", flush=True)
  
  if Qifdipole == 1:
    from src.JKelectrostatics import compute_energies_forces
    from src.JKdispersions import compute_d4_energy_forces as compute_dispersions
    #from src.JKdispersions import compute_d3bj_energy_forces as compute_dispersions
    if 1==0:
      #D_dipole = np.array([np.array(i)[0] for i in D_dipole])
      Qifdipole = 2
    else:
      #print(D_dipole)
      D_dipole = []
      for i in range(0,len(Q_charges)):
        Dx = 0
        Dy = 0
        Dz = 0
        #dipole moment
        tmpatoms = strs.values[i].copy()
        if 1==0:
          from tblite.ase import TBLite
          tmpatoms.calc = TBLite(method="GFN1-xTB", cache_api=True, charge=float(0), verbosity = 0, max_iterations = 30000, accuracy = 1.0)
          Q_charges_tmp = tmpatoms.get_charges()
        elif 1==1:
          from schnetpack.interfaces import SpkCalculator
          spk_calc_charges = SpkCalculator(
            model_file="/home/kubeckaj/TEST/Ivo/PN/S_100/charges.pkl",
            device="cpu",
            neighbor_list=trn.ASENeighborList(cutoff=10.0),
            position_unit="Ang",
          )
          tmpatoms.calc = spk_calc_charges
          print(tmpatoms.get_potential_energy())
          Q_charges_tmp = spk_calc_charges.model_results['partial_charges'].detach().numpy()
        else:
          Q_charges_tmp = Q_charges[i] 
        print(Q_charges_tmp)
        print(Q_charges[i])
        R_tmp = strs.values[i].get_positions()
        Ntmp = len(Q_charges_tmp)
        for k in range(0,Ntmp):
            Dx += Q_charges_tmp[k]*(R_tmp[k][0])
            Dy += Q_charges_tmp[k]*(R_tmp[k][1])
            Dz += Q_charges_tmp[k]*(R_tmp[k][2])
        D_dipole.append(list([Dx, Dy, Dz]))
        
        electrostatics_E, electrostatics_F = compute_energies_forces(strs.values[i].get_positions(), Q_charges_tmp)
        dispersions_E, dispersions_F = compute_dispersions(strs.values[i].get_positions(), symbols = np.array(strs.values[i].get_chemical_symbols()), totalcharge = 0)
        Y_train[i] -= electrostatics_E + dispersions_E
        print(f"EE:JKML(SchNetPack): {Y_train[i]+electrostatics_E+dispersions_E} {Y_train[i]} {electrostatics_E} {dispersions_E}")
        F_train[i] -= electrostatics_F + dispersions_F
        import numpy as np
        #print(F_train[i])
        #print(electrostatics_F)
        #print(np.sum((F_train[i])**2,axis=1))
        #print(np.sum((electrostatics_F)**2,axis=1))
        #print(np.sum((F_train[i]+electrostatics_F)**2,axis=1))
        #print(f"FF:(SchNetPack): {np.mean(np.sum((F_train[i]+electrostatics_F+dispersions_F)**2,axis=1)**0.5)} {np.mean(np.sum((F_train[i])**2,axis=1)**0.5)} {np.mean(np.sum((electrostatics_F)**2,axis=1)**0.5)} {np.mean(np.sum((dispersions_F)**2,axis=1)**0.5)}")
        print(f"FF:(SchNetPack): {np.sum(np.sum((F_train[i]+electrostatics_F+dispersions_F)**2,axis=1)**0.5)} {np.sum(np.sum((F_train[i])**2,axis=1)**0.5)} {np.sum(np.sum((electrostatics_F)**2,axis=1)**0.5)} {np.sum(np.sum((dispersions_F)**2,axis=1)**0.5)}")
        #print(f"FF:(SchNetPack): {np.sum(np.sum((F_train[i]+electrostatics_F+dispersions_F),axis=0)**2)**0.5} {np.sum(np.sum((F_train[i]),axis=0)**2)**0.5} {np.sum(np.sum((electrostatics_F),axis=0)**2)**0.5} {np.sum(np.sum((dispersions_F),axis=0)**2)**0.5}")
      Qifdipole = 0
      Qifcharges = 0  

  temperary_file_name = "training.db"
  if os.path.exists(temperary_file_name):
      os.remove(temperary_file_name)
  if Qforces == 0:
      new_dataset = ASEAtomsData.create(temperary_file_name,
                                        distance_unit='Ang',
                                        #NOTE that units of energy are not important in this part
                                        property_unit_dict={'energy': 'eV', 'total_charge': 'e'},
                                        atomrefs={'energy': [0] * 100}
                                        )
      properties = [{'energy': np.array([i]), 'total_charge': np.array([0], dtype=np.float32)} for i in
                    Y_train]
      target_properties = [spk.properties.energy]
      tradoffs = [1]
  elif Qifdipole == 1:
      new_dataset = ASEAtomsData.create(temperary_file_name,
                                        distance_unit='Ang',
                                        property_unit_dict={'energy': 'eV', 'forces': 'eV/Ang',
                                                            'total_charge': 'e', 'dipole_moment': 'e*Ang'},
                                        atomrefs={'energy': [0] * 100}
                                        )
      properties = [{'energy': 27.2107 * np.array([Y_train[i]]), 'forces': 27.2114 * np.array(F_train[i]),
                     'total_charge': np.array([0], dtype=np.float32), 'dipole_moment': np.array([D_dipole[i]], dtype=np.float32)} for i in range(len(Y_train))]
      target_properties = [spk.properties.energy, spk.properties.forces, spk.properties.dipole_moment]
      tradoffs = [Qenergytradoff, 1, Qenergytradoff]
  elif Qifdipole == 2: #DOES NOT WORK
      new_dataset = ASEAtomsData.create(temperary_file_name,
                                        distance_unit='Ang',
                                        property_unit_dict={'energy': 'eV', 'total_charge': 'e', 'dipole_moment': 'e*Ang'},
                                        atomrefs={'energy': [0] * 100}
                                        )
      properties = [{'energy': 27.2107 * np.array([0], dtype=np.float32), 'total_charge': np.array([0], dtype=np.float32), 'dipole_moment': np.array([D_dipole[i]], dtype=np.float32)} for i in range(len(Y_train))]
      target_properties = [spk.properties.energy, spk.properties.dipole_moment]
      tradoffs = [0, 1]
  elif Qifcharges == 1: #DOES NOT WORK
      new_dataset = ASEAtomsData.create(temperary_file_name,
                                        distance_unit='Ang',
                                        property_unit_dict={'energy': 'eV', 'forces': 'eV/Ang',
                                                            'total_charge': 'e', 'partial_charges': 'e'},
                                        atomrefs={'energy': [0] * 100}
                                        )
      properties = [{'energy': 27.2107 * np.array([Y_train[i]]), 'forces': 27.2114 * np.array(F_train[i]),
                     'total_charge': np.array([0], dtype=np.float32), 'partial_charges': np.array(Q_charges[i], dtype=np.float32)} for i in range(len(Y_train))]
      target_properties = [spk.properties.energy, spk.properties.forces, spk.properties.partial_charges]
      tradoffs = [Qenergytradoff, 1, 1]
  else:
      new_dataset = ASEAtomsData.create(temperary_file_name,
                                        distance_unit='Ang',
                                        property_unit_dict={'energy': 'eV', 'forces': 'eV/Ang',
                                                            'total_charge': 'e'},
                                        atomrefs={'energy': [0] * 100}
                                        )
      properties = [{'energy': 27.2107 * np.array([Y_train[i]]), 'forces': 27.2114 * np.array(F_train[i]),
      #properties = [{'energy': np.round(27.2107 * np.array([Y_train[i]]), 4), 'forces': np.round(27.2114 * np.array(F_train[i]), 4),
                     'total_charge': np.array([0], dtype=np.float32)} for i in range(len(Y_train))]
      target_properties = [spk.properties.energy, spk.properties.forces]
      tradoffs = [Qenergytradoff, 1]

  new_dataset.add_systems(properties, strs)

  n_train = int(np.round(nn_tvv * len(strs)))
  n_val = len(strs) - n_train

  if cuda.is_available():
      pin_memory = True
      device = 'gpu'
      logging.info(cuda.get_device_name(0))
  else:
      pin_memory = False
      device = 'cpu'

  dataset = AtomsDataModule(
      temperary_file_name,
      batch_size=Qbatch_size,
      num_train=n_train,
      num_val=n_val,
      transforms=[
          trn.ASENeighborList(cutoff=nn_cutoff),
          trn.RemoveOffsets(spk.properties.energy, remove_mean=True,
                            remove_atomrefs=False),
          trn.CastTo32()
      ],
      num_workers=nw,  # use 2-4
      pin_memory=pin_memory,
      data_workdir="./"
  )

  print("JKML(SchNetPack): Setting up training database for NN.", flush=True)
  dataset.prepare_data()
  dataset.setup()
  print("JKML(SchNetPack): Setting up NN.", flush=True)

  # PainNN representation
  pairwise_distance = spk.atomistic.PairwiseDistances()
  radial_basis = spk.nn.GaussianRBF(n_rbf=nn_rbf, cutoff=nn_cutoff)

  if Qrepresentation == "painn":
      X_train = spk.representation.PaiNN(
          n_atom_basis=nn_atom_basis,
          n_interactions=nn_interactions,
          radial_basis=radial_basis,
          cutoff_fn=spk.nn.CosineCutoff(nn_cutoff)
      )
  elif Qrepresentation == "schnet":
      X_train = spk.representation.SchNet(
          n_atom_basis=nn_atom_basis,
          n_interactions=nn_interactions,
          radial_basis=radial_basis,
          cutoff_fn=spk.nn.CosineCutoff(nn_cutoff)
      )
  elif Qrepresentation == "so3net":
      X_train = spk.representation.SO3net(
          n_atom_basis=nn_atom_basis,
          n_interactions=nn_interactions,
          radial_basis=radial_basis,
          cutoff_fn=spk.nn.CosineCutoff(nn_cutoff),
      )
  else:
      raise Exception("JKML(SchNetPack): You have probably enetered some weird molecular representation. [EXIT]")

  output_modules = []
  output_losses = []
  for p, w in zip(target_properties, tradoffs):
      if p == spk.properties.energy:
          pred = spk.atomistic.Atomwise(n_in=nn_atom_basis, output_key=p)
      elif p == spk.properties.forces:
          pred = spk.atomistic.Forces(energy_key=spk.properties.energy, force_key=spk.properties.forces)
      elif p == spk.properties.dipole_moment:
          #pred = spk.atomistic.DipoleMoment(n_in=nn_atom_basis, return_charges=True) #, use_vector_representation = True)
          pred = spk.atomistic.DipoleMoment(n_in=nn_atom_basis, return_charges=True, use_vector_representation = True)
      elif p == spk.properties.partial_charges:
          pred = spk.atomistic.Atomwise(n_in=nn_atom_basis, output_key=p)
      else:
          raise NotImplementedError(f'{p} property does not exist')

      loss = spk.task.ModelOutput(
          name=p,
          loss_fn=nn.MSELoss(),
          loss_weight=w,
          metrics={
              "MAE": torchmetrics.MeanAbsoluteError(),
              "RMSE": torchmetrics.MeanSquaredError(squared=False)
          }
      )
      output_modules.append(pred)
      output_losses.append(loss)

      # MODEL (this could be for instance also Atomistic Model)
  #print(output_modules)
  #print(output_losses)
  nnpot = spk.model.NeuralNetworkPotential(
      representation=X_train,
      input_modules=[pairwise_distance],
      output_modules=output_modules,
      postprocessors=[
          trn.CastTo64(),
          trn.AddOffsets(spk.properties.energy, add_mean=True, add_atomrefs=False)
      ]
  )

  task = spk.task.AtomisticTask(
      model=nnpot,
      outputs=output_losses,
      optimizer_cls=optim.AdamW,
      optimizer_args={"lr": Qlearningrate},
      scheduler_cls=spk.train.ReduceLROnPlateau,
      scheduler_args={'factor': 0.5, 'patience': 20, 'min_lr': 1e-7},
      scheduler_monitor='val_loss'
  )

  os.makedirs("lightning_logs")
  logger = pl.loggers.TensorBoardLogger(save_dir="./")

  class MyPrintingCallback(pl.callbacks.LearningRateMonitor):
      def __init__(self, *args, **kwargs):
          super().__init__()
          self.state_train = []
          self.state_val = []
          pl.seed_everything(seed)
          global start_time
          start_time = time.time()
          print("JKML(SchNetPack): NN training starting", flush=True)

      def on_fit_start(self, trainer, pl_module):
          super().on_fit_start(trainer, pl_module)
          # print(trainer.callbacks[3])
          if trainer.callbacks[3].best_model_score is not None:
              # trainer.callbacks[3].best_model_score  = None #torch.tensor(float("inf"))
              trainer.callbacks[3].best_k_models = {}
              trainer.callbacks[3].best_model_score = None

      def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx, unused=0):
          super().on_train_batch_end(trainer, pl_module, outputs, batch, batch_idx)
          #print(outputs)
          self.state_train.append(23.060541945329334 * outputs["loss"].item())
          print(".", end="", flush=True)

      def on_train_end(self, trainer, pl_module):
          super().on_train_end(trainer, pl_module)
          print("JKML(SchNetPack): Training is ending", flush=True)

      def on_validation_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
          super().on_validation_batch_end(trainer, pl_module, outputs, batch, batch_idx)
          self.state_val.append(23.060541945329334 * outputs["val_loss"].item())
          print("_", end="", flush=True)

      def on_validation_epoch_end(self, trainer, pl_module):
          super().on_validation_epoch_end(trainer, pl_module)
          # access output using state
          all_outputs_train = self.state_train
          all_outputs_val = self.state_val
          self.state_train = []
          self.state_val = []

          def mean_std(arr):
              if len(arr) == 0:
                  return np.nan, np.nan
              elif len(arr) == 1:
                  return arr[0], np.nan
              else:
                  return np.mean(arr), np.std(arr)

          tr_m, tr_s = mean_std(all_outputs_train)
          val_m, val_s = mean_std(all_outputs_val)
          ep = trainer.current_epoch
          run_time = time.time() - start_time
          print(f"\nEPOCH: %i TRAIN: %.5f +- %.5f [??] VAL: %.5f +- %.5f [??] LR: %.2e T: %.1f s " % (
              ep, tr_m, tr_s, val_m, val_s, task.lr, run_time))
          with open(parentdir + "/epoch_trE_trS_valE_valS_lr_t.txt", "a") as f:
              print(f"%i %.6f %.6f %.6f %.6f %e %.1f" % (ep, tr_m, tr_s, val_m, val_s, task.lr, run_time),
                    file=f)

  callbacks = [
      spk.train.ModelCheckpoint(
          model_path=os.path.join("./", varsoutfile),
          save_top_k=1,
          monitor="val_loss",
          save_last=True,
      ),
      MyPrintingCallback(
          monitor="val_loss",
          patience=Qearlystop,
      ),
  ]

  trainer = pl.Trainer(
      accelerator=device,
      # devices=1,
      enable_progress_bar=False,  # TODO testing
      log_every_n_steps=1,  # np.round(n_train/Qbatch_size/10),
      callbacks=callbacks,  # checks for early stopping and checkpoints
      logger=logger,  # This is the logger for tensorboard
      default_root_dir="./lightning_logs",
      max_epochs=nn_epochs,
      max_time=Qtime  # 90% e.g."00:12:00:00"
      # resume_from_checkpoint=model_checkpoint,
      # log_every_n_steps=1,
      # accumulate_grad_batches = 10,
      # TODO: gpus=2
      # TODO: precision=16
  )

  print("JKML(SchNetPack): Initiating training.", flush=True)

  if Qcheckpoint is None:
      trainer.fit(task, datamodule=dataset)
  else:
      trainer.fit(task, datamodule=dataset, ckpt_path=Qcheckpoint)
  print("JKML(SchNetPack): Training done", flush=True)

  #globals().update(locals())
  return 

#####################################################################################################
#####################################################################################################

def evaluate(Qforces,varsoutfile,nn_cutoff,clusters_df,method,Qmin,Qifcharges):

  from schnetpack.interfaces import SpkCalculator
  from torch import cuda
  import schnetpack.transform as trn
  from numpy import array

  if 1==1:
    import warnings
    warnings.filterwarnings(
        "ignore",
        ".*which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling.*"
    )

  if cuda.is_available():
      device = 'cuda'
  else:
      device = 'cpu'

  print("JKML(SchNetPack): Calculator loading", flush=True)
  if Qforces == 0:
      spk_calc = SpkCalculator(
          model_file=varsoutfile,
          device=device,
          neighbor_list=trn.ASENeighborList(cutoff=nn_cutoff),
          # neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
          # transforms=spk.transform.atomistic.SubtractCenterOfMass(),
          energy_key='energy',
          energy_unit="eV",  # YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
          position_unit="Ang",
      )
  elif Qifcharges == 1:
      spk_calc = SpkCalculator(
          model_file=varsoutfile,
          device=device,
          neighbor_list=trn.ASENeighborList(cutoff=nn_cutoff),
          # neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
          # transforms=spk.transform.atomistic.SubtractCenterOfMass(),
          energy_key='energy',
          energy_unit="eV",  # YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
          position_unit="Ang",
          #dipole_key="dipole_moment",
          #dipole_unit="e*Ang",
      )
  else:
      spk_calc = SpkCalculator(
          model_file=varsoutfile,
          device=device,
          neighbor_list=trn.ASENeighborList(cutoff=nn_cutoff),
          # neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
          # transforms=spk.transform.atomistic.SubtractCenterOfMass(),
          energy_key='energy',
          force_key='forces',
          energy_unit="eV",  # YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
          force_unit="eV/Ang",  # YEAH I have no idea what the output is :-D
          position_unit="Ang",
      )

  Y_predicted = []
  F_predicted = []
  Q_predicted = []
  for i in range(len(clusters_df)):
      atoms = clusters_df["xyz"]["structure"].values[i].copy()
      atoms.calc = spk_calc
      if Qifcharges == 1:
        from src.JKelectrostatics import compute_energies_forces
        from src.JKdispersions import compute_d4_energy_forces as compute_dispersions
        #from src.JKdispersions import compute_d3bj_energy_forces as compute_dispersions

        internal_E = atoms.get_potential_energy()
        internal_F = atoms.get_forces()
        #TODO
        if 1==0:
          from tblite.ase import TBLite
          tmpatoms = atoms.copy()
          tmpatoms.calc = TBLite(method="GFN1-xTB", cache_api=True, charge=float(0), verbosity = 0, max_iterations = 300, accuracy = 1.0)
          Q_charges_i = array([tmpatoms.get_charges()]) #.transpose()
        elif 1==1:
          from schnetpack.interfaces import SpkCalculator
          tmpatoms = atoms.copy()
          spk_calc_charges = SpkCalculator(
            model_file="/home/kubeckaj/TEST/Ivo/PN/S_100/charges.pkl",
            device="cpu",
            neighbor_list=trn.ASENeighborList(cutoff=10.0),
            position_unit="Ang",
          )
          tmpatoms.calc = spk_calc_charges
          print(tmpatoms.get_potential_energy())
          Q_charges_i = spk_calc_charges.model_results['partial_charges'].detach().numpy()
        else:
          Q_charges_i = spk_calc.model_results['partial_charges'].detach().numpy()
        #print(f"JKML(SchNetPack): {Q_charges_i}")
        #TODO
	#print(f"JKML(SchNetPack): {Q_charges_i}")
        #TODO add charge PERHAPS NOT NECESSARY STEP
        #Q_charges_i = Q_charges_i - (Q_charges_i - 0).mean()
        Q_predicted.append(Q_charges_i)
        electrostatics_E, electrostatics_F = compute_energies_forces(atoms.get_positions(), Q_charges_i)
        dispersions_E, dispersions_F = compute_dispersions(atoms.get_positions(), symbols = array(atoms.get_chemical_symbols()), totalcharge = 0)
        print(f"JKML(SchNetPack): {0.0367493 * internal_E} {electrostatics_E} {dispersions_E}")
        Y_predicted.append(0.0367493 * internal_E + electrostatics_E + dispersions_E)
        F_predicted.append(0.0367493 * internal_F + electrostatics_F + dispersions_F)
      else:
        Y_predicted.append(0.0367493 * atoms.get_potential_energy())
        if Qforces == 1:
          F_predicted.append(0.0367493 * atoms.get_forces())  # Hartree/Ang

  if Qforces == 1:
      Y_predicted = [array(Y_predicted)]  # Hartree
      F_predicted = [F_predicted]
  else:
      Y_predicted = [array(Y_predicted)]

  if method == "min":
      Y_predicted[0] += Qmin

  return Y_predicted, F_predicted, Q_predicted

#####################################################################################################
#####################################################################################################

def optimize(test_high_database,varsoutfile,nn_cutoff,Qopt,opt_maxstep,opt_dump,opt_steps,md_temperature,Qmd_timestep,md_thermostatfriction,md_dump,md_steps):

    from schnetpack.interfaces import SpkCalculator
    import torch
    import schnetpack as spk

    if 1==1:
      import warnings
      warnings.filterwarnings(
          "ignore",
          ".*which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling.*"
      )

    ### DATABASE LOADING ###
    ## The high level of theory
    clusters_df = test_high_database
    strs = clusters_df["xyz"]["structure"]

    if torch.cuda.is_available():
        device = 'cuda'
    else:
        device = 'cpu'

    spk_calc = SpkCalculator(
        model_file=varsoutfile,
        device=device,
        neighbor_list=spk.transform.ASENeighborList(cutoff=nn_cutoff),
        # neighbor_list=spk.transform.TorchNeighborList(cutoff=5.0),
        # transforms=spk.transform.atomistic.SubtractCenterOfMass(),
        energy_key='energy',
        force_key='forces',
        energy_unit="eV",  # YEAH BUT THE OUTPUT UNITS ARE ACTUALLY Hartree
        force_unit="eV/Ang",  # YEAH I have no idea what the output is :-D
        position_unit="Ang",
    )

    Y_predicted = []
    F_predicted = []
    from ase.optimize import BFGS
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.verlet import VelocityVerlet
    from ase.md.langevin import Langevin
    from ase import units
    from ase.io import read, write

    for i in range(len(clusters_df)):
        atoms = clusters_df["xyz"]["structure"].values[i].copy()
        print(atoms)
        atoms.calc = spk_calc
        if Qopt == 1:
            dyn = BFGS(atoms, maxstep=opt_maxstep)

            def printenergy(a=atoms):
                write("opt.xyz", a, append=True)

            dyn.attach(printenergy, interval=opt_dump)
            dyn.run(fmax=1e-6, steps=opt_steps)
        else:

            # Set the momenta corresponding to T
            MaxwellBoltzmannDistribution(atoms, temperature_K=md_temperature)
            # We want to run MD with constant energy using the VelocityVerlet algorithm.
            # dyn = VelocityVerlet(atoms, 1 * units.fs)  # 5 fs time step.
            dyn = Langevin(atoms, Qmd_timestep * units.fs, temperature_K=md_temperature,
                           friction=md_thermostatfriction)  # friction coeffitient 0.002

            def printenergy(a=atoms):  # store a reference to atoms in the definition.
                """Function to print the potential, kinetic and total energy."""
                epot = a.get_potential_energy() / len(a)
                ekin = a.get_kinetic_energy() / len(a)
                write("traj.xyz", a, append=True)
                print('Energy per atom: Epot = %.3f kcal/mol  Ekin = %.3f kcal/mol (T=%3.0f K)  '
                      'Etot = %.3f kcal/mol' % (
                          23.060541945329334 * epot, 23.060541945329334 * ekin, ekin / (1.5 * units.kB),
                          23.060541945329334 * (epot + ekin)))

            # Now run the dynamics
            dyn.attach(printenergy, interval=md_dump)
            printenergy()
            dyn.run(md_steps)

    globals().update(locals())
