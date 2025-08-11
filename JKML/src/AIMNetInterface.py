####################################################################################################
####################################################################################################

def training(Y_train,F_train,Z_atoms,N_atoms,Qenergytradoff,strs,nn_tvv,Qbatch_size,Qlearningrate,parentdir,varsoutfile,Qearlystop,nn_epochs,Q_charge,Q_charges,file_basenames,seed,Qmonomers):

    import numpy as np
    import h5py
    from collections import defaultdict
    import os
    import sys
    import subprocess
    import yaml
    import glob
    import shutil

    thepath=glob.glob(os.path.dirname(os.path.abspath(__file__))+'/../../JKQC/JKCS/AIMNET/lib/py*/site-packages/')[0]
    env = os.environ.copy()
    env["PYTHONPATH"] = thepath

    grouped_data = defaultdict(lambda: defaultdict(list))
    str_dtype = h5py.string_dtype(encoding="utf-8")

    print("Creating HDF5 file...")

    for i in range(len(Y_train)):
        natoms = N_atoms[i]
        group = f"{natoms:03d}"
        
        hartree_to_ev = 27.211386245988
        grouped_data[group]["coord"].append(strs.values[i].get_positions().astype(np.float32)) # (N_atoms, 3)
        grouped_data[group]["numbers"].append(Z_atoms[i].astype(np.int32))                     # (N_atoms,)
        grouped_data[group]["forces"].append(F_train[i].astype(np.float32) * hartree_to_ev)    # (N_atoms, 3)
        grouped_data[group]["energy"].append(np.float32(Y_train[i] * hartree_to_ev))           # scalar
        #print(np.float32(Y_train[i] * hartree_to_ev))
        grouped_data[group]["charges"].append(np.asarray(Q_charges[i], dtype=np.float32))      # (N_atoms,)
        grouped_data[group]["charge"].append(np.float32(Q_charge[i]))                          # scalar
        grouped_data[group]["basename"].append(file_basenames[i])

    with h5py.File("dataset.h5", "w") as f:
        for group, data in grouped_data.items():
            g = f.create_group(group)
            for key in ["coord", "numbers", "forces", "energy", "charges", "charge","basename"]:
                arr = np.array(data[key])
                if key == "numbers":
                    g.create_dataset(key, data=arr.astype(np.int32))
                elif key == "basename":
                    g.create_dataset(key, data=arr.astype(object), dtype=str_dtype)
                else:
                    g.create_dataset(key, data=arr.astype(np.float32))
    
    subprocess.run(["aimnet", "calc_sae", "dataset.h5", "dataset_sae.yaml"], env=env)
    
    if Qmonomers != 2:
        with open("dataset_sae.yaml", "r") as f:
            lines = f.readlines()

        with open("dataset_sae.yaml", "w") as f:
            for line in lines:
                if ':' in line:
                    key = line.split(':')[0].strip()
                    f.write(f"{key}: 0.00000000000000000\n")
                else:
                    f.write(line)
        
    run_name = os.path.basename(parentdir)
    
    # Build config dict
    conf = {
        "run_name": run_name,
        "data": {
            "train": "dataset.h5",
            "val": None,
            "sae": {
                "energy": {
                    "file": "dataset_sae.yaml",
                    "mode": "linreg"
                }
            },
            "val_fraction": 1.0 - nn_tvv,
            "separate_val": True,
            "ddp_load_full_dataset": False,
            "x": ["coord", "numbers", "charge"],
            "y": ["energy", "forces", "charges"],
            "datasets": {
                "train": {"class": "aimnet.data.SizeGroupedDataset", "kwargs": {}},
                "val": {"class": "aimnet.data.SizeGroupedDataset", "kwargs": {}}
            },
            "samplers": {
                "train": {
                    "class": "aimnet.data.SizeGroupedSampler",
                    "kwargs": {
                        "batch_size": Qbatch_size,
                        "batch_mode": "molecules",
                        "shuffle": True,
                        "batches_per_epoch": 1
                    }
                },
                "val": {
                    "class": "aimnet.data.SizeGroupedSampler",
                    "kwargs": {
                        "batch_size": Qbatch_size,
                        "batch_mode": "molecules",
                        "shuffle": False,
                        "batches_per_epoch": -1
                    }
                }
            },
            "loaders": {
                "train": {"num_workers": 0, "pin_memory": True},
                "val": {"num_workers": 0, "pin_memory": True}
            }
        },
        "loss": {
            "class": "aimnet.train.loss.MTLoss",
            "kwargs": {
                "components": {
                    "energy": {
                        "fn": "aimnet.train.loss.energy_loss_fn",
                        "weight": 1.0
                    },
                    "forces": {
                        "fn": "aimnet.train.loss.peratom_loss_fn",
                        "weight": 0.2,
                        "kwargs": {"key_true": "forces", "key_pred": "forces"}
                    },
                    "charges": {
                        "fn": "aimnet.train.loss.peratom_loss_fn",
                        "weight": 0.05,
                        "kwargs": {"key_true": "charges", "key_pred": "charges"}
                    }
                }
            }
        },
        "optimizer": {
            "force_no_train": [],
            "force_train": [],
            "class": "torch.optim.RAdam",
            "kwargs": {
                "lr": Qlearningrate,
                "weight_decay": 1e-8
            },
            "param_groups": {
                "shifts": {
                    "re": ".*.atomic_shift.shifts.weight$",
                    "weight_decay": 0.0
                }
            }
        },
        "scheduler": {
            "class": "ignite.handlers.param_scheduler.ReduceLROnPlateauScheduler",
            "kwargs": {
                "metric_name": "loss",
                "factor": 0.75,
                "patience": Qearlystop
            },
            "terminate_on_low_lr": 1e-5
        },
        "trainer": {
            "trainer": "aimnet.train.utils.default_trainer",
            "evaluator": "aimnet.train.utils.default_evaluator",
            "epochs": nn_epochs
        },
        "checkpoint": {
            "dirname": "checkpoints",
            "filename_prefix": run_name,
            "kwargs": {
                "n_saved": 1,
                "require_empty": False
            }
        },
        "wandb": {
            "init": {
                "name": run_name,
                "mode": "offline",
                "entity": None,
                "project": None,
                "notes": None
            },
            "watch_model": {
                "log": "all",
                "log_freq": 1000,
                "log_graph": True
            }
        },
        "metrics": {
            "class": "aimnet.train.metrics.RegMultiMetric",
            "kwargs": {
                "cfg": {
                    "energy": {"abbr": "E", "scale": 23.06},
                    "dipole": {"abbr": "D", "scale": 1.0, "mult": 3},
                    "quadrupole": {"abbr": "Q", "scale": 1.0, "mult": 6},
                    "charges": {"abbr": "q", "peratom": True},
                    "volumes": {"abbr": "V", "peratom": True},
                    "forces": {"abbr": "F", "peratom": True, "mult": 3, "scale": 23.06}
                }
            }
        }
    }

    # Write config file
    with open("extra_conf.yaml", "w") as f:
        yaml.dump(conf, f, sort_keys=False)

    # Run training
    subprocess.run([
        "aimnet", "train",
        "data.train=dataset.h5",
        "data.sae.energy.file=dataset_sae.yaml",
        "--config", "extra_conf.yaml"
    ], env=env, check=True)

    # Find the model file in the checkpoints directory
    model_pattern = os.path.join("checkpoints", f"{run_name}*.pt")
    matching_files = glob.glob(model_pattern)
    print(matching_files)
    pt_path = matching_files[0]
    subprocess.run(["aimnet", "jitcompile", pt_path, varsoutfile], env=env, check=True)

####################################################################################################
####################################################################################################

def evaluate(varsoutfile, clusters_df, method, Qmin, parentdir, Qmonomers, Z_atoms):
    import os, sys
    import yaml
    from torch import cuda
    from numpy import array

    # Suppress pickle warnings from torch
    if 1==1:
        import warnings
        warnings.filterwarnings(
            "ignore",
            ".*which uses the default pickle module implicitly. It is possible to construct malicious pickle data which will execute arbitrary code during unpickling.*"
        )
        warnings.filterwarnings(
            "ignore",
            "PySisiphus is not installed"
        )
    
    from aimnet.calculators.aimnet2ase import AIMNet2ASE

    Y_predicted = []
    F_predicted = []

    chargeREM = None
    for i in range(len(clusters_df)):
        atoms = clusters_df["xyz"]["structure"].values[i].copy()

        # Determine charge
        if ("log", "charge") in clusters_df.columns:
            charge = clusters_df["log"]["charge"].values[i]
            if charge != chargeREM:
                chargeREM = charge
                calc = AIMNet2ASE(varsoutfile, charge=charge, mult=1)
        else:
            charge = 0
            calc = AIMNet2ASE(varsoutfile, charge=charge, mult=1)

        # Set up AIMNet calculator
        atoms.calc = calc

        energy = atoms.get_potential_energy()  # eV

        Y_predicted.append(0.0367493 * energy) #energy[0])  # Convert to Hartree
        F_predicted.append(0.0367493 * atoms.get_forces())  # Hartree/Angstrom

    Y_predicted = [array(Y_predicted)]
    F_predicted = [F_predicted]
    #print(Y_predicted)

    if method == "min":
        Y_predicted[0] += Qmin

    return Y_predicted, F_predicted

