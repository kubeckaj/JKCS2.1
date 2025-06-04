import pandas as pd
import numpy as np
import h5py
import subprocess
from collections import defaultdict
from ase import Atoms
from ase.io import read

HARTREE_TO_EV = 27.211386245988  # CODATA 2018

def convert_pickle_to_aimnet_h5(df, output_file="database.h5"):
    grouped_data = defaultdict(lambda: defaultdict(list))

    print("Processing structures...")

    for i, row in df.iterrows():
        structure = row[("xyz", "structure")]
        atoms = structure if isinstance(structure, Atoms) else read(structure)

        natoms = len(atoms)
        group = f"{natoms:03d}"

        # Atomic info
        coords = atoms.get_positions()
        numbers = atoms.get_atomic_numbers()

        # Forces (in eV/Å)
        raw_forces = np.array(row[("extra", "forces")])
        dispersion_forces = np.array(row[("extra", "dispersion_forces")])
        forces = (raw_forces - dispersion_forces) * HARTREE_TO_EV

        # Energy (in eV)
        raw_energy = row[("log", "electronic_energy")]
        dispersion_energy = row[("extra", "dispersion_electronic_energy")]
        energy = (raw_energy - dispersion_energy) * HARTREE_TO_EV

        # Charges
        mulliken = np.array(row[("log", "mulliken_charges")])
        total_charge = row[("log", "charge")]

        # Basename
        basename = row[("info", "file_basename")]

        # Store in grouped structure
        grouped_data[group]["coord"].append(coords.astype(np.float32))
        grouped_data[group]["numbers"].append(numbers.astype(np.int32))
        grouped_data[group]["forces"].append(forces.astype(np.float32))
        grouped_data[group]["energy"].append(np.float32(energy))
        grouped_data[group]["charges"].append(mulliken.astype(np.float32))
        grouped_data[group]["charge"].append(np.float32(total_charge))
        grouped_data[group]["basename"].append(basename)

    print(f"Writing data to HDF5 file: {output_file}")

    str_dtype = h5py.string_dtype(encoding="utf-8")

    with h5py.File(output_file, "w") as f:
        for group, data in grouped_data.items():
            g = f.create_group(group)
            for key in ["coord", "numbers", "forces", "energy", "charges", "charge", "basename"]:
                arr = np.array(data[key])
                if key == "numbers":
                    g.create_dataset(key, data=arr.astype(np.int32))
                elif key == "basename":
                    arr = np.array(data[key], dtype=object)  # cast to object array of Python strings
                    g.create_dataset(key, data=arr, dtype=str_dtype)
                else:
                    g.create_dataset(key, data=arr.astype(np.float32))

    subprocess.run(["aimnet", "calc_sae", output_file, "dataset_sae.yaml"])

    print("Done.")

