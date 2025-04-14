import pandas as pd
from tblite.ase import TBLite
from joblib import Parallel, delayed
from ase import Atoms
from numpy import array
import gc
import sys
from os import environ

print("JKsp_grad started", flush=True)

# Number of cores
try:
    num_cores = int(environ['SLURM_JOB_CPUS_PER_NODE'])
except:
    from multiprocessing import cpu_count
    num_cores = cpu_count()
num_cores = min(num_cores, 8)

# Conversion constants
eV_to_Hartree = 1 / 27.2114
e_bohr_to_Debye = 2.541746

# Load the database
file = sys.argv[1]
db = pd.read_pickle(file).reset_index(drop=True)
structures = db.loc[:, ("xyz", "structure")].values

print("Database read", flush=True)
if len(structures) < num_cores:
    num_cores = len(structures)
print(f"Using {num_cores} cores", flush=True)

def JKsp_grad(num):
    struc0 = structures[num]
    struc = Atoms(symbols=struc0.get_chemical_symbols(), positions=struc0.get_positions())
    struc.calc = TBLite(method="GFN1-xTB", verbosity=0)

    try:
        energy_eV = struc.get_potential_energy()
        energy_Hartree = energy_eV * eV_to_Hartree
        forces_eV_per_A = struc.get_forces()
        forces_Hartree_per_A = [[x * eV_to_Hartree for x in f] for f in forces_eV_per_A]

        dipole_vector = struc.calc.results.get("dipole", None)
        if dipole_vector is not None:
            dipole_vector = array(dipole_vector)
            dipole_mag_Debye = (dipole_vector**2).sum()**0.5 * e_bohr_to_Debye
            dipole_vec_Debye = (dipole_vector * e_bohr_to_Debye).tolist()
        else:
            dipole_mag_Debye = None
            dipole_vec_Debye = None

        result = {
            "structure": struc.copy(),
            "energy": energy_Hartree,
            "forces": forces_Hartree_per_A,
            "dipole_moment": dipole_mag_Debye,
            "dipole_moments": dipole_vec_Debye,
        }
    except Exception as e:
        print(f"Error at index {num}: {e}", flush=True)
        result = {
            "structure": None,
            "energy": None,
            "forces": None,
            "dipole_moment": None,
            "dipole_moments": None,
        }

    del struc
    gc.collect()
    return result

# Batch execution
print("Start single point calculations", flush=True)
results = []
batch = 1000
for step in range(0, len(structures), batch):
    max_i = min(step + batch, len(structures))
    print(step, "-", max_i, flush=True)
    results_step = Parallel(n_jobs=num_cores, backend="loky")(
        delayed(JKsp_grad)(i) for i in range(step, max_i)
    )
    results.extend(results_step)

db_out = db.loc[:, pd.IndexSlice[['info'], :]].copy()

# Add missing columns if needed
for col in [
    ("xyz", "structure"),
    ("extra", "forces"),
    ("log", "program"),
    ("log", "method"),
    ("log", "dipole_moment"),
    ("log", "dipole_moments"),
    ("log", "electronic_energy"),
]:
    if col not in db_out.columns:
        db_out[col] = None

# Update database
for i, result in enumerate(results):
    if result["structure"] is not None:
        db_out.at[i, ("log", "program")] = "tblite"
        db_out.at[i, ("log", "method")] = "GFN1-xTB" 
        db_out.at[i, ("xyz", "structure")] = result["structure"]
        db_out.at[i, ("log", "electronic_energy")] = result["energy"]
        db_out.at[i, ("extra", "forces")] = result["forces"]
        db_out.at[i, ("log", "dipole_moment")] = result["dipole_moment"]
        db_out.at[i, ("log", "dipole_moments")] = result["dipole_moments"]

# Save the updated database
db_out.to_pickle("tblite_en_force_" + file)
print("Done", flush=True)

