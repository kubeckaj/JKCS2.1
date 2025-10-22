import os, sys
from time import time
from typing import Iterable

os.environ['TORCHINDUCTOR_CACHE_DIR'] = '/scratch/hvehkama/kaharaja/cache'
cache = '/scratch/hvehkama/kaharaja/cache/uma-s-1p1-cache.pkl'

import numpy as np
from ase import Atoms
#from ase.io import read, write
from ase.units import *
#from ase.calculators.calculator import all_changes
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

import argparse
parser = argparse.ArgumentParser(description='Optimize structure with NNP')

# Add the arguments
parser.add_argument('--model', type=str, required=False, default='uma',
                    help='Model name')
parser.add_argument('--input', type=str, required=True, default=None,
                    help='Path to input .pkl or .dat file')
parser.add_argument('--out', type=str, required=True, default=None,
                    help='Path to output pickle file')
parser.add_argument('--restart', action="store_true", default=False,
                    help='Continue optimizing or add free energies to output')
parser.add_argument('--redo', action="store_true", default=False,
                    help='Rerun optimization or add free energy calculation')

parser.add_argument('--fmax', type=float, required=False,  default=1e-3,
                    help='Maximum force tolerance for optimization.')
parser.add_argument('--maxiter', type=float, required=False,  default=500,
                    help='Maximum number iterations to optimize.')

parser.add_argument('--freq', action="store_true", default=False,
                    help='Compute frequencies and Gibbs free energies')
parser.add_argument('--delta', type=float, required=False,  default=0.01,
                    help='')
parser.add_argument('--cutr', type=float, required=False,  default=np.nan,
                    help='Skip free energy calculation above relative cutoff in kcal/mol')
parser.add_argument('--temp', type=float, required=False,  default=298.15,
                    help='Temperature in K')
parser.add_argument('--pres', type=float, required=False,  default=101325.0,
                    help='Pressure in Pa')

parser.add_argument('--turbo', action="store_true", default=False,
                    help='Compiles UMA for speed but requires fixed atomic composition')

args = parser.parse_args()

method = args.model.split('/')[-1]
if 'EGRET' in method:
    from mace.calculators import mace_off
    calculator = mace_off(model=args.model+".model", default_dtype="float64",  device='cuda')
    program = "MACE_calculator"
elif 'uma' in method:
    from fairchem.core import pretrained_mlip, FAIRChemCalculator
    from fairchem.core.units.mlip_unit import InferenceSettings
    
    method = 'uma-s-1p1.pt'
    uma_path = '/scratch/hvehkama/kaharaja/UMA/uma-s-1p1.pt'
    program = "FAIRChemCalculator"
    if args.turbo:
        turbo_settings = InferenceSettings(
            tf32=True,
            activation_checkpointing=True,
            merge_mole=True,
            compile=True,
            wigner_cuda=False,
            external_graph_gen=False,
            #internal_graph_gen_version=2,
        )  
        model = None # load model later
        calculator = None
    else:
        model = pretrained_mlip.load_predict_unit(uma_path, device='cuda')
        calculator = FAIRChemCalculator(model, 'omol')
else:
    print('Model not found!'); exit()

from data_readers import read_pickled_data
df = read_pickled_data(args.input)

if args.restart and os.path.isfile(args.out):
    df_new = read_pickled_data(args.out)
elif os.path.isfile(args.out):
    print('Error! '+args.out+' already exists. Use --restart option to continue.')
    exit()
else:
    args.restart = False
    df_new = df[['info']].copy()

    df_new[('log', 'termination')] = False

if args.freq:
    if ('log', 'vibrational_frequencies') not in df_new.columns:
        df_new[('log', 'vibrational_frequencies')] = np.nan

    T=args.temp
    P=args.pres
    df_new[('log', 'temperature')] = T
    df_new[('log', 'pressure')] = P

if not args.restart:
    df_new[('xyz','structure')] = object()
    df_new[('extra','forces')] = object()

fmax = args.fmax
maxiter = args.maxiter

unique_cluster_types = np.unique(df[('info','cluster_type')])
if not np.isnan(args.cutr):
    df_indeces = np.array([],dtype=int)
    cutr = args.cutr * kcal/mol/Ha
    for ct in np.unique(df[('info','cluster_type')]):
        subset = df[df[('info','cluster_type')]==ct]
        elim = cutr + np.min(subset[('log', 'electronic_energy')])
        passed = subset[subset[('log', 'electronic_energy')] <= elim].index.values
        df_indeces = np.append(df_indeces, passed)
else:
    df_indeces = df.index.values


# select optimizer
from ase.optimize import BFGS
def optimize(atoms, calc, fmax=1e-3, maxiter=500):
    results = []
    #for atoms in list_of_atoms
    try:
        atoms.calc = calc

        # Do geometry optimization
        opt = BFGS(atoms)
        opt.run(fmax, maxiter)

        converged =  opt.converged()
        el = atoms.get_potential_energy()
        forces = atoms.get_forces()
        atoms.calc = None
    except:
        converged = False
        el = np.nan
        forces = np.nan

    results.append([atoms, converged, el, forces])
    return results

def frequencies(atoms, calc, options):

    return

i = 0
save_every_n = 1
for ct in np.unique(df[('info','cluster_type')][df_indeces]):
    subset = df[df[('info','cluster_type')]==ct].index.values
    print(ct)

    if 'uma' in method and args.turbo:
        try:
            import pickle
            import torch
            import torch._dynamo
            # Load precompiled cache for faster startup
            with open(cache, 'rb') as f:
                artifact_bytes = pickle.load(f)
            torch.compiler.load_cache_artifacts(artifact_bytes)
            #print('Loaded cache artifacts', len(artifact_bytes))
        except:
            print('Cache loading failed')
            
        del model
        # model is compiled for each cluster type
        model = pretrained_mlip.load_predict_unit(uma_path, inference_settings=turbo_settings, device='cuda')
        calculator = FAIRChemCalculator(model, 'omol')

    for idx in subset:
        name = df.loc[idx, ('info','file_basename')]
        atoms = df.loc[idx,('xyz','structure')]

        charge = df.loc[idx, ('log', 'charge')] if ('log', 'charge') in df.columns else 0
        spin = df.loc[idx, ('log', 'multiplicity')] if ('log', 'multiplicity') in df.columns else 1
        
        idx = df_new[df_new[('info','file_basename')]==name].index[0]
        if args.restart and (('xyz','structure') in df_new.columns):
            if isinstance(df_new.loc[idx,('xyz','structure')], Atoms):
                atoms = df_new.loc[idx,('xyz','structure')]
                ####print('Continuing... ')

        if not isinstance(atoms, Atoms):
            print('Error! No atoms for '+name)
            continue

        df_new.loc[idx, ('log', 'program')] = program
        df_new.loc[idx, ('log', 'method')] = method
        df_new.loc[idx, ('log', 'NAtoms')] = len(atoms)
        
        charge = 0 if np.isnan(charge) else int(charge)
        df_new.loc[idx, ('log', 'charge')] = charge
        spin = 1 if np.isnan(spin) else int(spin)
        df_new.loc[idx, ('log', 'multiplicity')] = spin

        has_gibbs = (('log','gibbs_free_energy') in df_new.columns) and (df_new.loc[idx, ('log','gibbs_free_energy')] < 0)
        if has_gibbs:
            converged = True
        elif (('extra', 'forces') in df_new.columns) and isinstance(df_new.at[idx, ('extra', 'forces')], Iterable):
            forces = df_new.at[idx, ('extra', 'forces')]
            converged = (fmax > np.linalg.norm(forces, axis=1).max()*Ha)
        else: 
            converged = False

        if args.redo:
            converged = False
            has_gibbs = False
            df_new.loc[idx, ('log','zero_point_correction')] = np.nan
            df_new.loc[idx, ('log','enthalpy_energy')] = np.nan
            df_new.loc[idx, ('log','entropy_energy')] = np.nan
            df_new.loc[idx, ('log','gibbs_free_energy')] = np.nan

        df_new.loc[idx, ('log', 'termination')] = converged

        if has_gibbs or (converged and not args.freq):
            continue

        atoms.calc = calculator
        atoms.info = {'charge': charge, 'spin': spin}

        t1 = time()
        print('Optimizing '+name+':')

        # do a geometry relaxation
        opt = BFGS(atoms)
        opt.run(fmax, maxiter)
        
        converged = opt.converged()
        if not converged:
            print('Optimizer failed to converge in '+str(maxiter)+' iteration steps.')
        df_new.loc[idx, ('log', 'electronic_energy')] = atoms.get_potential_energy() / Ha
        forces = atoms.get_forces()

        df_new.at[idx, ('xyz', 'structure')] = atoms
        df_new.at[idx, ('extra', 'forces')] = object()
        df_new.at[idx, ('extra', 'forces')] = forces / Ha

        print('Final electronic energy = '+str(df_new.loc[idx, ('log', 'electronic_energy')])+' Ha')

        if converged and args.freq:
            tmpdir = 'TMP/vib_'+name
            if os.path.exists(tmpdir):
                for f in os.listdir(tmpdir):
                    os.remove(tmpdir+'/'+f)

            print('Running Frequency calculations')

            vib = Vibrations(atoms, delta=args.delta, name=tmpdir)
            vib.run()
            
            freq = np.real(vib.get_frequencies())
            df_new.at[idx, ('log','vibrational_frequencies')] = object()
            df_new.at[idx, ('log','vibrational_frequencies')] = freq
            vib_energies = vib.get_energies()
            vib_energies[:6] = 0.0
            try:
                thermo = IdealGasThermo(
                    vib_energies=vib_energies,
                    potentialenergy=atoms.get_potential_energy(),
                    atoms=atoms,
                    geometry='nonlinear',
                    symmetrynumber=1,
                    spin=0
                )
                G = thermo.get_gibbs_energy(temperature=T, pressure=P)
                
                df_new.loc[idx, ('log','zero_point_correction')] = thermo.get_ZPE_correction() / Ha
                df_new.loc[idx, ('log','enthalpy_energy')] = thermo.get_enthalpy(T, verbose=False) / Ha
                df_new.loc[idx, ('log','entropy_energy')] = thermo.get_entropy(T, P, verbose=False) / Ha
                df_new.loc[idx, ('log','gibbs_free_energy')] = G / Ha
            except:
                print('Error! Thermochemistry Failed.')
                print(vib_energies)

                df_new.loc[idx, ('log','zero_point_correction')] = np.nan
                df_new.loc[idx, ('log','enthalpy_energy')] = np.nan
                df_new.loc[idx, ('log','entropy_energy')] = np.nan
                df_new.loc[idx, ('log','gibbs_free_energy')] = np.nan

            # remove temporary files
            if os.path.exists(tmpdir):
                for f in os.listdir(tmpdir):
                    os.remove(tmpdir+'/'+f)
                os.rmdir(tmpdir)

        t2 = time() - t1
        print('TOTAL RUN TIME: '+str(round(t2,2))+' sec')
        df_new.loc[idx, ('log', 'termination')] = converged
        df_new.loc[idx, ('log', 'time')] = t2/60.
        
        atoms.info = {}
        atoms.calc = None
        
        if i % save_every_n == 0:
            df_new.to_pickle(args.out+'.tmp.pkl')
            os.rename(args.out+'.tmp.pkl',args.out)
        i+=1

        print('')
        print('###########################################################')
        print('')

if args.turbo and not os.path.isfile(cache):
    import pickle
    import torch
    artifacts = torch.compiler.save_cache_artifacts()

    assert artifacts is not None
    artifact_bytes, cache_info = artifacts
    with open(cache, 'wb') as f:
        pickle.dump(artifact_bytes, f)
              
df_new.to_pickle(args.out+'.tmp.pkl')
os.rename(args.out+'.tmp.pkl',args.out)

