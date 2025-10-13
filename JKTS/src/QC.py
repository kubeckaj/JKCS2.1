import os
import re
from types import SimpleNamespace
from .utils import is_aldehyde
from .config import config


def _default_args(molecule):
    """Return a namespace with safe default args for QC input generation."""
    # Defaults mirror cli.py
    return SimpleNamespace(
        dispersion=False,
        F12=False,
        cpu=4,
        mem=8000,
        hybrid=False,
        method=getattr(molecule, 'method', 'wB97XD'),
        basis_set='6-31++g(d,p)',
        CHO=None
    )


def _ensure_args(args, molecule):
    """Ensure an argparse-like object with required attributes exists.

    If args is None, try config.args, else create defaults. If provided args is
    missing attributes, fill them from defaults.
    """
    if args is None:
        args = config.args
    if args is None:
        args = _default_args(molecule)
    # Fill missing attributes
    defaults = _default_args(molecule)
    for key, value in defaults.__dict__.items():
        if not hasattr(args, key) or getattr(args, key) is None:
            setattr(args, key, value)
    return args


def QC_input(molecule, constrain, TS, args=None):
    # Allow callers to pass args explicitly or rely on the shared config with safe defaults
    args = _ensure_args(args, molecule)
    global DFT_method
    file_name = f"{molecule.name}{molecule.input}"
    file_path = os.path.join(molecule.directory, file_name)
    atoms = molecule.atoms
    coords = molecule.coordinates
    maxiter = 150 # Maximum iterations for geoemtry optimization
    freq = 'freq' if TS or not (molecule.product or 'opt_constrain' in molecule.current_step) else 'freq'
    SCF = '' #'NoTrah'
    disp = 'EmpiricalDispersion=GD3BJ' if args.dispersion else ''
    aldehyde = False

    # Normalize method casing and set defaults
    method = (args.method or molecule.method)
    mlow = str(method).lower() if method else ''
    if mlow == 'wb97xd' and molecule.program.lower() == 'orca':
        method = 'wb97x-d3bj'

    if 'ccsd' not in str(method).lower():
        DFT_method = method
    
    # Basis set handling with safe default
    basis_set = getattr(args, 'basis_set', '6-31++g(d,p)')
    methods_no_basis = {"b97-3c", "r2scan-3c", "pm3", "am1", "pm6", "pm7", 'g3mp2', 'g3'}
    if str(method).lower() in methods_no_basis:
        basis_set = ""

    if constrain or TS:
        if not molecule.constrained_indexes:
            molecule.find_active_site(indexes=args.CHO)
        if molecule.program.lower() == 'orca': # ORCA indexes from 0
            C_index = molecule.constrained_indexes['C']-1
            H_index = molecule.constrained_indexes['H']-1
            O_index = molecule.constrained_indexes['O']-1
        else:
            C_index = molecule.constrained_indexes['C']
            H_index = molecule.constrained_indexes['H']
            O_index = molecule.constrained_indexes['O']

        aldehyde, aldehyde_O = is_aldehyde(molecule, C_index, H_index)

    if molecule.program.lower() == "orca" and molecule.converged is False:
        disp = 'D3BJ' if args.dispersion else ''
        with open(file_path, "w") as f:
            if molecule.current_step == "DLPNO":
                if args.F12:
                    molecule.method = 'f12-ccsd(t)'
                    method = '! cc-pVTZ-F12 cc-pVTZ/C cc-pVTZ-F12-CABS DLPNO-CCSD(T)-F12 TightPNO TightSCF'
                else:
                    molecule.method = 'dlpno-ccsd(t)'
                    method = f'! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) VeryTightSCF RI-JK aug-cc-pVTZ/JK {SCF}'
                # Consider scaling of the memory with the molecule size
                f.write(f"{method}\n")
                f.write(f"%pal nprocs {args.cpu} end\n")
                if molecule.error_termination_count == 0:
                    f.write(f"%maxcore {round((args.mem+12000)/args.cpu)}\n") 
                else:
                    f.write(f"%maxcore {round((args.mem+18000)/args.cpu)}\n") 
            elif TS:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OptTS defgrid2 {freq}\n")
                f.write(f"%pal nprocs {args.cpu} end\n")
                f.write(f"%maxcore {round(args.mem/args.cpu)}\n")
            else:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OPT defgrid2 {freq}\n")
                f.write(f"%pal nprocs {args.cpu} end\n")
                f.write(f"%maxcore {round(args.mem/args.cpu)}\n")
            if constrain:
                f.write("%geom\n")
                f.write("Constraints\n")
                f.write(f"{{B {C_index} {H_index} C}}\n") # Bond being broken
                f.write(f"{{B {H_index} {O_index} C}}\n") # Bond being formed
                if aldehyde:
                    f.write(f"{{A {C_index} {H_index} {O_index} C}}\n") # CHO angle
                f.write("end\nend\n")
            elif TS:
                f.write("%geom\n")
                f.write(f"maxiter {maxiter}\n")
                f.write("Calc_Hess true\n")
                if molecule.error_termination_count == 0:
                    f.write("Recalc_Hess 5\n")
                else:
                    f.write("Recalc_Hess 1\n")
                    f.write("Trust -0.2\n") # Maybe -0.1?
                if args.hybrid:
                    f.write(f"Hybrid_Hess {{ {C_index} {H_index} {O_index} }} end\n")
                f.write(f"TS_Mode {{ B {H_index} {O_index} }} end\n")
                f.write(f"TS_Active_Atoms {{ {C_index} {H_index} {O_index} }} end\n")
                f.write("TS_Active_Atoms_Factor 3\n")
                f.write("end\n")
            f.write("\n")
            f.write(f"* xyz {molecule.charge} {molecule.mult}\n") 
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("*")
########################################################G16################################################################
    elif molecule.program.lower() == "g16" and molecule.converged is False:
        with open(file_path, "w") as f:
            f.write(f"%nprocshared={args.cpu}\n")
            f.write(f"%mem={args.mem}mb\n")
            if TS and constrain: # should only be in the case of hard convergence problems where some flexible parts should be constrained.
                f.write(f'# {method} {basis_set} opt=(tight,calcfc,ts,noeigen,modredundant,MaxCycles={maxiter}) {freq} {disp} int=UltraFine\n\n')
            elif TS and constrain is False:
                if molecule.error_termination_count == 0:
                    f.write(f'# {method} {basis_set} opt=(tight,calcfc,ts,noeigen,MaxCycles={maxiter},ReCalcFC=5) {freq} {disp} int=UltraFine\n\n') # RecalcFC=N also option, recalc Hessian every N iteration
                else:
                    f.write(f'# {method} {basis_set} opt=(tight,calcall,ts,noeigen,MaxCycles={maxiter}) {freq} {disp} int=UltraFine\n\n')
            elif TS is False and constrain:
                f.write(f'# {method} {basis_set} opt=(tight,modredundant) {freq} {disp} int=UltraFine\n\n')
            else:
                f.write(f"# {method} {basis_set} opt=tight {freq} {disp} int=UltraFine\n\n")
            f.write("Title\n\n")
            f.write(f"{molecule.charge} {molecule.mult}\n")
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("\n")
            if constrain:
                f.write(f"B {C_index} {H_index} F\n")
                f.write(f"B {H_index} {O_index} F\n")
                if aldehyde:
                    f.write(f"A {C_index} {H_index} {O_index} F\n") # CHO angle
                    if 'conf' in molecule.name:
                        f.write(f"D {aldehyde_O} {C_index} {O_index} {O_index+1} F\n") 
                        f.write(f"D {C_index} {H_index} {O_index} {O_index+1} F\n") 
                f.write("\n")
    else:
        print(f"QC_input was called but no program was specified for {molecule.name}")


def extract_normal_coordinates(molecule):
    '''Return format: [[-1.02, -0.02, 0.06], [0.07, -0.11, -0.02], ...]''' 
    if molecule.program.lower() == 'g15':
        with open(molecule.log_file_path, 'r') as file:
            lines = file.readlines()
        
        negative_freq_found = False
        read_xyz = False  # Flag to start reading XYZ coordinates
        xyz_coordinates = []
        for line in lines:
            if "Frequencies --" in line and '-' in line:
                negative_freq_found = True
                continue
            if negative_freq_found and not read_xyz:
                if "Atom  AN      X      Y      Z" in line:  # This indicates the start of the XYZ section
                    read_xyz = True  # Now start reading XYZ coordinates on subsequent lines
                continue  # Skip all lines until XYZ section is reached
            if read_xyz:
                # Assuming the coordinates section ends before a line not matching the XYZ format
                match = re.match(r'\s*\d+\s+\d+\s+(-?\d+\.\d+\s+)*', line)
                if match:
                    xyz = re.findall(r'-?\d+\.\d+', line)[:2]  # Extract first three floats
                    xyz = [float(coord) for coord in xyz]
                    xyz_coordinates.append(xyz)
                else:
                    break 
        return xyz_coordinates[:-2]
    else:
        print("extract_normal_coordinates() was called, but only works for G15 log files")