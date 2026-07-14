import os
import sys
from datetime import datetime

import runtime
from output import console
from ts_validation import is_aldehyde
from metadata import update_metadata


def crest_constrain(molecule, force_constant=1):
    if molecule.reactant or molecule.product:
        pass
    else:
        if not molecule.constrained_indexes:
            molecule.find_active_site(indexes=runtime.args.CHO)
        C_index = molecule.constrained_indexes['C']
        H_index = molecule.constrained_indexes['H']
        abstractor_index = molecule.constrained_indexes['X']

        aldehyde, aldehyde_O = is_aldehyde(molecule, C_index, H_index)
        CHO_angle = molecule.calculate_angle(molecule.coordinates[C_index-1], molecule.coordinates[H_index-1], molecule.coordinates[abstractor_index-1])

        with open(molecule.directory + "/constrain.inp", "w") as f:
            f.write("$constrain\n")
            f.write(f"  force constant={force_constant}\n")
            f.write(f"  distance: {C_index}, {H_index}, auto\n")
            f.write(f"  distance: {H_index}, {abstractor_index}, auto\n")
            f.write(f"  angle: {C_index}, {H_index}, {abstractor_index}, {CHO_angle}\n")
            if 'XH' in molecule.constrained_indexes and aldehyde:
                XH_index = molecule.constrained_indexes['XH']
                f.write(f"  dihedral: {aldehyde_O}, {C_index}, {abstractor_index}, {XH_index}, 0\n")
                f.write(f"  dihedral: {C_index}, {H_index}, {abstractor_index}, {XH_index}, 0\n")
            f.write("$end\n")


def mkdir(molecule, crest_constrain_flag=True):
    if not os.path.exists(molecule.directory):
        os.makedirs(molecule.directory, exist_ok=True)
        if not os.path.exists(os.path.join(molecule.directory, "input_files")):
            os.makedirs(os.path.join(molecule.directory, "input_files"))
        if not os.path.exists(os.path.join(molecule.directory, "log_files")):
            os.makedirs(os.path.join(molecule.directory, "log_files"))
        if not os.path.exists(os.path.join(molecule.directory, "slurm_output")):
            os.makedirs(os.path.join(molecule.directory, "slurm_output"))
    fields = {
        'molecule': molecule.name,
        'reaction': 'Cl' if runtime.args.Cl else 'NO3' if runtime.args.NO3 else 'OH',
        'method': runtime.args.method,
        'basis_set': runtime.args.basis_set,
        'program': runtime.QC_program,
        'slurm': {'par': runtime.args.par, 'cpu': runtime.args.cpu, 'mem': runtime.args.mem},
        'created': datetime.now().isoformat(timespec='seconds'),
        'jkts_command': ' '.join(sys.argv[1:]),
    }
    if crest_constrain_flag:
        # Active-site indices and σᵢ (conformers are rebuilt fresh downstream)
        fields['constrained_indexes'] = molecule.constrained_indexes
        fields['reaction_path_degeneracy'] = getattr(molecule, 'reaction_path_degeneracy', 1)
    update_metadata(molecule.directory, **fields)


def QC_input(molecule, constrain, TS, method=None, basis_set=None):
    file_name = f"{molecule.name}{molecule.input}"
    file_path = os.path.join(molecule.directory, file_name)
    atoms = molecule.atoms
    coords = molecule.coordinates
    maxiter = 150
    freq = 'freq' if TS or not (molecule.product or 'opt_constrain' in molecule.current_step) else 'freq'
    SCF = ''
    disp = 'EmpiricalDispersion=GD3BJ' if runtime.args.dispersion or runtime.args.method.upper() == 'B3LYP' else ''
    aldehyde = False

    method = method or molecule.method
    if method.lower() == 'wb97xd' and molecule.program.lower() == 'orca':
        method = 'wB97X-D3BJ'

    if not 'ccsd' in method:
        runtime.DFT_method = method

    methods_no_basis = {"b97-3c", "r2scan-3c", "pm3", "am1", "pm6", "pm7", 'g3mp2', 'g3'}
    if method.lower() in methods_no_basis:
        basis_set = ""
    else:
        basis_set = basis_set or runtime.args.basis_set

    if constrain or TS:
        if not molecule.constrained_indexes:
            molecule.find_active_site(indexes=runtime.args.CHO)
        if molecule.program.lower() == 'orca':
            C_index = molecule.constrained_indexes['C']-1
            H_index = molecule.constrained_indexes['H']-1
            abstractor_index = molecule.constrained_indexes['X']-1
        else:
            C_index = molecule.constrained_indexes['C']
            H_index = molecule.constrained_indexes['H']
            abstractor_index = molecule.constrained_indexes['X']

        aldehyde, aldehyde_O = is_aldehyde(molecule, C_index, H_index)

    if molecule.program.lower() == "orca" and molecule.converged is False:
        disp = 'D3BJ' if runtime.args.dispersion else ''
        with open(file_path, "w") as f:
            if molecule.current_step == "DLPNO":
                if runtime.args.F12:
                    molecule.method = 'f12-ccsd(t)'
                    method = f'! cc-pVTZ-F12 cc-pVTZ/C cc-pVTZ-F12-CABS DLPNO-CCSD(T)-F12 TightPNO TightSCF'
                else:
                    molecule.method = 'dlpno-ccsd(t)'
                    method = f'! aug-cc-pVTZ aug-cc-pVTZ/C DLPNO-CCSD(T) TightSCF RI-JK aug-cc-pVTZ/JK NoTRAH'
                f.write(f"{method}\n")
                f.write(f"%pal nprocs {runtime.args.cpu} end\n")
                if molecule.error_termination_count == 0:
                    f.write(f"%maxcore {round((runtime.args.mem+12000)/runtime.args.cpu)}\n")
                else:
                    f.write(f"%maxcore {round((runtime.args.mem+18000)/runtime.args.cpu)}\n")
            elif TS:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OptTS defgrid2 {freq}\n")
                f.write(f"%pal nprocs {runtime.args.cpu} end\n")
                f.write(f"%maxcore {round(runtime.args.mem/runtime.args.cpu)}\n")
            else:
                f.write(f"! {method} {basis_set} TightSCF SlowConv OPT defgrid2 {freq}\n")
                f.write(f"%pal nprocs {runtime.args.cpu} end\n")
                f.write(f"%maxcore {round(runtime.args.mem/runtime.args.cpu)}\n")
            if constrain:
                f.write("%geom\n")
                f.write("Constraints\n")
                f.write(f"{{B {C_index} {H_index} C}}\n")
                f.write(f"{{B {H_index} {abstractor_index} C}}\n")
                if aldehyde:
                    f.write(f"{{A {C_index} {H_index} {abstractor_index} C}}\n")
                f.write("end\nend\n")
            elif TS:
                f.write("%geom\n")
                f.write(f"maxiter {maxiter}\n")
                f.write("Calc_Hess true\n")
                if molecule.error_termination_count == 0:
                    f.write("Recalc_Hess 5\n")
                else:
                    f.write("Recalc_Hess 1\n")
                    f.write("Trust -0.2\n")
                if runtime.args.hybrid:
                    f.write(f"Hybrid_Hess {{ {C_index} {H_index} {abstractor_index} }} end\n")
                f.write(f"TS_Mode {{ B {H_index} {abstractor_index} }} end\n")
                f.write(f"TS_Active_Atoms {{ {C_index} {H_index} {abstractor_index} }} end\n")
                f.write("TS_Active_Atoms_Factor 3\n")
                f.write("end\n")
            f.write("\n")
            f.write(f"* xyz {molecule.charge} {molecule.mult}\n")
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("*")

    elif molecule.program.lower() == "g16" and molecule.converged is False:
        with open(file_path, "w") as f:
            f.write(f"%nprocshared={runtime.args.cpu}\n")
            f.write(f"%mem={runtime.args.mem}mb\n")
            if TS and constrain:
                f.write(f'# {method} {basis_set} opt=(tight,calcfc,ts,noeigen,modredundant,MaxCycles={maxiter}) {freq} {disp} int=UltraFine {SCF}\n\n')
            elif TS and constrain is False:
                if molecule.error_termination_count == 0:
                    f.write(f'# {method} {basis_set} opt=(tight,calcfc,ts,noeigen,MaxCycles={maxiter},ReCalcFC=5) {freq} {disp} int=UltraFine {SCF}\n\n')
                else:
                    f.write(f'# {method} {basis_set} opt=(tight,calcall,ts,noeigen,MaxCycles={maxiter}) {freq} {disp} int=UltraFine {SCF}\n\n')
            elif TS is False and constrain:
                f.write(f'# {method} {basis_set} opt=(tight,modredundant) {freq} {disp} int=UltraFine {SCF}\n\n')
            else:
                f.write(f"# {method} {basis_set} opt=tight {freq} {disp} int=UltraFine {SCF}\n\n")
            f.write("Title\n\n")
            f.write(f"{molecule.charge} {molecule.mult}\n")
            for atom, coord in zip(atoms, coords):
                f.write(f'{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n')
            f.write("\n")
            if constrain:
                f.write(f"B {C_index} {H_index} F\n")
                f.write(f"B {H_index} {abstractor_index} F\n")
                if aldehyde:
                    f.write(f"A {C_index} {H_index} {abstractor_index} F\n")
                    if 'conf' in molecule.name and 'XH' in molecule.constrained_indexes:
                        # XH_index is the H of the abstracting OH radical — only valid for OH abstraction
                        XH_index = molecule.constrained_indexes['XH']
                        f.write(f"D {aldehyde_O} {C_index} {abstractor_index} {XH_index} F\n")
                        f.write(f"D {C_index} {H_index} {abstractor_index} {XH_index} F\n")
                f.write("\n")
    else:
        console.warning(f"QC_input was called but no program was specified for {molecule.name}")
