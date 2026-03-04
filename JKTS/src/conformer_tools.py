import os
import re
from copy import deepcopy

from classes import Molecule
import runtime


def filter_molecules(molecules, logger=None, pickle=False, RMSD_threshold=0.34, Energy_threshold=1e-4, Dipole_threshold=1e-1):
    from ArbAlign import compare
    initial_len = len(molecules)
    unique_molecules = []
    energy_difference = 0
    dipole_difference = 0

    while molecules:
        reference = molecules.pop(0)
        unique_molecules.append(reference)

        non_similar_molecules = []

        for molecule in molecules:
            rmsd_check = False
            energy_check = False
            dipole_check = False

            rmsd_value = compare(reference, molecule)
            if rmsd_value <= RMSD_threshold:
                rmsd_check = True

            if reference.electronic_energy is not None and molecule.electronic_energy is not None:
                energy_difference = abs(reference.electronic_energy - molecule.electronic_energy)
                if energy_difference <= Energy_threshold:
                    energy_check = True
            else:
                energy_check = True

            if reference.dipole_moment is not None and molecule.dipole_moment is not None:
                dipole_difference = abs(reference.dipole_moment - molecule.dipole_moment)
                if dipole_difference <= Dipole_threshold:
                    dipole_check = True
            else:
                dipole_check = True

            if not logger:
                m_num = re.search(r'conf(\d+)', molecule.name).group(1)
                r_num = re.search(r'conf(\d+)', reference.name).group(1)
                print(f"R: conf{r_num:<3}  M: conf{m_num:<3} RMSD: {rmsd_value:.4f} E_diff: {energy_difference:.3e} D_diff: {dipole_difference:.3e} Identical: {all(i for i in [rmsd_check, energy_check, dipole_check])}")

            if rmsd_check and energy_check and dipole_check:
                if molecule.electronic_energy <= reference.electronic_energy:
                    unique_molecules.pop()
                    unique_molecules.append(molecule)
                    reference = molecule
            else:
                non_similar_molecules.append(molecule)

        molecules = non_similar_molecules

    if logger:
        logger.log(f"Filtered {initial_len} conformers to {len(unique_molecules)} conformers using RMSD threshold: {RMSD_threshold}, Energy difference threshold: {Energy_threshold} Hartree, and Dipole moment difference threshold: {Dipole_threshold} Debye")
    else:
        print(f"Filtered {initial_len} conformers to {len(unique_molecules)} conformers using RMSD threshold: {RMSD_threshold} Energy difference threshold: {Energy_threshold} Hartree Dipole moment difference threshold: {Dipole_threshold} Debye")

    if pickle:
        Molecule.molecules_to_pickle(unique_molecules, os.path.join(runtime.start_dir, "filtered_molecules.pkl"))
    return unique_molecules


def energy_cutoff(molecules):
    Htokcalmol = 627.509
    cutoff = runtime.args.energy_cutoff if runtime.args.energy_cutoff else 5
    lowest_energy = sorted(molecules, key=lambda molecule: molecule.electronic_energy)[0].electronic_energy * Htokcalmol
    filtered_molecules = [m for m in molecules if m.electronic_energy*Htokcalmol - lowest_energy <= cutoff]
    return filtered_molecules


def initiate_conformers(input_file=None):
    from pandas import read_pickle

    conformer_molecules = []

    if not input_file:
        print("No input file specified.")
        return []

    with open(input_file, 'rb') as f:
        df = read_pickle(f)

    for index, row in df.iterrows():
        el_energy = row[('log', 'electronic_energy')]
        xyz = row[('xyz', 'structure')]
        file_basename = row[('info', 'file_basename')]

        coords = xyz.get_positions()
        atoms = xyz.get_chemical_symbols()
        coordinates_list = [[atom] + list(coord) for atom, coord in zip(atoms, coords)]

        name = file_basename.replace("-str", "_conf")
        conformer_molecule = Molecule(name=name, directory=os.path.dirname(input_file), electronic_energy=el_energy, atoms=[], coordinates=[], program='CREST')

        if 'reactant' in input_file:
            conformer_molecule.reactant = True
            conformer_molecule.mult = 1
        elif 'product' in input_file:
            conformer_molecule.product = True
        else:
            conformer_molecule.reactant = False
            conformer_molecule.product = False

        conformer_molecule.electronic_energy = el_energy
        conformer_molecule.atoms = [atom[0] for atom in coordinates_list]
        conformer_molecule.coordinates = [atom[1:] for atom in coordinates_list]
        conformer_molecule.program = 'CREST'

        conformer_molecule.workflow = conformer_molecule.determine_workflow()
        conformer_molecule.set_current_step()
        conformer_molecule.method = conformer_molecule.log2method()

        conformer_molecules.append(conformer_molecule)

    return sorted(conformer_molecules, key=lambda m: m.electronic_energy)


def collect_DFT_and_DLPNO(molecules):
    collected_molecules = []

    for m in molecules:
        process_conditions = (m.current_step == 'optimization' and (m.reactant or m.product)) or ('TS_opt_conf' in m.current_step)

        if process_conditions:
            if m.name in ('OH', 'OH_DLPNO', 'H2O', 'H2O_DLPNO', 'Cl', 'Cl_DLPNO', 'HCl', 'HCl_DLPNO', 'NO3', 'NO3_DLPNO', 'HNO3', 'HNO3_DLPNO'):
                identifier = m.name.replace('_DLPNO', '')
            else:
                match = re.search(r'conf(\d+)', m.name)
                identifier = match.group() if match else None

            if identifier:
                pattern = re.compile(re.escape(identifier) + r'(?!\d)')
                for m_DLPNO in molecules:
                    if pattern.search(m_DLPNO.name) and m_DLPNO.current_step == 'DLPNO':
                        m.electronic_energy = m_DLPNO.electronic_energy
                        m.update_step()
                        m.name = m.name.replace('TS', 'DLPNO')
                        collected_molecules.append(m)
                        break

    return collected_molecules
