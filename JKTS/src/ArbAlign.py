import re

def energy_cutoff(molecules, energy_cutoff=5):
    Htokcalmol = 627.509

    lowest_energy = sorted(molecules, key=lambda molecule: molecule.electronic_energy)[0].electronic_energy * Htokcalmol
    filtered_molecules = [m for m in molecules if m.electronic_energy*Htokcalmol - lowest_energy <= energy_cutoff]

    return filtered_molecules


def filter_molecules(molecules, logger=None, RMSD_threshold=0.34, Energy_threshold=1e-4, Dipole_threshold=1e-1):
    # `compare` used to be imported here; the function isn't present as a
    # stand-alone symbol. Use the Molecule.compare_structures method instead.
    initial_len = len(molecules)
    unique_molecules = []
    energy_difference = 0
    dipole_difference = 0

    while molecules:
        reference = molecules.pop(0) # Take first molecule and assume its unique for now
        unique_molecules.append(reference)

        non_similar_molecules = []

        for molecule in molecules:
            rmsd_check = False
            energy_check = False
            dipole_check = False

            # Calculate RMSD using the Molecule.compare_structures helper
            rmsd_value = reference.compare_structures(molecule)
            if rmsd_value <= RMSD_threshold:
                rmsd_check = True

            # Calculate energy and dipole moment difference if applicable
            if reference.electronic_energy is not None and molecule.electronic_energy is not None:
                energy_difference = abs(reference.electronic_energy - molecule.electronic_energy)
                if energy_difference <= Energy_threshold:
                    energy_check = True
            else:
                energy_check = True # If we can't compare energy, consider it as passing the check

            if reference.dipole_moment is not None and molecule.dipole_moment is not None:
                dipole_difference = abs(reference.dipole_moment - molecule.dipole_moment)
                if dipole_difference <= Dipole_threshold:
                    dipole_check = True
            else:
                dipole_check = True # Same for dipole

            if not logger:
                m_num = re.search(r'conf(\d+)', molecule.name).group(1)
                r_num = re.search(r'conf(\d+)', reference.name).group(1)
                print(f"R: conf{r_num:<3}  M: conf{m_num:<3} RMSD: {rmsd_value:.4f} E_diff: {energy_difference:.3e} D_diff: {dipole_difference:.3e} Identical: {all(i for i in [rmsd_check, energy_check, dipole_check])}")

            if rmsd_check and energy_check and dipole_check:
                if molecule.electronic_energy <= reference.electronic_energy:
                    unique_molecules.pop()
                    unique_molecules.append(molecule)
                    reference = molecule  # Update reference to the new molecule
            else:
                non_similar_molecules.append(molecule)

        molecules = non_similar_molecules

    if logger:
        logger.log(f"Filtered {initial_len} conformers to {len(unique_molecules)} conformers using RMSD threshold: {RMSD_threshold}, Energy difference threshold: {Energy_threshold} Hartree, and Dipole moment difference threshold: {Dipole_threshold} Debye")
    else:
        print(f"Filtered {initial_len} conformers to {len(unique_molecules)} conformers using RMSD threshold: {RMSD_threshold} Energy difference threshold: {Energy_threshold} Hartree Dipole moment difference threshold: {Dipole_threshold} Debye")

    return unique_molecules

