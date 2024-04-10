import matplotlib.pyplot as plt
import numpy as np
import re
import pickle
import os

def plot_relative_energy(molecules_level1, molecules_level2):
    def extract_conf_number(molecule):
        match = re.search(r'conf(\d+)', molecule.name)
        return int(match.group(1)) if match else -1

    # Ensure the molecules are in the same order by configuration number
    molecules_level1.sort(key=extract_conf_number)
    molecules_level2.sort(key=extract_conf_number)

    # Determine the shortest length between the two sets to ensure matching pairs
    max_length = min(len(molecules_level1), len(molecules_level2))

    # Extract methods as labels
    xlabel = molecules_level1[0].method
    ylabel = molecules_level2[0].method
    Htokcalmol = 627.51

    # Vectorized operations to calculate energies
    energies_level1 = np.array([mol.electronic_energy for mol in molecules_level1[:max_length]]) * Htokcalmol
    energies_level2 = np.array([mol.electronic_energy for mol in molecules_level2[:max_length]]) * Htokcalmol

    # Calculating relative energies
    relative_energies_level1 = energies_level1 - np.min(energies_level1)
    relative_energies_level2 = energies_level2 - np.min(energies_level2)

    # Linear fit of the data
    fit = np.polyfit(relative_energies_level1, relative_energies_level2, 1)
    fit_fn = np.poly1d(fit)

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.scatter(relative_energies_level1, relative_energies_level2, color='blue', label='Data points')
    x_fit = np.linspace(0, max(relative_energies_level1), num=100)  # Ensuring smooth line
    plt.plot(x_fit, fit_fn(x_fit), color='red', label=f'Linear fit')  # Linear fit line

    plt.xlabel(f"Relative {xlabel} (kcal/mol)")
    plt.ylabel(f"Relative {ylabel} (kcal/mol)")
    plt.xlim([-0.1, max(relative_energies_level1) + 0.1])
    plt.ylim([-0.1, max(relative_energies_level2) + 0.1])

    # Calculate Pearson correlation coefficient
    r = np.corrcoef(relative_energies_level1, relative_energies_level2)[0, 1]
    plt.legend([f'Linear fit, Pearson r: {r:.2f}'])

    plt.savefig(f"relative_plot_{xlabel}_{ylabel}.png", dpi=500)
    plt.show()



def plot_rmsd(molecules, cache_file='.rmsd_cache.pkl'):
    rmsd_values = []
    
    # Check for cached RMSD values first
    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as f:
            rmsd_values = pickle.load(f)
    else:
        # Calculate RMSD values and store them in cache
        from ArbAlign import compare
        for i, mol1 in enumerate(molecules):
            for j, mol2 in enumerate(molecules):
                if i != j:
                    rmsd = compare(mol1, mol2)
                    rmsd_values.append(rmsd)
        
        with open(cache_file, 'wb') as f:
            pickle.dump(rmsd_values, f)
    
    sorted_rmsd_values = sorted(rmsd_values)
    x_values = list(range(1, len(sorted_rmsd_values) + 1))

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, sorted_rmsd_values, color='blue')
    plt.yscale('log')
    plt.xlabel('Pair number')
    plt.ylabel('RMSD (Å)')

    plt.show()
    plt.savefig(f"rmsd_values.png", dpi=500)


def plot_energy_dipole_differences(molecules):
    # Extract electronic energies and dipole moments into NumPy arrays
    electronic_energies = np.array([molecule.electronic_energy for molecule in molecules])
    dipole_moments = np.array([molecule.dipole_moment for molecule in molecules])
    
    # Calculate differences using broadcasting, avoiding redundant calculations
    energy_diff_matrix = np.abs(electronic_energies[:, np.newaxis] - electronic_energies)
    dipole_diff_matrix = np.abs(dipole_moments[:, np.newaxis] - dipole_moments)
    
    # Mask to remove diagonal and duplicate elements
    mask = np.triu(np.ones_like(energy_diff_matrix, dtype=bool), k=1)
    energy_diff_values = energy_diff_matrix[mask]
    dipole_diff_values = dipole_diff_matrix[mask]

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.scatter(energy_diff_values, dipole_diff_values, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy Difference (kcal/mol)')
    plt.ylabel('Dipole Moment Difference (Debye)')
    plt.show()
    plt.savefig(f"energy_dipole_difference.png", dpi=500)





