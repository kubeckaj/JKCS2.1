import numpy as np

def compute_energies_forces(positions, charges, sr_cut=5.0):
    """
    Compute Coulomb energies and forces between atoms using vectorized operations.

    Parameters:
    positions : (N, 3) numpy array
        Cartesian coordinates of N atoms in Angstrom.
    charges : (N) numpy array
        Atomic charges of N atoms in elementary charge.
    sr_cut : float, optional
        Cut-off distance for switching function in Angstrom.

    Returns:
    energy : float
        Total electrostatic energy in Hartree.
    forces : (N, 3) numpy array
        Forces on each atom due to Coulomb interactions in Hartree per Angstrom.
    """
    conv_hartree = 0.529177  # Conversion factor to Hartree

    # Compute pairwise distance vectors and distances
    pos_diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    r2 = np.sum(pos_diff**2, axis=-1)
    r = np.sqrt(r2)

    # Compute pairwise charge products
    qiqj = np.outer(charges, charges)

    # Define switch function
    def switch(d):
        cut = sr_cut / 2
        x = d / cut
        x3 = x**3
        x4 = x3 * x
        x5 = x4 * x
        return np.where(d < cut, 6 * x5 - 15 * x4 + 10 * x3, 1.0)

    # Compute switching function
    switch_value = switch(r)
    cswitch = 1.0 - switch_value

    #Avoid self-interactions
    np.fill_diagonal(r2, np.inf) 
    r = np.sqrt(r2)
    
    # Compute electrostatic energy
    E_ordinary = 1.0 / r
    E_shielded = 1.0 / np.sqrt(r + 1.0)

    energy = np.sum(np.triu(qiqj * (cswitch * E_shielded + switch_value * E_ordinary), k=1)) * conv_hartree

    # Compute forces
    F_ordinary = qiqj / r2
    F_shielded = qiqj / (r2 + 2*r + 1.0) 
    force_magnitude = (cswitch * F_shielded + switch_value * F_ordinary)
    forces = - np.sum(force_magnitude[:, :, np.newaxis] * pos_diff / r[:, :, np.newaxis], axis=1) * conv_hartree**2

    return energy, forces
