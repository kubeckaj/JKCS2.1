import numpy as np


def compute_energies_forces(positions, charges, epsilon=1e-12):
    """
    Compute Coulomb energies and forces between atoms using vectorized operations.

    Parameters:
    positions : (N, 3) numpy array
        Cartesian coordinates of N atoms in Angstrom.
    charges : (N) numpy array
        Atomic charges of N atoms in elementary charge.
    epsilon : float, optional
        Small number to avoid division by zero.

    Returns:
    energy : float
        Total electrostatic energy in Hartree.
    forces : (N, 3) numpy array
        Forces on each atom due to Coulomb interactions in Hartree per Angstrom.

    Conversion:
    q1q2/r where the unit of q is 1.602177e-19 Coulomb and r is 1e-10 m: 2.566971e-28 C^2/m.
    Divide by 4pi*eps_0 (eps_0=vacuum permittivity with unit C^2*s^2/kg/m^3): 2.307078e-18 J.
    Hartree is 4.359745e-18 J. The conversion factor is then 0.529177.
    """

    conv_hartree = 0.529177

    # Create pairwise distance vectors
    pos_diff = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]  # (N, N, 3)

    # Compute pairwise distances
    r_sq = np.sum(pos_diff ** 2, axis=-1)  # (N, N), squared distances
    r_sq[r_sq < 5**2] = np.inf
    np.fill_diagonal(r_sq, np.inf)  # Remove self-interactions
    r = np.sqrt(r_sq)  # (N, N)

    # Compute pairwise charge products
    qiqj = np.outer(charges, charges)  # (N, N)

    # Compute energy (only upper triangle to avoid double counting)
    energy = np.sum(np.triu(qiqj / r, k=1)) * conv_hartree

    # Compute force magnitudes
    force_magnitude = qiqj / (r_sq + epsilon)  # Adding epsilon to avoid dividing by zero

    # Compute force vectors
    forces = np.sum(force_magnitude[:, :, np.newaxis] * pos_diff / r[:, :, np.newaxis], axis=1) * conv_hartree

    return energy, forces
