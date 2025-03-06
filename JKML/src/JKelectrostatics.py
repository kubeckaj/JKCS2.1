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
    #def switch(d):
    #    cut = sr_cut / 2
    #    x = d / cut
    #    x3 = x**3
    #    x4 = x3 * x
    #    x5 = x4 * x
    #    return np.where(d < cut, (6 * x5 - 15 * x4 + 10 * x3), 1.0)

    #def Dswitch(d):
    #    cut = sr_cut / 2
    #    x = d / cut
    #    x2 = x**2
    #    x3 = x2 * x
    #    x4 = x2**2
    #    return np.where(d < cut, (30 * x4 - 60 * x3 + 30 * x2), 0.0)

    #def switch(d):
    #    cut = sr_cut 
    #    x =  ( d / cut - 3.0 / 5.0  ) * 5.0 / 2.0
    #    return np.where( x - 1e-2 <= 0.0, 0.0, np.where(x + 1e-2 >= 1.0, 1.0, 0.5 * (np.tanh(5 * (2 * x - 1) / (2 * np.sqrt((1 - x) * x))) + 1)))

    #def Dswitch(d):
    #    cut = sr_cut 
    #    x = (d / cut - 3.0 / 5.0) * 5.0 / 2.0
    #    return np.where( x - 1e-2 <= 0.0, 0.0, np.where(x + 1e-2 >= 1.0, 0.0, 0.5 * (5 / np.sqrt((1 - x) * x) - (5 * (1 - 2 * x) * (-1 + 2 * x)) / (4 * ((1 - x) * x)**(3/2))) * (1 / np.cosh((5 * (-1 + 2 * x)) / (2 * np.sqrt((1 - x) * x))))**2))


    ## Compute switching function
    #switch_value = switch(r)
    #cswitch = 1.0 - switch_value
    #dswitch = Dswitch(r)

    ##Avoid self-interactions
    #np.fill_diagonal(r2, np.inf) 
    #np.fill_diagonal(r, np.inf) 
    #
    ## Compute electrostatic energy
    #E_ordinary = qiqj/r #np.where(r < sr_cut, 1.0 / sr_cut, 1.0 / r)
    #E_shielded = 0     #1.0 / (3/5*sr_cut) #1.0 / np.sqrt(r2 + 1.0)

    #energy = np.sum(np.triu(cswitch * E_shielded + switch_value * E_ordinary, k=1)) * conv_hartree

    ## Compute forces
    #F_ordinary = qiqj/r2  #np.where(r < sr_cut, 0, qiqj / r2 )
    #F_shielded = 0        #qiqj * np.where(r2 > 100000, 0.0, np.sqrt(r2 / (r2 + 1.0)**3)) 
    #force_magnitude = (cswitch * F_shielded + switch_value * F_ordinary + dswitch * E_ordinary - dswitch * E_shielded)

    def switch(x):
      x = x / sr_cut
      return 3*x**2 - 2*x**3

    def intsw(x):
      x = x / sr_cut
      return - x**3 + 3*x**4/4   

    np.fill_diagonal(r2, np.inf) 
    np.fill_diagonal(r, np.inf) 

    c = -sr_cut + 1/4    

    energy = np.sum(np.triu(np.where(r < sr_cut, -qiqj / sr_cut**2 * (intsw(r) + c), qiqj / r), k=1)) * conv_hartree    

    force_magnitude = np.where(r < sr_cut, switch(r) * qiqj / sr_cut**2, qiqj / r2 )
    forces = np.sum(force_magnitude[:, :, np.newaxis] * pos_diff / r[:, :, np.newaxis], axis=1) * conv_hartree**2

    return energy, forces
