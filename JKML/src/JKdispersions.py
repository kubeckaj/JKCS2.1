import torch
import tad_dftd4 as d4
import tad_dftd3 as d3
import tad_mctc as mctc
import numpy as np

# Conversion factors
ang_to_bohr = 1.0 / 0.52917721092


def compute_d4_energy_forces(positions: np.ndarray, symbols: np.ndarray, totalcharge: float = 0.0, 
                          s6: float = 1.0, s8: float = 1.85897750, s9: float = 1.0, 
                          a1: float = 0.44286966, a2: float = 4.60230534):
    """
    Compute Grimme D4 dispersion energy and forces for a given molecular structure.
    
    Parameters:
    positions : (N, 3) numpy array
        Cartesian coordinates of N atoms in Angstrom.
    symbols : (N) numpy array
        Symbols of N atoms as strings.
    totalcharge : float, optional
        Total charge of the system (default: 0.0).
    s6, s8, s9, a1, a2 : float, optional
        D4 dispersion parameters
    
    Returns:
    total_energy : float
        Total energy in Hartree.
    forces : (N, 3) numpy array
        Forces in Hartree/Ångström.
    """
    # Convert coordinates from Angstrom to Bohr
    pos = torch.tensor(positions * ang_to_bohr, dtype=torch.float64, requires_grad=True)
    
    # Convert symbols to atomic numbers
    numbers = mctc.convert.symbol_to_number(symbols=symbols.tolist())
    
    # Convert charge to torch tensor
    charge = torch.tensor(totalcharge, dtype=torch.float64)
    
    # TPSSh-D4-ATM parameters
    param = {
        "s6": pos.new_tensor(s6),
        "s8": pos.new_tensor(s8),
        "s9": pos.new_tensor(s9),
        "a1": pos.new_tensor(a1),
        "a2": pos.new_tensor(a2),
    }
    
    # Compute energy
    energy = d4.dftd4(numbers, pos, charge, param)
    total_energy = energy.sum().item()  # Convert to Python float
    
    # Compute forces via autograd
    (grad,) = torch.autograd.grad(energy.sum(), pos)
    forces = -grad.detach().numpy() / ang_to_bohr  # Convert from Hartree/Bohr to Hartree/Å
    
    return total_energy, forces


def compute_d3bj_energy_forces(positions: np.ndarray, symbols: np.ndarray, totalcharge: float = 0.0,
                          a1: float = 0.49484001, s8: float = 0.78981345, a2: float = 5.73083694):
    """
    Compute Grimme D3BJ dispersion energy and forces for a given molecular structure.

    Parameters:
    positions : (N, 3) numpy array
        Cartesian coordinates of N atoms in Angstrom.
    symbols : (N) numpy array
        Symbols of N atoms as strings.
    totalcharge : float, optional #ACTUALLY NOT USED, JUST KEPT FOR CONSISTENCY
    a1, s8, a2 : float, optional
        D3BJ dispersion parameters

    Returns:
    total_energy : float
        Total energy in Hartree.
    forces : (N, 3) numpy array
        Forces in Hartree/Ångström.
    """
    # Convert coordinates from Angstrom to Bohr
    pos = torch.tensor(positions * ang_to_bohr, dtype=torch.float64, requires_grad=True)

    # Convert symbols to atomic numbers
    numbers = mctc.convert.symbol_to_number(symbols=symbols.tolist())

    # D3BJ parameters
    param = {
        "a1": pos.new_tensor(a1),
        "s8": pos.new_tensor(s8),
        "a2": pos.new_tensor(a2),
    }

    # Compute energy
    energy = d3.dftd3(numbers, pos, param)
    total_energy = energy.sum().item()  # Convert to Python float

    # Compute forces via autograd
    (grad,) = torch.autograd.grad(energy.sum(), pos)
    forces = -grad.detach().numpy() / ang_to_bohr  # Convert from Hartree/Bohr to Hartree/Å

    return total_energy, forces

