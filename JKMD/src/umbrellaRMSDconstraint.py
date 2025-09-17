def RMSD(atoms):
    import numpy as np
    dist1 = np.linalg.norm(atoms.get_positions()[0]-atoms.get_positions()[1])
    dist2 = np.linalg.norm(atoms.get_positions()[1]-atoms.get_positions()[2])
    dist3 = np.linalg.norm(atoms.get_positions()[0]-atoms.get_positions()[2])
    #rmsd1 = ((dist1 - 1.089)**2 + (dist2 - 3.750)**2 + (dist3 - 1.089 - 3.75)**2)**0.5
    #rmsd2 = ((dist1 - 3.750)**2 + (dist2 - 1.930)**2 + (dist3 - 1.930 - 3.75)**2)**0.5
    rmsd = (dist1 - 1.38)/(3.75 - 1.38) + (dist2 - 3.06)/(1.93 - 3.06)   
    return rmsd
    #return rmsd1 - rmsd2
    #return (rmsd1 - rmsd2)/((3.750-1.089)**2+(3.75-1.930)**2+(1.089-1.930)**2)**0.5

def RMSD3(atoms):
    import numpy as np
    dist1 = np.linalg.norm(atoms.get_positions()[0]-atoms.get_positions()[1])
    dist2 = np.linalg.norm(atoms.get_positions()[1]-atoms.get_positions()[2])
    dist3 = np.linalg.norm(atoms.get_positions()[0]-atoms.get_positions()[2])
    rmsd = (dist1 - 1.38)/(3.75 - 1.38) + (dist2 - 3.06)/(1.93 - 3.06)   
    rmsd1 = ((dist1 - 1.38)**2 + (dist2 - 3.06)**2)**0.5
    rmsd2 = ((dist1 - 3.75)**2 + (dist2 - 1.93)**2)**0.5
    return rmsd, rmsd1, rmsd2
    #return rmsd1 - rmsd2, rmsd1, rmsd2
    #return (rmsd1 - rmsd2)/((3.750-1.089)**2+(3.75-1.930)**2+(1.089-1.930)**2)**0.5, rmsd1, rmsd2

def penalty(atoms):
    import numpy as np
    dist1 = np.linalg.norm(atoms.get_positions()[0]-atoms.get_positions()[1])
    dist2 = np.linalg.norm(atoms.get_positions()[1]-atoms.get_positions()[2])
    dist3 = np.linalg.norm(atoms.get_positions()[0]-atoms.get_positions()[2])
    rmsd = ((dist1 - 1.38)/(3.75 - 1.38) - (dist2 - 3.06)/(1.93 - 3.06))
    return rmsd
   
def numerical_derivative(atoms, REF, k, epsilon=1e-4):
    """
    Compute numerical derivatives using finite differences.
    :param func: Function that takes atoms and returns a scalar (e.g., RMSD)
    :param positions: np.ndarray of shape (N, 3), atomic positions
    :param epsilon: Small step for finite differences
    :return: np.ndarray of shape (N, 3), numerical gradient
    """
    import numpy as np
    gradient = np.zeros_like(atoms.get_positions())
    here = 0.5*k*(RMSD(atoms)-REF)**2
    here = 0.5*k*(RMSD(atoms)-REF)**2+0.5*100*penalty(atoms)**2

    for i in range(3): #range(len(gradient)):  # Loop over atoms
        for j in range(3):  # Loop over x, y, z coordinates
            atoms_forward = atoms.copy()
            positions_forward = atoms.get_positions()
            positions_forward[i, j] += epsilon
            atoms_forward.set_positions(positions_forward)
            #gradient[i, j] = ((0.5*k*(RMSD(atoms_forward)-REF)**2) - here) / epsilon
            gradient[i, j] = ((0.5*k*(RMSD(atoms_forward)-REF)**2 + 0.5*100*penalty(atoms_forward)**2) - here) / epsilon
    return gradient

class UmbrellaConstraint:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, a, k, r, adjuststeps = 0):
        self.a = a
        self.k = k*0.043364115308770496
        self.r = r
        self.step = 0
        self.adjusthalf = 0
        self.adjuststeps = adjuststeps
        if self.adjuststeps > 0:
          self.k_to_be = self.k
          self.k = self.k_to_be * self.step / self.adjuststeps
    def adjust_positions(self, atoms, newpositions):
        pass
    def adjust_potential_energy(self, atoms):
        import numpy as np
        CS = atoms.constraints
        #print(CS)
        del atoms.constraints
        #rmsd=RMSD(atoms)
        bias_en = 0.5*self.k*(RMSD(atoms)-self.r)**2
        #bias_en = 0.5*self.k*RMSD2(atoms,self.r)**2
        #print(bias_en)
        atoms.set_constraint(CS)
        return bias_en
    def adjust_forces(self, atoms, forces):
        import numpy as np
        CS = atoms.constraints
        del atoms.constraints
        grad = numerical_derivative(atoms, self.r, self.k) 
        #adjustment = -self.k*(RMSD(atoms)-self.r)*grad/np.linalg.norm(grad)
        forces += -grad
        #TODO make the forces applied as weighted distribution on atoms
        #forces[0:self.splitpoint] -= np.array(atoms[0:self.splitpoint].get_masses())[:,np.newaxis]/np.sum(atoms[0:self.splitpoint].get_masses())*adjustment
        #forces[self.splitpoint:] += np.array(atoms[self.splitpoint:].get_masses())[:,np.newaxis]/np.sum(atoms[self.splitpoint:].get_masses())*adjustment
        atoms.set_constraint(CS)
        if self.adjuststeps > self.step:
          if self.adjusthalf == 0:
            self.adjusthalf = 1 
          else:
            self.adjusthalf = 0
            self.step += 1
            self.k = self.k_to_be * self.step / self.adjuststeps
    def index_shuffle(self, atoms, ind):
        pass
