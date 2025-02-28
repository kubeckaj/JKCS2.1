def compare_pair(arg, ArbAlign = 0, Qreturn_geometry = 0, REF = None):
    from ArbAlign import compare, kabsch, sorted_xyz
    if REF is None:
      exit()
    if ArbAlign:
      return compare(arg, REF, Qreturn_geometry = Qreturn_geometry)
    else:
      _, arg_pos, _, _ = sorted_xyz(arg, 0)
      #_, REF_pos, _, _ = sorted_xyz(REF, 0)
      REF_pos = REF.get_positions()
      return kabsch(arg_pos, REF_pos)

def numerical_derivative(func, atoms, REF, epsilon=1e-2):
    """
    Compute numerical derivatives using finite differences.
    :param func: Function that takes atoms and returns a scalar (e.g., RMSD)
    :param positions: np.ndarray of shape (N, 3), atomic positions
    :param epsilon: Small step for finite differences
    :return: np.ndarray of shape (N, 3), numerical gradient
    """
    import numpy as np
    gradient = np.zeros_like(atoms.get_positions())
    newref, here = func(atoms, 1, 1, REF = REF)
    
    bias_en = 0.5*100*(here-0)**2
    print("RMSD: " + str(here) + " Bias energy: " + str(bias_en * 23.0609) + " kcal/mol")
    #print(here)
    #print(func(atoms,1,REF = newref))
    #print(func(atoms,REF = newref))
    for i in range(len(gradient)):  # Loop over atoms
        for j in range(3):  # Loop over x, y, z coordinates
            atoms_forward = atoms.copy()
            positions_forward = atoms.get_positions()
            positions_forward[i, j] += epsilon
            atoms_forward.set_positions(positions_forward)
            #print(func(atoms_forward,REF = newref))
            gradient[i, j] = (func(atoms_forward, REF = newref) - here) / epsilon
    #print(gradient)
    return gradient

class RMSDConstraint:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, a, k, r, rmsdfile, adjuststeps = 0):
        from ase.io import read
        print("RMSD: " + rmsdfile)
        self.REF = read(rmsdfile)
        self.delay_step = 0
        self.a = a
        self.k = k*0.043364115308770496
        self.r = r
        self.remembered_forces = None
        self.step = 0
        self.adjusthalf = 0
        self.adjuststeps = adjuststeps
        if self.adjuststeps > 0:
          self.k_to_be = self.k
          self.k = self.k_to_be * self.step / self.adjuststeps
    def adjust_positions(self, atoms, newpositions):
        pass
    def adjust_potential_energy(self, atoms):
        #import numpy as np
        #CS = atoms.constraints
        #print(CS)
        #del atoms.constraints
        #rmsd = compare_pair(atoms, ArbAlign = 1)
        #rmsd = 0
        #TODO: testing some stupid adjustment
        #print(bias_en)
        #atoms.set_constraint(CS)
        return 0
    def adjust_forces(self, atoms, forces):
        import numpy as np
        if self.delay_step == 0:
          CS = atoms.constraints
          del atoms.constraints
          adjustment = numerical_derivative(compare_pair, atoms, self.REF)
          adjustment *= self.k
          #print(adjustment)
          maximumF = np.max(np.abs(forces))
          maximumA = np.max(np.abs(adjustment))
          factor = 1.1
          if maximumA < maximumF:
            adjustment *= factor*maximumF/maximumA
          else:
            adjustment[adjustment > maximumF] = factor*maximumF
            adjustment[adjustment < -maximumF] = -factor*maximumF
          forces -= adjustment
          forces[forces > maximumF] = maximumF
          forces[forces < -maximumF] = -maximumF
          atoms.set_constraint(CS)
          #print(forces[0])
          #print("-----")
          if self.adjuststeps > self.step:
            if self.adjusthalf == 0:
              self.adjusthalf = 1 
            else:
              self.adjusthalf = 0
              self.step += 1
              self.k = self.k_to_be * self.step / self.adjuststeps
          self.remembered_forces = adjustment
        else:
          forces -= self.remembered_forces
        if np.mod(self.delay_step, 1) == 0 and self.delay_step > 0:
          self.delay_step = 0
        else:
          self.delay_step += 1
    def index_shuffle(self, atoms, ind):
        pass
