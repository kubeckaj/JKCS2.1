class UmbrellaConstraint:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, a, k, splitpoint, r):
        self.a = a
        self.k = k*0.043364115308770496
        self.splitpoint = splitpoint
        self.r = r
    def adjust_positions(self, atoms, newpositions):
        pass
    def adjust_potential_energy(self, atoms):
        import numpy as np
        CS = atoms.constraints
        #print(CS)
        del atoms.constraints
        vec_COM = atoms[0:self.splitpoint].get_center_of_mass()-atoms[self.splitpoint:].get_center_of_mass()
        norm_vec_COM = vec_COM/np.sqrt(np.sum(vec_COM**2))
        bias_en = 0.5*self.k*(np.sqrt(np.sum(vec_COM**2))-self.r)**2
        #print(bias_en)
        atoms.set_constraint(CS)
        return bias_en
    def adjust_forces(self, atoms, forces):
        import numpy as np
        CS = atoms.constraints
        del atoms.constraints
        vec_COM = atoms[0:self.splitpoint].get_center_of_mass()-atoms[self.splitpoint:].get_center_of_mass()
        norm_vec_COM = vec_COM/np.sqrt(np.sum(vec_COM**2))
        adjustment = self.k*(-self.r*norm_vec_COM+vec_COM)
        #print(adjustment)
        #print(forces[0])
        #TODO make the forces applied as weighted distribution on atoms
        forces[0:self.splitpoint] -= np.array(atoms[0:self.splitpoint].get_masses())[:,np.newaxis]/np.sum(atoms[0:self.splitpoint].get_masses())*adjustment
        forces[self.splitpoint:] += np.array(atoms[self.splitpoint:].get_masses())[:,np.newaxis]/np.sum(atoms[self.splitpoint:].get_masses())*adjustment
        atoms.set_constraint(CS)
        #print(forces[0])
        #print("-----")
    def index_shuffle(self, atoms, ind):
        pass
