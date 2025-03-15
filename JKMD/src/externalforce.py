def heaviside(x, x0):
  from numpy.linalg import norm
  return norm(x, axis = 1) >= x0

class ExternalForce:

    def __init__(self, QEF, QEF_par, QEF_systems):
        ### HARMONIC POTENTIAL ###
        if QEF == 'h_A' or QEF == 'h_A_xyz':
          self.qef = QEF
          self.k_ext = QEF_par[1]*0.0433634 #eV
          self.end = QEF_par[0]
          self.mfrom = QEF_systems[0]
          self.mto = QEF_systems[1]
        ### FLAT BOTTOM HARMONIC POTENTIAL ###
        if QEF == 'fbh_A' or QEF == 'fbh_A_xyz':
          self.qef = QEF
          self.end = QEF_par[0]             #target position in Angstrom
          self.k_ext = QEF_par[2]*0.0433634 #eV
          self.r0 = QEF_par[1]
          self.mfrom = QEF_systems[0]
          self.mto = QEF_systems[1] 
        ### CONSTANT EXT FORCE ###
        if QEF == 'c_COM':
          self.qef = QEF
          self.force = np.array(QEF_par)*0.0433634 #eV
          self.mfrom = QEF_systems[0]
          self.mto = QEF_systems[1]

    def adjust_positions(self, atoms, new):
        pass

    def index_shuffle(self, atoms, ind):
        pass

    def adjust_forces(self, atoms, forces):
        import numpy as np
        CS = atoms.constraints
        del atoms.constraints
        ### HARMONIC POTENTIAL ###
        if self.qef == 'h_A' or self.qef == 'h_A_xyz':
          if self.qef == 'h_A':
            centrum = atoms[self.mfrom:self.mto].get_center_of_mass()
          else:
            centrum = self.end
          vec = atoms[self.mfrom:self.mto].get_positions() - centrum
          external_force = -self.k_ext*vec
        ### FLAT BOTTOM HARMONIC POTENTIAL ###
        elif self.qef == 'fbh_A' or self.qef == 'fbh_A_xyz':
          if self.qef == 'fbh_A':
            centrum = atoms[self.mfrom:self.mto].get_center_of_mass()
          else:
            centrum = self.end
          vec = atoms[self.mfrom:self.mto].get_positions() - centrum
          external_force = - self.k_ext * heaviside(vec, self.r0).reshape(-1,1) * (vec - self.r0 * vec/np.linalg.norm(vec,axis=1).reshape(-1,1))
        ### CONSTANT EXT FORCE ###
        elif self.qef == 'c_COM':
          masses = atoms[self.mfrom:self.mto].get_masses()
          external_force = np.array(masses)[:,np.newaxis]/np.sum(masses) * self.force
        forces[self.mfrom:self.mto] += external_force
        atoms.set_constraint(CS)

    def adjust_potential_energy(self, atoms):
        import numpy as np
        CS = atoms.constraints
        del atoms.constraints
        ### HARMONIC POTENTIAL ###
        if self.qef == 'h_A' or self.qef == 'h_A_xyz':
          if self.qef == 'h_A':
            centrum = atoms[self.mfrom:self.mto].get_center_of_mass()
          else:
            centrum = self.end
          vec = atoms[self.mfrom:self.mto].get_positions() - centrum
          external_energy = 0.5 * self.k_ext * np.sum(vec**2)
        ### FLAT BOTTOM HARMONIC POTENTIAL ###
        if self.qef == 'fbh_A' or self.qef == 'fbh_A_xyz':
          if self.qef == 'fbh_A':
            centrum = atoms[self.mfrom:self.mto].get_center_of_mass()
          else:
            centrum = self.end
          vec = atoms[self.mfrom:self.mto].get_positions() - centrum
          external_energy = 0.5 * self.k_ext * np.sum(heaviside(vec, self.r0).reshape(-1,1) * (vec - self.r0 * vec/np.linalg.norm(vec,axis=1).reshape(-1,1))**2)
        ### CONSTANT EXT FORCE ###
        elif self.qef == 'c_COM':
          external_energy = 0
        atoms.set_constraint(CS)
        return external_energy

