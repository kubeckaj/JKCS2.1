class EnergyConstraint:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, energy, ref_energy):
        self.energy = energy
        self.ref_energy = ref_energy
    def adjust_positions(self, atoms, newpositions):
        pass
    def adjust_potential_energy(self, atoms):
        return 0
    def adjust_forces(self, atoms, forces):
        from numpy import max,abs
        adjustmet = (self.energy - self.ref_energy)**(-3)*(self.energy - self.ref_energy)*forces 
        forces -= 0.5*max(abs(forces))/max(abs(adjustmet))*adjustmet 
    def index_shuffle(self, atoms, ind):
        pass
