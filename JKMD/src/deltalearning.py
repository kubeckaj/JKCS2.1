class DeltaLearning:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, QEF = "deltalearning", QEF_par = "GFN1-xTB", QEF_systems = "all"):
        from src.calculator import calculator
        self.calls = 0
        self.pos = None     
        self.calc = calculator(Qcalculator = "XTB", Qcalculator_input = QEF_par, Qcharge = 0, Qout = 0)
        self.potential_energy = None
        self.forces = None  

    def adjust_positions(self, atoms, newpositions):
        pass

    def get_low_level(self, atoms):
        atoms.calc = self.calc
        self.potential_energy = atoms.get_potential_energy()
        self.forces = atoms.get_forces()
        self.pos = atoms.get_positions()
        self.calls += 1

    def adjust_potential_energy(self, atoms):
        CS = atoms.constraints 
        del atoms.constraints
        from numpy import allclose
        if self.calls == 0:
          self.get_low_level(atoms.copy())
        if not allclose(self.pos, atoms.get_positions()):
          self.get_low_level(atoms.copy())
        atoms.set_constraint(CS)
        #print(self.calls, flush = True)
        return self.potential_energy

    def adjust_forces(self, atoms, forces):
        CS = atoms.constraints
        del atoms.constraints
        from numpy import allclose
        if self.calls == 0:
          self.get_low_level(atoms.copy())
        if not allclose(self.pos, atoms.get_positions()):
          self.get_low_level(atoms.copy())
        atoms.set_constraint(CS)
        forces += self.forces

    def index_shuffle(self, atoms, ind):
        pass
