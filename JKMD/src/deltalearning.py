class DeltaLearning:
    """Constrain an atom to move along a given direction only."""
    def __init__(self, QEF = "deltalearning", QEF_par = "GFN1-xTB", QEF_systems = "all"):
        from src.calculator import calculator
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

    def adjust_potential_energy(self, atoms):
        CS = atoms.constraints 
        del atoms.constraints
        if self.pos != atoms.get_positions():
          get_low_level(atoms.copy())
        atoms.set_constraint(CS)
        return self.potential_energy

    def adjust_forces(self, atoms, forces):
        import numpy as np
        CS = atoms.constraints
        del atoms.constraints
        if self.pos != atoms.get_positions():
          get_low_level(atoms.copy())
        atoms.set_constraint(CS)
        forces += self.forces

    def index_shuffle(self, atoms, ind):
        pass
