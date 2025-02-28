class QMMM:
    """QM/MM"""
    def __init__(self, QMMMspecies):
        from src.calculator import calculator
        self.calc_low = calculator(Qcalculator = "XTB", Qcalculator_input = "GFNFF", Qcharge = 0, Qout = 0)
        self.calc_high = calculator(Qcalculator = "XTB", Qcalculator_input = "GFN1-xTB", Qcharge = 0, Qout = 0)
        self.FROM=QMMMspecies[0]
        self.TO=QMMMspecies[1]
    def adjust_positions(self, atoms, newpositions):
        pass
    def adjust_potential_energy(self, atoms):
        CS = atoms.constraints
        del atoms.constraints
        tmp_atoms = atoms[self.FROM:self.TO].copy()
        atoms.set_constraint(CS)
        tmp_atoms.calc = self.calc_low
        Elow = tmp_atoms.get_potential_energy()
        tmp_atoms.calc = self.calc_high
        Ehigh = tmp_atoms.get_potential_energy()
        return Ehigh-Elow
    def adjust_forces(self, atoms, forces):
        CS = atoms.constraints
        del atoms.constraints
        tmp_atoms = atoms[self.FROM:self.TO].copy()
        atoms.set_constraint(CS)
        tmp_atoms.calc = self.calc_low
        Flow = tmp_atoms.get_forces()
        tmp_atoms.calc = self.calc_high
        Fhigh = tmp_atoms.get_forces()
        forces[self.FROM:self.TO] += Fhigh-Flow
    def index_shuffle(self, atoms, ind):
        pass
