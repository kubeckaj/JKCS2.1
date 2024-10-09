class ExternalForce(FixConstraint):

    def __init__(self, from, a2, towards_point=[0,0,0], f_ext, measure):
        self.atom_from = a1
	self.atom_to = a2
        self.external_force = f_ext
        self.end = towards_point
	self.measure = measure

    def adjust_positions(self, atoms, new):
        pass

    def adjust_forces(self, atoms, forces):
        dist = np.subtract.reduce(atoms.positions[self.atom_from:self.atom_to])
        force = self.external_force * dist / np.linalg.norm(dist)
        forces[self.indices] += (force, -force)

    def adjust_potential_energy(self, atoms):
        dist = np.subtract.reduce(atoms.positions[self.indices])
        return -np.linalg.norm(dist) * self.external_force

    #def todict(self):
    #    return {'name': 'ExternalForce',
    #            'kwargs': {'a1': self.indices[0], 'a2': self.indices[1],
    #                       'f_ext': self.external_force}}
