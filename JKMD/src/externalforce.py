class ExternalForce(FixConstraint):

    def __init__(self, k_ext, towards_point=[0,0,0]):
        self.k_ext = k_ext
        self.end = towards_point

    def adjust_positions(self, atoms, new):
        pass

    def adjust_forces(self, atoms, forces):
        external_force = 2*self.k_ext*np.array([np.linalg.norm(i)*(np.array([0,0,0])-i)/(np.linalg.norm(np.array([0,0,0])-i)+1e-8) for i in pos])
        forces[self.indices] += force

    def adjust_potential_energy(self, atoms):
        external_energy = np.sum(self.k_ext*np.array([np.linalg.norm(np.array([0,0,0])-i)**2 for i in pos]))
        return external_energy

    #def todict(self):
    #    return {'name': 'ExternalForce',
    #            'kwargs': {'a1': self.indices[0], 'a2': self.indices[1],
    #                       'f_ext': self.external_force}}
