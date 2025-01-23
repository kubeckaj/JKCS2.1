# Standard imports
import torch
import numpy as np
import argparse
#ASE importations
import ase
from ase.calculators.calculator import Calculator
from ase.neighborlist import neighbor_list
#Neural network imports
from Neural_Net_evid import PhysNet
from layers.activation_fn import  *
''' 
Calculator for the atomic simulation environment (ASE) 
that evaluates energies and forces using a neural network.
'''

class PhysNetCalculator(Calculator):

    def __init__(self,
                 # ASE atoms object
                 atoms,
                 # ckpt file to restore the model (can also be a list for ensembles)
                 checkpoint,
                 # Respective config file for PhysNet architecture
                 config,
                 # System charge
                 charge=0,
                 # Cutoff distance for long range interactions (default: no cutoff)
                 lr_cut = None,
                 # Activation function
                 activation_fn="shift_softplus",
                 # Single or double precision
                 dtype=torch.float64):
        # Read config file to ensure same PhysNet architecture as during fit
        # Initiate parser
        parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

        # Add arguments
        parser.add_argument("--restart", type=str, default='No',
                            help="Restart training from a specific folder")
        parser.add_argument("--num_features", default=128, type=int)
        parser.add_argument("--num_basis", default=64, type=int)
        parser.add_argument("--num_blocks", default=5, type=int)
        parser.add_argument("--num_residual_atomic", default=2, type=int)
        parser.add_argument("--num_residual_interaction", default=3, type=int)
        parser.add_argument("--num_residual_output", default=1, type=int)
        parser.add_argument("--cutoff", default=10.0, type=float)
        parser.add_argument("--use_electrostatic", default=1, type=int)
        parser.add_argument("--use_dispersion", default=1, type=int)
        parser.add_argument("--grimme_s6", default=None, type=float)
        parser.add_argument("--grimme_s8", default=None, type=float)
        parser.add_argument("--grimme_a1", default=None, type=float)
        parser.add_argument("--grimme_a2", default=None, type=float)
        parser.add_argument("--dataset", type=str)
        parser.add_argument("--num_train", type=int)
        parser.add_argument("--num_valid", type=int)
        parser.add_argument("--batch_size", type=int)
        parser.add_argument("--valid_batch_size", type=int)
        parser.add_argument("--seed", default=None, type=int)
        parser.add_argument("--max_steps", default=10000, type=int)
        parser.add_argument("--learning_rate", default=0.001, type=float)
        parser.add_argument("--decay_steps", default=1000, type=int)
        parser.add_argument("--decay_rate", default=0.1, type=float)
        parser.add_argument("--max_norm", default=1000.0, type=float)
        parser.add_argument("--ema_decay", default=0.999, type=float)
        parser.add_argument("--rate", default=0.0, type=float)
        parser.add_argument("--l2lambda", default=0.0, type=float)
        parser.add_argument("--nhlambda", default=0.1, type=float)
        parser.add_argument("--lambda_conf",default=0.2,type=float)
        parser.add_argument("--summary_interval", default=5, type=int)
        parser.add_argument("--validation_interval", default=5, type=int)
        parser.add_argument("--show_progress", default=True, type=bool)
        parser.add_argument("--save_interval", default=5, type=int)
        parser.add_argument("--record_run_metadata", default=0, type=int)
        parser.add_argument('--device', default='cuda', type=str)

        # Read config file
        args = parser.parse_args(["@" + config])
        # Create neighborlist
        if lr_cut is None:
            self._sr_cutoff = args.cutoff
            self._lr_cutoff = None
            self._use_neighborlist = False
        else:
            self._sr_cutoff = args.cutoff
            self._lr_cutoff = lr_cut
            self._use_neighborlist = True

        # Periodic boundary conditions
        self.pbc = atoms.pbc
        self.cell = atoms.cell.diagonal()

        # Set up device
        self.device = args.device
        # Initiate calculator
        Calculator.__init__(self)

        # Set checkpoint file(s)
        self._checkpoint = checkpoint

        # Create PhysNet model
        self._model = PhysNet(
            F=args.num_features,
            K=args.num_basis,
            sr_cut=args.cutoff,
            num_blocks=args.num_blocks,
            num_residual_atomic=args.num_residual_atomic,
            num_residual_interaction=args.num_residual_interaction,
            num_residual_output=args.num_residual_output,
            use_electrostatic=(args.use_electrostatic == 1),
            use_dispersion=(args.use_dispersion == 1),
            s6=args.grimme_s6,
            s8=args.grimme_s8,
            a1=args.grimme_a1,
            a2=args.grimme_a2,
            writer=False,
            activation_fn=shifted_softplus,
            device=args.device)

        self._Z = torch.tensor(atoms.get_atomic_numbers(), dtype=torch.int32,device=self.device)
        self._R = torch.tensor(atoms.get_positions(), dtype=torch.float32,requires_grad=False,device=self.device)
        self._Q_tot = torch.tensor([charge],dtype=dtype,device=self.device)
        self._idx_i, self._idx_j = self.get_indices(atoms,device=self.device)



        def load_checkpoint(path):
            if path is not None:
                checkpoint = torch.load(path, map_location=self.device)
                return checkpoint
        # Load neural network parameter
        latest_ckpt = load_checkpoint(self.checkpoint)
        self._model.load_state_dict(latest_ckpt['model_state_dict'])
        self._model.eval()
        # Calculate properties once to initialize everything
        self._calculate_all_properties(atoms)
        # Set last_atoms to None as pcpot get enabled later and recalculation
        # becomes necessary again
        self._last_atoms = None

    def get_indices(self, atoms,device='cpu'):
        # Number of atoms
        N = len(atoms)
        # Indices pointing to atom at each batch image
        idx = torch.arange(end=N,dtype=torch.int32).to(device)
        # Indices for atom pairs ij - Atom i
        idx_i = idx.repeat(int(N) - 1)
        # Indices for atom pairs ij - Atom j
        idx_j = torch.roll(idx, -1, dims=0)

        if N>=2:
            for Na in torch.arange(2, N):
                Na_tmp = Na.cpu()
                idx_j = torch.concat(
                    [idx_j, torch.roll(idx, int(-Na_tmp.numpy()), dims=0)],
                    dim=0)
        idx_i = torch.sort(idx_i)[0]
        return idx_i.to(device), idx_j.to(device)

    def calculation_required(self, atoms):

        # Check positions, atomic numbers, unit cell and pbc
        return atoms != self.last_atoms

    def _calculate_all_properties(self, atoms):
        # find neighbors and offsets

        if self.use_neighborlist or any(atoms.get_pbc()):

            idx_i, idx_j, S = neighbor_list('ijS', atoms, self.lr_cutoff)
            offsets = np.dot(S, atoms.get_cell())
            sr_idx_i, sr_idx_j, sr_S = neighbor_list(
                'ijS', atoms, self.sr_cutoff)
            sr_offsets = np.dot(sr_S, atoms.get_cell())

        else:
            idx_i = self.idx_i
            idx_j = self.idx_j
            offsets = None
            sr_idx_i = None
            sr_idx_j = None
            sr_offsets = None

        # Calculate energy
        # (in case multiple NNs are used as ensemble, take the average)
        # Only one NN
        #with torch.no_grad():
        #    self.model.eval()
        #    self._last_energy, lambdas, alpha, beta = self.model.energy_evidential(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None,
        #                                    offsets=offsets, sr_idx_i=sr_idx_i, sr_idx_j=sr_idx_j,
        #                                    sr_offsets=sr_offsets)

        #    # self._last_energy, lambdas, alpha, beta = torch.split(out1,
        #    #                                                              out1.shape[1]//4,
        #    #                                                              dim=1)
        #    self._sigma2 = beta.detach().cpu().numpy()/(alpha.detach().cpu().numpy()-1)
        #    self._var = (1/lambdas.detach().cpu().numpy())*self.sigma2
        self.model.eval()
        self.R.requires_grad = True
        self._last_energy,forces = self.model.energy_and_forces(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None)
        #print(forces)
        #self.R.requires_grad = True
        #self._last_energy, lambdas, alpha, beta = self.model.energy_evidential(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None,offsets=offsets, sr_idx_i=sr_idx_i, sr_idx_j=sr_idx_j, sr_offsets=sr_offsets)
        #reduced_energy = torch.sum(self._last_energy.clone())
        #forces = -torch.autograd.grad([reduced_energy], [self.R], create_graph=True)[0]

        #print("HERE")
        #print(self._last_energy)
        #self._last_energy = energy_now
        #self._last_energy, self._last_forces, thrash_charges = self.model.energy_and_forces_and_charges(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None,offsets=offsets, sr_idx_i=sr_idx_i, sr_idx_j=sr_idx_j,sr_offsets=sr_offsets)
        #self.Z.requires_grad = True
        #thrash, self._last_forces = self.model.energy_and_forces(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None,offsets=offsets, sr_idx_i=sr_idx_i, sr_idx_j=sr_idx_j,sr_offsets=sr_offsets)

        #if self.last_forces is None:
        #self._last_energy = self.model.JKenergy(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None,offsets=offsets, sr_idx_i=sr_idx_i, sr_idx_j=sr_idx_j, sr_offsets=sr_offsets)
        #print(self._last_energy)
        #self._last_forces = torch.zeros_like(self.R)
 
        #print("self._last_energy:", self._last_energy)
        #print("self.R.requires_grad:", self.R.requires_grad)


        #if idx_i.numel() > 0:  # autograd will fail if there are no distances
        #    grad = torch.autograd.grad([self._last_energy], [self.R], create_graph=True)[0]
        #    if grad is not None:  # necessary for torch.jit compatibility
        #        forces = -grad
        #    else:
        #        forces = torch.zeros_like(self.R)
        #else:  # if there are no distances, the forces are zero
        #    forces = torch.zeros_like(self.R)
        #self._last_forces = forces
        #print(forces)

        #calculate atomic charges, energy and force evaluation nodes
        #Ea, Qa, Dij, nhloss = self.model.atomic_properties(self.Z, self.R, idx_i, idx_j, offsets, sr_idx_i, sr_idx_j, sr_offsets)
        #self._charges = self.model.scaled_charges(self.Z, Qa, self.Q_tot)
        #self._energy, self._forces = self.model.energy_and_forces_from_scaled_atomic_properties(Ea, self._charges, Dij, self.Z, self.R, idx_i, idx_j)
        #self._last_forces = self._forces

        # Convert results to numpy array
        self._last_energy = self.last_energy.detach().cpu().numpy()
        self._last_forces = forces.detach().cpu().numpy()

        # prevents some problems... but not for me, it actually does one
        # self._last_energy = np.array(1*[self.last_energy])

        # Store a copy of the atoms object
        self._last_atoms = atoms.copy()

    def get_potential_energy(self, atoms, force_consistent=False):

        if self.calculation_required(atoms):
            self._calculate_all_properties(atoms)

        return self.last_energy #, self.variance, self.sigma2

    def get_forces(self, atoms):
 
        if self.calculation_required(atoms):
            self._calculate_all_properties(atoms)
        return self.last_forces

    @property
    def last_atoms(self):
        return self._last_atoms

    @property
    def last_energy(self):
        return self._last_energy

    @property
    def last_forces(self):
        return self._last_forces

    @property
    def variance(self):
        return self._var

    @property
    def sigma2(self):
        return self._sigma2


    @property
    def sr_cutoff(self):
        return self._sr_cutoff

    @property
    def lr_cutoff(self):
        return self._lr_cutoff

    @property
    def use_neighborlist(self):
        return self._use_neighborlist

    @property
    def model(self):
        return self._model

    @property
    def checkpoint(self):
        return self._checkpoint

    @property
    def Z(self):
        return self._Z

    @property
    def Q_tot(self):
        return self._Q_tot

    @property
    def R(self):
        return self._R

    @property
    def idx_i(self):
        return self._idx_i

    @property
    def idx_j(self):
        return self._idx_j

    @property
    def energy(self):
        return self._energy

    @property
    def forces(self):
        return self._forces
