import torch
import torch.nn as nn
import torch.nn.functional as fnc
from tensorboardX import SummaryWriter
import numpy as np
from layers.utils import segment_sum
from layers.RBFLayer import RBFLayer
from layers.InteractionBlock import InteractionBlock
from layers.OutputBlock import OutputBlock
from layers.activation_fn import *
from grimme_d3.grimme_d3 import *

def softplus_inverse(x):
    '''numerically stable inverse of softplus transform'''
    return x + np.log(-np.expm1(-x))


def gather_nd(params, indices):
    '''
    the input indices must be a 2d tensor in the form of [[a,b,..,c],...],
    which represents the location of the elements.

    This function comes from:
        https://discuss.pytorch.org/t/implement-tf-gather-nd-in-pytorch/37502/6
    '''
    # Normalize indices values
    params_size = list(params.size())

    assert len(indices.size()) == 2
    assert len(params_size) >= indices.size(1)

    # Generate indices
    indices = indices.t().long()
    ndim = indices.size(0)
    idx = torch.zeros_like(indices[0]).long()
    m = 1

    for i in range(ndim)[::-1]:
        idx = idx + indices[i] * m
        m *= params.size(i)

    params = params.reshape((-1, *tuple(torch.tensor(params.size()[ndim:]))))
    return params[idx]

class PhysNet(nn.Module):

    def __init__(self,
                 # Dimensionality of feature vector
                 F=128,
                 # Number of radial basis functions
                 K=64,
                 # Cutoff distance for short range interactions
                 sr_cut=10.0,
                 # Cutoff distance for long range interactions
                 # (default: no cutoff)
                 lr_cut=None,
                 # Number of building blocks to be stacked
                 num_blocks=5,
                 # Number of residual layers for atomic refinements of
                 # feature vector
                 num_residual_atomic=2,
                 # Number of residual layers for refinement of message vector
                 num_residual_interaction=3,
                 # Number of residual layers for the output blocks
                 num_residual_output=1,
                 # Adds electrostatic contributions to atomic energy
                 use_electrostatic=True,
                 # Adds dispersion contributions to atomic energy
                 use_dispersion=True,
                 # s6 coefficient for d3 dispersion, by default is learned
                 s6=None,
                 # s8 coefficient for d3 dispersion, by default is learned
                 s8=None,
                 # a1 coefficient for d3 dispersion, by default is learned
                 a1=None,
                 # a2 coefficient for d3 dispersion, by default is learned
                 a2=None,
                 # Initial value for output energy shift
                 # (makes convergence faster)
                 Eshift=0.0,
                 # Initial value for output energy scale
                 # (makes convergence faster)
                 Escale=1.0,
                 # Initial value for output charge shift
                 Qshift=0.0,
                 # Initial value for output charge scale
                 Qscale=1.0,
                 # Half (else double counting) of the Coulomb constant
                 # (default is in units e=1, eV=1, A=1)
                 kehalf=7.199822675975274,
                 # Activation function
                 activation_fn=shifted_softplus,
                 # Single or double precision
                 dtype=torch.float32,
                 # Rate for dropout,
                 rate=0.0,
                 # Device to use
                 device="cuda",
                 #Summary writter
                 writer = None
                 ):
        super(PhysNet, self).__init__()

        assert (num_blocks > 0)
        self.num_blocks = num_blocks
        self.dtype = dtype
        self.kehalf = kehalf
        self.F = F
        self.K = K
        self.sr_cut = sr_cut  # cutoff for neural network interactions
        self.lr_cut = lr_cut  # cutoff for long-range interactions
        self.use_electrostatic = use_electrostatic
        self.use_dispersion = use_dispersion
        self.activation_fn = activation_fn
        self.rate = rate

        #Check if your model can be passed to cuda.
        if device=="cuda":
            # assert torch.cuda.is_available()
            cuda_device = torch.device("cuda")
            self.device = cuda_device
        else:
            #print('You do not select cuda, code will be run in cpu')
            #print('Your calculations might be slow.')
            self.device ="cpu"

        # Atom embeddings (we go up to Pu(94): 95 - 1 ( for index 0))
        self.embeddings = nn.Parameter(torch.empty(95, self.F,device=self.device).uniform_(-np.sqrt(3), np.sqrt(3)).requires_grad_(True))

        # Initialize the radial basis functions
        self.rbf_layer = RBFLayer(K, sr_cut,device=self.device)
        # Initialize variables for d3 dispersion (the way this is done,
        # positive values are guaranteed)
        if s6 is None:
            self.s6 = nn.Parameter(fnc.softplus(
                torch.tensor(softplus_inverse(d3_s6), requires_grad=True, dtype=dtype, device=self.device)))
        else:
            self.s6 = torch.tensor(s6, requires_grad=False, dtype=dtype, device=self.device)


        if s8 is None:
            self.s8 = nn.Parameter(fnc.softplus(
                torch.tensor(softplus_inverse(d3_s8), requires_grad=True, dtype=dtype, device=self.device)))
        else:
            self.s8 = torch.tensor(s8, requires_grad=False, dtype=dtype, device=self.device)


        if a1 is None:
            self.a1 = nn.Parameter(fnc.softplus(
                torch.tensor(softplus_inverse(d3_a1), requires_grad=True, dtype=dtype, device=self.device)))
        else:
            self.a1 = torch.tensor(a1, requires_grad=False, dtype=dtype, device=self.device)


        if a2 is None:
            self.a2 = nn.Parameter(fnc.softplus(
                torch.tensor(softplus_inverse(d3_a2), requires_grad=True, dtype=dtype, device=self.device)))
        else:
            self.a2 = torch.tensor(a2, requires_grad=False, dtype=dtype, device=self.device)


        if writer is None:
            self.writer = SummaryWriter()
        elif writer == False:
            pass
        else:
            self.writer = writer
            self.writer.add_histogram("embeddings", self.embeddings, 0)
            self.writer.add_scalar("d3-s6", self.s6)
            self.writer.add_scalar("d3-s8", self.s8)
            self.writer.add_scalar("d3-a1", self.a1)
            self.writer.add_scalar("d3-a2", self.a2)

        # # Initialize output scale/shift variables
        self.Eshift = nn.Parameter(torch.empty(95,device=self.device).new_full((95,), Eshift).type(dtype))
        self.Escale = nn.Parameter(torch.empty(95,device=self.device).new_full((95,), Escale).type(dtype))
        self.Qshift = nn.Parameter(torch.empty(95,device=self.device).new_full((95,), Qshift).type(dtype))
        self.Qscale = nn.Parameter(torch.empty(95,device=self.device).new_full((95,), Qscale).type(dtype))

        # Output scale for extra variables
        self.ascale = nn.Parameter(torch.ones(95, device=self.device, dtype=dtype))
        self.bscale = nn.Parameter(torch.ones(95, device=self.device, dtype=dtype))
        self.lscale = nn.Parameter(torch.ones(95, device=self.device, dtype=dtype))


        self.interaction_block = nn.ModuleList([InteractionBlock(
            K, F, num_residual_atomic, num_residual_interaction,
            activation_fn=self.activation_fn, rate=self.rate,device=self.device)
            for _ in range(self.num_blocks)])

        self.output_block = nn.ModuleList([OutputBlock(
            F, num_residual_output, n_output=1, activation_fn=self.activation_fn, rate=self.rate,device=self.device)
            for _ in range(self.num_blocks)])

        self.output_block_evid = nn.ModuleList([OutputBlock(
            F, num_residual_output, n_output=4, activation_fn=self.activation_fn, rate=self.rate,device=self.device)
            for _ in range(self.num_blocks)])


        self.build_requires_grad_dict()
    def train(self, mode=True):
        """ Turn on training mode. """
        super(PhysNet, self).train(mode=mode)
        for name, param in self.named_parameters():
            param.requires_grad = self.requires_grad_dict[name]

    def eval(self):
        super(PhysNet,self).eval()
        for name, param in self.named_parameters():
            param.requires_grad = False

    def build_requires_grad_dict(self):
        """
        Build a dictionary of which parameters require gradient information (are
        trained). Can be manually edited to freeze certain parameters)
        """
        self.requires_grad_dict = {}
        for name, param in self.named_parameters():
            self.requires_grad_dict[name] = param.requires_grad

    # ------------------------------------
    # Evidential layer
    # ------------------------------------

    def evidential_layer(self, loglambdas, logalphas, logbetas):
        min_val = 1e-6
        lambdas = torch.nn.Softplus()(loglambdas) + min_val
        alphas = torch.nn.Softplus()(logalphas) + min_val + 1  # add 1 for numerical contraints of Gamma function
        betas = torch.nn.Softplus()(logbetas) + min_val

        return lambdas, alphas, betas


    def calculate_interatomic_distances(self, R, idx_i, idx_j, offsets=None):
        ''' Calculate interatomic distances '''

        Ri = torch.gather(R, 0, idx_i.type(torch.int64).view(-1, 1).repeat(1, 3))
        Rj = torch.gather(R, 0, idx_j.type(torch.int64).view(-1, 1).repeat(1, 3))
        if offsets is not None:
            Rj = Rj + offsets
        p = nn.ReLU(inplace=True)
        m = p(torch.sum((Ri - Rj) ** 2, dim=-1))
        Dij = torch.sqrt(m)
        # ReLU: y = max(0, x), prevent negative sqrt
        return Dij
    @torch.jit.export
    def evidential_atomic_properties(self, Z, R, idx_i, idx_j, offsets=None, sr_idx_i=None, sr_idx_j=None,
                    sr_offsets=None):
        ''' Calculate evidential atomic properties '''

        # Calculate distances (for long range interaction)
        Dij_lr = self.calculate_interatomic_distances(R, idx_i, idx_j, offsets=offsets)
        # Optionally, it is possible to calculate separate distances
        # for short range interactions (computational efficiency)
        if sr_idx_i is not None and sr_idx_j is not None:
            Dij_sr = self.calculate_interatomic_distances(R, sr_idx_i, sr_idx_j, offsets=sr_offsets)
        else:
            sr_idx_i = idx_i
            sr_idx_j = idx_j
            Dij_sr = Dij_lr

        # Calculate radial basis function expansion
        rbf = self.rbf_layer(Dij_sr.to(self.device))

        # Initialize feature vectors according to embeddings for
        # nuclear charges
        z_pros = Z.view(-1, 1).expand(-1, self.F).type(torch.int64)
        x = torch.gather(self.embeddings, 0, z_pros)
        # Apply blocks
        Ea = 0  # atomic energy
        Qa = 0  # atomic charge
        lambdas, alpha, beta = 0, 0, 0
        nhloss = 0  # non-hierarchicality loss
        for i in range(self.num_blocks):
            x = self.interaction_block[i](x, rbf, sr_idx_i, sr_idx_j)
            out = self.output_block_evid[i](x)
            Ea = Ea + out[:, 0]
            Qa = Qa + out[:, 4]
            lambdas = lambdas + out[:, 1]
            alpha = alpha + out[:, 2]
            beta = beta + out[:, 3]
            # Compute non-hierarchicality loss
            out2 = out ** 2
            if i > 0:
                nhloss = nhloss + torch.mean(out2 / (out2 + lastout2 + 1e-7))
            lastout2 = out2
            # Apply scaling/shifting
        Ea = self.Escale[Z.type(torch.int64)] * Ea \
            + self.Eshift[Z.type(torch.int64)]

        Qa = self.Qscale[Z.type(torch.int64)] * Qa \
            + self.Qshift[Z.type(torch.int64)]

        lambdas = self.lscale[Z.type(torch.int64)] * lambdas
        alpha = self.ascale[Z.type(torch.int64)]* alpha
        beta = self.bscale[Z.type(torch.int64)] * beta

        return Ea,lambdas, alpha, beta, Qa, Dij_lr, nhloss

    @torch.jit.export
    def atomic_properties(self, Z, R, idx_i, idx_j, offsets=None, sr_idx_i=None, sr_idx_j=None,
                    sr_offsets=None):
        ''' Calculates the atomic energies, charges and distances
                    (needed if unscaled charges are wanted e.g. for loss function) '''

        # Calculate distances (for long range interaction)
        Dij_lr = self.calculate_interatomic_distances(R, idx_i, idx_j,offsets=offsets)
        # Optionally, it is possible to calculate separate distances
        # for short range interactions (computational efficiency)
        if sr_idx_i is not None and sr_idx_j is not None:
            Dij_sr = self.calculate_interatomic_distances(R, sr_idx_i, sr_idx_j, offsets=sr_offsets)
        else:
            sr_idx_i = idx_i
            sr_idx_j = idx_j
            Dij_sr = Dij_lr

        # Calculate radial basis function expansion
        rbf = self.rbf_layer(Dij_sr.to(self.device)).to(self.device)

        # Initialize feature vectors according to embeddings for
        # nuclear charges
        z_pros = Z.view(-1, 1).expand(-1, self.F).type(torch.int64)
        x = torch.gather(self.embeddings, 0, z_pros)
        # Apply blocks
        Ea = 0  # atomic energy
        Qa = 0  # atomic charge
        nhloss = 0  # non-hierarchicality loss
        for i in range(self.num_blocks):
            x = self.interaction_block[i](x, rbf, sr_idx_i, sr_idx_j)
            out = self.output_block[i](x)
            Ea = Ea + out[:, 0]
            Qa = Qa + out[:, 1]
            # Compute non-hierarchicality loss
            out2 = out ** 2
            if i > 0:
                nhloss = nhloss + torch.mean(out2 / (out2 + lastout2 + 1e-7))
            lastout2 = out2

        # Apply scaling/shifting
        Ea = self.Escale[Z.type(torch.int64)] * Ea \
             + self.Eshift[Z.type(torch.int64)]

        # Last term necessary to guarantee no "None" in force evaluation
        Qa = self.Qscale[Z.type(torch.int64)] * Qa \
                 + self.Qshift[Z.type(torch.int64)]


        return Ea, Qa, Dij_lr, nhloss


    @torch.jit.export
    def energy_evidential_from_scaled_atomic_properties(
            self, Ea,lambdas, alpha, beta, Qa, Dij, Z, idx_i, idx_j, batch_seg=None):
        ''' Calculates the energy given the scaled atomic properties (in order
            to prevent recomputation if atomic properties are calculated) '''
        if batch_seg is None:
            batch_seg = torch.zeros_like(Z).type(torch.int64)

        # Add electrostatic and dispersion contribution to atomic energy
        if self.use_electrostatic:
            Ea = Ea + self.electrostatic_energy_per_atom(Dij, Qa, idx_i, idx_j)
        if self.use_dispersion:
            if self.lr_cut is not None:
                Ea = Ea + d3_autoev * edisp(Z, Dij / d3_autoang, idx_i, idx_j,
                                            s6=self.s6, s8=self.s8, a1=self.a1, a2=self.a2,
                                            cutoff=self.lr_cut / d3_autoang, device=self.device)
            else:
                Ea = Ea + d3_autoev * edisp(Z, Dij / d3_autoang, idx_i, idx_j,
                                            s6=self.s6, s8=self.s8, a1=self.a1, a2=self.a2,device=self.device)

        Ea = torch.squeeze(segment_sum(Ea,batch_seg,device=self.device))
        lambdas = torch.squeeze(segment_sum(lambdas,batch_seg,device=self.device))
        alpha = torch.squeeze(segment_sum(alpha,batch_seg,device=self.device))
        beta = torch.squeeze(segment_sum(beta,batch_seg,device=self.device))

        lambdas, alpha, beta = self.evidential_layer(lambdas, alpha, beta)

        return Ea,lambdas,alpha,beta

    @torch.jit.export
    def energy_from_scaled_atomic_properties(
            self, Ea, Qa, Dij, Z, idx_i, idx_j, batch_seg=None):
        ''' Calculates the energy given the scaled atomic properties (in order
            to prevent recomputation if atomic properties are calculated) '''
        if batch_seg is None:
            batch_seg = torch.zeros_like(Z)

        # Add electrostatic and dispersion contribution to atomic energy
        if self.use_electrostatic:
            Ea = Ea + self.electrostatic_energy_per_atom(Dij, Qa, idx_i, idx_j)
        if self.use_dispersion:
            if self.lr_cut is not None:
                Ea = Ea + d3_autoev * edisp(Z, Dij / d3_autoang, idx_i, idx_j,
                                            s6=self.s6, s8=self.s8, a1=self.a1, a2=self.a2,
                                            cutoff=self.lr_cut / d3_autoang, device=self.device)
            else:
                Ea = Ea + d3_autoev * edisp(Z, Dij / d3_autoang, idx_i, idx_j,
                                            s6=self.s6, s8=self.s8, a1=self.a1, a2=self.a2,device=self.device)

        Ea = torch.squeeze(segment_sum(Ea,batch_seg,device=self.device))
        return Ea

    @torch.jit.export
    def energy_and_forces_from_scaled_atomic_properties(
            self, Ea, Qa, Dij, Z, R, idx_i, idx_j, batch_seg=None, create_graph=True):
        ''' Calculates the energy and forces given the scaled atomic atomic
            properties (in order to prevent recomputation if atomic properties
            are calculated .
            Calculation of the forces was done following the implementation of
            spookynet. '''

        Dij = self.calculate_interatomic_distances(
            R, idx_i, idx_j)

        energy = self.energy_from_scaled_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, batch_seg)

        reduced_energy = torch.sum(energy.clone())
        if idx_i.numel() > 0:  # autograd will fail if there are no distances
            grad = torch.autograd.grad([reduced_energy], [R], create_graph=create_graph)[0]

            if grad is not None:  # necessary for torch.jit compatibility
                forces = -grad
            else:
                forces = torch.zeros_like(R)
        else:  # if there are no distances, the forces are zero
            forces = torch.zeros_like(R)

        return energy, forces

    @torch.jit.export
    def energy_evidential_from_atomic_properties(
            self, Ea, lambdas, alpha, beta, Qa, Dij, Z, idx_i, idx_j, Q_tot=None, batch_seg=None):
        ''' Calculates the energy given the atomic properties (in order to
            prevent recomputation if atomic properties are calculated) '''

        if batch_seg is None:
            batch_seg = torch.zeros_like(Z).type(torch.int64)

            # Scale charges such that they have the desired total charge
        Qa = self.scaled_charges(Z, Qa, Q_tot, batch_seg)

        return self.energy_evidential_from_scaled_atomic_properties(
            Ea, lambdas, alpha, beta, Qa, Dij, Z, idx_i, idx_j, batch_seg)


    @torch.jit.export
    def energy_from_atomic_properties(
            self, Ea, Qa, Dij, Z, idx_i, idx_j, Q_tot=None, batch_seg=None):
        ''' Calculates the energy given the atomic properties (in order to
            prevent recomputation if atomic properties are calculated) '''

        if batch_seg is None:
            batch_seg = torch.zeros_like(Z).type(torch.int64)

        # Scale charges such that they have the desired total charge
        Qa = self.scaled_charges(Z, Qa, Q_tot, batch_seg)

        return self.energy_from_scaled_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, batch_seg)

    @torch.jit.export
    def energy_and_forces_from_atomic_properties(
            self, Ea, Qa, Dij, Z, R, idx_i, idx_j, Q_tot=None, batch_seg=None, create_graph=True):
        ''' Calculates the energy and force given the atomic properties
            (in order to prevent recomputation if atomic properties are
            calculated) '''

        Dij = self.calculate_interatomic_distances(
            R, idx_i, idx_j)

        energy = self.energy_from_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, Q_tot, batch_seg)

        if idx_i.numel() > 0:  # autograd will fail if there are no distances
            grad = torch.autograd.grad([energy.sum()], [R], create_graph=create_graph)[0]
            if grad is not None:  # necessary for torch.jit compatibility
                forces = -grad
            else:
                forces = torch.zeros_like(R)
        else:  # if there are no distances, the forces are zero
            forces = torch.zeros_like(R)

        return energy, forces

    @torch.jit.export
    def energy_evidential(self, Z, R, idx_i, idx_j, Q_tot=None,batch_seg=None, offsets=None,
            sr_idx_i=None, sr_idx_j=None, sr_offsets=None):
        ''' Calculates the total energy (including electrostatic
            interactions) '''
        Ea, lambdas, alpha, beta, Qa, Dij, _ = self.evidential_atomic_properties(
            Z, R, idx_i, idx_j, offsets, sr_idx_i, sr_idx_j, sr_offsets)

        energy, lambdas, alpha, beta = self.energy_evidential_from_atomic_properties(
            Ea, lambdas, alpha, beta, Qa, Dij, Z, idx_i, idx_j, Q_tot, batch_seg)

        return energy, lambdas, alpha, beta


    @torch.jit.export
    def energy(self, Z, R, idx_i, idx_j, Q_tot=None, batch_seg=None, offsets=None,
            sr_idx_i=None, sr_idx_j=None, sr_offsets=None):
        ''' Calculates the total energy (including electrostatic
            interactions) '''

        Ea, Qa, Dij, _ = self.atomic_properties(
            Z, R, idx_i, idx_j, offsets, sr_idx_i, sr_idx_j, sr_offsets)

        energy = self.energy_from_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, Q_tot, batch_seg)

        return energy

    @torch.jit.export
    def energy_and_forces(
            self, Z, R, idx_i, idx_j, Q_tot=None, batch_seg=None, offsets=None,
            sr_idx_i=None, sr_idx_j=None, sr_offsets=None, create_graph=True):
        ''' Calculates the total energy and forces (including electrostatic
            interactions)'''
        Ea, Qa, Dij, _ = self.atomic_properties(Z, R, idx_i, idx_j, offsets,
                                                sr_idx_i, sr_idx_j, sr_offsets)

        energy = self.energy_from_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, Q_tot, batch_seg)

        reduced_energy = torch.sum(energy)


        if idx_i.numel() > 0:  # autograd will fail if there are no distances
            grad = torch.autograd.grad([reduced_energy], [R], create_graph=create_graph)[0]
            if grad is not None:  # necessary for torch.jit compatibility
                forces = -grad
            else:
                forces = torch.zeros_like(R)
        else:  # if there are no distances, the forces are zero
            forces = torch.zeros_like(R)

        return energy, forces

    @torch.jit.export
    def energy_and_forces_and_atomic_properties(
            self, Z, R, idx_i, idx_j, Q_tot=None, batch_seg=None, offsets=None,
            sr_idx_i=None, sr_idx_j=None, sr_offsets=None, create_graph=True):

        ''' Calculates the total energy and forces (including electrostatic
            interactions)'''
        Ea, Qa, Dij, nhloss = self.atomic_properties(Z, R, idx_i, idx_j, offsets,
                                                     sr_idx_i, sr_idx_j, sr_offsets)

        energy = self.energy_from_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, Q_tot, batch_seg)

        if idx_i.numel() > 0:  # autograd will fail if there are no distances
            grad = torch.autograd.grad(
                [torch.sum(energy)], [R], create_graph=create_graph)[0]
            if grad is not None:  # necessary for torch.jit compatibility
                forces = -grad
            else:
                forces = torch.zeros_like(R)
        else:  # if there are no distances, the forces are zero
            forces = torch.zeros_like(R)

        return energy, forces, Ea, Qa, nhloss

    @torch.jit.export
    def energy_and_forces_and_charges(
            self, Z, R, idx_i, idx_j, Q_tot=None, batch_seg=None, offsets=None,
            sr_idx_i=None, sr_idx_j=None, sr_offsets=None, create_graph=True):
        ''' Calculates the total energy and forces (including electrostatic
            interactions)'''
        Ea, Qa, Dij, nhloss = self.atomic_properties(
            Z, R, idx_i, idx_j, offsets,
            sr_idx_i, sr_idx_j, sr_offsets)

        energy = self.energy_from_atomic_properties(
            Ea, Qa, Dij, Z, idx_i, idx_j, Q_tot, batch_seg)

        reduced_energy = torch.sum(energy)

        if idx_i.numel() > 0:  # autograd will fail if there are no distances
            grad = torch.autograd.grad(
                [reduced_energy], [R], create_graph=create_graph)[0]
            if grad is not None:  # necessary for torch.jit compatibility
                forces = -grad
            else:
                forces = torch.zeros_like(R)
        else:  # if there are no distances, the forces are zero
            forces = torch.zeros_like(R)

        return energy, forces, Qa

    def scaled_charges(self, Z, Qa, Q_tot=None, batch_seg=None):
        ''' Returns scaled charges such that the sum of the partial atomic
            charges equals Q_tot (defaults to 0) '''

        if batch_seg is None:
            batch_seg = torch.zeros_like(Z)

        # Number of atoms per batch (needed for charge scaling)
        Na_helper = torch.ones_like(batch_seg, dtype=self.dtype)
        Na_per_batch = segment_sum(Na_helper,batch_seg.type(torch.int64),device=self.device)


        if Q_tot is None:  # Assume desired total charge zero if not given
            Q_tot = torch.zeros_like(Na_per_batch, dtype=self.dtype)

        # Return scaled charges (such that they have the desired total charge)
        Q_correct = Q_tot - segment_sum(Qa,batch_seg.type(torch.int64),device=self.device)
        Q_scaled = Qa + torch.gather((Q_correct / Na_per_batch), 0, batch_seg.type(torch.int64))

        return Q_scaled

    def _switch(self, Dij):
        ''' Switch function for electrostatic interaction (switches between
            shielded and unshielded electrostatic interaction) '''

        cut = self.sr_cut / 2
        x = Dij / cut
        x3 = x * x * x
        x4 = x3 * x
        x5 = x4 * x

        return torch.where(Dij < cut, 6 * x5 - 15 * x4 + 10 * x3, torch.ones_like(Dij))

    def electrostatic_energy_per_atom(self, Dij, Qa, idx_i, idx_j):
        ''' Calculates the electrostatic energy per atom for very small
            distances, the 1/r law is shielded to avoid singularities '''

        # Gather charges
        Qi = torch.gather(Qa, 0, idx_i.type(torch.int64))
        Qj = torch.gather(Qa, 0, idx_j.type(torch.int64))

        # Calculate variants of Dij which we need to calculate
        # the various shielded/non-shielded potentials
        DijS = torch.sqrt(Dij * Dij + 1.0)  # shielded distance

        # Calculate value of switching function
        switch = self._switch(Dij)  # normal switch
        cswitch = 1.0 - switch  # complementary switch

        # Calculate shielded/non-shielded potentials
        if self.lr_cut is None:  # no non-bonded cutoff

            Eele_ordinary = 1.0 / Dij  # ordinary electrostatic energy
            Eele_shielded = 1.0 / DijS  # shielded electrostatic energy

            # Combine shielded and ordinary interactions and apply prefactors
            Eele = self.kehalf * Qi * Qj * (
                    cswitch * Eele_shielded + switch * Eele_ordinary)

        else:  # with non-bonded cutoff

            cut = self.lr_cut
            cut2 = self.lr_cut * self.lr_cut

            Eele_ordinary = 1.0 / Dij + Dij / cut2 - 2.0 / cut
            Eele_shielded = 1.0 / DijS + DijS / cut2 - 2.0 / cut

            # Combine shielded and ordinary interactions and apply prefactors
            Eele = self.kehalf * Qi * Qj * (
                    cswitch * Eele_shielded + switch * Eele_ordinary)
            Eele = torch.where(Dij <= cut, Eele, torch.zeros_like(Eele))

        Eele_f = segment_sum(Eele,idx_i,device=self.device)

        return Eele_f


