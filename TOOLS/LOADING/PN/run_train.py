# Standard importations
import torch
import torch.nn as nn
from tensorboardX import SummaryWriter
import os
import sys
import numpy as np
import argparse
import logging
import string
import random
from torch_ema import ExponentialMovingAverage
# For Time measurement
from datetime import datetime
from time import time
# Neural Network importations
from Neural_Net_evid import PhysNet
from layers.utils import segment_sum
from layers.activation_fn import *
from Neural_Net_evid import gather_nd
from DataContainer import DataContainer
#Other importations
import functools
# Configure logging environment
logging.basicConfig(level=logging.DEBUG)

# ------------------------------------------------------------------------------
# Command line arguments
# ------------------------------------------------------------------------------

# Initiate parser
parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

# Add arguments
parser.add_argument("--restart", type=str, default='No',
                    help="Restart training from a specific folder")
parser.add_argument("--checkpoint-file", type=str, default=None,
                    help="File to be loaded if model is restarted")
parser.add_argument("--num_features", default=128, type=int,
                    help="Dimensionality of feature vectors")
parser.add_argument("--num_basis", default=64, type=int,
                    help="Number of radial basis functions")
parser.add_argument("--num_blocks", default=5, type=int,
                    help="Number of interaction blocks")
parser.add_argument("--num_residual_atomic", default=2, type=int,
                    help="Number of residual layers for atomic refinements")
parser.add_argument("--num_residual_interaction", default=3, type=int,
                    help="Number of residual layers for the message phase")
parser.add_argument("--num_residual_output", default=1, type=int,
                    help="Number of residual layers for the output blocks")
parser.add_argument("--cutoff", default=10.0, type=float,
                    help="Cutoff distance for short range interactions")
parser.add_argument("--use_electrostatic", default=1, type=int,
                    help="Use electrostatics in energy prediction (0/1)")
parser.add_argument("--use_dispersion", default=1, type=int,
                    help="Use dispersion in energy prediction (0/1)")
parser.add_argument("--grimme_s6", default=None, type=float,
                    help="Grimme s6 dispersion coefficient")
parser.add_argument("--grimme_s8", default=None, type=float,
                    help="Grimme s8 dispersion coefficient")
parser.add_argument("--grimme_a1", default=None, type=float,
                    help="Grimme a1 dispersion coefficient")
parser.add_argument("--grimme_a2", default=None, type=float,
                    help="Grimme a2 dispersion coefficient")
parser.add_argument("--dataset", type=str,
                    help="File path to dataset")
# This number is configured for the size of the QM9 dataset
parser.add_argument("--num_train", default=103130, type=int,
                    help="Number of training samples")
# This number is configured for the size of the QM9 dataset
parser.add_argument("--num_valid", default=12891, type=int,
                    help="Number of validation samples")
parser.add_argument("--batch_size", default=100, type=int,
                    help="Batch size used per training step")
parser.add_argument("--valid_batch_size", default=20, type=int,
                    help="Batch size used for going through validation_set")
parser.add_argument("--seed", default=np.random.randint(1000000), type=int,
                    help="Seed for splitting dataset into " + \
                         "training/validation/test")
parser.add_argument("--max_steps", default=10000, type=int,
                    help="Maximum number of training steps")
parser.add_argument("--learning_rate", default=0.001, type=float,
                    help="Learning rate used by the optimizer")
parser.add_argument("--decay_steps", default=1000, type=int,
                    help="Decay the learning rate every N steps by decay_rate")
parser.add_argument("--decay_rate", default=0.1, type=float,
                    help="Factor with which the learning rate gets " + \
                         "multiplied by every decay_steps steps")
parser.add_argument("--max_norm", default=1000.0, type=float,
                    help="Max norm for gradient clipping")
parser.add_argument("--ema_decay", default=0.999, type=float,
                    help="Exponential moving average decay used by the " + \
                         "trainer")
parser.add_argument("--rate", default=0.0, type=float,
                    help="Rate probability for dropout regularization of " + \
                         "rbf layer")
parser.add_argument("--l2lambda", default=0.0, type=float,
                    help="Lambda multiplier for l2 loss (regularization)")
#Note: This parameter is setup to 0.2 as it was the best value on the paper of Amini...
parser.add_argument("--lambda_conf", default=0.2, type=float,
                    help="Lambda value of the confidence of the prediction")
parser.add_argument('--summary_interval', default=5, type=int,
                    help="Write a summary every N steps")
parser.add_argument('--validation_interval', default=5, type=int,
                    help="Check performance on validation set every N steps")
parser.add_argument('--show_progress', default=True, type=bool,
                    help="Show progress of the epoch")
parser.add_argument('--save_interval', default=5, type=int,
                    help="Save progress every N steps")
parser.add_argument('--record_run_metadata', default=0, type=int,
                    help="Records metadata like memory consumption etc.")
parser.add_argument('--device',default='cuda',type=str,
                    help='Selects the device that will be used for training')

# ------------------------------------------------------------------------------
# Read Parameters and define output files
# ------------------------------------------------------------------------------


# Generate an (almost) unique id for the training session
def id_generator(size=8,
                 chars=(string.ascii_uppercase
                        + string.ascii_lowercase
                        + string.digits)):
    return ''.join(random.SystemRandom().choice(chars) for _ in range(size))


# Read config file if no arguments are given
config_file = 'config.txt'
if len(sys.argv) == 1:
    if os.path.isfile(config_file):
        args = parser.parse_args(["@" + config_file])
    else:
        args = parser.parse_args(["--help"])
else:
    args = parser.parse_args()

# Create output directory for training session and
# load config file arguments if restart
if args.restart == 'No':
    directory = (
            datetime.utcnow().strftime("%Y%m%d%H%M%S")
            + "_" + id_generator() + "_F" + str(args.num_features)
            + "K" + str(args.num_basis) + "b" + str(args.num_blocks)
            + "a" + str(args.num_residual_atomic)
            + "i" + str(args.num_residual_interaction)
            + "o" + str(args.num_residual_output) + "cut" + str(args.cutoff)
            + "e" + str(args.use_electrostatic) + "d" + str(args.use_dispersion)
            + "rate" + str(args.rate))
    checkpoint_file = args.checkpoint_file
else:
    directory = args.restart
    args = parser.parse_args(["@" + os.path.join(args.restart, config_file)])
    checkpoint_file = os.path.join(args.restart, args.checkpoint_file)

# Create sub directories
logging.info("Creating directories...")

if not os.path.exists(directory):
    os.makedirs(directory)

best_dir = os.path.join(directory, 'best')
if not os.path.exists(best_dir):
    os.makedirs(best_dir)

log_dir = os.path.join(directory, 'logs')
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

ckpt_dir = os.path.join(directory, 'ckpt')
if not os.path.exists(ckpt_dir):
    os.makedirs(ckpt_dir)

# Define output files
best_loss_file = os.path.join(best_dir, 'best_loss.npz')


# Write config file of current training session
logging.info("Writing args to file...")

with open(os.path.join(directory, config_file), 'w') as f:
    for arg in vars(args):
        f.write('--' + arg + '=' + str(getattr(args, arg)) + "\n")
args.validation_interval = 1
args.summary_interval = 1

logging.info("device: {}".format(args.device))

# ------------------------------------------------------------------------------
# Define utility functions
# ------------------------------------------------------------------------------
def save_checkpoint(model, epoch, optimizer, name_of_ckpt=None,best=False):
    state = {'model_state_dict': model.state_dict(),
             'epoch': epoch, 'optimizer_state_dict': optimizer.state_dict()}
    if best:
        path = os.path.join(best_dir, 'best_model.pt')
    else:
        name = 'model' + str(name_of_ckpt) + '.pt'
        path = os.path.join(ckpt_dir, name)
    torch.save(state, path)


def load_checkpoint(path):
    if path is not None:
        checkpoint = torch.load(path)
        return checkpoint
    else:
        return None


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='#', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(
        100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print("\r{0} |{1}| {2}% {3}".format(
        prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def reset_averages(type,device='cpu'):
    ''' Reset counter and average values '''
    null_float = torch.tensor(0.0, dtype=torch.float32,device=device)
    if type == "train":
        return null_float, null_float, null_float, null_float, null_float, null_float
    elif type == "valid":
        return null_float, null_float, null_float, null_float

def l2_regularizer(model,l2_lambda=args.l2lambda):
    l2_norm = sum(p.pow(2.0).sum() for p in model.parameters())/sum(1 for p in model.parameters())
    return l2_lambda*l2_norm

#====================================
# Some functions
#====================================
def compute_pnorm(model):
    """Computes the norm of the parameters of a model."""
    return np.sqrt(sum([p.norm().item() ** 2 for p in model.parameters()]))


def compute_gnorm(model):
    """Computes the norm of the gradients of a model."""
    return np.sqrt(sum([p.grad.norm().item() ** 2 for p in model.parameters() if p.grad is not None]))
# ------------------------------------------------------------------------------
# Load data and initiate PhysNet model
# ------------------------------------------------------------------------------

# Load dataset
logging.info("Loading dataset...")

data = DataContainer(
    args.dataset, args.num_train, args.num_valid,
    args.batch_size, args.valid_batch_size, seed=args.seed)

# Initiate PhysNet model
logging.info("Creating PhysNet model...")

# Initiate summary writer
summary_writer = SummaryWriter(log_dir)
model = PhysNet(
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
    Eshift=data.EperA_m_n,
    Escale=data.EperA_s_n,
    activation_fn=shifted_softplus,
    device=args.device,
    writer=summary_writer)


if os.path.isfile(best_loss_file):
    loss_file = np.load(best_loss_file)
    best_loss = loss_file["loss"].item()
    best_emae = loss_file["emae"].item()
    best_fmae = loss_file["fmae"].item()
else:
    best_loss = np.inf
    best_emae = np.inf
    best_fmae = np.inf
    best_epoch = 0.
    np.savez(
        best_loss_file, loss=best_loss, emae=best_emae, fmae=best_fmae,
        epoch=best_epoch)

#------------------------------------
# Loss function
#------------------------------------

def evidential_loss_new(mu, v, alpha, beta, targets, lam=args.lambda_conf, epsilon=1e-4):
    """
    Use Deep Evidential Regression negative log likelihood loss + evidential
        regularizer
    We will use the new version on the paper..
    :mu: pred mean parameter for NIG
    :v: pred lam parameter for NIG
    :alpha: predicted parameter for NIG
    :beta: Predicted parmaeter for NIG
    :targets: Outputs to predict
    :return: Loss
    """
    # Calculate NLL loss
    twoBlambda = 2*beta*(1+v)
    nll = 0.5*torch.log(np.pi/v) \
        - alpha*torch.log(twoBlambda) \
        + (alpha+0.5) * torch.log(v*(targets-mu)**2 + twoBlambda) \
        + torch.lgamma(alpha) \
        - torch.lgamma(alpha+0.5)

    L_NLL = nll #torch.mean(nll, dim=-1)

    # Calculate regularizer based on absolute error of prediction
    error = torch.abs((targets - mu))
    reg = error * (2 * v + alpha)
    L_REG = reg #torch.mean(reg, dim=-1)

    # Loss = L_NLL + L_REG
    # TODO If we want to optimize the dual- of the objective use the line below:
    loss = L_NLL + lam * (L_REG - epsilon) + l2_regularizer(model)

    return loss

def gauss_loss(mu,sigma,targets):
    """
    This defines a simple loss function for learning the log-likelihood of a gaussian distribution.
    In the future, we should use a regularizer.
    """
    loss = 0.5*np.log(2*np.pi) + 0.5*torch.log(sigma**2) + ((targets-mu)**2)/(2*sigma**2) + l2_regularizer(model)

    return loss

def evidential_loss(mu, v, alpha, beta, targets):
    """
    Use Deep Evidential Regression Sum of Squared Error loss
    :mu: Pred mean parameter for NIG
    :v: Pred lambda parameter for NIG
    :alpha: predicted parameter for NIG
    :beta: Predicted parmaeter for NIG
    :targets: Outputs to predict
    :return: Loss
    """

    # Calculate SOS
    # Calculate gamma terms in front
    def Gamma(x):
        return torch.exp(torch.lgamma(x))

    coeff_denom = 4 * Gamma(alpha) * v * torch.sqrt(beta)
    coeff_num = Gamma(alpha - 0.5)
    coeff = coeff_num / coeff_denom

    # Calculate target dependent loss
    second_term = 2 * beta * (1 + v)
    second_term += (2 * alpha - 1) * v * torch.pow((targets - mu), 2)
    L_SOS = coeff * second_term

    # Calculate regularizer
    L_REG = torch.pow((targets - mu), 2) * (2 * alpha + v)

    loss_val = L_SOS + L_REG  + l2_regularizer(model)

    return loss_val
# ------------------------------------------------------------------------------
# Define training step
# ------------------------------------------------------------------------------
def get_indices(Nref,device='cpu'):
    # Get indices pointing to batch image
    # For some reason torch does not make repetition for float
    batch_seg = torch.arange(0, Nref.size()[0],device=device).repeat_interleave(Nref.type(torch.int64))
    # Initiate auxiliary parameter
    Nref_tot = torch.tensor(0, dtype=torch.int32).to(device)

    # Indices pointing to atom at each batch image
    idx = torch.arange(end=Nref[0], dtype=torch.int32).to(device)
    # Indices for atom pairs ij - Atom i
    Ntmp = Nref.cpu()
    idx_i = idx.repeat(int(Ntmp.numpy()[0]) - 1) + Nref_tot
    # Indices for atom pairs ij - Atom j
    idx_j = torch.roll(idx, -1, dims=0) + Nref_tot
    for Na in torch.arange(2, Nref[0]):
        Na_tmp = Na.cpu()
        idx_j = torch.concat(
            [idx_j, torch.roll(idx, int(-Na_tmp.numpy()), dims=0) + Nref_tot],
            dim=0)

    # Increment auxiliary parameter
    Nref_tot = Nref_tot + Nref[0]

    # Complete indices arrays
    for Nref_a in Nref[1:]:

        rng_a = torch.arange(end=Nref_a).to(device)
        Nref_a_tmp = Nref_a.cpu()
        idx = torch.concat([idx, rng_a], axis=0)
        idx_i = torch.concat(
            [idx_i, rng_a.repeat(int(Nref_a_tmp.numpy()) - 1) + Nref_tot],
            dim=0)
        for Na in torch.arange(1, Nref_a):
            Na_tmp = Na.cpu()
            idx_j = torch.concat(
                [idx_j, torch.roll(rng_a, int(-Na_tmp.numpy()), dims=0) + Nref_tot],
                dim=0)

        # Increment auxiliary parameter
        Nref_tot = Nref_tot + Nref_a
    #Reorder the idx in i
    idx_i = torch.sort(idx_i)[0]
    # Combine indices for batch image and respective atoms
    idx = torch.stack([batch_seg, idx], dim=1)
    return idx.type(torch.int64), idx_i.type(torch.int64), idx_j.type(torch.int64), batch_seg.type(torch.int64)

def evid_train_step(batch,num_t,loss_avg_t, fmae_avg_t, emae_avg_t,pnorm,gnorm,device,maxnorm=1000):
    model.train()
    # lr_schedule = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.8, patience=2,
    #                                                          min_lr=1e-4,verbose=True)
    batch = [i.to(device) for i in batch]
    N_t, Z_t, R_t, Eref_t, Earef_t, Fref_t, Qref_t, Qaref_t, Dref_t = batch

    # Get indices
    idx_t, idx_i_t, idx_j_t, batch_seg_t = get_indices(N_t,device=device)

    # Gather data
    Z_t = gather_nd(Z_t, idx_t)
    R_t = gather_nd(R_t, idx_t)

    if torch.count_nonzero(Earef_t) != 0:
        Earef_t = gather_nd(Earef_t, idx_t)
    if torch.count_nonzero(Fref_t) != 0:
        Fref_t = gather_nd(Fref_t, idx_t)
    if torch.count_nonzero(Qaref_t) != 0:
        Qaref_t = gather_nd(Qaref_t, idx_t)

    if 1==0:
      energy_t, lambdas_t, alpha_t, beta_t = \
              model.energy_evidential(Z_t, R_t, idx_i_t, idx_j_t, Qref_t, batch_seg=batch_seg_t)
    else: 
      #energy_t, forces_t = model.energy_and_forces(Z_t, R_t, idx_i_t, idx_j_t, Qref_t, batch_seg=batch_seg_t)
      #energy_t, forces_t, charges_t = model.energy_and_forces_and_charges(Z_t, R_t, idx_i_t, idx_j_t, Qref_t, batch_seg=batch_seg_t)
      energy_t, forces_t, Ea_t, Qa_t, nhloss_t = \
            model.energy_and_forces_and_atomic_properties(
                Z_t, R_t, idx_i_t, idx_j_t, Qref_t, batch_seg_t)
      # Get total charge
      #Qtot_t = tf.math.segment_sum(Qa_t, batch_seg_t)
      segment_sum = torch.zeros(batch_seg_t.max().item() + 1, dtype=Qa_t.dtype)
      device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
      segment_sum = segment_sum.to(device)
      segment_sum = segment_sum.scatter_add(0, batch_seg_t, Qa_t)
      #print(segment_sum)
      # Get dipole moment vector
      #QR_t = tf.stack([Qa_t*R_t[:,0], Qa_t*R_t[:,1], Qa_t*R_t[:,2]], 1)
      #D_t = tf.math.segment_sum(QR_t, batch_seg_t)
      QR_t = torch.stack([Qa_t * R_t[:, 0], Qa_t * R_t[:, 1], Qa_t * R_t[:, 2]], dim=1)
      num_segments = batch_seg_t.max().item() + 1
      D_t = torch.zeros((num_segments, QR_t.size(1)), dtype=QR_t.dtype)
      D_t = D_t.to(device)
      D_t = D_t.scatter_add(0, batch_seg_t.unsqueeze(1).expand(-1, QR_t.size(1)), QR_t)
      #print(D_t)
      #with torch.no_grad():
      #    energy_t, lambdas_t, alpha_t, beta_t = \
      #        model.energy_evidential(Z_t, R_t, idx_i_t, idx_j_t, Qref_t, batch_seg=batch_seg_t)

    # loss
    if 1==0:
      loss_t = evidential_loss_new(energy_t, lambdas_t, alpha_t, beta_t, Eref_t).sum()
    else:
      loss_t = 50 * torch.mean(torch.abs(forces_t - Fref_t))
      if not torch.isnan(Dref_t).any(): 
        loss_t = loss_t + 27 * torch.mean(torch.abs(D_t - Dref_t)) 
      loss_t = loss_t + torch.mean(torch.abs(energy_t - Eref_t)) 
      loss_t = loss_t + l2_regularizer(model).sum() 
      loss_t = loss_t + 14 * torch.mean(torch.abs(segment_sum - Qref_t)) 
      loss_t = loss_t + 0.01 * nhloss_t
      #loss_t = 0.99 * torch.mean(torch.abs(forces_t - Fref_t)) + 0.01 * torch.sum(torch.abs(charges_t - Qaref_t)) + 0.01 * torch.sum(torch.abs(energy_t - Eref_t)) + l2_regularizer(model).sum()
      #loss_t = 0.99 * torch.mean((forces_t - Fref_t) ** 2)**0.5 + 0.01 * torch.mean((charges_t - Qaref_t) ** 2) ** 0.5 + 0.01 * torch.mean((energy_t - Eref_t) ** 2) ** 0.5 + l2_regularizer(model).sum()
      #loss_t = 0.99 * torch.mean(((forces_t - Fref_t) ** 2).sum())**0.5 + 0.01 * evidential_loss_new(energy_t, lambdas_t, alpha_t, beta_t, Eref_t).sum()

    mae_energy = torch.mean(torch.abs(energy_t - Eref_t))
    loss_t.backward(retain_graph=True)

    # #Gradient clip
    nn.utils.clip_grad_norm_(model.parameters(),maxnorm)
    pnorm = pnorm + compute_pnorm(model)
    gnorm = gnorm + compute_gnorm(model)
    f = num_t /(num_t + N_t.dim())
    loss_avg_t = f * loss_avg_t + (1.0 -f)*float(loss_t)
    #ivo_output = 0.99 * torch.mean(torch.abs(forces_t - Fref_t)) + 0.01 torch.mean(torch.abs(energy_t - Eref_t))
    fmae_avg_t = f * fmae_avg_t + (1.0 - f) * float(torch.mean(torch.abs(forces_t - Fref_t)))
    emae_avg_t = f * emae_avg_t + (1.0 - f) * float(mae_energy)
    num_t = num_t + N_t.dim()
    return num_t, loss_avg_t, fmae_avg_t, emae_avg_t, pnorm,gnorm

def evid_valid_step(batch,num_v,loss_avg_v, fmae_avg_v, emae_avg_v,device):
    model.eval()
    batch = [i.to(device) for i in batch]
    N_v, Z_v, R_v, Eref_v, Earef_v, Fref_v, Qref_v, Qaref_v, Dref_v = batch
    # Get indices
    idx_v, idx_i_v, idx_j_v, batch_seg_v = get_indices(N_v,device=device)
    Z_v = gather_nd(Z_v, idx_v)
    R_v = gather_nd(R_v, idx_v)

    if torch.count_nonzero(Earef_v) != 0:
        Earef_v = gather_nd(Earef_v, idx_v)
    if torch.count_nonzero(Fref_v) != 0:
        Fref_v = gather_nd(Fref_v, idx_v)
    if torch.count_nonzero(Qaref_v) != 0:
        Qaref_v = gather_nd(Qaref_v, idx_v)

    if 1==0:
      with torch.no_grad():
          energy_v, lambdas_v, alpha_v, beta_v = \
              model.energy_evidential(Z_v, R_v, idx_i_v, idx_j_v, Qref_v, batch_seg=batch_seg_v)
    else: 
      #energy_v, forces_v = model.energy_and_forces(Z_v, R_v, idx_i_v, idx_j_v, Qref_v, batch_seg=batch_seg_v)
      #energy_v, forces_v, charges_v = model.energy_and_forces_and_charges(Z_v, R_v, idx_i_v, idx_j_v, Qref_v, batch_seg=batch_seg_v)
      energy_v, forces_v, Ea_v, Qa_v, nhloss_v = \
            model.energy_and_forces_and_atomic_properties(
                Z_v, R_v, idx_i_v, idx_j_v, Qref_v, batch_seg_v)
      # Get total charge
      #Qtot_t = tf.math.segment_sum(Qa_t, batch_seg_t)
      segment_sum = torch.zeros(batch_seg_v.max().item() + 1, dtype=Qa_v.dtype)
      device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
      segment_sum = segment_sum.to(device)
      segment_sum = segment_sum.scatter_add(0, batch_seg_v, Qa_v)
      #print(segment_sum)
      # Get dipole moment vector
      #QR_t = tf.stack([Qa_t*R_t[:,0], Qa_t*R_t[:,1], Qa_t*R_t[:,2]], 1)
      #D_t = tf.math.segment_sum(QR_t, batch_seg_t)
      QR_v = torch.stack([Qa_v * R_v[:, 0], Qa_v * R_v[:, 1], Qa_v * R_v[:, 2]], dim=1)
      num_segments = batch_seg_v.max().item() + 1
      D_v = torch.zeros((num_segments, QR_v.size(1)), dtype=QR_v.dtype)
      D_v = D_v.to(device)
      D_v = D_v.scatter_add(0, batch_seg_v.unsqueeze(1).expand(-1, QR_v.size(1)), QR_v)
      #print(D_v)
      #with torch.no_grad():
      #    energy_v, lambdas_v, alpha_v, beta_v = \
      #        model.energy_evidential(Z_v, R_v, idx_i_v, idx_j_v, Qref_v, batch_seg=batch_seg_v)

    # loss
    if 1==0:
      loss_v = evidential_loss_new(energy_v, lambdas_v, alpha_v, beta_v, Eref_v).sum()
    else:
      #loss_v = 0.99 * torch.mean((forces_v - Fref_v) ** 2)**0.5 + 0.01 * torch.mean((charges_v - Qaref_v) ** 2) ** 0.5 + 0.01 * torch.mean((energy_v - Eref_v) ** 2) ** 0.5 + l2_regularizer(model).sum()
      #loss_v = 0.99 * torch.mean(torch.abs(forces_v - Fref_v)) + 0.01 * torch.sum(torch.abs(charges_v - Qaref_v)) + 0.01 * torch.sum(torch.abs(energy_v - Eref_v)) + l2_regularizer(model).sum()
      loss_v = 50 * torch.mean(torch.abs(forces_v - Fref_v)) 
      if not torch.isnan(Dref_v).any(): 
        loss_v = loss_v + 27 * torch.mean(torch.abs(D_v - Dref_v)) 
      loss_v = loss_v + torch.mean(torch.abs(energy_v - Eref_v)) 
      loss_v = loss_v + l2_regularizer(model).sum() 
      loss_v = loss_v + 14 * torch.mean(torch.abs(segment_sum - Qref_v)) 
      loss_v = loss_v + 0.01 * nhloss_v
      #print(0.99 * torch.mean((forces_v - Fref_v) ** 2)**0.5)
      #print(0.01 * torch.mean((energy_v - Eref_v) ** 2) ** 0.5)
      #print(0.01 * torch.mean((charges_v - Qaref_v) ** 2) ** 0.5)
      print(torch.mean(torch.abs(forces_v - Fref_v)))
      print(torch.mean(torch.abs(D_v - Dref_v)))
      print(torch.mean(torch.abs(energy_v - Eref_v)))
      print(torch.mean(torch.abs(segment_sum - Qref_v)))
      print(0.01 * nhloss_v)
      print(l2_regularizer(model).sum())
      #loss_v = 0.99 * torch.mean(((forces_v - Fref_v) ** 2).sum())**0.5 + 0.01 * evidential_loss_new(energy_v, lambdas_v, alpha_v, beta_v, Eref_v).sum()     
 
    mae_energy = torch.mean(torch.abs(energy_v - Eref_v))
    f = num_v / (num_v + N_v.dim())
    loss_avg_v = f * loss_avg_v + (1.0 - f) * float(loss_v)
    fmae_avg_v = f * fmae_avg_v + (1.0 - f) * float(torch.mean(torch.abs(forces_v - Fref_v)))
    emae_avg_v = f * emae_avg_v + (1.0 - f) * float(mae_energy)
    # c = 1 /((preds[2]-1)*preds[1])
    # var = preds[3]*c
    num_v = num_v + N_v.dim()
    return num_v,loss_avg_v, fmae_avg_v, emae_avg_v

# ------------------------------------------------------------------------------
# Train PhysNet model
# ------------------------------------------------------------------------------
logging.info("starting training...")

# Define Optimizer
optimizer = torch.optim.Adam(params=model.parameters(), lr=args.learning_rate,
                             weight_decay=args.l2lambda,amsgrad=True)

lr_schedule = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=np.power(args.decay_rate,1/args.decay_steps))

# Define Exponential Moving Average
ema = ExponentialMovingAverage(model.parameters(),decay=args.ema_decay)


def model_summary(model):
  print("model_summary")
  print()
  print("Layer_name"+"\t"*7+"Number of Parameters")
  print("="*100)
  model_parameters = [layer for layer in model.parameters() if layer.requires_grad]
  layer_name = [child for child in model.children()]
  j = 0
  total_params = 0
  print("\t"*10)
  for i in layer_name:
    print()
    param = 0
    try:
      bias = (i.bias is not None)
    except:
      bias = False  
    if not bias:
      param =model_parameters[j].numel()+model_parameters[j+1].numel()
      j = j+2
    else:
      param =model_parameters[j].numel()
      j = j+1
    print(str(i)+"\t"*3+str(param))
    total_params+=param
  print("="*100)
  print(f"Total Params:{total_params}")       

model_summary(model)

# Initiate epoch and step counter
epoch = torch.tensor(1, requires_grad=False, dtype=torch.int64)
step = torch.tensor(1, requires_grad=False, dtype=torch.int64)

# Initiate checkpoints and load last checkpoint
latest_ckpt = load_checkpoint(checkpoint_file)
if latest_ckpt is not None:
    model.load_state_dict(latest_ckpt['model_state_dict'])
    optimizer.load_state_dict(latest_ckpt['optimizer_state_dict'])
    epoch = latest_ckpt['epoch']

# Create validation batches
valid_batches = data.get_valid_batches()

# Initialize counter for estimated time per epoch
time_train_estimation = np.nan
time_train = 0.0
best_loss = np.inf

# Training loop
# Terminate training when maximum number of iterations is reached
while epoch <= args.max_steps:

    # Reset error averages
    num_t, loss_avg_t, fmae_avg_t, emae_avg_t,\
    pnorm_t, gnorm_t = reset_averages("train",device=args.device)

    # Create train batches
    train_batches, N_train_batches = data.get_train_batches()

    # Start train timer
    train_start = time()

    # Iterate over batches
    for ib, batch in enumerate(train_batches):
        optimizer.zero_grad()
        # Start batch timer
        batch_start = time()

        # Show progress bar
        if args.show_progress:
            printProgressBar(
                ib, N_train_batches, prefix="Epoch {0: 5d}".format(
                    epoch.numpy()),
                suffix=("Complete - Remaining Epoch Time: "
                        + "{0: 4.1f} s     ".format(time_train_estimation)),
                length=42)

        # Training step

        num_t, loss_avg_t, fmae_avg_t, emae_avg_t,\
        pnorm_t, gnorm_t = evid_train_step(batch,num_t,loss_avg_t,
                                                 fmae_avg_t, emae_avg_t,pnorm_t,gnorm_t,args.device)
        optimizer.step()

        ema.update()
        # Stop batch timer
        batch_end = time()

        # Actualize time estimation
        if args.show_progress:
            if ib == 0:
                time_train_estimation = (
                        (batch_end - batch_start) * (N_train_batches - 1))
            else:
                time_train_estimation = (
                        0.5 * (time_train_estimation - (batch_end - batch_start))
                        + 0.5 * (batch_end - batch_start) * (N_train_batches - ib - 1))

        # Increment step number
        step = step + 1

    # Stop train timer
    train_end = time()
    time_train = train_end - train_start

    # Show final progress bar and time
    if args.show_progress:
        loss_ev_t_temp = loss_avg_t.detach().cpu()
        lat = float(loss_ev_t_temp.numpy())
        printProgressBar(
            N_train_batches, N_train_batches, prefix="Epoch {0: 5d}".format(
                epoch.numpy()),
            suffix=("Done - Epoch Time: "
                    + "{0: 4.1f} s, Average Loss: {1: 4.4f}   ".format(
                        time_train, lat)))  # length=42))

    # Save progress
    if (epoch % args.save_interval == 0):
        number_of_ckpt = int(epoch / args.save_interval)
        save_checkpoint(model=model, epoch=epoch, optimizer=optimizer, name_of_ckpt=number_of_ckpt)


    # Check performance on the validation set

    if (epoch % args.validation_interval) == 0:
        # Update training results
        results_t = {}
        loss_ev_t_temp = loss_avg_t.detach().cpu()
        results_t["loss_train"] = loss_ev_t_temp.numpy()
        results_t["time_train"] = time_train
        results_t["norm_parm"] = pnorm_t
        results_t["norm_grad"] = gnorm_t
        if data.include_E:
            emae_t_temp = emae_avg_t.detach().cpu()
            fmae_t_temp = fmae_avg_t.detach().cpu()
            results_t["energy_mae_train"] = emae_t_temp.numpy()
            results_t["force_mae_train"] = fmae_t_temp.numpy()

        # Write Results to tensorboard
        for key, value in results_t.items():
            summary_writer.add_scalar(key, value, global_step=epoch)

        # # Backup variables and assign EMA variables
        # backup_vars = [tf.identity(var) for var in model.trainable_variables]
        # for var in model.trainable_variables:
        #     var.assign(ema.average(var))

        # Reset error averages
        num_v,loss_avg_v, fmae_avg_v, emae_avg_v = reset_averages('valid',device=args.device)
        # Start valid timer
        valid_start = time()

        for ib, batch in enumerate(valid_batches):
            num_v,loss_avg_v, fmae_avg_v, emae_avg_v =\
                evid_valid_step(batch, num_v,loss_avg_v, fmae_avg_v, emae_avg_v,
                                device=args.device)

        # Stop valid timer
        valid_end = time()
        time_valid = valid_end - valid_start

        # Update validation results
        results_v = {}
        loss_avg_v_temp = loss_avg_v.detach().cpu()
        results_v["loss_valid"] = loss_avg_v_temp.numpy()
        results_v["time_valid"] = time_valid
        if data.include_E:
            emae_v_temp = emae_avg_v.detach().cpu()
            fmae_v_temp = fmae_avg_v.detach().cpu()
            results_v["energy_mae_valid"] = emae_v_temp.numpy()
            results_v["force_mae_valid"] = fmae_v_temp.numpy()

        for key, value in results_v.items():
            summary_writer.add_scalar(key, value, global_step=epoch)

        if results_v["loss_valid"] < best_loss:

            # Assign results of best validation
            best_loss = results_v["loss_valid"]
            if data.include_E:
                best_emae = results_v["energy_mae_valid"]
                best_fmae = results_v["force_mae_valid"]
            else:
                best_emae = np.inf
                best_fmae = np.inf

            best_epoch = epoch.numpy()

            # Save best results
            np.savez(
                best_loss_file, loss=best_loss,
                emae=best_emae, fmae=best_fmae,
                epoch=best_epoch)

            # Save best model variables
            save_checkpoint(model=model, epoch=epoch, optimizer=optimizer,best=True)

            # Update best results
            results_b = {}
            results_b["loss_best"] = best_loss
            if data.include_E:
                results_b["energy_mae_best"] = best_emae
                results_b["force_mae_best"] = best_fmae

             # Write the results to tensorboard
            for key, value in results_b.items():
                summary_writer.add_scalar(key, value, global_step=epoch)


        # for var, bck in zip(model.trainable_variables, backup_vars):
        #     var.assign(bck)

    #Generate summaries
    if ((epoch % args.summary_interval == 0)
            and (epoch >= args.validation_interval)):

        if data.include_E:
            print(
                "Summary Epoch: " + \
                str(epoch.numpy()) + '/' + str(args.max_steps),
                "\n  ", str(epoch.numpy()),"  Loss   train/valid: {0: 1.3e} / {1: 1.3e} ".format(
                    results_t["loss_train"],
                    results_v["loss_valid"]),
                " Best valid loss:   {0: 1.3e} ".format(
                    results_b["loss_best"]),
                "\n  ", str(epoch.numpy()),"  MAE(E) train/valid: {0: 1.3e} / {1: 1.3e} ".format(
                    results_t["energy_mae_train"],
                    results_v["energy_mae_valid"]),
                " Best valid MAE(E): {0: 1.3e} ".format(
                    results_b["energy_mae_best"]),
                "\n  ", str(epoch.numpy()),"  MAE(F) train/valid: {0: 1.3e} / {1: 1.3e} ".format(
                    results_t["force_mae_train"],
                    results_v["force_mae_valid"]),
                " Best valid MAE(F): {0: 1.3e} ".format(
                    results_b["force_mae_best"]))

    # Increment epoch number
    lr_schedule.step()
    epoch += 1

#print("Jakub's stupid test")
##self._last_energy, lambdas, alpha, beta = self.model.energy_evidential(self.Z, self.R, idx_i, idx_j, Q_tot=self.Q_tot, batch_seg=None,offsets=offsets, sr_idx_i=sr_idx_i, sr_idx_j=sr_idx_j, sr_offsets=sr_offsets)
#data = np.load('database.npz')
#E, lambdas, alpha, beta = model.energy_evidential(data['R'][0], data['R'][0], Q_tot=0, batch_seg=None)
#import pandas as pd
#
#clusters_df = pd.read_pickle("SA_SA_100_WB97X3C_FORCES.pkl")[0]
#for i in range(len(clusters_df)):
#    atoms = clusters_df["xyz"]["structure"].values[i].copy()
#    #if ("log","charge") in clusters_df.columns:
#    #  charge = clusters_df["log"]["charge"].values[i]
#    #else:
#    #  charge = 0
#
#    #calc = PhysNetCalculator(
#    #  checkpoint=varsoutfile,
#    #  atoms=atoms,
#    #  charge=charge,
#    #  config='input.inp')
#
#    #atoms.set_calculator(calc)
#    atoms.calc = model
#    Y_predicted.append(0.0367493 * atoms.get_potential_energy())
#    F_predicted.append(0.0367493 * atoms.get_forces())  # Hartree/Ang
##energy_t, lambdas_t, alpha_t, beta_t = \
##        model.energy_evidential(Z_t, R_t, idx_i_t, idx_j_t, Qref_t, batch_seg=batch_seg_t)
