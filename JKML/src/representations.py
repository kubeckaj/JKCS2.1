"""This module contains methods related to creating and manipulating chemical descriptors."""

from ase import Atoms
from typing import Iterable, Dict, Union
import numpy as np
from collections import defaultdict


def generate_global_fchl19(
    strs: Iterable[Atoms], max_atoms=None, elements=None, rcut=8.0, acut=8.0, **kwargs
) -> np.ndarray:
    X = generate_fchl19(strs, max_atoms, elements, rcut, acut, **kwargs)
    X = np.sum(X, axis=1)
    return X


def generate_fchl19(
    strs: Iterable[Atoms], max_atoms=None, elements=None, rcut=8.0, acut=8.0, **kwargs
) -> np.ndarray:
    from qmllib.representations import generate_fchl19 as generate_representation

    if elements is None:
        elements = [1, 6, 7, 8, 16]
    if max_atoms is None:
        max_atoms = max([len(s.get_atomic_numbers()) for s in strs])
    n = len(strs)
    representation = generate_representation(
        strs[0].get_atomic_numbers(),
        strs[0].get_positions(),
        elements=elements,
        rcut=rcut,
        acut=acut,
        pad=max_atoms,
    )
    X = np.zeros((n, max_atoms, representation.shape[-1]))
    X[0, :, :] = representation
    for i in range(1, n):
        X[i, :] = generate_representation(
            strs[i].get_atomic_numbers(),
            strs[i].get_positions(),
            elements=elements,
            rcut=rcut,
            acut=acut,
            pad=max_atoms,
        )
    if np.isnan(X).any():
        raise ValueError("NaNs in FCHL representation!")
    return X


def generate_global_mbdf(
    strs: Iterable[Atoms], max_atoms=None, cutoff: float = 8.0, **kwargs
) -> np.ndarray:
    from MBDF import generate_mbdf as generate_representation

    if max_atoms is None:
        max_atoms = max([len(s.get_atomic_numbers()) for s in strs])
    n = len(strs)
    ragged_atomic_numbers = np.empty(n, dtype=object)
    ragged_atomic_numbers[:] = [i.get_atomic_numbers() for i in strs]
    ragged_positions = np.empty(n, dtype=object)
    ragged_positions[:] = [i.get_positions() for i in strs]
    X = generate_representation(
        ragged_atomic_numbers,
        ragged_positions,
        cutoff_r=cutoff,
        normalized=False,
        local=False,
        pad=max_atoms,
    )
    return X


def generate_mbdf(strs: Iterable[Atoms], max_atoms=None, neighbors=None, cutoff=8.0):
    from MBDF import generate_mbdf as generate_representation

    if max_atoms is None:
        max_atoms = max([len(s.get_atomic_numbers()) for s in strs])
    if neighbors is None:
        neighbors = max_atoms
    X = generate_representation(
        np.array([i.get_atomic_numbers() for i in strs]),
        np.array([i.get_positions() for i in strs]),
        cutoff_r=cutoff,
        normalized=False,
    )
    return X


def generate_bob(
    strs: Iterable[Atoms],
    max_atoms: int = None,
    asize: Dict[str, Union[np.int64, int]] = None,
    **kwargs,
):
    from qmllib.representations import generate_bob as generate_representation

    n = len(strs)
    if max_atoms is None:
        max_atoms = max([len(x) for x in strs])

    if asize is None:
        asize = defaultdict(int)
        for struct in strs:
            elements, counts = np.unique(
                struct.get_chemical_symbols(), return_counts=True
            )
            for e, c in zip(elements, counts):
                if c > asize[e]:
                    asize[e] = c

    # get the first representation to find out the length
    repres = generate_representation(
        strs[0].get_atomic_numbers(),
        strs[0].get_positions(),
        atomtypes=None,
        size=max_atoms,
        asize=asize,
    )
    X = np.zeros((n, repres.shape[0]))
    X[0, :] = repres
    for i in range(1, n):
        X[i, :] = generate_representation(
            strs[i].get_atomic_numbers(),
            strs[i].get_positions(),
            # this argument is not used for anything, but it's mandatory :) thanks QML!
            atomtypes=None,
            size=max_atoms,
            asize=asize,
        )
    return X


def generate_coulomb(strs: Iterable[Atoms], max_atoms: int = None, **kwargs):
    from qmllib.representations import (
        generate_coulomb_matrix as generate_representation,
    )

    n = len(strs)
    if max_atoms is None:
        max_atoms = max([len(x) for x in strs])

    # get the first representation to find out the length
    X = np.zeros((n, max_atoms * (max_atoms + 1) // 2))
    for i, struct in enumerate(strs):
        X[i, :] = generate_representation(
            struct.get_atomic_numbers(),
            struct.get_positions(),
            size=max_atoms,
        )
    return X


def generate_global_mbtr(
    strs: Iterable[Atoms],
    k2_width=0.1,
    k2_weight=0.5,
    k3_width=0.1,
    k3_weight=0.5,
    **kwargs,
):
    from dscribe.descriptors import MBTR

    unique_atoms = set()
    for structure in strs:
        atom_set = set(structure.get_chemical_symbols())
        unique_atoms = unique_atoms | atom_set
    unique_atoms = list(unique_atoms)
    mbtr = MBTR(
        species=unique_atoms,
        k2={
            "geometry": {"function": "inverse_distance"},
            "grid": {"min": 0, "max": 1, "n": 200, "sigma": k2_width},
            "weighting": {"function": "exp", "scale": k2_weight, "threshold": 3e-3},
        },
        k3={
            "geometry": {"function": "cosine"},
            "grid": {"min": -1, "max": 1, "n": 200, "sigma": k3_width},
            "weighting": {"function": "exp", "scale": k3_weight, "threshold": 3e-3},
        },
        periodic=False,
        normalization="none",
    )
    X = np.zeros((len(strs), mbtr.get_number_of_features()))
    for i, struct in enumerate(strs):
        X[i, :] = mbtr.create(struct)
    return X


def generate_fchl18(strs: Iterable[Atoms], max_atoms=None, cutoff=8.0):
    from qmllib.representations import generate_fchl18 as generate_representation

    if max_atoms is None:
        max_atoms = max([len(s.get_atomic_numbers()) for s in strs])

    representations = []
    for struct in strs:
        representations.append(
            generate_representation(
                struct.get_atomic_numbers(),
                struct.get_positions(),
                max_size=max_atoms,
                neighbors=max_atoms,
                cut_distance=cutoff,
            )
        )
    return np.array(representations)


def correct_fchl18_kernel_size(X_test: np.ndarray, X_train: np.ndarray):
    if X_train.shape[1] != X_test.shape[1]:
        if X_train.shape[1] > X_test.shape[1]:
            small = X_test
            large = X_train
        else:
            small = X_train
            large = X_test
        newmatrix = np.zeros([small.shape[0], large.shape[1], 5, large.shape[3]])
        newmatrix[:, :, 0, :] = 1e100
        newmatrix[0 : small.shape[0], 0 : small.shape[1], 0:5, 0 : small.shape[3]] = (
            small
        )
        if X_train.shape[1] > X_test.shape[1]:
            X_test = newmatrix
        else:
            X_train = newmatrix
    return X_test, X_train
