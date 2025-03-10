###############################################################################
###############################################################################
###############################################################################

import numpy as np
import pickle
from sklearn.neighbors import KNeighborsRegressor
from metric_learn import MLKR
from ase.atoms import Atoms
from typing import List, Tuple, Union, Dict
import os
from collections import defaultdict
import time


def _generate_fchl19(strs: List[Atoms], krr_cutoff: float = 8) -> np.ndarray:
    from qmllib.representations import generate_fchl19 as generate_representation

    n = len(strs)
    representation = generate_representation(
        strs[0].get_atomic_numbers(),
        strs[1].get_positions(),
        rcut=krr_cutoff,
        acut=krr_cutoff,
    )
    X = np.zeros((n, representation.shape[1]))
    X[0, :] = np.sum(representation, axis=0)
    for i in range(1, n):
        X[i, :] = generate_representation(
            strs.iloc[i].get_atomic_numbers(),
            strs.iloc[i].get_positions(),
            rcut=krr_cutoff,
            acut=krr_cutoff,
        ).sum(axis=0)
    if np.isnan(X).any():
        raise ValueError("NaNs in FCHL representation!")
    return X


def _generate_mbdf(strs: List[Atoms], cutoff_r: float) -> np.ndarray:
    from MBDF import generate_mbdf as generate_representation

    X = generate_representation(
        np.array([i.get_atomic_numbers() for i in strs]),
        np.array([i.get_positions() for i in strs]),
        cutoff_r=cutoff_r,
        normalized=False,
        local=False,
    )
    return X


def _generate_bob(
    strs: List[Atoms], max_atoms: int, asize: Dict[str, Union[np.int64, int]]
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


def _generate_coulomb(strs: List[Atoms], max_atoms: int):
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


def calculate_representation(
    Qrepresentation, strs, krr_cutoff=8, max_atoms=None, asize=None
):
    if Qrepresentation == "fchl":
        return _generate_fchl19(strs, krr_cutoff)
    elif Qrepresentation == "mbdf":
        return _generate_mbdf(strs, krr_cutoff)
    elif Qrepresentation == "bob":
        return _generate_bob(strs, max_atoms, asize)
    elif Qrepresentation == "coulomb":
        return _generate_coulomb(strs, max_atoms)
    else:
        raise NotImplementedError(
            f"Representation 'f{Qrepresentation}' not supported with the k-NN model!"
        )


def training(
    Qrepresentation: str,
    strs: List[Atoms],
    Y_train: np.ndarray,
    krr_cutoff: float,
    varsoutfile: Union[str, os.PathLike],
    max_atoms: int = None,
    asize: Dict[str, Union[np.int64, int]] = None,
    no_metric=False,
):

    ### REPRESENTATION CALCULATION ###
    repr_wall_start = time.perf_counter()
    repr_cpu_start = time.process_time()
    print(
        f"JKML(Q-kNN): Calculating {Qrepresentation.upper()} representation.",
        flush=True,
    )
    X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
    X_train = calculate_representation(
        Qrepresentation, strs, krr_cutoff, max_atoms, asize
    )
    repr_train_wall = time.perf_counter() - repr_wall_start
    repr_train_cpu = time.process_time() - repr_cpu_start

    # some info about the full representation
    print(
        "JKML(Q-kNN): Shape of the training representation: " + str(X_train.shape),
        flush=True,
    )
    # save train input files for off-site debugging
    with open(varsoutfile, "wb") as f:
        pickle.dump([X_train, Y_train], f)
        print(f"Saved pretrain vars to {str(f)}.", flush=True)
    train_wall_start = time.perf_counter()
    train_cpu_start = time.process_time()
    if not no_metric:
        print("JKML(Q-kNN): Training MLKR metric.", flush=True)
        mlkr = MLKR()
        mlkr.fit(X_train, Y_train)
        A = mlkr.get_mahalanobis_matrix()
        print("JKML(Q-kNN): Training k-NN regressor with MLKR metric.")
        knn = KNeighborsRegressor(
            metric=mlkr.get_metric(), n_jobs=-1, algorithm="ball_tree"
        )
    else:
        knn = KNeighborsRegressor(n_jobs=-1, algorithm="auto")
    knn.fit(X_train, Y_train)
    train_wall = time.perf_counter() - train_wall_start
    train_cpu = time.process_time() - train_cpu_start
    n_train, d_train = X_train.shape
    print("JKML(Q-kNN): Training completed.", flush=True)
    knn_params = knn.get_params()
    knn_params["metric"] = "MLKR_placeholder"
    with open(varsoutfile, "wb") as f:
        if not no_metric:
            pickle.dump([X_train, Y_train, X_atoms, A, mlkr, knn_params], f)
        else:
            pickle.dump([X_train, Y_train, X_atoms, knn_params], f)
    return {
        key: value
        for key, value in locals().items()
        if key
        in [
            "X_train",
            "X_atoms",
            "A",
            "knn",
            "repr_train_wall",
            "repr_train_cpu",
            "train_wall",
            "train_cpu",
            "n_train",
            "d_train",
        ]
    }


###############################################################################
###############################################################################
###############################################################################


def evaluate(Qrepresentation, krr_cutoff, X_train, strs, knn_model):

    import numpy as np

    ### REPRESENTATION CALCULATION ###
    X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
    repr_wall_start = time.perf_counter()
    repr_cpu_start = time.process_time()
    X_test = calculate_representation(Qrepresentation, strs, krr_cutoff)
    repr_test_wall = time.perf_counter() - repr_wall_start
    repr_test_cpu = time.process_time() - repr_cpu_start

    # some info about the full representation
    print(
        "JKML(QML): Shape of the testing representation: " + str(X_test.shape),
        flush=True,
    )

    ### CORRECTING THE FCHL MATRIX SIZES
    # IF YOU ARE EXTENDING THIS WILL MAKE THE MATRIXES OF THE SAME SIZE
    if Qrepresentation == "fchl":
        if X_train.shape[1] != X_test.shape[1]:
            if X_train.shape[1] > X_test.shape[1]:
                small = X_test
                large = X_train
            else:
                small = X_train
                large = X_test
            newmatrix = np.zeros([small.shape[0], large.shape[1], 5, large.shape[3]])
            newmatrix[:, :, 0, :] = 1e100
            newmatrix[
                0 : small.shape[0], 0 : small.shape[1], 0:5, 0 : small.shape[3]
            ] = small
            if X_train.shape[1] > X_test.shape[1]:
                X_test = newmatrix
            else:
                X_train = newmatrix

    ### THE EVALUATION
    test_wall_start = time.perf_counter()
    test_cpu_start = time.process_time()
    Y_predicted = knn_model.predict(X_test)
    test_wall = time.perf_counter() - test_wall_start
    test_cpu = time.process_time() - test_cpu_start
    Y_predicted = Y_predicted[None, :]
    d_test = X_test.shape[1]
    return Y_predicted, repr_test_wall, repr_test_cpu, test_wall, test_cpu, d_test
