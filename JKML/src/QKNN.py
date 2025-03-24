###############################################################################
###############################################################################
###############################################################################

import numpy as np
import pandas as pd
import pickle
from sklearn.neighbors import KNeighborsRegressor
from metric_learn import MLKR
from ase.atoms import Atoms
from typing import List, Tuple, Union, Dict
import os
from collections import defaultdict
import time
import warnings

# ignore sklearn futurewarning
warnings.filterwarnings(
    "ignore",
    "'force_all_finite' was renamed to 'ensure_all_finite' in 1.6 and will be removed in 1.8.",
)


def _generate_fchl19(strs: List[Atoms], rcut=8.0, acut=8.0, **kwargs) -> np.ndarray:
    from qmllib.representations import generate_fchl19 as generate_representation

    n = len(strs)
    representation = generate_representation(
        strs[0].get_atomic_numbers(),
        strs[0].get_positions(),
        rcut=rcut,
        acut=acut,
    )
    X = np.zeros((n, representation.shape[1]))
    X[0, :] = np.sum(representation, axis=0)
    for i in range(1, n):
        X[i, :] = generate_representation(
            strs[i].get_atomic_numbers(),
            strs[i].get_positions(),
            rcut=rcut,
            acut=acut,
        ).sum(axis=0)
    if np.isnan(X).any():
        raise ValueError("NaNs in FCHL representation!")
    return X


def _generate_mbdf(strs: List[Atoms], cutoff: float = 8.0, **kwargs) -> np.ndarray:
    from MBDF import generate_mbdf as generate_representation

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
    )
    return X


def _generate_bob(
    strs: List[Atoms],
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


def _generate_coulomb(strs: List[Atoms], max_atoms: int = None, **kwargs):
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


def _generate_mbtr(
    strs: List[Atoms],
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


def calculate_representation(Qrepresentation, strs, **repr_kwargs):
    if Qrepresentation == "fchl":
        return _generate_fchl19(strs, **repr_kwargs)
    elif Qrepresentation == "mbdf":
        return _generate_mbdf(strs, **repr_kwargs)
    elif Qrepresentation == "bob":
        return _generate_bob(strs, **repr_kwargs)
    elif Qrepresentation == "coulomb":
        return _generate_coulomb(strs, **repr_kwargs)
    elif Qrepresentation == "mbtr":
        return _generate_mbtr(strs, **repr_kwargs)
    else:
        raise NotImplementedError(
            f"Representation 'f{Qrepresentation}' not supported with the k-NN model!"
        )


def load_hyperparams(hyper_cache: str):
    if hyper_cache is not None:
        with open(hyper_cache, "rb") as f:
            hyperparams = pickle.load(f)
        print(f"JKML(Q-kNN): Loaded hyperparameters from {hyper_cache}:", flush=True)
        print(hyperparams, flush=True)
    else:
        # use defaults
        hyperparams = {
            "knn": {"n_neighbors": 5, "weights": "uniform"},
            "representation": {"cutoff": 8.0},
        }
        print(f"JKML(Q-kNN): Using default hyperparams {hyperparams}", flush=True)
    return hyperparams


def training(
    Qrepresentation: str,
    strs: pd.Series,
    Y_train: np.ndarray,
    varsoutfile: Union[str, os.PathLike],
    no_metric=False,
    hyper_cache=None,
):

    hyperparams = load_hyperparams(hyper_cache)
    # Convert structures to list to avoid indexing issues
    strs = strs.values
    ### REPRESENTATION CALCULATION ###
    repr_wall_start = time.perf_counter()
    repr_cpu_start = time.process_time()
    print(
        f"JKML(Q-kNN): Calculating {Qrepresentation.upper()} representation.",
        flush=True,
    )
    X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
    X_train = calculate_representation(
        Qrepresentation, strs, **hyperparams["representation"]
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
        # Limit the number of MLKR components for faster training
        mlkr = MLKR(n_components=50)
        mlkr.fit(X_train, Y_train)
        A = mlkr.get_mahalanobis_matrix()
        print("JKML(Q-kNN): Training k-NN regressor with MLKR metric.")
        knn = KNeighborsRegressor(
            metric=mlkr.get_metric(),
            n_jobs=-1,
            algorithm="ball_tree",
            **hyperparams["knn"],
        )
    else:
        knn = KNeighborsRegressor(n_jobs=-1, algorithm="auto", **hyperparams["knn"])
    knn.fit(X_train, Y_train)
    train_wall = time.perf_counter() - train_wall_start
    train_cpu = time.process_time() - train_cpu_start
    n_train, d_train = X_train.shape
    train_metadata = {
        "repr_train_wall": repr_train_wall,
        "repr_train_cpu": repr_train_cpu,
        "train_wall": train_wall,
        "train_cpu": train_cpu,
        "n_train": n_train,
        "d_train": d_train,
    }
    print("JKML(Q-kNN): Training completed.", flush=True)
    knn_params = knn.get_params()
    if not no_metric:
        knn_params["metric"] = "MLKR_placeholder"
    with open(varsoutfile, "wb") as f:
        print(f"JKML(Q-kNN): Saving training data to {varsoutfile}")
        if not no_metric:
            pickle.dump(
                [X_train, Y_train, X_atoms, A, mlkr, knn_params, train_metadata], f
            )
        else:
            pickle.dump([X_train, Y_train, X_atoms, knn_params, train_metadata], f)
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


def evaluate(Qrepresentation, X_train, strs, knn_model, hyper_cache=None):

    import numpy as np

    hyperparams = load_hyperparams(hyper_cache)
    ### REPRESENTATION CALCULATION ###
    X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
    repr_wall_start = time.perf_counter()
    repr_cpu_start = time.process_time()
    X_test = calculate_representation(
        Qrepresentation, strs, **hyperparams["representation"]
    )
    repr_test_wall = time.perf_counter() - repr_wall_start
    repr_test_cpu = time.process_time() - repr_cpu_start

    # some info about the full representation
    print(
        "JKML(QML): Shape of the testing representation: " + str(X_test.shape),
        flush=True,
    )

    ### CORRECTING THE FCHL MATRIX SIZES
    # IF YOU ARE EXTENDING THIS WILL MAKE THE MATRIXES OF THE SAME SIZE
    # THIS IS COMMENTED OUT AS IT DOES NOT APPLY TO FCHL'19 (global)
    #    if Qrepresentation == "fchl":
    #        if X_train.shape[1] != X_test.shape[1]:
    #            if X_train.shape[1] > X_test.shape[1]:
    #                small = X_test
    #                large = X_train
    #            else:
    #                small = X_train
    #                large = X_test
    #            newmatrix = np.zeros([small.shape[0], large.shape[1], 5, large.shape[3]])
    #            newmatrix[:, :, 0, :] = 1e100
    #            newmatrix[
    #                0 : small.shape[0], 0 : small.shape[1], 0:5, 0 : small.shape[3]
    #            ] = small
    #            if X_train.shape[1] > X_test.shape[1]:
    #                X_test = newmatrix
    #            else:
    #                X_train = newmatrix
    #
    ### THE EVALUATION
    test_wall_start = time.perf_counter()
    test_cpu_start = time.process_time()
    Y_predicted = knn_model.predict(X_test)
    test_wall = time.perf_counter() - test_wall_start
    test_cpu = time.process_time() - test_cpu_start
    Y_predicted = Y_predicted[None, :]
    d_test = X_test.shape[1]
    return Y_predicted, repr_test_wall, repr_test_cpu, test_wall, test_cpu, d_test


def hyperopt(
    Qrepresentation,
    strs,
    Y_train,
    hyperparamfile,
    no_metric=False,
    cv_folds=5,
    verbose=True,
    optimise_representation=False,
):

    import skopt
    from skopt.space import Real, Integer, Categorical
    from functools import lru_cache
    from sklearn.model_selection import cross_val_score

    # hard-coded search spaces (for now)
    space = []
    if optimise_representation:
        if Qrepresentation == "fchl":
            space += [
                Real(name="rcut", low=1.0, high=20.0, prior="uniform"),
                Real(name="acut", low=1.0, high=20.0, prior="uniform"),
            ]
        elif Qrepresentation == "mbdf":
            space += [
                Real(name="cutoff", low=1.0, high=20.0, prior="uniform"),
            ]
        elif Qrepresentation == "mbtr":
            # these values are based on stuke2021efficient
            raise NotImplementedError("MBTR is still under construction!")
        else:
            raise NotImplementedError(
                "Hyperparameter tuning not yet implement for representation {Qrepresentation}!"
            )

        print(
            f"JKML(Q-kNN): Begin hyperparameter optimisation with {Qrepresentation.upper()} representation.",
            flush=True,
        )
    else:
        print(
            f"JKML(Q-kNN): Begin hyperparameter optimisation with {Qrepresentation.upper()} representation (only k-NN).",
            flush=True,
        )
        global X
        X = calculate_representation(Qrepresentation, strs)

    # add k-nn specific hyperparameters
    space.append(skopt.space.Integer(3, 15, name="n_neighbors"))
    space.append(skopt.space.Categorical(["uniform", "distance"], name="weights"))

    knn_param_names = ["n_neighbors", "weights"]

    # need to have two different optimisation functions to use precalculate X
    # in the case where we don't optimise the representation hyperparams
    if optimise_representation:

        @skopt.utils.use_named_args(space)
        @lru_cache
        def objective(n_neighbors, weights, **repr_params):
            X = calculate_representation(Qrepresentation, strs, **repr_params)
            if not no_metric:
                mlkr = MLKR(n_components=50)
                mlkr.fit(X, Y_train)
                knn = KNeighborsRegressor(
                    metric=mlkr.get_metric(),
                    n_neighbors=n_neighbors,
                    weights=weights,
                    algorithm="ball_tree",
                )
            else:
                knn = KNeighborsRegressor(
                    n_jobs=-1,
                    n_neighbors=n_neighbors,
                    weights=weights,
                    algorithm="auto",
                )
            return -np.mean(
                cross_val_score(
                    knn,
                    X,
                    Y_train,
                    cv=cv_folds,
                    n_jobs=-1,
                    scoring="neg_mean_absolute_error",
                )
            )

    else:

        @skopt.utils.use_named_args(space)
        @lru_cache
        def objective(n_neighbors, weights, **repr_params):
            if not no_metric:
                mlkr = MLKR(n_components=50)
                mlkr.fit(X, Y_train)
                knn = KNeighborsRegressor(
                    metric=mlkr.get_metric(),
                    n_neighbors=n_neighbors,
                    weights=weights,
                    algorithm="ball_tree",
                )
            else:
                knn = KNeighborsRegressor(
                    n_jobs=-1,
                    n_neighbors=n_neighbors,
                    weights=weights,
                    algorithm="auto",
                )
            return -np.mean(
                cross_val_score(
                    knn,
                    X,
                    Y_train,
                    cv=cv_folds,
                    n_jobs=-1,
                    scoring="neg_mean_absolute_error",
                )
            )

    start_time = time.perf_counter()
    res = skopt.gp_minimize(
        objective,
        space,
        # can afford more calls if not optimising representation
        n_calls=10 if optimise_representation else 50,
        random_state=42,
        verbose=verbose,
        n_jobs=-1,
    )
    elapsed = time.perf_counter() - start_time
    print(f"JKML: Hyperparameter tuning done, took {elapsed:.2f} s.", flush=True)
    params = {"knn": {}, "representation": {}}
    for s, v in zip(space, res.x):
        if s.name in knn_param_names:
            params["knn"][s.name] = v
        else:
            params["representation"][s.name] = v

    with open(hyperparamfile, "wb") as f:
        pickle.dump(params, f)
        print(f"JKML(Q-kNN): Saved hyperparams to {hyperparamfile}")

    return params
