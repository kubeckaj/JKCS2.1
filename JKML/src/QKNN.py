###############################################################################
###############################################################################
###############################################################################

import numpy as np
import pandas as pd
import pickle
from sklearn.neighbors import KNeighborsRegressor
from src.batched_mlkr import BatchedMLKR as MLKR
from ase.atoms import Atoms
from typing import Iterable, Tuple, Union, Dict, Literal, Callable, List, Optional, Any
import os
from collections import defaultdict
import time
import warnings
import sys
import os
from src.representations import *

sys.path.append(os.path.join(os.path.dirname(__file__), "fortran"))
from qmllib.utils.alchemy import get_alchemy
from qmllib.representations.fchl.fchl_kernel_functions import get_kernel_parameters
from src.fortran.ffchl19_vp_tree import fchl19_vp_tree
from src.fortran.ffchl18_vp_tree import fchl18_vp_tree

# ignore sklearn futurewarning
warnings.filterwarnings(
    "ignore",
    "'force_all_finite' was renamed to 'ensure_all_finite' in 1.6 and will be removed in 1.8.",
)


def calculate_representation(Qrepresentation, strs, **repr_kwargs):
    if Qrepresentation == "fchl":
        return generate_global_fchl19(strs, **repr_kwargs)
    elif Qrepresentation == "mbdf":
        return generate_global_mbdf(strs, **repr_kwargs)
    elif Qrepresentation == "bob":
        return generate_bob(strs, **repr_kwargs)
    elif Qrepresentation == "coulomb":
        return generate_coulomb(strs, **repr_kwargs)
    elif Qrepresentation == "mbtr":
        return generate_global_mbtr(strs, **repr_kwargs)
    elif Qrepresentation == "fchl-kernel":
        return generate_fchl18(strs, **repr_kwargs)
    elif Qrepresentation == "fchl19-kernel":
        return generate_fchl19(strs, **repr_kwargs)
    else:
        raise NotImplementedError(
            f"Representation '{Qrepresentation}' not supported with the k-NN model!"
        )


class VPTreeKNN19:

    def __init__(
        self,
        n_neighbors: int = 5,
        weights: Literal["uniform", "distance"] = "uniform",
        sigma=20.0,
        verbose=False,
    ):
        """
        Class for fast k-NN with FCHL19 representation.

        Parameters
        ----------
        n_neighbors: int, optional
            Number of neighbours (k). Defaults to 5.
        weights: "uniform" or "distance", optional
            Weighing scheme; either uniform (equal weights) or distance (reciprocal).
            Defaults to uniform.
        sigma: float, optional
            Standard deviation for the FCHL19 kernel. Defaults to 20.
        verbose: boolean, optional
            Currently unused.
        """
        # knn params
        self.k = n_neighbors
        self.weights = weights
        # kernel params
        self.sigma = sigma
        # data containers
        self.X_train: Optional[np.ndarray] = None
        self.Y_train: Optional[np.ndarray] = None
        self.X_test: Optional[np.ndarray] = None
        self.verbose: Optional[bool] = verbose

    def fit(
        self,
        X: np.ndarray,
        X_atoms: np.ndarray,
        Y: np.ndarray,
    ):
        """
        Fit k-NN model.

        Parameters
        ----------
        X: np.ndarray
            FCHL19 representation of the training data
        X_atoms: np.ndarray
            Nuclear charges (element numbers) of the molecules in the training data
        Y: np.ndarray
            Labels to be learned.
        """

        na = np.array([len(x) for x in X_atoms])
        Q_tr = np.zeros((max(na), X.shape[0]), dtype=np.int32)
        for i, q in enumerate(X_atoms):
            Q_tr[: len(q), i] = q

        nm1 = X.shape[0]

        self.X_train = X
        self.Y_train = Y

        fchl19_vp_tree.train(
            x_in=X,
            q_in=Q_tr,
            y_in=Y,
            verbose_in=self.verbose,
            n1=na,
            nm1=nm1,
            sigma_in=self.sigma,
        )
        return

    def _check_test(self, X_test):
        if self.X_test is not None:
            f_X_test = self.X_test
            # the test data is already in the fortran module, no need to reset
            if np.array_equal(f_X_test, X_test):
                return True
        return False

    def _set_test(self, X_test: np.ndarray, X_atoms_test: np.ndarray):
        if self._check_test(X_test):
            return

        na_test = np.array([len(x) for x in X_atoms_test])
        Q_te = np.zeros((max(na_test), X_test.shape[0]), dtype=np.int32)
        for i, q in enumerate(X_atoms_test):
            Q_te[: len(q), i] = q

        nm2 = X_test.shape[0]
        fchl19_vp_tree.set_up_test(
            X_test,
            Q_te,
            na_test,
            nm2,
        )
        self.X_test = X_test

    def kneighbours(
        self,
        X: np.ndarray,
        X_atoms_test: np.ndarray,
        n_neighbors: Optional[int] = None,
        return_distance=True,
    ):
        """Find closest k neighbors in the tree."""
        if n_neighbors is None:
            n_neighbors = self.k

        self._set_test(X, X_atoms_test)
        if return_distance:
            (k_neighbors, distances) = fchl19_vp_tree.kneighbors(
                n_neighbors, X.shape[0], return_distances=True
            )
            # fortran indexing is 1-based
            k_neighbors = k_neighbors - 1
            # transpose required for sk-learn compatibility as fortran is column-based
            return distances.T, k_neighbors.T
        else:
            k_neighbors = fchl19_vp_tree.kneighbors(n_neighbors, X.shape[0])
            # fortran indexing is 1-based
            k_neighbors = k_neighbors - 1
            # transpose required for sk-learn compatibility as fortran is column-based
            return k_neighbors.T

    def predict(
        self, X: np.ndarray, X_atoms_test: np.ndarray, n_neighbors: Optional[int] = None
    ):
        """Wrapper to replicate sklearn k-NN model behaviour."""
        self._set_test(X, X_atoms_test)
        if n_neighbors is None:
            n_neighbors = self.k
        if self.weights == "uniform":
            return fchl19_vp_tree.predict(n_neighbors, X.shape[0])
        elif self.weights == "distance":
            return fchl19_vp_tree.predict(
                n_neighbors, X.shape[0], weight_by_distance=True
            )

    def get_params(self, deep=True):
        """
        Get parameters for this estimator.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        out = dict()
        out["n_neighbors"] = self.k
        out["weights"] = self.weights
        return out

    def get_tree_params(self) -> Dict[str, np.ndarray]:
        """Get the vp-tree collections, which can be used to rebuild the tree later."""
        out = dict()
        out["vp_index"] = fchl19_vp_tree.vp_index
        out["vp_left"] = fchl19_vp_tree.vp_left
        out["vp_right"] = fchl19_vp_tree.vp_right
        out["vp_threshold"] = fchl19_vp_tree.vp_threshold
        return out


class VPTreeKNN18:

    def __init__(
        self,
        n_neighbors: int = 5,
        weights: Literal["uniform", "distance"] = "uniform",
        alchemy: str = "periodic-table",
        alchemy_period_width: float = 1.6,
        alchemy_group_width: float = 1.6,
        verbose: bool = False,
        two_body_scaling: float = np.sqrt(8),
        three_body_scaling: float = 1.6,
        two_body_width: float = 0.2,
        three_body_width: float = np.pi,
        two_body_power: float = 4.0,
        three_body_power: float = 2.0,
        cut_start: float = 1.0,
        cut_distance: float = 5.0,
        fourier_order: int = 1,
        kernel: str = "gaussian",
        kernel_args=None,
    ):
        """
        :param two_body_scaling: Weight for 2-body terms.
        :type two_body_scaling: float
        :param three_body_scaling: Weight for 3-body terms.
        :type three_body_scaling: float

        :param two_body_width: Gaussian width for 2-body terms
        :type two_body_width: float
        :param three_body_width: Gaussian width for 3-body terms.
        :type three_body_width: float

        :param two_body_power: Powerlaw for :math:`r^{-n}` 2-body terms.
        :type two_body_power: float
        :param three_body_power: Powerlaw for Axilrod-Teller-Muto 3-body term
        :type three_body_power: float

        :param cut_start: The fraction of the cut-off radius at which cut-off damping start.
        :type cut_start: float
        :param cut_distance: Cut-off radius. (default=5 angstrom)
        :type cut_distance: float

        :param fourier_order: 3-body Fourier-expansion truncation order.
        :type fourier_order: integer
        :param alchemy: Type of alchemical interpolation ``"periodic-table"`` or ``"off"`` are possible options. Disabling alchemical interpolation can yield dramatic speedups.
        :type alchemy: string

        :param alchemy_period_width: Gaussian width along periods (columns) in the periodic table.
        :type alchemy_period_width: float
        :param alchemy_group_width: Gaussian width along groups (rows) in the periodic table.
        :type alchemy_group_width: float
        """
        # knn params
        self.k = n_neighbors
        self.weights = weights
        # kernel params
        self.alchemy = alchemy
        self.alchemy_period_width = alchemy_period_width
        self.alchemy_group_width = alchemy_group_width
        self.verbose = verbose
        self.two_body_scaling = two_body_scaling
        self.three_body_scaling = three_body_scaling
        self.two_body_width = two_body_width
        self.three_body_width = three_body_width
        self.two_body_power = two_body_power
        self.three_body_power = three_body_power
        self.cut_start = cut_start
        self.cut_distance = cut_distance
        self.fourier_order = fourier_order
        self.kernel = kernel
        self.doalchemy, self.pd = get_alchemy(
            alchemy, emax=100, r_width=alchemy_group_width, c_width=alchemy_period_width
        )
        self.kernel_args = kernel_args
        # data containers
        self.X_train: Optional[np.ndarray] = None
        self.Y_train: Optional[np.ndarray] = None
        self.X_test: Optional[np.ndarray] = None

    def fit(
        self,
        X: np.ndarray,
        Y: np.ndarray,
    ):
        """Calculates the Gaussian kernel matrix K elements, where :math:`K_{ij}`:

            :math:`K_{ij} = \\exp \\big( -\\frac{\\|A_i - B_j\\|_2^2}{2\\sigma^2} \\big)`

        Where :math:`A_{i}` and :math:`B_{j}` are FCHL representation vectors.
        K is calculated analytically using an OpenMP parallel Fortran routine.
        Note, that this kernel will ONLY work with FCHL representations as input.

        :param A: Array of FCHL representation - shape=(N, maxsize, 5, maxneighbors).
        :type A: numpy array
        :param B: Array of FCHL representation - shape=(M, maxsize, 5, maxneighbors).
        :type B: numpy array

        """

        atoms_max = X.shape[1]

        nm1 = X.shape[0]

        N1 = np.zeros((nm1), dtype=np.int32)

        for a in range(nm1):
            N1[a] = len(np.where(X[a, :, 1, 0] > 0.0001)[0])

        neighbors1 = np.zeros((nm1, atoms_max), dtype=np.int32)

        for a, representation in enumerate(X):
            ni = N1[a]
            for i, x in enumerate(representation[:ni]):
                neighbors1[a, i] = len(np.where(x[0] < self.cut_distance)[0])

        kernel_idx, kernel_parameters, n_kernels = get_kernel_parameters(
            self.kernel, self.kernel_args
        )
        self.X_train = X
        self.Y_train = Y

        fchl18_vp_tree.train(
            x_in=X,
            y_in=Y,
            verbose_in=self.verbose,
            n1=N1,
            nneigh1_in=neighbors1,
            nm1=nm1,
            nsigmas_in=n_kernels,
            t_width_in=self.three_body_width,
            d_width_in=self.two_body_width,
            cut_start_in=self.cut_start,
            cut_distance_in=self.cut_distance,
            order_in=self.fourier_order,
            pd_in=self.pd,
            distance_scale_in=self.two_body_scaling,
            angular_scale_in=self.three_body_scaling,
            alchemy_in=self.doalchemy,
            two_body_power_in=self.two_body_power,
            three_body_power_in=self.three_body_power,
            kernel_idx_in=kernel_idx,
            parameters_in=kernel_parameters,
        )
        return

    def _check_test(self, X_test):
        if self.X_test is not None:
            f_X_test = self.X_test
            # the test data is already in the fortran module, no need to reset
            if np.array_equal(f_X_test, X_test):
                return True
        return False

    def _set_test(self, X_test: np.ndarray, cut_distance=5.0):
        if self._check_test(X_test):
            return

        atoms_max = self.X_train.shape[1]
        neighbors_max = self.X_train.shape[3]
        if not X_test.shape[1] == atoms_max:
            raise ValueError("Check FCHL representation sizes")
        if not X_test.shape[3] == neighbors_max:
            raise ValueError("Check FCHL representation sizes")
        nm2 = X_test.shape[0]
        neighbors2 = np.zeros((nm2, atoms_max), dtype=np.int32)
        N2 = np.zeros((nm2), dtype=np.int32)
        for a in range(nm2):
            N2[a] = len(np.where(X_test[a, :, 1, 0] > 0.0001)[0])
        for a, representation in enumerate(X_test):
            ni = N2[a]
            for i, x in enumerate(representation[:ni]):
                neighbors2[a, i] = len(np.where(x[0] < cut_distance)[0])
        fchl18_vp_tree.set_up_test(
            X_test,
            N2,
            neighbors2,
            nm2,
        )
        self.X_test = X_test

    def kneighbours(
        self, X: np.ndarray, n_neighbors: Optional[int] = None, return_distance=True
    ):
        """Find closest k neighbors in the tree."""
        if n_neighbors is None:
            n_neighbors = self.k

        self._set_test(X)
        if return_distance:
            (k_neighbors, distances) = fchl18_vp_tree.kneighbors(
                n_neighbors, X.shape[0], return_distances=True
            )
            # fortran indexing is 1-based
            k_neighbors = k_neighbors - 1
            # transpose required for sk-learn compatibility as fortran is column-based
            return distances.T, k_neighbors.T
        else:
            k_neighbors = fchl18_vp_tree.kneighbors(n_neighbors, X.shape[0])
            # fortran indexing is 1-based
            k_neighbors = k_neighbors - 1
            # transpose required for sk-learn compatibility as fortran is column-based
            return k_neighbors.T

    def predict(self, X: np.ndarray, n_neighbors: Optional[int] = None):
        """Wrapper to replicate sklearn k-NN model behaviour."""
        self._set_test(X)
        if n_neighbors is None:
            n_neighbors = self.k
        if self.weights == "uniform":
            return fchl18_vp_tree.predict(n_neighbors, X.shape[0])
        elif self.weights == "distance":
            return fchl18_vp_tree.predict(
                n_neighbors, X.shape[0], weight_by_distance=True
            )

    def get_params(self, deep=True):
        """
        Get parameters for this estimator.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        out = dict()
        out["n_neighbors"] = self.k
        out["weights"] = self.weights
        return out

    def get_tree_params(self) -> Dict[str, np.ndarray]:
        """Get the vp-tree collections, which can be used to rebuild the tree later."""
        out = dict()
        out["vp_index"] = fchl18_vp_tree.vp_index
        out["vp_left"] = fchl18_vp_tree.vp_left
        out["vp_right"] = fchl18_vp_tree.vp_right
        out["vp_threshold"] = fchl18_vp_tree.vp_threshold
        return out


def load_fchl18_vp_knn(X_train, Y_train, vp_params, **knn_params):
    knn = VPTreeKNN18(**knn_params)
    atoms_max = X_train.shape[1]

    nm1 = X_train.shape[0]

    N1 = np.zeros((nm1), dtype=np.int32)

    for a in range(nm1):
        N1[a] = len(np.where(X[a, :, 1, 0] > 0.0001)[0])

    neighbors1 = np.zeros((nm1, atoms_max), dtype=np.int32)

    for a, representation in enumerate(X):
        ni = N1[a]
        for i, x in enumerate(representation[:ni]):
            neighbors1[a, i] = len(np.where(x[0] < knn.cut_distance)[0])

    kernel_idx, kernel_parameters, n_kernels = get_kernel_parameters(
        knn.kernel, knn.kernel_args
    )
    knn.X_train = X_train
    knn.Y_train = Y_train
    fchl18_vp_tree.load(
        X_train,
        Y_train,
        vp_params["index"],
        vp_params["left"],
        vp_params["right"],
        vp_params["threshold"],
        knn.verbose,
        N1,
        neighbors1,
        nm1,
        n_kernels,
        knn.three_body_width,
        knn.two_body_width,
        knn.cut_start,
        knn.cut_distance,
        knn.fourier_order,
        knn.pd,
        knn.two_body_scaling,
        knn.three_body_scaling,
        knn.doalchemy,
        knn.two_body_power,
        knn.three_body_power,
        kernel_idx,
        kernel_parameters,
    )
    return knn


def load_fchl19_vp_knn(X_train, X_atoms, Y_train, vp_params, **knn_params):
    knn = VPTreeKNN19(**knn_params)
    na = np.array([len(x) for x in X_atoms])
    Q_tr = np.zeros((max(na), X.shape[0]), dtype=np.int32)
    for i, q in enumerate(X_atoms):
        Q_tr[: len(q), i] = q

    nm1 = X.shape[0]

    atoms_max = X_train.shape[1]

    nm1 = X_train.shape[0]

    knn.X_train = X_train
    knn.Y_train = Y_train
    fchl19_vp_tree.load(
        x_in=X_train,
        q_in=Q_tr,
        y_in=Y_train,
        vp_index_in=vp_params["index"],
        vp_left_in=vp_params["left"],
        vp_right_in=vp_params["right"],
        vp_threshold_in=vp_params["threshold"],
        verbose_in=knn.verbose,
        n1=na,
        nm1=nm1,
        sigma_in=knn.sigma,
    )
    return knn


def load_hyperparams(hyper_cache: Optional[Union[str, os.PathLike]]):
    """Load hyperparameters from hyper_cache or use defaults."""
    if hyper_cache is not None:
        with open(hyper_cache, "rb") as f:
            hyperparams: dict = pickle.load(f)
        print(f"JKML(Q-kNN): Loaded hyperparameters from {hyper_cache}:", flush=True)
        if "representation" not in hyperparams:
            hyperparams["representation"] = {"cutoff": 8.0}
            print(
                f"JKML (Q-kNN): Representation params not defined in hyperparams, use default value."
            )
        if "knn" not in hyperparams:
            hyperparams["knn"] = {"n_neighbors": 5, "weights": "uniform"}
            print(
                f"JKML (Q-kNN): knn params not defined in hyperparams, use default values."
            )
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
    sigmas,
    no_metric: bool = False,
    hyper_cache: Optional[Union[str, os.PathLike]] = None,
    subsample_mlkr: bool = False,
):

    hyperparams = load_hyperparams(hyper_cache)
    # Convert structures to list to avoid indexing issues
    for i, struct in enumerate(strs.values):
        assert struct == strs.iloc[i]
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
        pickle.dump([X_train, Y_train, strs], f)
        print(f"Saved pretrain vars to {varsoutfile}.", flush=True)
    train_wall_start = time.perf_counter()
    train_cpu_start = time.process_time()
    if not no_metric and "-kernel" not in Qrepresentation:
        print("JKML(Q-kNN): Training MLKR metric.", flush=True)
        # Limit the number of MLKR components for faster training
        mlkr = MLKR(n_components=50)
        subsample_size = 25_000
        if subsample_mlkr and (X_train.shape[0] > subsample_size):
            subsample_indices = np.random.permutation(X_train.shape[0])[:subsample_size]
            X_mlkr, Y_mlkr = X_train[subsample_indices, :], Y_train[subsample_indices]
        else:
            X_mlkr, Y_mlkr = X_train, Y_train
        if X_mlkr.shape[0] > 40_000:
            # cast to float 32 to save memory
            X_mlkr = X_mlkr.astype(np.float32)
            Y_mlkr = Y_mlkr.astype(np.float32)
        mlkr.fit(X_mlkr, Y_mlkr)
        A = mlkr.get_mahalanobis_matrix()
        print("JKML(Q-kNN): Training k-NN regressor with MLKR metric.")
        knn = KNeighborsRegressor(
            metric=mlkr.get_metric(),
            n_jobs=-1,
            algorithm="ball_tree",
            **hyperparams["knn"],
        )
    elif Qrepresentation == "fchl-kernel":

        print("JKML(Q-kNN): Learn VP-tree of kernel distances.")
        knn = VPTreeKNN18(kernel_args={"sigma": [1.0]}, **hyperparams["knn"])
    elif Qrepresentation == "fchl19-kernel":

        print("JKML(Q-kNN): Learn VP-tree of kernel distances.")
        knn = VPTreeKNN19(sigma=sigmas[0], **hyperparams["knn"])
    else:
        # "vanilla" k-NN
        knn = KNeighborsRegressor(n_jobs=-1, algorithm="auto", **hyperparams["knn"])

    if Qrepresentation == "fchl19-kernel":
        knn.fit(X_train, X_atoms, Y_train)
    else:
        knn.fit(X_train, Y_train)
    train_wall = time.perf_counter() - train_wall_start
    train_cpu = time.process_time() - train_cpu_start
    n_train, d_train = X_train.shape[0], np.sum(X_train.shape[1:])
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
        if no_metric:
            pickle.dump([X_train, Y_train, X_atoms, knn_params, train_metadata], f)
        elif "-kernel" in Qrepresentation:
            vp_params = knn.get_tree_params()
            pickle.dump(
                [
                    X_train,
                    Y_train,
                    X_atoms,
                    knn_params,
                    vp_params,
                    train_metadata,
                ],
                f,
            )
        else:
            pickle.dump(
                [
                    X_train,
                    Y_train,
                    X_atoms,
                    A,
                    mlkr,
                    knn_params,
                    train_metadata,
                ],
                f,
            )
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
    # Convert structures to list to avoid indexing issues
    for i, struct in enumerate(strs.values):
        assert struct == strs.iloc[i]
    strs = strs.values
    X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]
    repr_wall_start = time.perf_counter()
    repr_cpu_start = time.process_time()
    X_test = calculate_representation(
        Qrepresentation, strs, **hyperparams["representation"]
    )
    if Qrepresentation == "fchl-kernel":
        X_test, X_train = correct_fchl18_kernel_size(X_test, X_train)
    print("JKML(k-NN): Calculate test kernel(s).", flush=True)
    repr_test_wall = time.perf_counter() - repr_wall_start
    repr_test_cpu = time.process_time() - repr_cpu_start

    # some info about the full representation
    print(
        "JKML(k-NN): Shape of the testing representation: " + str(X_test.shape),
        flush=True,
    )

    test_wall_start = time.perf_counter()
    test_cpu_start = time.process_time()
    if Qrepresentation == "fchl19-kernel":
        Y_predicted = knn_model.predict(X_test, X_atoms)
    else:
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
    sigmas,
    no_metric=False,
    cv_folds=5,
    verbose=True,
    optimise_representation=False,
    hyper_cache=None,
):

    import skopt
    from skopt.space import Real
    from functools import lru_cache
    from sklearn.model_selection import cross_val_score, KFold
    from sklearn.metrics import pairwise_distances

    # Convert structures to list to avoid indexing issues
    for i, struct in enumerate(strs.values):
        assert struct == strs.iloc[i]
    strs = strs.values
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
        if hyper_cache is not None:
            print("JKML(Q-kNN): Loading representation hyperparams.", flush=True)
            hyperparams = load_hyperparams(hyper_cache)
            repr_params = hyperparams["representation"]
        else:
            repr_params = {}
            hyperparams = {"knn": {}, "representation": {}}
        global X
        global X_atoms
        X = calculate_representation(Qrepresentation, strs, **repr_params)
        X_atoms = [strs[i].get_atomic_numbers() for i in range(len(strs))]

    # add k-nn specific hyperparameters
    max_k = 15
    space.append(skopt.space.Integer(3, max_k, name="n_neighbors"))
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

        print("JKML(Q-kNN): Precalculating distance matrices for hyperopt.", flush=True)
        precalc_start = time.perf_counter()
        # precalculate distances and sorted Y matrices for SPEED
        kf = KFold(cv_folds)
        # could preallocate, but won't bother >:)
        distance_matrices = []
        sorted_Ys = []
        if "-kernel" in Qrepresentation:
            knns = []
        if not no_metric:
            mlkrs = []

        for i, (train_index, test_index) in enumerate(kf.split(X)):
            fold_start = time.perf_counter()
            X_fold, Y_fold = (
                X[train_index],
                Y_train[train_index],
            )
            X_test = X[test_index]
            X_atoms_fold = [X_atoms[j] for j in train_index]
            X_atoms_test = [X_atoms[j] for j in test_index]
            if no_metric:
                D = pairwise_distances(X_test, X_fold, n_jobs=-1)
            elif Qrepresentation == "fchl-kernel":
                knn = VPTreeKNN18(kernel_args={"sigma": [1.0]})
                knn.fit(X_fold, Y_fold)
                D, neighbors = knn.kneighbours(X_test, n_neighbors=max_k)
                Y_fold = Y_fold[neighbors]
            elif Qrepresentation == "fchl19-kernel":
                knn = VPTreeKNN19(sigma=sigmas[0], **hyperparams["knn"])
                knn.fit(X_fold, X_atoms_fold, Y_fold)
                D, neighbors = knn.kneighbours(X_test, X_atoms_test, n_neighbors=max_k)
                Y_fold = Y_fold[neighbors]
            else:
                mlkr = MLKR(n_components=50)
                mlkr.fit(X_fold, Y_fold)
                mlkrs.append(mlkr)
                D = pairwise_distances(
                    X_test, X_fold, metric=mlkr.get_metric(), n_jobs=-1
                )
            sorted_indices = np.argsort(D, axis=1)
            # presort distance matrix
            D = np.take_along_axis(D, sorted_indices, axis=1)
            distance_matrices.append(D)
            # also store sorted y_train for prediction
            if "-kernel" in Qrepresentation:
                # Y_fold is already made into a matrix; hence need to take along axis
                Y_sorted = np.take_along_axis(Y_fold, indices=sorted_indices, axis=1)
            else:
                # Y_fold is still a vector; can index directly
                Y_sorted = Y_fold[sorted_indices]
            sorted_Ys.append(Y_sorted)
            print(
                f"\tFold {i+1}/{cv_folds} done, took {time.perf_counter() - fold_start:.1f} s.",
                flush=True,
            )

        print(
            f"JKML(k-NN): Precalculation done, took {time.perf_counter() - precalc_start:.1f} s.",
            flush=True,
        )

        @skopt.utils.use_named_args(space)
        @lru_cache
        def objective(n_neighbors, weights, **repr_params):
            maes = np.zeros(cv_folds)
            # manual k-NN implementation
            for i, (_, test_index) in enumerate(kf.split(X)):
                Y_test = Y_train[test_index]
                if weights == "uniform":
                    yhat = np.mean(sorted_Ys[i][:, :n_neighbors], axis=1)
                elif weights == "distance":
                    w = 1 / distance_matrices[i][:, :n_neighbors]
                    yhat = np.average(sorted_Ys[i][:, :n_neighbors], axis=1, weights=w)
                maes[i] = np.mean(np.abs(Y_test - yhat))
            return np.mean(maes)

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
    for key, value in hyperparams["knn"].items():
        if key not in params["knn"]:
            # copy over value to final hyperdict (such as FCHL19 sigma)
            params["knn"][key] = value

    if (hyper_cache is not None) and (not optimise_representation):
        # use provided representation hyperparams
        params["representation"] = repr_params
    with open(hyperparamfile, "wb") as f:
        pickle.dump(params, f)
        print(f"JKML(Q-kNN): Saved hyperparams to {hyperparamfile}")

    return params
