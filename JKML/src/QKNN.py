###############################################################################
###############################################################################
###############################################################################

import numpy as np
import pandas as pd
import pickle
from sklearn.neighbors import KNeighborsRegressor
from metric_learn import MLKR
from ase.atoms import Atoms
from typing import Iterable, Tuple, Union, Dict, Literal, Callable, List, Optional
import os
from collections import defaultdict
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, Future
import heapq

# ignore sklearn futurewarning
warnings.filterwarnings(
    "ignore",
    "'force_all_finite' was renamed to 'ensure_all_finite' in 1.6 and will be removed in 1.8.",
)


def _generate_fchl19(
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
    X = np.zeros((n, representation.shape[1]))
    X[0, :] = np.sum(representation, axis=0)
    for i in range(1, n):
        X[i, :] = generate_representation(
            strs[i].get_atomic_numbers(),
            strs[i].get_positions(),
            elements=elements,
            rcut=rcut,
            acut=acut,
            pad=max_atoms,
        ).sum(axis=0)
    if np.isnan(X).any():
        raise ValueError("NaNs in FCHL representation!")
    return X


def _generate_mbdf(
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


def _generate_bob(
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


def _generate_coulomb(strs: Iterable[Atoms], max_atoms: int = None, **kwargs):
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


def _generate_fchl18(strs: Iterable[Atoms], max_atoms=None, cutoff=8.0):
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
    elif Qrepresentation == "fchl-kernel":
        return _generate_fchl18(strs, **repr_kwargs)
    else:
        raise NotImplementedError(
            f"Representation 'f{Qrepresentation}' not supported with the k-NN model!"
        )


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


def calculate_kernel(X: np.ndarray, X_other: np.ndarray = None, **kernel_kwargs):
    from qmllib.representations.fchl import (
        get_local_symmetric_kernels as JKML_sym_kernel,
    )
    from qmllib.representations.fchl import get_local_kernels as JKML_kernel

    if X_other is not None:
        return JKML_kernel(X, X_other, **kernel_kwargs)
    else:
        return JKML_sym_kernel(X, **kernel_kwargs)


def induced_kernel_distance(
    K: np.ndarray, K_train: np.ndarray = None, K_test: np.ndarray = None
):

    def _remove_extra_dim(K: np.ndarray, kernel_name: str):
        # FCHL kernels return a matrix with an extra leading dimension; this function drops that.
        if K.ndim > 2:
            if K.shape[0] == 1:
                K = K[0]
            else:
                raise ValueError(
                    f"Incompatible {kernel_name} kernel shape! (got shape {K.shape})"
                )
        return K

    K = _remove_extra_dim(K, "K")
    if K_train is None:
        diag_1 = diag_2 = np.diag(K)
    else:
        assert K_test is not None, "Need both self-similarity kernels!"
        K_train = _remove_extra_dim(K_train, "train")
        K_test = _remove_extra_dim(K_test, "test")
        diag_1 = np.diag(K_test)
        diag_2 = np.diag(K_train)
    # first two terms broadcast the self-similarities to a [n x m] matrix
    D_squared = diag_1[:, None] + diag_2[None, :] - 2 * K
    # get rid of possible numerical problems
    D_squared = np.maximum(D_squared, 0.0)
    D = np.sqrt(D_squared)
    return D


class VPTreeNode:
    def __init__(self, idx: int, threshold: float, left=None, right=None):
        self.idx = idx
        self.threshold = threshold
        self.left = left
        self.right = right


class VPTreeKNN:
    def __init__(
        self,
        kernel_fun: Callable[[int, int, np.ndarray, np.ndarray], float],
        n_neighbors: int = 5,
        n_jobs: int = -1,
        weights: Literal["uniform", "distance"] = "uniform",
    ):
        self.kernel = kernel_fun
        self.k = n_neighbors
        self.max_workers = int(
            os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())
            if n_jobs == -1
            else n_jobs
        )
        self.weights = weights

    def _build_vptree(
        self,
        indices: List[int],
        dist_fn: Callable[[int, int], float],
        executor: Optional[ThreadPoolExecutor] = None,
        depth: int = 0,
        max_parallel_depth: int = 4,
    ) -> Optional[VPTreeNode]:
        if not indices:
            return None
        if len(indices) == 1:
            return VPTreeNode(idx=indices[0], threshold=0.0)

        vp = indices[-1]
        others = indices[:-1]
        dists = [dist_fn(vp, i) for i in others]
        median = float(np.median(np.array(dists)))

        left = [i for i, d in zip(others, dists) if d <= median]
        right = [i for i, d in zip(others, dists) if d > median]

        if executor and depth < max_parallel_depth:
            left_future: Future = executor.submit(
                self._build_vptree,
                left,
                dist_fn,
                executor,
                depth + 1,
                max_parallel_depth,
            )
            right_future: Future = executor.submit(
                self._build_vptree,
                right,
                dist_fn,
                executor,
                depth + 1,
                max_parallel_depth,
            )
            left_node = left_future.result()
            right_node = right_future.result()
        else:
            left_node = self._build_vptree(
                left, dist_fn, executor, depth + 1, max_parallel_depth
            )
            right_node = self._build_vptree(
                right, dist_fn, executor, depth + 1, max_parallel_depth
            )
        return VPTreeNode(idx=vp, threshold=median, left=left_node, right=right_node)

    def _search_vptree(
        self,
        node: VPTreeNode,
        query_idx: int,
        dist_fn: Callable[[int, int], float],
        heap: List[Tuple[float, int]],
        k: int = None,
    ):
        if k is None:
            k = self.k
        if node is None:
            return

        d = dist_fn(query_idx, node.idx)
        heapq.heappush(heap, (-d, node.idx))
        if len(heap) > k:
            heapq.heappop(heap)

        if node.left is None and node.right is None:
            return

        if d < node.threshold:
            near, far = node.left, node.right
        else:
            near, far = node.right, node.left

        self._search_vptree(near, query_idx, dist_fn, heap, k=k)

        # Prune using triangle inequality
        if len(heap) < k or abs(d - node.threshold) < -heap[0][0]:
            self._search_vptree(far, query_idx, dist_fn, heap, k=k)

    def fit(
        self,
        X: np.ndarray,
        Y: np.ndarray,
    ):
        """Build VP-tree on training data."""

        def _train_dist_fn(i: int, j: int):
            d_squared = (
                self.kernel(i, i, X, X)
                + self.kernel(j, j, X, X)
                - 2 * self.kernel(i, j, X, X)
            )
            d_squared = np.maximum(d_squared, 0.0)
            return np.sqrt(d_squared)

        self.vp_tree = self._build_vptree(list(range(X.shape[0])), _train_dist_fn)
        self.X_train = X
        self.Y_train = Y

    def kneighbours(self, X: np.ndarray, n_neighbors: int = None, return_distance=True):
        """Find closest k neighbors in the tree."""
        if n_neighbors is None:
            n_neighbors = self.k

        def _test_dist_fn(i: int, j: int):
            d_squared = (
                self.kernel(i, i, X, X)
                + self.kernel(j, j, self.X_train, self.X_train)
                - 2 * self.kernel(i, j, X, self.X_train)
            ).sum()
            d_squared = np.maximum(d_squared, 0.0)
            return np.sqrt(d_squared)

        def _find_single(test_idx: int):
            heap = []
            self._search_vptree(
                self.vp_tree, test_idx, _test_dist_fn, heap, k=n_neighbors
            )
            neighbors = np.array([idx for (_, idx) in heap])
            d = np.array([-d for (d, _) in heap])
            return neighbors, d

        m_test = X.shape[0]
        if self.max_workers is not None:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                n_d = list(executor.map(_find_single, range(m_test)))
            neighbors = np.stack([x[0] for x in n_d])
            D = np.stack([x[1] for x in n_d])
        else:
            neighbors = np.zeros((m_test, n_neighbors))
            D = np.zeros_like(neighbors)
            for i in range(m_test):
                n, d = _find_single(i)
                neighbors[i, :] = n
                D[i, :] = d
        neighbors = neighbors.astype(np.uint32)
        if return_distance:
            return D, neighbors
        else:
            return neighbors

    def predict(self, X: np.ndarray, k: int = None):
        """Wrapper to replicate sklearn k-NN model behaviour."""
        D, neighbors = self.kneighbours(X, k)
        Y = self.Y_train[neighbors]
        if self.weights == "uniform":
            return np.mean(Y, axis=1)
        else:
            w = 1 / D
            return np.average(Y, weights=w, axis=1)

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
        out["n_jobs"] = self.max_workers
        return out


def fast_kernel(
    alchemy: str = "periodic-table",
    alchemy_period_width: float = 1.6,
    alchemy_group_width: float = 1.6,
):

    from qmllib.utils.alchemy import get_alchemy

    doalchemy, pd = get_alchemy(
        alchemy, emax=100, r_width=alchemy_group_width, c_width=alchemy_period_width
    )
    return lambda X, Y, k_args: get_local_kernels(
        X, Y, doalchemy, pd, kernel_args=k_args
    )


def get_local_kernels(
    A: np.ndarray,
    B: np.ndarray,
    doalchemy,
    pd,
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
    kernel_args: Optional[Dict[str, List[float]]] = None,
) -> np.ndarray:
    """Calculates the Gaussian kernel matrix K, where :math:`K_{ij}`:

        :math:`K_{ij} = \\exp \\big( -\\frac{\\|A_i - B_j\\|_2^2}{2\\sigma^2} \\big)`

    Where :math:`A_{i}` and :math:`B_{j}` are FCHL representation vectors.
    K is calculated analytically using an OpenMP parallel Fortran routine.
    Note, that this kernel will ONLY work with FCHL representations as input.

    :param A: Array of FCHL representation - shape=(N, maxsize, 5, maxneighbors).
    :type A: numpy array
    :param B: Array of FCHL representation - shape=(M, maxsize, 5, maxneighbors).
    :type B: numpy array

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

    :return: Array of FCHL kernel matrices matrix - shape=(n_sigmas, N, M),
    :rtype: numpy array
    """

    from qmllib.representations.fchl.ffchl_module import fget_kernels_fchl
    from qmllib.representations.fchl.fchl_kernel_functions import get_kernel_parameters

    atoms_max = A.shape[1]
    neighbors_max = A.shape[3]

    if not B.shape[1] == atoms_max:
        raise ValueError("Check FCHL representation sizes")
    if not B.shape[3] == neighbors_max:
        raise ValueError("Check FCHL representation sizes")

    nm1 = A.shape[0]
    nm2 = B.shape[0]

    N1 = np.zeros((nm1), dtype=np.int32)
    N2 = np.zeros((nm2), dtype=np.int32)

    for a in range(nm1):
        N1[a] = len(np.where(A[a, :, 1, 0] > 0.0001)[0])

    for a in range(nm2):
        N2[a] = len(np.where(B[a, :, 1, 0] > 0.0001)[0])

    neighbors1 = np.zeros((nm1, atoms_max), dtype=np.int32)
    neighbors2 = np.zeros((nm2, atoms_max), dtype=np.int32)

    for a, representation in enumerate(A):
        ni = N1[a]
        for i, x in enumerate(representation[:ni]):
            neighbors1[a, i] = len(np.where(x[0] < cut_distance)[0])

    for a, representation in enumerate(B):
        ni = N2[a]
        for i, x in enumerate(representation[:ni]):
            neighbors2[a, i] = len(np.where(x[0] < cut_distance)[0])

    kernel_idx, kernel_parameters, n_kernels = get_kernel_parameters(
        kernel, kernel_args
    )

    return fget_kernels_fchl(
        A,
        B,
        verbose,
        N1,
        N2,
        neighbors1,
        neighbors2,
        nm1,
        nm2,
        n_kernels,
        three_body_width,
        two_body_width,
        cut_start,
        cut_distance,
        fourier_order,
        pd,
        two_body_scaling,
        three_body_scaling,
        doalchemy,
        two_body_power,
        three_body_power,
        kernel_idx,
        kernel_parameters,
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
    for i, struct in enumerate(strs.values):
        assert struct == strs.iloc[i]
    print("Values works!")
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
    if not no_metric and Qrepresentation != "fchl-kernel":
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
    elif Qrepresentation == "fchl-kernel":
        kernel_fun = fast_kernel()

        def kernel(i, j, X1, X2):
            return kernel_fun(X1[None, i], X2[None, j], {"sigma": [1.0]})

        print("JKML(Q-kNN): Learn VP-tree of kernel distances.")
        knn = VPTreeKNN(kernel_fun=kernel, n_jobs=-1, **hyperparams["knn"])
    else:
        # "vanilla" k-NN
        knn = KNeighborsRegressor(n_jobs=-1, algorithm="auto", **hyperparams["knn"])

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
        elif Qrepresentation == "fchl-kernel":
            pickle.dump(
                [
                    X_train,
                    Y_train,
                    X_atoms,
                    knn_params,
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
        global X
        X = calculate_representation(Qrepresentation, strs)
        if Qrepresentation == "fchl-kernel":
            kernel_fun = fast_kernel()

            def kernel(i, j, X1, X2):
                return kernel_fun(X1[None, i], X2[None, j], {"sigma": [1.0]})

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

        print("JKML(k-NN): Precalculating distance matrices for hyperopt.", flush=True)
        precalc_start = time.perf_counter()
        # precalculate distances and sorted Y matrices for SPEED
        kf = KFold(cv_folds)
        # could preallocate, but won't bother >:)
        distance_matrices = []
        sorted_Ys = []
        if Qrepresentation == "fchl-kernel":
            knns = []
        if not no_metric:
            mlkrs = []

        for i, (train_index, test_index) in enumerate(kf.split(X)):
            fold_start = time.perf_counter()
            X_fold, Y_fold = X[train_index], Y_train[train_index]
            X_test = X[test_index]
            if no_metric:
                D = pairwise_distances(X_test, X_fold, n_jobs=-1)
            elif Qrepresentation == "fchl-kernel":
                knn = VPTreeKNN(kernel, n_neighbors=15, weights="uniform")
                knn.fit(X_fold, Y_fold)
                D, neighbors = knn.kneighbours(X_test, n_neighbors=max_k)
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
            if Qrepresentation == "fchl-kernel":
                # Y_fold is already made into a matrix; hence need to take along axis
                Y_sorted = np.take_along_axis(Y_fold, indices=sorted_indices, axis=1)
            else:
                # Y_fold is still a vector; can index directly
                Y_sorted = Y_fold[sorted_indices]
            sorted_Ys.append(Y_sorted)
            print(f"\tFold {i+1}/{cv_folds} done, took {time.perf_counter() - fold_start:.1f} s.", flush=True)


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

    with open(hyperparamfile, "wb") as f:
        pickle.dump(params, f)
        print(f"JKML(Q-kNN): Saved hyperparams to {hyperparamfile}")

    return params
