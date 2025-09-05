"""
Batched variant of Metric Learning for Kernel Regression (MLKR)
"""

import time
import sys
import os
import warnings
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from sklearn.base import TransformerMixin
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import pairwise_distances_chunked, pairwise_distances

from metric_learn.base_metric import MahalanobisMixin
from metric_learn._util import _initialize_components, _check_n_components

EPS = np.finfo(float).eps


class BatchedMLKR(MahalanobisMixin, TransformerMixin):
    """Batched variant Metric Learning for Kernel Regression (MLKR)

  Adapted from the original MLKR implementation to support batched
  computation of the loss and gradient. This allows for more efficient
  training on large datasets by processing the data in chunks, which is
  particularly useful when the dataset is too large to fit into memory
  all at once. The original MLKR algorithm is designed by Weinberger and Tesauro [1].

  MLKR is an algorithm for supervised metric learning, which learns a
  distance function by directly minimizing the leave-one-out regression error.
  This algorithm can also be viewed as a supervised variation of PCA and can be
  used for dimensionality reduction and high dimensional data visualization.

  Parameters
  ----------
  n_components : int or None, optional (default=None)
    Dimensionality of reduced space (if None, defaults to dimension of X).

  init : string or numpy array, optional (default='auto')
    Initialization of the linear transformation. Possible options are
    'auto', 'pca', 'identity', 'random', and a numpy array of shape
    (n_features_a, n_features_b).

    'auto'
      Depending on ``n_components``, the most reasonable initialization
      will be chosen. If ``n_components < min(n_features, n_samples)``,
      we use 'pca', as it projects data in meaningful directions (those
      of higher variance). Otherwise, we just use 'identity'.

    'pca'
      ``n_components`` principal components of the inputs passed
      to :meth:`fit` will be used to initialize the transformation.
      (See `sklearn.decomposition.PCA`)

    'identity'
      If ``n_components`` is strictly smaller than the
      dimensionality of the inputs passed to :meth:`fit`, the identity
      matrix will be truncated to the first ``n_components`` rows.

    'random'
      The initial transformation will be a random array of shape
      `(n_components, n_features)`. Each value is sampled from the
      standard normal distribution.

    numpy array
      n_features_b must match the dimensionality of the inputs passed to
      :meth:`fit` and n_features_a must be less than or equal to that.
      If ``n_components`` is not None, n_features_a must match it.

  tol : float, optional (default=None)
    Convergence tolerance for the optimization.

  max_iter : int, optional (default=1000)
    Cap on number of conjugate gradient iterations.

  verbose : bool, optional (default=False)
    Whether to print progress messages or not.

  preprocessor : array-like, shape=(n_samples, n_features) or callable
    The preprocessor to call to get tuples from indices. If array-like,
    tuples will be formed like this: X[indices].

  random_state : int or numpy.RandomState or None, optional (default=None)
    A pseudo random number generator object or a seed for it if int. If
    ``init='random'``, ``random_state`` is used to initialize the random
    transformation. If ``init='pca'``, ``random_state`` is passed as an
    argument to PCA when initializing the transformation.
  
  dtype : numpy.dtype, optional (default=np.float64)
    The data type of the input data and the learned transformation.
  
  max_non_batched : int, optional (default=20_000)
    The maximum number of samples for which the non-batched loss computation
    will be used. If the number of samples exceeds this value, the batched
    loss computation will be used instead. This is useful for large datasets
    where the batched computation can reduce the memory footprint.

  Attributes
  ----------
  n_iter_ : `int`
    The number of iterations the solver has run.

  components_ : `numpy.ndarray`, shape=(n_components, n_features)
    The learned linear transformation ``L``.

  References
  ----------
  .. [1] K.Q. Weinberger and G. Tesauto. `Metric Learning for Kernel
         Regression <http://proceedings.mlr.press/v2/weinberger07a\
         /weinberger07a.pdf>`_. AISTATS 2007.
  """

    def __init__(
        self,
        n_components=None,
        init="auto",
        tol=None,
        max_iter=1000,
        verbose=False,
        preprocessor=None,
        random_state=None,
        dtype=np.float64,
        max_non_batched=25_000,
    ):
        self.n_components = n_components
        self.init = init
        self.tol = tol
        self.max_iter = max_iter
        self.verbose = verbose
        self.random_state = random_state
        self.dtype: np.dtype = dtype
        self.max_non_batched = max_non_batched
        # get working memory;
        # JKML sets SLURM_MEM_PER_CPU rather than SLURM_MEM_PER_NODE
        self.dist_memory = os.environ.get("SLURM_MEM_PER_CPU")
        if self.dist_memory is not None:
            self.dist_memory = float(self.dist_memory)
        super(BatchedMLKR, self).__init__(preprocessor)

    def fit(self, X, y):
        """
        Fit MLKR model

        Parameters
        ----------
        X : (n x d) array of samples
        y : (n) data labels
        """
        X, y = self._prepare_inputs(
            X, y, y_numeric=True, ensure_min_samples=2, **{"dtype": self.dtype}
        )
        n, d = X.shape
        if y.shape[0] != n:
            raise ValueError(
                "Data and label lengths mismatch: %d != %d" % (n, y.shape[0])
            )

        m = _check_n_components(d, self.n_components)
        m = self.n_components
        if m is None:
            m = d
        # if the init is the default (None), we raise a warning
        A = _initialize_components(
            m,
            X,
            y,
            init=self.init,
            random_state=self.random_state,
            # MLKR works on regression targets:
            has_classes=False,
        )

        # Measure the total training time
        train_time = time.time()

        self.n_iter_ = 0
        if n > self.max_non_batched:
            print("Using batched loss computation for MLKR.", flush=True)
            # enforce verbose
            verbose_temp = self.verbose
            self.verbose = True
            res = minimize(
                self._batched_loss,
                A.ravel(),
                (X, y),
                method="L-BFGS-B",
                jac=True,
                tol=self.tol,
                options=dict(maxiter=self.max_iter),
            )
            self.verbose = verbose_temp
        else:
            res = minimize(
                self._loss,
                A.ravel(),
                (X, y),
                method="L-BFGS-B",
                jac=True,
                tol=self.tol,
                options=dict(maxiter=self.max_iter),
            )
        self.components_ = res.x.reshape(A.shape)

        # Stop timer
        train_time = time.time() - train_time
        if self.verbose:
            cls_name = self.__class__.__name__
            # Warn the user if the algorithm did not converge
            if not res.success:
                warnings.warn(
                    "[{}] MLKR did not converge: {}".format(cls_name, res.message),
                    ConvergenceWarning,
                )
            print("[{}] Training took {:8.2f}s.".format(cls_name, train_time))

        return self

    def _batched_loss(self, flatA, X, y):

        if self.n_iter_ == 0 and self.verbose:
            header_fields = ["Iteration", "Objective Value", "Time(s)"]
            header_fmt = "{:>10} {:>20} {:>10}"
            header = header_fmt.format(*header_fields)
            cls_name = self.__class__.__name__
            print("[{cls}]".format(cls=cls_name))
            print(
                "[{cls}] {header}\n[{cls}] {sep}".format(
                    cls=cls_name, header=header, sep="-" * len(header)
                )
            )

        start_time = time.time()

        A = flatA.reshape((-1, X.shape[1]))
        N, D = X.shape
        X_embedded = X @ A.T
        E = X_embedded.shape[1]

        # Accumulators
        total_cost = 0.0
        C1 = np.zeros((E, D), dtype=self.dtype)  # for Z^T W X
        C2 = np.zeros((E, D), dtype=self.dtype)  # for (W Z)^T X
        colsum = np.zeros(N, dtype=self.dtype)  # for diag term: col sums of W

        row_start = 0
        chunk_size = -1

        for dist_chunk in pairwise_distances_chunked(
            X_embedded,
            X_embedded,
            metric="sqeuclidean",
            n_jobs=-1,
            working_memory=self.dist_memory,
        ):
            B = dist_chunk.shape[0]
            if chunk_size == -1:
                chunk_size = B
            i0, i1 = row_start, row_start + B

            # Exclude self-edges in this block: set diagonal entries to +inf
            rows = np.arange(B)
            dist_chunk[rows, rows + row_start] = np.inf

            # Stable softmax over s = -dist
            s = -dist_chunk
            s_max = np.max(s, axis=1, keepdims=True)  # (B,1)
            w = np.exp(s - s_max)  # unnormalized weights (B,N)
            denom = w.sum(axis=1, keepdims=True)  # (B,1)
            zero_mask = denom == 0
            if np.any(zero_mask):
                # handle degenerate case (e.g., N==1)
                denom[zero_mask] = 1.0

            soft = w / denom  # softmax (B,N)

            # Predictions and residuals for this block
            yhat_block = soft @ y  # (B,)
            ydiff_block = yhat_block - y[i0:i1]  # (B,)
            total_cost += float(np.dot(ydiff_block, ydiff_block))

            # Build W_block in-place to save memory:
            # W_ij = soft_ij * ydiff_i * (y_j - yhat_i)
            soft *= ydiff_block[:, None]  # (B,N)
            soft *= y[None, :] - yhat_block[:, None]  # (B,N)  now 'soft' == W_bl

            # Accumulate gradient pieces:
            X_block = X[i0:i1, :]  # (B,D)
            X_emb_block = X_embedded[i0:i1, :]  # (B,E)

            # Term 1: Z^T W X  -> sum_i Z_i^T * (sum_j W_ij X_j)  == Z_block^T @ (W_block @ X)
            WX = soft @ X  # (B,D)
            C1 += X_emb_block.T @ WX  # (E,D)

            # Term 2: (W Z)^T X -> sum_i (sum_j W_ij Z_j)^T * X_i == (W_block @ Z)^T @ X_block
            WZ = soft @ X_embedded  # (B,E)
            C2 += WZ.T @ X_block  # (E,D)

            # Column sum for diagonal correction: diag(W_sym) = -colsum(W)
            colsum += soft.sum(axis=0)  # (N,)

            # Next block
            row_start = i1

            # help GC
            del dist_chunk, s, s_max, w, denom, yhat_block, ydiff_block, WX, WZ, soft

        # Term 3: - ((colsum[:,None] * Z).T @ X)  (diagonal overwrite)
        C3 = -((colsum[:, None] * X_embedded).T @ X)  # (E,D)

        grad = 4.0 * (C1 + C2 + C3)  # (E,D)

        if self.verbose:
            if self.n_iter_ == 0:
                print(f"----------Batch size = {chunk_size}-----------")
            start_time = time.time() - start_time
            values_fmt = "[{cls}] {n_iter:>10} {loss:>20.6e} {start_time:>10.2f}"
            print(
                values_fmt.format(
                    cls=self.__class__.__name__,
                    n_iter=self.n_iter_,
                    loss=total_cost,
                    start_time=start_time,
                )
            )
            sys.stdout.flush()

        self.n_iter_ += 1

        return total_cost, grad.ravel()

    def _loss(self, flatA, X, y):

        if self.n_iter_ == 0 and self.verbose:
            header_fields = ["Iteration", "Objective Value", "Time(s)"]
            header_fmt = "{:>10} {:>20} {:>10}"
            header = header_fmt.format(*header_fields)
            cls_name = self.__class__.__name__
            print("[{cls}]".format(cls=cls_name))
            print(
                "[{cls}] {header}\n[{cls}] {sep}".format(
                    cls=cls_name, header=header, sep="-" * len(header)
                )
            )

        start_time = time.time()

        A = flatA.reshape((-1, X.shape[1]))
        X_embedded = np.dot(X, A.T)
        dist = pairwise_distances(X_embedded, squared=True)
        np.fill_diagonal(dist, np.inf)
        softmax = np.exp(-dist - logsumexp(-dist, axis=1)[:, np.newaxis])
        yhat = softmax.dot(y)
        ydiff = yhat - y
        cost = (ydiff**2).sum()

        # also compute the gradient
        W = softmax * ydiff[:, np.newaxis] * (y - yhat[:, np.newaxis])
        W_sym = W + W.T
        np.fill_diagonal(W_sym, -W.sum(axis=0))
        grad = 4 * (X_embedded.T.dot(W_sym)).dot(X)

        if self.verbose:
            start_time = time.time() - start_time
            values_fmt = "[{cls}] {n_iter:>10} {loss:>20.6e} {start_time:>10.2f}"
            print(
                values_fmt.format(
                    cls=self.__class__.__name__,
                    n_iter=self.n_iter_,
                    loss=cost,
                    start_time=start_time,
                )
            )
            sys.stdout.flush()

        self.n_iter_ += 1

        return cost, grad.ravel()
