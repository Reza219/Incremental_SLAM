from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp


@dataclass(frozen=True)
class InformationGain:
    """Log-determinant information quantities used by IGG.

    eta is the paper's information surrogate
        eta_t = sum_i log |rho_{t,i}| = 0.5 log det(H_t), H_t = J_t.T J_t.

    delta_eta is the detrended increment
        Delta eta_t = eta_t - (N_{t-1}/N_t) eta_{t-1}.
    """

    eta: float
    delta_eta: float
    threshold: float

    @property
    def triggers_global_update(self) -> bool:
        return bool(self.delta_eta >= self.threshold)


def _dense_eta_from_spd(H: sp.spmatrix | npt.ArrayLike) -> float:
    if sp.issparse(H):
        A = H.toarray()
    else:
        A = np.asarray(H, dtype=float)
    A = 0.5 * (A + A.T)
    try:
        R = np.linalg.cholesky(A).T
        diag = np.diag(R)
        return float(np.sum(np.log(np.maximum(np.abs(diag), np.finfo(float).tiny))))
    except np.linalg.LinAlgError:
        sign, logdet = np.linalg.slogdet(A)
        if sign <= 0 or not math.isfinite(float(logdet)):
            return float("-inf")
        return 0.5 * float(logdet)


def eta_from_R(R: npt.ArrayLike) -> float:
    R_arr = np.asarray(R, dtype=float)
    if R_arr.ndim != 2 or R_arr.shape[0] != R_arr.shape[1]:
        raise ValueError("R must be a square triangular factor")
    diag = np.diag(R_arr)
    return float(np.sum(np.log(np.maximum(np.abs(diag), np.finfo(float).tiny))))


def eta_from_factor(factor: object) -> float:
    """Best-effort eta extraction from a scikit-sparse CHOLMOD factor.

    For LDL factors, det(H)=prod(D). For LL factors, CHOLMOD exposes a logdet
    method in some versions. The fallbacks deliberately keep the public API small
    and robust across scikit-sparse versions.
    """
    if hasattr(factor, "logdet"):
        try:
            return 0.5 * float(factor.logdet())
        except Exception:
            pass
    if hasattr(factor, "D"):
        try:
            D = np.asarray(factor.D(), dtype=float)
            if D.ndim == 2:
                D = np.diag(D)
            D = D.reshape(-1)
            return 0.5 * float(np.sum(np.log(np.maximum(np.abs(D), np.finfo(float).tiny))))
        except Exception:
            pass
    if hasattr(factor, "L"):
        try:
            L = factor.L()
            diag = L.diagonal() if sp.issparse(L) else np.diag(np.asarray(L, dtype=float))
            return float(np.sum(np.log(np.maximum(np.abs(diag), np.finfo(float).tiny))))
        except Exception:
            pass
    raise ValueError("could not extract log-determinant information from factor")


def eta_from_system(*, factor: object | None = None, R: npt.ArrayLike | None = None, H: sp.spmatrix | npt.ArrayLike | None = None) -> float:
    if factor is not None:
        try:
            return eta_from_factor(factor)
        except Exception:
            if H is None:
                raise
    if R is not None:
        return eta_from_R(R)
    if H is not None:
        return _dense_eta_from_spd(H)
    raise ValueError("provide factor, R, or H")


def detrended_delta_eta(eta_t: float, eta_prev: float, n_prev: int, n_t: int) -> float:
    if n_t <= 0:
        return float(eta_t)
    if n_prev <= 0:
        return float(eta_t)
    return float(eta_t - (float(n_prev) / float(n_t)) * float(eta_prev))


def compute_information_gain(eta_t: float, eta_prev: float, n_prev: int, n_t: int, threshold: float) -> InformationGain:
    return InformationGain(eta=float(eta_t), delta_eta=detrended_delta_eta(eta_t, eta_prev, n_prev, n_t), threshold=float(threshold))
