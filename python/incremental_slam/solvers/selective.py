from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import solve_triangular

from incremental_slam.solvers.full import solve_full


@dataclass
class SelectiveSolveResult:
    dx: npt.NDArray[np.float64]
    reach_perm: npt.NDArray[np.int64]
    active_perm: npt.NDArray[np.int64]
    used_full_fallback: bool


def _to_int_vec(x: npt.ArrayLike | None) -> npt.NDArray[np.int64]:
    if x is None:
        return np.empty(0, dtype=np.int64)
    return np.asarray(x, dtype=np.int64).reshape(-1)


def _to_float_vec(x: npt.ArrayLike) -> npt.NDArray[np.float64]:
    return np.asarray(x, dtype=float).reshape(-1)


def _symmetric_etree(A: sp.spmatrix | npt.NDArray[np.float64]) -> npt.NDArray[np.int64]:
    if sp.issparse(A):
        Ac = A.tocsc()
    else:
        Ac = sp.csc_matrix(np.asarray(A, dtype=float))

    n = Ac.shape[0]
    parent = np.full(n, -1, dtype=np.int64)
    ancestor = np.full(n, -1, dtype=np.int64)

    for k in range(n):
        col_rows = Ac.indices[Ac.indptr[k] : Ac.indptr[k + 1]]
        for ii in col_rows:
            i = int(ii)
            if i >= k:
                continue
            while i != -1 and i < k:
                inext = ancestor[i]
                ancestor[i] = k
                if inext == -1:
                    parent[i] = k
                i = int(inext)
    return parent


def _parents_of(nodes: npt.NDArray[np.int64], parent: npt.NDArray[np.int64] | None) -> npt.NDArray[np.int64]:
    if nodes.size == 0:
        return np.empty(0, dtype=np.int64)
    if parent is None or parent.size == 0:
        return np.unique(nodes)
    out: set[int] = set()
    for node in np.asarray(nodes, dtype=np.int64).reshape(-1):
        j = int(node)
        while j >= 0 and j not in out:
            out.add(j)
            j = int(parent[j]) if j < parent.size else -1
    return np.asarray(sorted(out), dtype=np.int64)


def _solve_with_R_subsystem(R: npt.NDArray[np.float64], b: npt.NDArray[np.float64], active_perm: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
    A = R.T @ R
    A_sub = A[np.ix_(active_perm, active_perm)]
    b_sub = b[active_perm]
    try:
        R_sub = np.linalg.cholesky(A_sub).T
        y = solve_triangular(R_sub.T, b_sub, lower=True, check_finite=False)
        return solve_triangular(R_sub, y, lower=False, check_finite=False)
    except np.linalg.LinAlgError:
        return np.linalg.solve(A_sub, b_sub)


def _solve_reduced_normal(H_perm: sp.spmatrix | npt.NDArray[np.float64], b_perm: npt.NDArray[np.float64], active_perm: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
    if sp.issparse(H_perm):
        H = H_perm.tocsc()
    else:
        H = sp.csc_matrix(np.asarray(H_perm, dtype=float))
    H_sub = H[active_perm, :][:, active_perm].tocsc()
    b_sub = b_perm[active_perm]
    try:
        x_sub = spla.spsolve(H_sub, b_sub)
    except Exception:
        x_sub = np.linalg.solve(H_sub.toarray(), b_sub)
    return np.asarray(x_sub, dtype=float).reshape(-1)


def solve_affected(
    *,
    b_perm: npt.ArrayLike,
    p: npt.ArrayLike,
    affected_vars: npt.ArrayLike,
    parent: npt.ArrayLike | None = None,
    factor: object | None = None,
    R: npt.NDArray[np.float64] | None = None,
    H_perm: sp.spmatrix | npt.NDArray[np.float64] | None = None,
    alpha: float = 0.3,
) -> SelectiveSolveResult:
    b = _to_float_vec(b_perm)
    p_vec = _to_int_vec(p)
    affected = _to_int_vec(affected_vars)
    parent_vec = None if parent is None else _to_int_vec(parent)
    n = p_vec.size

    if affected.size == 0:
        return SelectiveSolveResult(dx=np.zeros(n, dtype=float), reach_perm=np.empty(0, dtype=np.int64), active_perm=np.empty(0, dtype=np.int64), used_full_fallback=False)

    invp = np.empty(n, dtype=np.int64)
    invp[p_vec] = np.arange(n, dtype=np.int64)
    affected_perm = np.unique(invp[affected])

    if parent_vec is None and H_perm is not None:
        parent_vec = _symmetric_etree(H_perm)
    reach = _parents_of(affected_perm, parent_vec)

    if reach.size / max(1, n) > float(alpha):
        dx = solve_full(b_perm=b, p=p_vec, factor=factor, R=R)
        dx_masked = np.zeros_like(dx)
        dx_masked[affected] = dx[affected]
        return SelectiveSolveResult(dx=dx_masked, reach_perm=reach, active_perm=np.arange(n, dtype=np.int64), used_full_fallback=True)

    if H_perm is None:
        if R is not None:
            x_sub = _solve_with_R_subsystem(np.asarray(R, dtype=float), b, reach)
        elif factor is not None:
            dx = solve_full(b_perm=b, p=p_vec, factor=factor)
            dx_masked = np.zeros_like(dx)
            dx_masked[affected] = dx[affected]
            return SelectiveSolveResult(dx=dx_masked, reach_perm=reach, active_perm=reach, used_full_fallback=True)
        else:
            raise ValueError("Need at least one of H_perm, R, or factor")
    else:
        x_sub = _solve_reduced_normal(H_perm, b, reach)

    dx_perm = np.zeros(n, dtype=float)
    dx_perm[reach] = x_sub
    dx = np.zeros(n, dtype=float)
    dx[p_vec] = dx_perm
    dx_masked = np.zeros_like(dx)
    dx_masked[affected] = dx[affected]
    return SelectiveSolveResult(dx=dx_masked, reach_perm=reach, active_perm=reach, used_full_fallback=False)
