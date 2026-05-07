from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from scipy.linalg import cholesky as dense_cholesky
from scipy.sparse.csgraph import reverse_cuthill_mckee
from incremental_slam.linearization.edge_jacobian import jacobian_edge_jr
from incremental_slam.types import Graph
@dataclass
class LinearizeAffectedResult:
    factor: object | None
    R: npt.NDArray[np.float64] | None
    A_perm: sp.csc_matrix
    b: npt.NDArray[np.float64]
    J: sp.csc_matrix
    r: npt.NDArray[np.float64]
    p: npt.NDArray[np.int64]
def _ensure_csc(J): return J.tocsc() if sp.issparse(J) else sp.csc_matrix(np.asarray(J, dtype=float))
def _ensure_vec(r): return np.asarray(r, dtype=float).reshape(-1)
def _colamd_like(J: sp.csc_matrix) -> npt.NDArray[np.int64]:
    try:
        from sksparse.ccolamd import ccolamd
        return np.asarray(ccolamd(J), dtype=np.int64).reshape(-1)
    except Exception:
        try:
            from sksparse.colamd import colamd
            return np.asarray(colamd(J), dtype=np.int64).reshape(-1)
        except Exception:
            return np.asarray(reverse_cuthill_mckee((J.T @ J).tocsc(), symmetric_mode=True), dtype=np.int64).reshape(-1)
def _factorize_spd(A_perm):
    try:
        from sksparse.cholmod import cholesky
        return cholesky(A_perm, ordering_method='natural'), None
    except Exception:
        return None, np.asarray(dense_cholesky(A_perm.toarray(), lower=False, check_finite=False), dtype=float)
def _starts(edge_j_index): return sorted(int(v) for v in edge_j_index.values())
def _infer(row_start, edge_j_index, total_rows):
    starts = _starts(edge_j_index); next_candidates = [s for s in starts if s > row_start]; next_start = next_candidates[0] if next_candidates else total_rows; stored = next_start - row_start
    if stored <= 0: raise ValueError('invalid stored row layout')
    return stored
def linearize_affected(J, r, edge_j_index, p, g: Graph, affected_edge_ids):
    J_csc = _ensure_csc(J); r_vec = _ensure_vec(r); affected = np.asarray(affected_edge_ids, dtype=np.int64).reshape(-1) if affected_edge_ids is not None else np.empty(0, dtype=np.int64); J_lil = J_csc.tolil(copy=True)
    for eid in affected:
        edge = g.edges[int(eid)]; Je, re = jacobian_edge_jr(edge, g, g.state_size); m = Je.shape[0]; row_start = int(edge_j_index[int(eid)]); row_stop = row_start + m
        if _infer(row_start, edge_j_index, J_lil.shape[0]) != m: raise ValueError('row block size mismatch')
        J_lil[row_start:row_stop, :] = Je.tolil(); r_vec[row_start:row_stop] = np.asarray(re, dtype=float).reshape(-1)
    J_csc = J_lil.tocsc(); p_vec = _colamd_like(J_csc) if p is None or np.asarray(p).size == 0 else np.asarray(p, dtype=np.int64).reshape(-1); Jp = J_csc[:, p_vec].tocsc(); A_perm = (Jp.T @ Jp).tocsc(); b = np.asarray(Jp.T @ r_vec, dtype=float).reshape(-1); factor, R = _factorize_spd(A_perm)
    return LinearizeAffectedResult(factor=factor, R=R, A_perm=A_perm, b=b, J=J_csc, r=r_vec, p=p_vec)
