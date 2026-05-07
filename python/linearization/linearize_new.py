from __future__ import annotations
from dataclasses import dataclass
from typing import MutableMapping
import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from scipy.linalg import cholesky as dense_cholesky
from scipy.sparse.csgraph import reverse_cuthill_mckee
from incremental_slam.linearization.edge_jacobian import jacobian_edge_jr
from incremental_slam.types import Graph
@dataclass
class LinearizeNewResult:
    factor: object | None
    R: npt.NDArray[np.float64] | None
    A_perm: sp.csc_matrix
    b: npt.NDArray[np.float64]
    J: sp.csc_matrix
    r: npt.NDArray[np.float64]
    edge_j_index: MutableMapping[int, int]
    p: npt.NDArray[np.int64]
    parent: npt.NDArray[np.int64] | None

def _as_1d_int_array(edge_ids: npt.ArrayLike | None) -> npt.NDArray[np.int64]: return np.empty(0, dtype=np.int64) if edge_ids is None else np.asarray(edge_ids, dtype=np.int64).reshape(-1)
def _ensure_csc(J): return J.tocsc() if sp.issparse(J) else sp.csc_matrix(np.asarray(J, dtype=float))
def _ensure_vec(r): return np.asarray(r, dtype=float).reshape(-1)
def _expand_jacobian_width(J: sp.csc_matrix, target_cols: int) -> sp.csc_matrix: return J if J.shape[1] >= target_cols else sp.hstack([J, sp.csc_matrix((J.shape[0], target_cols - J.shape[1]), dtype=J.dtype)], format='csc')
def _stack_edge_blocks(J, r, edge_j_index, g, edge_ids, g_full=None):
    """Append newly added *local* edge blocks.

    edge_ids must refer to indices in g.edges, not to the original full-graph edge
    array. This matters whenever the current graph's state layout differs from the
    source graph's original offsets.
    """
    next_row = J.shape[0]; J_blocks = []; r_blocks = []
    for eid in edge_ids:
        local_eid = int(eid)
        Je, re = jacobian_edge_jr(g.edges[local_eid], g, state_size=g.state_size); edge_j_index[local_eid] = next_row; next_row += Je.shape[0]; J_blocks.append(Je); r_blocks.append(np.asarray(re, dtype=float).reshape(-1))
    if J_blocks:
        J_new = sp.vstack(J_blocks, format='csc'); r_new = np.concatenate(r_blocks) if r_blocks else np.empty(0, dtype=float); J = sp.vstack([J, J_new], format='csc'); r = np.concatenate([r, r_new]) if r.size else r_new
    return J, r, edge_j_index
def _symamd_like(A: sp.csc_matrix) -> npt.NDArray[np.int64]:
    try:
        from sksparse.ccolamd import csymamd
        p = np.asarray(csymamd(A), dtype=np.int64).reshape(-1)
    except Exception:
        try:
            from sksparse.colamd import symamd
            p = np.asarray(symamd(A), dtype=np.int64).reshape(-1)
        except Exception:
            p = np.asarray(reverse_cuthill_mckee(A, symmetric_mode=True), dtype=np.int64).reshape(-1)
    return p
def _etree_parent_dense_from_R(R, tol: float = 1e-14):
    n = R.shape[0]; parent = -np.ones(n, dtype=np.int64)
    for j in range(n):
        rows = np.flatnonzero(np.abs(R[:j, j]) > tol)
        if rows.size: parent[j] = int(rows[-1])
    return parent
def _factorize_spd(A_perm):
    try:
        from sksparse.cholmod import cholesky
        return cholesky(A_perm, ordering_method='natural'), None, None
    except Exception:
        R = dense_cholesky(A_perm.toarray(), lower=False, check_finite=False)
        return None, np.asarray(R, dtype=float), _etree_parent_dense_from_R(R)
def linearize_new(J, r, edge_j_index, g: Graph, edge_ids, g_full):
    J_csc = _ensure_csc(J); r_vec = _ensure_vec(r); edge_ids_vec = _as_1d_int_array(edge_ids); J_csc = _expand_jacobian_width(J_csc, g.state_size)
    if edge_ids_vec.size: J_csc, r_vec, edge_j_index = _stack_edge_blocks(J_csc, r_vec, edge_j_index, g, edge_ids_vec, g_full)
    A = (J_csc.T @ J_csc).tocsc(); p = _symamd_like(A); Jp = J_csc[:, p].tocsc(); A_perm = (Jp.T @ Jp).tocsc(); b = np.asarray(Jp.T @ r_vec, dtype=float).reshape(-1); factor, R, parent = _factorize_spd(A_perm)
    return LinearizeNewResult(factor=factor, R=R, A_perm=A_perm, b=b, J=J_csc, r=r_vec, edge_j_index=edge_j_index, p=p, parent=parent)
