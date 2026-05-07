from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from incremental_slam.metrics.information import eta_from_factor, eta_from_system

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

try:
    from sksparse.cholmod import ldl_factor  # type: ignore
    _HAVE_CHOLMOD = True
except Exception:  # pragma: no cover - optional dependency
    ldl_factor = None  # type: ignore
    _HAVE_CHOLMOD = False

CSC = sp.csc_array | sp.csc_matrix


def cholmod_available() -> bool:
    return bool(_HAVE_CHOLMOD)


@dataclass
class EdgeContribution:
    C: CSC
    r: npt.NDArray[np.float64]
    touched_vars: npt.NDArray[np.int64] = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    touched_nodes: npt.NDArray[np.int64] = field(default_factory=lambda: np.empty(0, dtype=np.int64))


def _to_csc(A: sp.spmatrix | npt.ArrayLike) -> CSC:
    if sp.issparse(A):
        return A.tocsc()
    return sp.csc_array(np.asarray(A, dtype=float))


def _to_vec(x: npt.ArrayLike) -> npt.NDArray[np.float64]:
    return np.asarray(x, dtype=float).reshape(-1)


def _empty_A(n_state: int) -> CSC:
    return sp.csc_array((n_state, 0), dtype=float)


def _hstack(blocks: list[CSC], n_state: int) -> CSC:
    if not blocks:
        return _empty_A(n_state)
    return sp.hstack(blocks, format="csc")


def _factor_spd(H: CSC, order: str):
    """Compatibility wrapper across scikit-sparse versions."""
    if not _HAVE_CHOLMOD:
        raise RuntimeError("CHOLMOD backend is unavailable. Install scikit-sparse and SuiteSparse.")
    H = _to_csc(H)
    for kwargs in ({"order": order}, {"ordering_method": order}, {}):
        try:
            return ldl_factor(H, **kwargs)
        except TypeError:
            continue
    return ldl_factor(H)


@dataclass
class CholmodNormalManager:
    """Maintains an LDL factorization of H = J.T @ J via A A.T with A = J.T."""

    state_size: int
    factor: object | None = None
    g: npt.NDArray[np.float64] = field(default_factory=lambda: np.empty(0, dtype=float))
    edge_cache: dict[int, EdgeContribution] = field(default_factory=dict)
    _order: str = "amd"

    def __post_init__(self) -> None:
        self.g = np.zeros(self.state_size, dtype=float)

    @property
    def num_cached_edges(self) -> int:
        return len(self.edge_cache)

    def clear(self) -> None:
        self.factor = None
        self.g = np.zeros(self.state_size, dtype=float)
        self.edge_cache.clear()

    def rebuild(self, state_size: int, edge_blocks: Iterable[tuple[int, CSC, npt.ArrayLike, npt.ArrayLike | None, npt.ArrayLike | None]], order: str = "amd") -> None:
        if not _HAVE_CHOLMOD:
            raise RuntimeError("CHOLMOD backend is unavailable. Install scikit-sparse and SuiteSparse.")
        self.state_size = int(state_size)
        self._order = str(order)
        self.g = np.zeros(self.state_size, dtype=float)
        self.edge_cache.clear()
        C_blocks: list[CSC] = []
        for eid, C, r, touched_vars, touched_nodes in edge_blocks:
            C = _to_csc(C)
            r = _to_vec(r)
            if C.shape[0] != self.state_size:
                raise ValueError(f"edge {eid}: C has {C.shape[0]} rows, expected {self.state_size}")
            if C.shape[1] != r.size:
                raise ValueError(f"edge {eid}: C columns {C.shape[1]} != residual size {r.size}")
            self.edge_cache[int(eid)] = EdgeContribution(C=C, r=r, touched_vars=np.asarray([] if touched_vars is None else touched_vars, dtype=np.int64).reshape(-1), touched_nodes=np.asarray([] if touched_nodes is None else touched_nodes, dtype=np.int64).reshape(-1))
            self.g += np.asarray(C @ r, dtype=float).reshape(-1)
            C_blocks.append(C)
        A = _hstack(C_blocks, self.state_size)
        H = (A @ A.T).tocsc()
        self.factor = _factor_spd(H, self._order)

    def add_edge(self, eid: int, C: CSC, r: npt.ArrayLike, touched_vars: npt.ArrayLike | None = None, touched_nodes: npt.ArrayLike | None = None) -> None:
        if not _HAVE_CHOLMOD:
            raise RuntimeError("CHOLMOD backend is unavailable. Install scikit-sparse and SuiteSparse.")
        if self.factor is None:
            raise RuntimeError("factor is not initialized; call rebuild() first")
        if eid in self.edge_cache:
            raise KeyError(f"edge {eid} already exists in cache")
        C = _to_csc(C)
        r = _to_vec(r)
        if C.shape[0] != self.state_size:
            raise ValueError(f"edge {eid}: C has {C.shape[0]} rows, expected {self.state_size}")
        if C.shape[1] != r.size:
            raise ValueError(f"edge {eid}: C columns {C.shape[1]} != residual size {r.size}")
        self.factor.update(C)
        self.g += np.asarray(C @ r, dtype=float).reshape(-1)
        self.edge_cache[int(eid)] = EdgeContribution(C=C, r=r, touched_vars=np.asarray([] if touched_vars is None else touched_vars, dtype=np.int64).reshape(-1), touched_nodes=np.asarray([] if touched_nodes is None else touched_nodes, dtype=np.int64).reshape(-1))

    def remove_edge(self, eid: int) -> None:
        if self.factor is None:
            raise RuntimeError("factor is not initialized; call rebuild() first")
        if eid not in self.edge_cache:
            raise KeyError(f"edge {eid} not found in cache")
        old = self.edge_cache.pop(int(eid))
        self.factor.downdate(old.C)
        self.g -= np.asarray(old.C @ old.r, dtype=float).reshape(-1)

    def replace_edge(self, eid: int, C_new: CSC, r_new: npt.ArrayLike, touched_vars: npt.ArrayLike | None = None, touched_nodes: npt.ArrayLike | None = None) -> None:
        if self.factor is None:
            raise RuntimeError("factor is not initialized; call rebuild() first")
        if eid not in self.edge_cache:
            raise KeyError(f"edge {eid} not found in cache")
        old = self.edge_cache[int(eid)]
        self.factor.downdate(old.C)
        self.g -= np.asarray(old.C @ old.r, dtype=float).reshape(-1)
        C_new = _to_csc(C_new)
        r_new = _to_vec(r_new)
        if C_new.shape[0] != self.state_size:
            raise ValueError(f"edge {eid}: C_new has {C_new.shape[0]} rows, expected {self.state_size}")
        if C_new.shape[1] != r_new.size:
            raise ValueError(f"edge {eid}: C_new columns {C_new.shape[1]} != residual size {r_new.size}")
        self.factor.update(C_new)
        self.g += np.asarray(C_new @ r_new, dtype=float).reshape(-1)
        self.edge_cache[int(eid)] = EdgeContribution(C=C_new, r=r_new, touched_vars=np.asarray([] if touched_vars is None else touched_vars, dtype=np.int64).reshape(-1), touched_nodes=np.asarray([] if touched_nodes is None else touched_nodes, dtype=np.int64).reshape(-1))

    def solve(self, rhs: npt.ArrayLike | None = None) -> npt.NDArray[np.float64]:
        if self.factor is None:
            raise RuntimeError("factor is not initialized; call rebuild() first")
        b = self.g if rhs is None else _to_vec(rhs)
        if b.size != self.state_size:
            raise ValueError(f"rhs length {b.size} does not match state size {self.state_size}")
        return np.asarray(self.factor.solve(b), dtype=float).reshape(-1)

    def get_perm(self) -> npt.NDArray[np.int64]:
        if self.factor is None:
            raise RuntimeError("factor is not initialized; call rebuild() first")
        return np.asarray(self.factor.get_perm(), dtype=np.int64).reshape(-1)

    def build_A(self) -> CSC:
        blocks = [self.edge_cache[eid].C for eid in sorted(self.edge_cache)]
        return _hstack(blocks, self.state_size)

    def build_normal_matrix(self) -> CSC:
        A = self.build_A()
        return (A @ A.T).tocsc()

    def information_eta(self) -> float:
        """Return eta = 0.5 log det(J.T J) for the current normal system."""
        if self.factor is None:
            raise RuntimeError("factor is not initialized; call rebuild() first")
        try:
            return eta_from_factor(self.factor)
        except Exception:
            return eta_from_system(H=self.build_normal_matrix())

    def build_permuted_system(self) -> tuple[CSC, npt.NDArray[np.float64], npt.NDArray[np.int64]]:
        H = self.build_normal_matrix()
        p = self.get_perm()
        H_perm = H[p, :][:, p].tocsc()
        b_perm = self.g[p].copy()
        return H_perm, b_perm, p
