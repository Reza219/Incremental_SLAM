import numpy as np
import scipy.sparse as sp

from incremental_slam.solvers.full import solve_full
from incremental_slam.solvers.selective import solve_affected


def _build_spd_system(n: int = 6):
    H = np.zeros((n, n), dtype=float)
    for i in range(n):
        H[i, i] = 4.0
        if i + 1 < n:
            H[i, i + 1] = H[i + 1, i] = -1.0
    H[0, 0] += 1.0
    b = np.linspace(1.0, float(n), n)
    p = np.array([2, 0, 1, 5, 3, 4], dtype=np.int64)[:n]
    H_perm = H[np.ix_(p, p)]
    b_perm = b[p]
    R = np.linalg.cholesky(H_perm).T
    return H, H_perm, b, b_perm, p, R


def test_selective_solver_matches_reduced_subsystem_solution():
    H, H_perm, b, b_perm, p, _ = _build_spd_system()
    affected_vars = np.array([0, 1], dtype=np.int64)
    res = solve_affected(b_perm=b_perm, p=p, affected_vars=affected_vars, H_perm=sp.csc_matrix(H_perm), alpha=1.0)

    active = res.active_perm
    H_sub = H_perm[np.ix_(active, active)]
    x_sub = np.linalg.solve(H_sub, b_perm[active])
    x_perm = np.zeros(len(p), dtype=float)
    x_perm[active] = x_sub
    x = np.zeros(len(p), dtype=float)
    x[p] = x_perm
    x_masked = np.zeros_like(x)
    x_masked[affected_vars] = x[affected_vars]

    assert not res.used_full_fallback
    assert np.allclose(res.dx, x_masked)


def test_selective_solver_matches_legacy_R_subsystem_solution():
    _, H_perm, _, b_perm, p, R = _build_spd_system()
    affected_vars = np.array([3, 4], dtype=np.int64)
    res = solve_affected(b_perm=b_perm, p=p, affected_vars=affected_vars, R=R, alpha=1.0)
    active = res.active_perm
    H_sub = H_perm[np.ix_(active, active)]
    x_sub = np.linalg.solve(H_sub, b_perm[active])
    x_perm = np.zeros(len(p), dtype=float)
    x_perm[active] = x_sub
    x = np.zeros(len(p), dtype=float)
    x[p] = x_perm
    x_masked = np.zeros_like(x)
    x_masked[affected_vars] = x[affected_vars]

    assert not res.used_full_fallback
    assert np.allclose(res.dx, x_masked)


def test_selective_solver_full_fallback_matches_masked_full_solution():
    _, H_perm, _, b_perm, p, R = _build_spd_system()
    affected_vars = np.array([0, 5], dtype=np.int64)
    res = solve_affected(b_perm=b_perm, p=p, affected_vars=affected_vars, H_perm=H_perm, R=R, alpha=0.01)
    full = solve_full(b_perm=b_perm, p=p, R=R)
    masked = np.zeros_like(full)
    masked[affected_vars] = full[affected_vars]

    assert res.used_full_fallback
    assert np.allclose(res.dx, masked)
