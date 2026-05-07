from __future__ import annotations
import numpy as np
import numpy.typing as npt
from scipy.linalg import solve_triangular

def solve_full(*, b_perm: npt.ArrayLike, p: npt.ArrayLike, factor: object | None = None, R: npt.NDArray[np.float64] | None = None) -> npt.NDArray[np.float64]:
    b = np.asarray(b_perm, dtype=float).reshape(-1); p = np.asarray(p, dtype=np.int64).reshape(-1); n = p.size
    if (factor is None) == (R is None): raise ValueError('provide exactly one of factor or R')
    if factor is not None:
        dx_p = np.asarray(factor.solve_A(b), dtype=float).reshape(-1)
    else:
        y = solve_triangular(np.asarray(R).T, b, lower=True, check_finite=False); dx_p = solve_triangular(np.asarray(R), y, lower=False, check_finite=False)
    dx = np.zeros(n, dtype=float); dx[p] = dx_p; return dx
def converged(dx: npt.ArrayLike, dx_th: float) -> bool:
    dx = np.asarray(dx, dtype=float).reshape(-1); return True if dx.size == 0 else float(np.max(np.abs(dx))) <= float(dx_th)
