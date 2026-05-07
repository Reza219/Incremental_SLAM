from __future__ import annotations
import numpy as np
import numpy.typing as npt

def v2t(v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    v = np.asarray(v, dtype=float).reshape(-1)
    x, y, theta = v
    c = np.cos(theta); s = np.sin(theta)
    return np.array([[c, -s, x], [s, c, y], [0.0, 0.0, 1.0]], dtype=float)
def t2v(A: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    A = np.asarray(A, dtype=float)
    if A.ndim == 2:
        return np.array([A[0,2], A[1,2], np.arctan2(A[1,0], A[0,0])], dtype=float)
    out = np.zeros((3, A.shape[2]), dtype=float)
    for k in range(A.shape[2]):
        out[:,k] = np.array([A[0,2,k], A[1,2,k], np.arctan2(A[1,0,k], A[0,0,k])], dtype=float)
    return out
def invt(m: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    m = np.asarray(m, dtype=float); R = m[0:2,0:2]; t = m[0:2,2]; Rt = R.T; out = np.eye(3, dtype=float); out[0:2,0:2] = Rt; out[0:2,2] = -Rt @ t; return out
