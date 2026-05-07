from __future__ import annotations
import numpy as np
import numpy.typing as npt

def linearize_pose_landmark(x: npt.NDArray[np.float64], l: npt.NDArray[np.float64], z: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    x = np.asarray(x, dtype=float).reshape(-1); l = np.asarray(l, dtype=float).reshape(-1); z = np.asarray(z, dtype=float).reshape(-1)
    xi, yi, th = x; lx, ly = l; c = np.cos(th); s = np.sin(th)
    Rt = np.array([[c, s], [-s, c]], dtype=float)
    e = Rt @ (np.array([lx, ly], dtype=float) - np.array([xi, yi], dtype=float)) - z
    A = np.array([[-c, -s, (xi - lx) * s + (ly - yi) * c], [s, -c, (xi - lx) * c + (yi - ly) * s]], dtype=float)
    B = np.array([[c, s], [-s, c]], dtype=float)
    return e, A, B
