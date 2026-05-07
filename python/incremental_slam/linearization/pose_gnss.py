from __future__ import annotations
import numpy as np
import numpy.typing as npt

def linearize_pose_gnss(xf: npt.NDArray[np.float64], z: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    xf = np.asarray(xf, dtype=float).reshape(-1); z = np.asarray(z, dtype=float).reshape(-1)
    e = np.array([xf[0] - z[0], xf[1] - z[1]], dtype=float)
    A = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float)
    return e, A
