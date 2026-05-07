from __future__ import annotations
import numpy as np
import numpy.typing as npt
from incremental_slam.geometry.angles import normalize_angle

def linearize_pose_pose(x1: npt.NDArray[np.float64], x2: npt.NDArray[np.float64], z: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    x1 = np.asarray(x1, dtype=float).reshape(-1); x2 = np.asarray(x2, dtype=float).reshape(-1); z = np.asarray(z, dtype=float).reshape(-1)
    x_1, y_1, th1 = x1; x_2, y_2, th2 = x2; cz = np.cos(z[2]); sz = np.sin(z[2]); c1 = np.cos(th1); s1 = np.sin(th1)
    R1 = np.array([[c1, -s1], [s1, c1]], dtype=float); Rz = np.array([[cz, -sz], [sz, cz]], dtype=float)
    delta_t = np.array([x_2 - x_1, y_2 - y_1], dtype=float); pred_xy = R1.T @ delta_t; e_xy = Rz.T @ (pred_xy - z[0:2]); e_th = normalize_angle(th2 - th1 - z[2])
    e = np.array([e_xy[0], e_xy[1], e_th], dtype=float)
    Q1 = np.array([[-c1, -s1, (x_1 - x_2) * s1 + (y_2 - y_1) * c1], [s1, -c1, (x_1 - x_2) * c1 + (y_1 - y_2) * s1]], dtype=float)
    Q2 = np.array([[c1, s1, 0.0], [-s1, c1, 0.0]], dtype=float)
    A = np.vstack([Rz.T @ Q1, np.array([[0.0, 0.0, -1.0]], dtype=float)])
    B = np.vstack([Rz.T @ Q2, np.array([[0.0, 0.0, 1.0]], dtype=float)])
    return e, A, B
