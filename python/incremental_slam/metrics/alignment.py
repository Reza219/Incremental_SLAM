from __future__ import annotations
import numpy as np
import numpy.typing as npt

def kabsch2d_local(est_xy: npt.NDArray[np.float64], gt_xy: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    est = np.asarray(est_xy, dtype=float); gt = np.asarray(gt_xy, dtype=float); est_mu = est.mean(axis=0); gt_mu = gt.mean(axis=0); X = est - est_mu; Y = gt - gt_mu; H = X.T @ Y; U, _, Vt = np.linalg.svd(H); R = Vt.T @ U.T
    if np.linalg.det(R) < 0: Vt[-1,:] *= -1; R = Vt.T @ U.T
    t = gt_mu - (R @ est_mu)
    return R, t
