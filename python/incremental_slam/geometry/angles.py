from __future__ import annotations
import numpy as np
import numpy.typing as npt

def normalize_angle(phi: float | npt.NDArray[np.float64]) -> float | npt.NDArray[np.float64]:
    arr = np.asarray(phi, dtype=float)
    wrapped = np.mod(arr + np.pi, 2.0 * np.pi) - np.pi
    return float(wrapped) if np.isscalar(phi) else wrapped
