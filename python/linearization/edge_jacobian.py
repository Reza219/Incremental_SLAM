from __future__ import annotations
import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from scipy.linalg import sqrtm
from incremental_slam.linearization.pose_gnss import linearize_pose_gnss
from incremental_slam.linearization.pose_landmark import linearize_pose_landmark
from incremental_slam.linearization.pose_pose import linearize_pose_pose
from incremental_slam.types import Edge, Graph
DIM_POSE = 3; DIM_LANDMARK = 2
def _symmetrize(mat: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]: return 0.5 * (mat + mat.T)
def _whitening_factor(info: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]: return np.asarray(np.real_if_close(sqrtm(_symmetrize(np.asarray(info, dtype=float))), tol=1000), dtype=float)
def jacobian_edge_jr(edge: Edge, g: Graph, state_size: int | None = None) -> tuple[sp.csc_matrix, npt.NDArray[np.float64]]:
    state_size = g.state_size if state_size is None else state_size
    if edge.type == 'P':
        xf = g.block(edge.from_idx, DIM_POSE); xt = g.block(edge.to_idx, DIM_POSE); e, A, B = linearize_pose_pose(xf, xt, edge.measurement); dim_from, dim_to = DIM_POSE, DIM_POSE
    elif edge.type == 'L':
        xf = g.block(edge.from_idx, DIM_POSE); xt = g.block(edge.to_idx, DIM_LANDMARK); e, A, B = linearize_pose_landmark(xf, xt, edge.measurement); dim_from, dim_to = DIM_POSE, DIM_LANDMARK
    elif edge.type == 'G':
        xf = g.block(edge.from_idx, DIM_POSE); e, A = linearize_pose_gnss(xf, edge.measurement); B = None; dim_from, dim_to = DIM_POSE, 0
    else:
        raise ValueError(f'unsupported edge type {edge.type!r}')
    L = _whitening_factor(edge.information); m = int(np.asarray(e).size); Je = sp.lil_matrix((m, state_size), dtype=float); Je[:, edge.from_idx: edge.from_idx + dim_from] = L @ np.asarray(A, dtype=float)
    if dim_to > 0: Je[:, edge.to_idx: edge.to_idx + dim_to] = L @ np.asarray(B, dtype=float)
    return Je.tocsc(), L @ np.asarray(e, dtype=float).reshape(-1)
