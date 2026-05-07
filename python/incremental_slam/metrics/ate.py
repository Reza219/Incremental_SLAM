from __future__ import annotations
import numpy as np
from incremental_slam.metrics.alignment import kabsch2d_local
from incremental_slam.types import Graph

def _pose_xy_from_graph(g: Graph):
    pts = []; ids = []
    for node_id, ref in sorted(g.id_lookup.items(), key=lambda kv: kv[1].offset):
        if ref.dimension == 3: pts.append(g.x[ref.offset: ref.offset + 2]); ids.append(node_id)
    return (np.vstack(pts), ids) if pts else (np.empty((0,2), dtype=float), [])
def compute_ate_rmse(g_est: Graph, g_gt: Graph):
    est_xy, est_ids = _pose_xy_from_graph(g_est); gt_xy, gt_ids = _pose_xy_from_graph(g_gt); gt_map = {nid: gt_xy[i] for i, nid in enumerate(gt_ids)}; common_ids = [nid for nid in est_ids if nid in gt_map]
    if not common_ids: return float('nan'), np.empty(0, dtype=float)
    est = np.vstack([est_xy[est_ids.index(nid)] for nid in common_ids]); gt = np.vstack([gt_map[nid] for nid in common_ids]); R, t = kabsch2d_local(est, gt); est_aligned = (R @ est.T).T + t; err = np.linalg.norm(est_aligned - gt, axis=1); return float(np.sqrt(np.mean(err ** 2))), err
