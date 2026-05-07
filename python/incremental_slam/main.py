from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.io import loadmat

from incremental_slam.config import RunConfig, DatasetSpec, dataset_spec
from incremental_slam.graph.affected import (
    all_vars,
    find_affected,
    initial_active_vars_from_new_measurements,
)
from incremental_slam.graph.loop_closure import LoopClosureState
from incremental_slam.graph.reorder import reorder_edges
from incremental_slam.graph.update import update_graph
from incremental_slam.io.g2o import read_graph_g2o
from incremental_slam.io.matlab_graph import read_graph_mat
from incremental_slam.io.toro import read_graph_toro
from incremental_slam.linearization.edge_blocks import linearize_all_current_blocks, linearize_edge_block
from incremental_slam.linearization.linearize_affected import linearize_affected
from incremental_slam.linearization.linearize_new import LinearizeNewResult, linearize_new
from incremental_slam.metrics.ate import compute_ate_rmse
from incremental_slam.metrics.error import compute_global_error
from incremental_slam.metrics.information import compute_information_gain, eta_from_system
from incremental_slam.solvers.cholmod_manager import CholmodNormalManager, cholmod_available
from incremental_slam.solvers.full import converged, solve_full
from incremental_slam.solvers.selective import solve_affected
from incremental_slam.types import Graph


@dataclass
class RunHistory:
    gn_iters: list[int] = field(default_factory=list)
    loop_closure_flags: list[bool] = field(default_factory=list)
    global_update_flags: list[bool] = field(default_factory=list)
    num_nodes: list[int] = field(default_factory=list)
    num_edges: list[int] = field(default_factory=list)
    active_set_size: list[int] = field(default_factory=list)
    dx_inf_norm: list[float] = field(default_factory=list)
    global_error: list[float] = field(default_factory=list)
    ate_rmse: list[float] = field(default_factory=list)
    information_eta: list[float] = field(default_factory=list)
    information_delta_eta: list[float] = field(default_factory=list)


def load_graph_auto(path: Path) -> Graph:
    ext = path.suffix.lower()
    if ext == ".g2o":
        return read_graph_g2o(path)
    if ext == ".graph":
        return read_graph_toro(path)
    if ext == ".mat":
        return read_graph_mat(path)
    raise ValueError(f"unsupported dataset format: {ext}")


def load_ground_truth_mat(path: Path) -> Graph | None:
    if not path.exists():
        return None
    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    return read_graph_mat(path) if "g" in data else None


def _batch_edge_ids(num_edges: int, batch_size: int):
    for start in range(0, num_edges, batch_size):
        yield np.arange(start, min(start + batch_size, num_edges), dtype=np.int64)


def _local_edge_ids_for_latest_batch(g_current: Graph, edge_ids: np.ndarray) -> np.ndarray:
    """Return local edge ids just appended to g_current for this batch."""
    n_new = int(np.asarray(edge_ids, dtype=np.int64).size)
    start = len(g_current.edges) - n_new
    return np.arange(start, len(g_current.edges), dtype=np.int64)

def _seed_anchor():
    return sp.eye(3, format="csc", dtype=float), np.zeros(3, dtype=float)


def _finalize_hist(
    hist: RunHistory,
    g_current: Graph,
    loop_closure: bool,
    gn_used: int,
    last_dx_norm: float,
    cur_gerr: float,
    gt: Graph | None,
    *,
    eta: float = float("nan"),
    delta_eta: float = float("nan"),
    global_update: bool = False,
    active_set_size: int = 0,
) -> None:
    hist.gn_iters.append(int(gn_used))
    hist.loop_closure_flags.append(bool(loop_closure))
    hist.global_update_flags.append(bool(global_update))
    hist.num_nodes.append(len(g_current.id_lookup))
    hist.num_edges.append(len(g_current.edges))
    hist.active_set_size.append(int(active_set_size))
    hist.dx_inf_norm.append(float(last_dx_norm))
    hist.global_error.append(float(cur_gerr))
    hist.ate_rmse.append(float("nan") if gt is None else compute_ate_rmse(g_current, gt)[0])
    hist.information_eta.append(float(eta))
    hist.information_delta_eta.append(float(delta_eta))


def _prefer_cholmod(cfg: RunConfig) -> bool:
    backend = getattr(cfg, "linear_backend", "auto")
    if backend == "cholmod":
        return True
    if backend == "legacy":
        return False
    return cholmod_available()


def _build_or_extend_nm(
    nm: CholmodNormalManager | None,
    g_current: Graph,
    edge_ids: np.ndarray,
    new_vars: np.ndarray,
    order: str,
) -> CholmodNormalManager:
    if nm is None or new_vars.size > 0 or nm.state_size != g_current.state_size:
        nm = CholmodNormalManager(state_size=g_current.state_size)
        nm.rebuild(g_current.state_size, linearize_all_current_blocks(g_current), order=order)
        return nm

    start_local = len(g_current.edges) - len(edge_ids)
    for local_eid in range(start_local, len(g_current.edges)):
        edge = g_current.edges[local_eid]
        C, re, touched_vars, touched_nodes = linearize_edge_block(edge, g_current)
        if local_eid in nm.edge_cache:
            nm.replace_edge(local_eid, C, re, touched_vars, touched_nodes)
        else:
            nm.add_edge(local_eid, C, re, touched_vars, touched_nodes)
    return nm


def _refresh_nm_edges(nm: CholmodNormalManager, g_current: Graph, affected_edge_ids: np.ndarray | None = None) -> None:
    if affected_edge_ids is None:
        edge_iter = range(len(g_current.edges))
    else:
        edge_iter = [int(eid) for eid in np.asarray(affected_edge_ids, dtype=np.int64).reshape(-1)]
    for local_eid in edge_iter:
        edge = g_current.edges[local_eid]
        C, re, touched_vars, touched_nodes = linearize_edge_block(edge, g_current)
        if local_eid in nm.edge_cache:
            nm.replace_edge(local_eid, C, re, touched_vars, touched_nodes)
        else:
            nm.add_edge(local_eid, C, re, touched_vars, touched_nodes)


def _eta_from_linearization(lin: LinearizeNewResult) -> float:
    return eta_from_system(factor=lin.factor, R=lin.R, H=lin.A_perm)


def _active_set_from_gate(
    *,
    cfg: RunConfig,
    spec: DatasetSpec,
    g_current: Graph,
    edge_nodes: np.ndarray,
    loop_closure: bool,
    eta: float,
    prev_eta: float,
    prev_state_size: int,
) -> tuple[np.ndarray, bool, float]:
    threshold = float(spec.ent_th if spec.ent_th is not None else 0.0)
    info = compute_information_gain(eta, prev_eta, prev_state_size, g_current.state_size, threshold)
    mode = getattr(cfg, "gating_mode", "igg")

    if mode == "always":
        global_update = True
    elif mode == "lcg":
        global_update = bool(loop_closure)
    elif mode == "igg":
        global_update = info.triggers_global_update
    else:
        raise ValueError(f"unsupported gating_mode={mode!r}")

    active_vars = all_vars(g_current) if global_update else initial_active_vars_from_new_measurements(g_current, edge_nodes)
    return active_vars, global_update, info.delta_eta


def run_full_demo(cfg: RunConfig):
    if _prefer_cholmod(cfg):
        return run_full_demo_cholmod(cfg)
    spec = dataset_spec(cfg)
    g_full = load_graph_auto(spec.data_file)
    g_full = reorder_edges(g_full) if cfg.use_reorder_edges else g_full
    gt = load_ground_truth_mat(spec.gt_file) if spec.gt_file is not None else None
    g_current = None
    lc_state = LoopClosureState()
    J, r = _seed_anchor()
    edge_j_index = {}
    hist = RunHistory()
    for edge_ids in _batch_edge_ids(len(g_full.edges), cfg.batch_size):
        g_current, _, _, _, loop_closure, lc_state = update_graph(edge_ids, g_full, g_current, spec.lc_gap, lc_state)
        local_edge_ids = _local_edge_ids_for_latest_batch(g_current, edge_ids)
        lin = linearize_new(J, r, edge_j_index, g_current, local_edge_ids, g_full)
        J, r = lin.J, lin.r
        edge_j_index = dict(lin.edge_j_index)
        gn_used = 0
        last_dx_norm = 0.0
        for k in range(1, cfg.max_gn_iter + 1):
            dx = solve_full(b_perm=lin.b, p=lin.p, factor=lin.factor, R=lin.R)
            last_dx_norm = float(np.max(np.abs(dx))) if dx.size else 0.0
            if converged(dx, spec.dx_th):
                gn_used = k - 1
                break
            g_current.x = g_current.x - dx
            all_edge_ids = np.arange(len(g_current.edges), dtype=np.int64)
            aff = linearize_affected(J, r, edge_j_index, lin.p, g_current, all_edge_ids)
            J, r = aff.J, aff.r
            lin = LinearizeNewResult(
                factor=aff.factor,
                R=aff.R,
                A_perm=aff.A_perm,
                b=aff.b,
                J=aff.J,
                r=aff.r,
                edge_j_index=edge_j_index,
                p=aff.p,
                parent=None,
            )
            gn_used = k
        _finalize_hist(hist, g_current, loop_closure, gn_used, last_dx_norm, compute_global_error(g_current), gt, global_update=True, active_set_size=g_current.state_size)
    return g_current, hist


def run_full_demo_cholmod(cfg: RunConfig):
    if not cholmod_available():
        raise RuntimeError("CHOLMOD backend requested, but scikit-sparse/CHOLMOD is unavailable.")
    spec = dataset_spec(cfg)
    g_full = load_graph_auto(spec.data_file)
    g_full = reorder_edges(g_full) if cfg.use_reorder_edges else g_full
    gt = load_ground_truth_mat(spec.gt_file) if spec.gt_file is not None else None
    g_current = None
    lc_state = LoopClosureState()
    hist = RunHistory()
    nm: CholmodNormalManager | None = None
    for edge_ids in _batch_edge_ids(len(g_full.edges), cfg.batch_size):
        g_current, _, _, new_vars, loop_closure, lc_state = update_graph(edge_ids, g_full, g_current, spec.lc_gap, lc_state)
        nm = _build_or_extend_nm(nm, g_current, edge_ids, new_vars, getattr(cfg, "cholmod_order", "amd"))
        gn_used = 0
        last_dx_norm = 0.0
        for k in range(1, cfg.max_gn_iter + 1):
            dx = nm.solve()
            last_dx_norm = float(np.max(np.abs(dx))) if dx.size else 0.0
            if converged(dx, spec.dx_th):
                gn_used = k - 1
                break
            g_current.x = g_current.x - dx
            _refresh_nm_edges(nm, g_current)
            gn_used = k
        _finalize_hist(hist, g_current, loop_closure, gn_used, last_dx_norm, compute_global_error(g_current), gt, global_update=True, active_set_size=g_current.state_size)
    return g_current, hist


def run_selective_demo(cfg: RunConfig):
    if _prefer_cholmod(cfg):
        return run_selective_demo_hybrid(cfg)
    spec = dataset_spec(cfg)
    g_full = load_graph_auto(spec.data_file)
    g_full = reorder_edges(g_full) if cfg.use_reorder_edges else g_full
    gt = load_ground_truth_mat(spec.gt_file) if spec.gt_file is not None else None
    g_current = None
    lc_state = LoopClosureState()
    J, r = _seed_anchor()
    edge_j_index = {}
    hist = RunHistory()
    prev_eta = 0.0
    prev_state_size = 0

    for edge_ids in _batch_edge_ids(len(g_full.edges), cfg.batch_size):
        g_current, edge_nodes, _, _, loop_closure, lc_state = update_graph(edge_ids, g_full, g_current, spec.lc_gap, lc_state)
        local_edge_ids = _local_edge_ids_for_latest_batch(g_current, edge_ids)
        lin = linearize_new(J, r, edge_j_index, g_current, local_edge_ids, g_full)
        J, r = lin.J, lin.r
        edge_j_index = dict(lin.edge_j_index)

        eta = _eta_from_linearization(lin)
        active_vars, global_update, delta_eta = _active_set_from_gate(
            cfg=cfg,
            spec=spec,
            g_current=g_current,
            edge_nodes=edge_nodes,
            loop_closure=loop_closure,
            eta=eta,
            prev_eta=prev_eta,
            prev_state_size=prev_state_size,
        )
        initial_active_size = int(active_vars.size)

        gn_used = 0
        last_dx_norm = 0.0
        for k in range(1, cfg.max_gn_iter + 1):
            if active_vars.size == 0:
                break
            sres = solve_affected(
                b_perm=lin.b,
                p=lin.p,
                affected_vars=active_vars,
                parent=lin.parent,
                factor=lin.factor,
                R=lin.R,
                H_perm=lin.A_perm,
                alpha=getattr(cfg, "selective_alpha", 0.3),
            )
            dx = sres.dx
            last_dx_norm = float(np.max(np.abs(dx))) if dx.size else 0.0
            if converged(dx, spec.dx_th):
                gn_used = k - 1
                break

            g_current.x = g_current.x - dx
            active_vars, affected_edge_ids = find_affected(dx, g_current, spec.dx_th)
            if affected_edge_ids.size == 0 or active_vars.size == 0:
                gn_used = k
                break

            aff = linearize_affected(J, r, edge_j_index, lin.p, g_current, affected_edge_ids)
            J, r = aff.J, aff.r
            lin = LinearizeNewResult(
                factor=aff.factor,
                R=aff.R,
                A_perm=aff.A_perm,
                b=aff.b,
                J=aff.J,
                r=aff.r,
                edge_j_index=edge_j_index,
                p=aff.p,
                parent=None,
            )
            gn_used = k

        _finalize_hist(
            hist,
            g_current,
            loop_closure,
            gn_used,
            last_dx_norm,
            compute_global_error(g_current),
            gt,
            eta=eta,
            delta_eta=delta_eta,
            global_update=global_update,
            active_set_size=initial_active_size,
        )
        prev_eta = eta
        prev_state_size = g_current.state_size
    return g_current, hist


def run_selective_demo_hybrid(cfg: RunConfig):
    if not cholmod_available():
        raise RuntimeError("Hybrid CHOLMOD selective backend requested, but scikit-sparse/CHOLMOD is unavailable.")
    spec = dataset_spec(cfg)
    g_full = load_graph_auto(spec.data_file)
    g_full = reorder_edges(g_full) if cfg.use_reorder_edges else g_full
    gt = load_ground_truth_mat(spec.gt_file) if spec.gt_file is not None else None
    g_current = None
    lc_state = LoopClosureState()
    hist = RunHistory()
    nm: CholmodNormalManager | None = None
    prev_eta = 0.0
    prev_state_size = 0

    for edge_ids in _batch_edge_ids(len(g_full.edges), cfg.batch_size):
        g_current, edge_nodes, _, new_vars, loop_closure, lc_state = update_graph(edge_ids, g_full, g_current, spec.lc_gap, lc_state)
        nm = _build_or_extend_nm(nm, g_current, edge_ids, new_vars, getattr(cfg, "cholmod_order", "amd"))

        eta = nm.information_eta()
        active_vars, global_update, delta_eta = _active_set_from_gate(
            cfg=cfg,
            spec=spec,
            g_current=g_current,
            edge_nodes=edge_nodes,
            loop_closure=loop_closure,
            eta=eta,
            prev_eta=prev_eta,
            prev_state_size=prev_state_size,
        )
        initial_active_size = int(active_vars.size)

        gn_used = 0
        last_dx_norm = 0.0
        for k in range(1, cfg.max_gn_iter + 1):
            if active_vars.size == 0:
                break
            H_perm, b_perm, p = nm.build_permuted_system()
            sres = solve_affected(
                b_perm=b_perm,
                p=p,
                affected_vars=active_vars,
                H_perm=H_perm,
                alpha=getattr(cfg, "selective_alpha", 0.3),
            )
            dx = sres.dx
            last_dx_norm = float(np.max(np.abs(dx))) if dx.size else 0.0
            if converged(dx, spec.dx_th):
                gn_used = k - 1
                break

            g_current.x = g_current.x - dx
            active_vars, affected_edge_ids = find_affected(dx, g_current, spec.dx_th)
            if affected_edge_ids.size == 0 or active_vars.size == 0:
                gn_used = k
                break

            _refresh_nm_edges(nm, g_current, affected_edge_ids)
            gn_used = k

        _finalize_hist(
            hist,
            g_current,
            loop_closure,
            gn_used,
            last_dx_norm,
            compute_global_error(g_current),
            gt,
            eta=eta,
            delta_eta=delta_eta,
            global_update=global_update,
            active_set_size=initial_active_size,
        )
        prev_eta = eta
        prev_state_size = g_current.state_size
    return g_current, hist
