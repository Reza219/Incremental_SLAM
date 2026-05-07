import numpy as np

from incremental_slam.config import DatasetSpec, RunConfig
from incremental_slam.graph.affected import prune_expand_active_set
from incremental_slam.main import _active_set_from_gate
from incremental_slam.metrics.error import compute_global_error
from incremental_slam.metrics.information import detrended_delta_eta, eta_from_R
from incremental_slam.types import Edge, Graph, NodeRef


def _three_pose_chain():
    g = Graph(
        x=np.zeros(9),
        edges=[],
        id_lookup={
            0: NodeRef(0, 3),
            1: NodeRef(3, 3),
            2: NodeRef(6, 3),
        },
        var2node=np.array([0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=np.int64),
    )
    g.edges = [
        Edge("P", 0, 1, 0, 3, np.zeros(3), np.eye(3)),
        Edge("P", 1, 2, 3, 6, np.zeros(3), np.eye(3)),
    ]
    g.node_edges = {0: [0], 1: [0, 1], 2: [1]}
    return g


def test_eta_and_detrended_delta_eta_match_paper_formula():
    R = np.diag([2.0, 4.0, 8.0])
    assert np.isclose(eta_from_R(R), np.log(2.0) + np.log(4.0) + np.log(8.0))
    assert np.isclose(detrended_delta_eta(10.0, 6.0, n_prev=3, n_t=6), 7.0)


def test_igg_gate_uses_information_gain_not_objective_change():
    g = _three_pose_chain()
    spec = DatasetSpec(data_file=None, gt_file=None, dx_th=1e-3, lc_gap=4, ent_th=1.0)
    cfg = RunConfig(gating_mode="igg")

    local_active, is_global, delta = _active_set_from_gate(
        cfg=cfg,
        spec=spec,
        g_current=g,
        edge_nodes=np.array([1], dtype=np.int64),
        loop_closure=True,  # IGG mode should not rely on this flag.
        eta=3.0,
        prev_eta=2.9,
        prev_state_size=9,
    )
    assert not is_global
    assert np.array_equal(local_active, np.array([3, 4, 5]))
    assert delta < spec.ent_th

    all_active, is_global, delta = _active_set_from_gate(
        cfg=cfg,
        spec=spec,
        g_current=g,
        edge_nodes=np.array([1], dtype=np.int64),
        loop_closure=False,
        eta=8.0,
        prev_eta=2.9,
        prev_state_size=9,
    )
    assert is_global
    assert np.array_equal(all_active, np.arange(g.state_size))
    assert delta >= spec.ent_th


def test_spo_prune_then_expand_keeps_full_node_blocks_and_neighbors():
    g = _three_pose_chain()
    dx = np.zeros(9)
    dx[4] = 0.1  # one coordinate of node 1 remains active
    active_vars, affected_edges = prune_expand_active_set(dx, g, dx_th=1e-3)
    assert np.array_equal(active_vars, np.arange(9))
    assert np.array_equal(affected_edges, np.array([0, 1]))


def test_normalized_chi_squared_divides_by_scalar_measurements():
    g = Graph(
        x=np.array([3.0, 4.0, 0.0]),
        id_lookup={0: NodeRef(0, 3)},
        var2node=np.array([0, 0, 0], dtype=np.int64),
    )
    g.edges = [Edge("G", 0, None, 0, None, np.array([0.0, 0.0]), np.eye(2))]
    assert np.isclose(compute_global_error(g), 25.0 / 2.0)


def test_linearize_new_uses_current_graph_local_edge_indices_not_full_offsets():
    # Source graph has node ids/offsets that are not the same as the current
    # incremental state layout. The new edge must be linearized from
    # g_current.edges, whose from_idx/to_idx are local offsets.
    from incremental_slam.linearization.linearize_new import linearize_new
    import scipy.sparse as sp

    g_full = Graph(
        x=np.array([100.0, 100.0, 0.0, 0.0, 0.0, 0.0]),
        edges=[],
        id_lookup={10: NodeRef(0, 3), 0: NodeRef(3, 3)},
        var2node=np.array([10, 10, 10, 0, 0, 0], dtype=np.int64),
    )
    # Current graph stores the same nodes in a different incremental order.
    g = Graph(
        x=np.array([0.0, 0.0, 0.0, 100.0, 100.0, 0.0]),
        edges=[],
        id_lookup={0: NodeRef(0, 3), 10: NodeRef(3, 3)},
        var2node=np.array([0, 0, 0, 10, 10, 10], dtype=np.int64),
    )
    g.edges = [Edge("P", 0, 10, 0, 3, np.array([100.0, 100.0, 0.0]), np.eye(3))]
    g_full.edges = [Edge("P", 0, 10, 3, 0, np.array([100.0, 100.0, 0.0]), np.eye(3))]

    J0 = sp.eye(3, format="csc")
    r0 = np.zeros(3)
    lin = linearize_new(J0, r0, {}, g, np.array([0], dtype=np.int64), g_full)
    assert lin.J.shape[1] == 6
    assert lin.r.shape[0] == 6
    assert np.linalg.norm(lin.r[-3:]) < 1e-12
