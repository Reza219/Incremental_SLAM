from __future__ import annotations

import numpy as np
import numpy.typing as npt

from incremental_slam.types import Graph


def _vars_for_node(g: Graph, node_id: int) -> list[int]:
    ref = g.id_lookup[int(node_id)]
    return list(range(ref.offset, ref.offset + ref.dimension))


def vars_for_nodes(g: Graph, node_ids: npt.ArrayLike | None) -> npt.NDArray[np.int64]:
    if node_ids is None:
        return np.empty(0, dtype=np.int64)
    nodes = np.asarray(node_ids, dtype=np.int64).reshape(-1)
    out: set[int] = set()
    for nid in nodes:
        if int(nid) in g.id_lookup:
            out.update(_vars_for_node(g, int(nid)))
    return np.asarray(sorted(out), dtype=np.int64)


def all_vars(g: Graph) -> npt.NDArray[np.int64]:
    return np.arange(g.state_size, dtype=np.int64)


def edge_ids_for_nodes(g: Graph, node_ids: npt.ArrayLike | None) -> npt.NDArray[np.int64]:
    if node_ids is None:
        return np.empty(0, dtype=np.int64)
    node_edges = getattr(g, "node_edges", {})
    out: set[int] = set()
    for nid in np.asarray(node_ids, dtype=np.int64).reshape(-1):
        out.update(int(eid) for eid in node_edges.get(int(nid), []))
    return np.asarray(sorted(out), dtype=np.int64)


def nodes_for_vars(g: Graph, vars_: npt.ArrayLike | None) -> npt.NDArray[np.int64]:
    if vars_ is None or g.var2node is None:
        return np.empty(0, dtype=np.int64)
    vv = np.asarray(vars_, dtype=np.int64).reshape(-1)
    vv = vv[(vv >= 0) & (vv < g.var2node.size)]
    if vv.size == 0:
        return np.empty(0, dtype=np.int64)
    return np.asarray(sorted(set(int(g.var2node[v]) for v in vv)), dtype=np.int64)


def _edge_endpoint_nodes(g: Graph, edge_ids: npt.ArrayLike | None) -> set[int]:
    nodes: set[int] = set()
    if edge_ids is None:
        return nodes
    for eid in np.asarray(edge_ids, dtype=np.int64).reshape(-1):
        edge = g.edges[int(eid)]
        nodes.add(int(edge.from_id))
        if edge.to_id is not None:
            nodes.add(int(edge.to_id))
    return nodes


def initial_active_vars_from_new_measurements(g: Graph, edge_nodes: npt.ArrayLike | None) -> npt.NDArray[np.int64]:
    """Paper Algorithm 1 line 3: variables involved in new measurements."""
    return vars_for_nodes(g, edge_nodes)


def prune_expand_active_set(
    dx: npt.ArrayLike,
    g: Graph,
    dx_th: float,
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """SPO active-set update from the paper.

    1. Prune at node/block level: keep nodes with at least one increment above tau_d.
    2. Expand through pose-graph connectivity: include all nodes that share an incident
       edge with the retained nodes.
    3. Return the corresponding variable blocks and the incident edges to relinearize.
    """
    dx_vec = np.asarray(dx, dtype=float).reshape(-1)
    raw_vars = np.flatnonzero(np.abs(dx_vec) > float(dx_th)).astype(np.int64)
    if raw_vars.size == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)

    retained_nodes = set(int(n) for n in nodes_for_vars(g, raw_vars))
    incident_edges = edge_ids_for_nodes(g, np.asarray(sorted(retained_nodes), dtype=np.int64))
    expanded_nodes = set(retained_nodes)
    expanded_nodes.update(_edge_endpoint_nodes(g, incident_edges))

    active_vars = vars_for_nodes(g, np.asarray(sorted(expanded_nodes), dtype=np.int64))
    return active_vars, incident_edges


def find_affected(dx: npt.ArrayLike, g: Graph, dx_th: float) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
    """Backward-compatible alias for the paper-faithful prune/expand rule."""
    return prune_expand_active_set(dx, g, dx_th)
