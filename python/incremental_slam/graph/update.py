from __future__ import annotations
import numpy as np
import numpy.typing as npt
from incremental_slam.graph.loop_closure import LoopClosureState, detect_loop_closure_unordered
from incremental_slam.types import Edge, Graph, NodeRef

def _as_1d_int_array(x: npt.ArrayLike | None) -> npt.NDArray[np.int64]:
    return np.empty(0, dtype=np.int64) if x is None else np.asarray(x, dtype=np.int64).reshape(-1)
def _empty_graph_like() -> Graph:
    g_current = Graph(x=np.empty(0, dtype=float), edges=[], id_lookup={}, var2node=np.empty(0, dtype=np.int64)); setattr(g_current, 'node_edges', {}); return g_current
def _ensure_current_graph(g_current: Graph | None) -> Graph:
    if g_current is None: return _empty_graph_like()
    if not hasattr(g_current, 'node_edges'): setattr(g_current, 'node_edges', {})
    if g_current.var2node is None: g_current.var2node = np.empty(0, dtype=np.int64)
    return g_current
def _copy_edge_with_current_indices(edge: Edge, g_current: Graph) -> Edge:
    from_ref = g_current.id_lookup[edge.from_id]; from_idx = from_ref.offset
    to_idx = None if edge.type == 'G' else g_current.id_lookup[edge.to_id].offset
    return Edge(edge.type, edge.from_id, edge.to_id, from_idx, to_idx, np.array(edge.measurement, copy=True), np.array(edge.information, copy=True))
def update_graph(edge_ids: npt.ArrayLike, g_full: Graph, g_current: Graph | None = None, lc_gap: int = 5, lc_state: LoopClosureState | None = None) -> tuple[Graph, npt.NDArray[np.int64], npt.NDArray[np.int64], npt.NDArray[np.int64], bool, LoopClosureState]:
    g_current = _ensure_current_graph(g_current); edge_ids_vec = _as_1d_int_array(edge_ids); edges_global = [g_full.edges[int(eid)] for eid in edge_ids_vec]
    touched_nodes = []
    for edge in edges_global:
        touched_nodes.append(int(edge.from_id));
        if edge.to_id is not None: touched_nodes.append(int(edge.to_id))
    edge_nodes = np.array(sorted(set(touched_nodes)), dtype=np.int64) if touched_nodes else np.empty(0, dtype=np.int64)
    seen = set(g_current.id_lookup.keys()); new_nodes_list = [nid for nid in edge_nodes.tolist() if nid not in seen]; new_nodes = np.array(new_nodes_list, dtype=np.int64)
    loop_closure, lc_state = detect_loop_closure_unordered(edges_global, lc_state, order_gap=lc_gap)
    new_vars_list = []
    if new_nodes.size > 0:
        current_size = g_current.state_size; x_blocks = []; v_blocks = []
        for nid in new_nodes:
            ref_full = g_full.id_lookup[int(nid)]; dim = ref_full.dimension; src = g_full.x[ref_full.offset: ref_full.offset + dim]
            g_current.id_lookup[int(nid)] = NodeRef(offset=current_size, dimension=dim); x_blocks.append(np.array(src, dtype=float, copy=True)); v_blocks.append(np.full(dim, int(nid), dtype=np.int64)); new_vars_list.extend(range(current_size, current_size + dim)); current_size += dim
        g_current.x = np.concatenate([g_current.x, *x_blocks]) if x_blocks else g_current.x; g_current.var2node = np.concatenate([g_current.var2node, *v_blocks]) if v_blocks else g_current.var2node
    new_vars = np.array(sorted(set(new_vars_list)), dtype=np.int64) if new_vars_list else np.empty(0, dtype=np.int64)
    new_edges_local = [_copy_edge_with_current_indices(edge, g_current) for edge in edges_global]
    node_edges = getattr(g_current, 'node_edges'); first_new_local_edge_idx = len(g_current.edges)
    for k, edge in enumerate(new_edges_local):
        local_edge_idx = first_new_local_edge_idx + k; node_edges.setdefault(int(edge.from_id), []).append(local_edge_idx)
        if edge.type != 'G' and edge.to_id is not None: node_edges.setdefault(int(edge.to_id), []).append(local_edge_idx)
    g_current.edges.extend(new_edges_local)
    return g_current, edge_nodes, new_nodes, new_vars, loop_closure, lc_state
