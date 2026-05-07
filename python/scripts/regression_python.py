from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from incremental_slam.graph.loop_closure import LoopClosureState
from incremental_slam.graph.reorder import reorder_edges
from incremental_slam.graph.update import update_graph
from incremental_slam.io.g2o import read_graph_g2o
from incremental_slam.io.toro import read_graph_toro
from incremental_slam.linearization.linearize_new import linearize_new

def load_graph_auto(path: Path): return read_graph_g2o(path) if path.suffix.lower() == '.g2o' else read_graph_toro(path)
def dataset_path(name: str) -> tuple[Path, int]:
    if name == 'MIT': return Path('data') / 'MITb_g2o.g2o', 4
    if name == 'Intel': return Path('data') / 'INTEL_g2o.g2o', 4
    if name == 'CSAIL': return Path('data') / 'CSAIL_P_toro.graph', 4
    raise ValueError('Unsupported dataset')
def edge_summary(edge):
    return {'type': edge.type, 'from_id': int(edge.from_id), 'to_id': None if edge.to_id is None else int(edge.to_id), 'from_idx': int(edge.from_idx), 'to_idx': None if edge.to_idx is None else int(edge.to_idx), 'meas_dim': int(edge.measurement.size), 'info_shape': list(edge.information.shape)}
def main(dataset_name: str = 'MIT', batch_size: int = 1):
    data_file, lc_gap = dataset_path(dataset_name); g = reorder_edges(load_graph_auto(data_file)); summary = {'dataset': dataset_name, 'state_size': int(g.x.size), 'num_edges': int(len(g.edges)), 'first_edges': [edge_summary(e) for e in g.edges[:5]]}; edge_ids = np.arange(min(batch_size, len(g.edges)), dtype=np.int64); J = sp.eye(3, format='csc', dtype=float); r = np.zeros(3, dtype=float); edge_j_index = {}; g_current = None; lc_state = LoopClosureState(); g_current, edge_nodes, new_nodes, new_vars, loop_closure, lc_state = update_graph(edge_ids, g, g_current, lc_gap, lc_state); lin = linearize_new(J, r, edge_j_index, g_current, edge_ids, g); summary['first_batch'] = {'edge_ids': edge_ids.tolist(), 'edge_nodes': edge_nodes.tolist(), 'new_nodes': new_nodes.tolist(), 'new_vars': new_vars.tolist(), 'loop_closure': bool(loop_closure), 'J_size': list(lin.J.shape), 'r_size': [int(lin.r.size), 1], 'J_nnz': int(lin.J.nnz), 'b_size': [int(lin.b.size), 1], 'p_len': int(lin.p.size), 'parent_len': None if lin.parent is None else int(lin.parent.size)}; print(json.dumps(summary, indent=2))
if __name__ == '__main__': main('MIT', 1)
