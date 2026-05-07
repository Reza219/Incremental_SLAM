from __future__ import annotations

import numpy as np
import numpy.typing as npt

from incremental_slam.linearization.edge_jacobian import jacobian_edge_jr
from incremental_slam.types import Edge, Graph


def edge_touched_vars(edge: Edge) -> npt.NDArray[np.int64]:
    touched: list[int] = []
    if edge.type in ("P", "L", "G"):
        touched.extend(range(edge.from_idx, edge.from_idx + 3))
    if edge.to_idx is not None:
        dim_to = 3 if edge.type == "P" else 2
        touched.extend(range(edge.to_idx, edge.to_idx + dim_to))
    return np.asarray(sorted(set(touched)), dtype=np.int64)


def edge_touched_nodes(edge: Edge) -> npt.NDArray[np.int64]:
    nodes = [int(edge.from_id)]
    if edge.to_id is not None:
        nodes.append(int(edge.to_id))
    return np.asarray(sorted(set(nodes)), dtype=np.int64)


def linearize_edge_block(edge: Edge, g: Graph):
    Je, re = jacobian_edge_jr(edge, g, state_size=g.state_size)
    return Je.T, re, edge_touched_vars(edge), edge_touched_nodes(edge)


def linearize_all_current_blocks(g: Graph):
    blocks = []
    for eid, edge in enumerate(g.edges):
        C, r, touched_vars, touched_nodes = linearize_edge_block(edge, g)
        blocks.append((int(eid), C, r, touched_vars, touched_nodes))
    return blocks
