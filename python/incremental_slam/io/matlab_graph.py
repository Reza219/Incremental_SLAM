from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
from scipy.io import loadmat
from incremental_slam.types import Edge, Graph, NodeRef

def _mat_struct_to_graph(g: Any) -> Graph:
    x = np.asarray(getattr(g, 'x')).reshape(-1).astype(float)
    edges_in = getattr(g, 'edges')
    id_lookup = {}
    edges = []
    flat_edges = np.ravel(edges_in)
    for e in flat_edges:
        etype = str(np.asarray(getattr(e, 'type')).reshape(-1)[0])
        from_idx = int(np.asarray(getattr(e, 'fromIdx')).reshape(-1)[0])
        to_attr = getattr(e, 'toIdx', None)
        to_idx = None if to_attr is None or np.size(to_attr) == 0 else int(np.asarray(to_attr).reshape(-1)[0])
        z = np.asarray(getattr(e, 'measurement')).reshape(-1).astype(float)
        info = np.asarray(getattr(e, 'information')).astype(float)
        from_id = int(np.asarray(getattr(e, 'fromId', from_idx)).reshape(-1)[0])
        to_id = None
        if hasattr(e, 'toId'):
            tmp = getattr(e, 'toId'); to_id = None if tmp is None or np.size(tmp) == 0 else int(np.asarray(tmp).reshape(-1)[0])
        edges.append(Edge(etype, from_id, to_id, from_idx, to_idx, z, info))
    var2node = np.asarray(getattr(g, 'var2node', np.empty(0))).reshape(-1).astype(np.int64) if hasattr(g, 'var2node') else None
    return Graph(x=x, edges=edges, id_lookup=id_lookup, var2node=var2node)
def read_graph_mat(path: str | Path) -> Graph:
    path = Path(path)
    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    if 'g' in data:
        return _mat_struct_to_graph(data['g'])
    if 'S' in data and hasattr(data['S'], 'g'):
        return _mat_struct_to_graph(data['S'].g)
    raise KeyError("Could not find graph struct 'g' or 'S.g' in MAT file")
