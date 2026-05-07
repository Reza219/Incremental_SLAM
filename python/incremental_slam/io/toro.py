from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import numpy.typing as npt
from incremental_slam.types import Edge, Graph, NodeRef
@dataclass(frozen=True)
class _RawToroEdge:
    from_id: int
    to_id: int
    measurement: npt.NDArray[np.float64]
    information: npt.NDArray[np.float64]
def _strip_comment(line: str) -> str:
    hash_pos = line.find("#")
    if hash_pos >= 0:
        line = line[:hash_pos]
    return line.strip()
def _symmetrize(mat: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    return 0.5 * (mat + mat.T)
def read_graph_toro(path: str | Path) -> Graph:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    x_blocks = []; id_lookup = {}; raw_edges = []; current_offset = 0
    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = _strip_comment(raw_line)
            if not line:
                continue
            tokens = line.split(); tag = tokens[0]
            if tag == "VERTEX2" and len(tokens) >= 5:
                node_id = int(float(tokens[1])); vals = np.array([float(tokens[2]), float(tokens[3]), float(tokens[4])], dtype=float)
                id_lookup[node_id] = NodeRef(offset=current_offset, dimension=3); x_blocks.append(vals); current_offset += 3
            elif tag == "EDGE2" and len(tokens) >= 12:
                from_id = int(float(tokens[1])); to_id = int(float(tokens[2]))
                meas = np.array([float(tokens[3]), float(tokens[4]), float(tokens[5])], dtype=float)
                info = np.zeros((3,3), dtype=float)
                info[0,0] = float(tokens[6]); info[0,1] = float(tokens[7]); info[1,0] = info[0,1]; info[1,1] = float(tokens[8]); info[2,2] = float(tokens[9]); info[0,2] = float(tokens[10]); info[2,0] = info[0,2]; info[1,2] = float(tokens[11]); info[2,1] = info[1,2]
                raw_edges.append(_RawToroEdge(from_id, to_id, meas, _symmetrize(info)))
    x = np.concatenate(x_blocks) if x_blocks else np.empty(0, dtype=float)
    edges = []
    for e in raw_edges:
        from_ref = id_lookup[e.from_id]; to_ref = id_lookup[e.to_id]
        edges.append(Edge("P", e.from_id, e.to_id, from_ref.offset, to_ref.offset, e.measurement, e.information))
    var2node = np.zeros(x.size, dtype=np.int64)
    for node_id, ref in id_lookup.items():
        var2node[ref.offset: ref.offset + ref.dimension] = node_id
    return Graph(x=x, edges=edges, id_lookup=id_lookup, var2node=var2node)
