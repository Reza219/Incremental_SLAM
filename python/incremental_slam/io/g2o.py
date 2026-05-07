from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import numpy.typing as npt
from incremental_slam.types import Edge, Graph, NodeRef
@dataclass(frozen=True)
class _RawEdge:
    type: str
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
def _parse_edge_se2(vals: list[float]) -> _RawEdge:
    if len(vals) < 11:
        raise ValueError("EDGE_SE2 requires 11 numeric values")
    from_id = int(vals[0]); to_id = int(vals[1])
    measurement = np.asarray(vals[2:5], dtype=float)
    ut = vals[5:11]
    info = np.array([[ut[0], ut[1], ut[2]], [ut[1], ut[3], ut[4]], [ut[2], ut[4], ut[5]]], dtype=float)
    return _RawEdge("P", from_id, to_id, measurement, _symmetrize(info))
def _parse_edge_se2_xy(vals: list[float]) -> _RawEdge:
    if len(vals) < 7:
        raise ValueError("EDGE_SE2_XY requires 7 numeric values")
    from_id = int(vals[0]); to_id = int(vals[1])
    measurement = np.asarray(vals[2:4], dtype=float)
    ut = vals[4:7]
    info = np.array([[ut[0], ut[1]], [ut[1], ut[2]]], dtype=float)
    return _RawEdge("L", from_id, to_id, measurement, _symmetrize(info))
def read_graph_g2o(path: str | Path) -> Graph:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    x_blocks = []
    id_lookup = {}
    raw_edges = []
    current_offset = 0
    with path.open("r", encoding="utf-8") as f:
        for line_num, raw_line in enumerate(f, start=1):
            line = _strip_comment(raw_line)
            if not line:
                continue
            tokens = line.split(); tag = tokens[0]
            try:
                vals = [float(tok) for tok in tokens[1:]]
            except ValueError as exc:
                raise ValueError(f"line {line_num}: failed numeric parse") from exc
            if tag == "VERTEX_SE2":
                node_id = int(vals[0]); state = np.asarray(vals[1:4], dtype=float)
                id_lookup[node_id] = NodeRef(offset=current_offset, dimension=3); x_blocks.append(state); current_offset += 3
            elif tag == "VERTEX_XY":
                node_id = int(vals[0]); state = np.asarray(vals[1:3], dtype=float)
                id_lookup[node_id] = NodeRef(offset=current_offset, dimension=2); x_blocks.append(state); current_offset += 2
            elif tag == "EDGE_SE2":
                raw_edges.append(_parse_edge_se2(vals))
            elif tag == "EDGE_SE2_XY":
                raw_edges.append(_parse_edge_se2_xy(vals))
    x = np.concatenate(x_blocks) if x_blocks else np.empty(0, dtype=float)
    edges = []
    for raw in raw_edges:
        from_ref = id_lookup[raw.from_id]; to_ref = id_lookup[raw.to_id]
        edges.append(Edge(raw.type, raw.from_id, raw.to_id, from_ref.offset, to_ref.offset, raw.measurement, raw.information))
    var2node = np.zeros(x.size, dtype=np.int64)
    for node_id, ref in id_lookup.items():
        var2node[ref.offset: ref.offset + ref.dimension] = node_id
    return Graph(x=x, edges=edges, id_lookup=id_lookup, var2node=var2node)
