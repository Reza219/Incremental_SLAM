from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal
import numpy as np
import numpy.typing as npt
EdgeType = Literal["P", "L", "G"]
@dataclass(frozen=True)
class NodeRef:
    offset: int
    dimension: int
@dataclass
class Edge:
    type: EdgeType
    from_id: int
    to_id: int | None
    from_idx: int
    to_idx: int | None
    measurement: npt.NDArray[np.float64]
    information: npt.NDArray[np.float64]
    def __post_init__(self) -> None:
        self.measurement = np.asarray(self.measurement, dtype=float).reshape(-1)
        self.information = np.asarray(self.information, dtype=float)
        if self.information.ndim != 2 or self.information.shape[0] != self.information.shape[1]:
            raise ValueError("information must be square")
        if self.information.shape[0] != self.measurement.size:
            raise ValueError("measurement dimension must match information matrix size")
@dataclass
class Graph:
    x: npt.NDArray[np.float64]
    edges: list[Edge] = field(default_factory=list)
    id_lookup: dict[int, NodeRef] = field(default_factory=dict)
    var2node: npt.NDArray[np.int64] | None = None
    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=float).reshape(-1)
        if self.var2node is not None:
            self.var2node = np.asarray(self.var2node, dtype=np.int64).reshape(-1)
            if self.var2node.shape[0] != self.x.shape[0]:
                raise ValueError("var2node must have same length as x")
    @property
    def state_size(self) -> int:
        return int(self.x.size)
    def block(self, start: int, dim: int) -> npt.NDArray[np.float64]:
        stop = start + dim
        if start < 0 or stop > self.state_size:
            raise IndexError(f"state block [{start}:{stop}] out of range")
        return self.x[start:stop]
    def node_block(self, node_id: int) -> npt.NDArray[np.float64]:
        ref = self.id_lookup[node_id]
        return self.block(ref.offset, ref.dimension)
