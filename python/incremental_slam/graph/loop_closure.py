from __future__ import annotations
from dataclasses import dataclass, field
from incremental_slam.types import Edge
@dataclass
class LoopClosureState:
    pose_order: dict[int, int] = field(default_factory=dict)
    lm2orders: dict[int, list[int]] = field(default_factory=dict)
    next_order: int = 1

def detect_loop_closure_unordered(edges: list[Edge], state: LoopClosureState | None = None, order_gap: int = 5) -> tuple[bool, LoopClosureState]:
    if state is None:
        state = LoopClosureState()
    for edge in edges:
        if edge.type == 'P' and edge.to_id is not None:
            i = int(edge.from_id); j = int(edge.to_id)
            i_in = i in state.pose_order; j_in = j in state.pose_order
            if i_in and j_in and abs(state.pose_order[i] - state.pose_order[j]) > order_gap:
                return True, state
            if not i_in: state.pose_order[i] = state.next_order; state.next_order += 1
            if not j_in: state.pose_order[j] = state.next_order; state.next_order += 1
        elif edge.type == 'L' and edge.to_id is not None:
            pid = int(edge.from_id); lid = int(edge.to_id)
            if pid not in state.pose_order: state.pose_order[pid] = state.next_order; state.next_order += 1
            po = state.pose_order[pid]
            if lid in state.lm2orders:
                prev = state.lm2orders[lid]
                if any(abs(old_po - po) > order_gap for old_po in prev): return True, state
                prev.append(po)
            else:
                state.lm2orders[lid] = [po]
        elif edge.type == 'G':
            pid = int(edge.from_id)
            if pid not in state.pose_order: state.pose_order[pid] = state.next_order; state.next_order += 1
    return False, state
