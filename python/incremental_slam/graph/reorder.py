from __future__ import annotations
from incremental_slam.types import Edge, Graph

def reorder_edges(g: Graph) -> Graph:
    if 0 not in g.id_lookup:
        raise ValueError("Vertex id0 not found in g.id_lookup")
    E = g.edges
    initialized = {0: True}; inserted = [False] * len(E); ordered = []
    while True:
        changed = False; new_batch = []
        for i, ei in enumerate(E):
            if inserted[i] or ei.type == 'G' or ei.to_id is None:
                continue
            from_ok = ei.from_id in initialized; to_ok = ei.to_id in initialized
            if from_ok ^ to_ok:
                new_batch.append(ei); inserted[i] = True; initialized[ei.from_id] = True; initialized[ei.to_id] = True; changed = True; break
        for i, ei in enumerate(E):
            if inserted[i]: continue
            if ei.type == 'G':
                if ei.from_id in initialized: new_batch.append(ei); inserted[i] = True; changed = True
                continue
            if ei.to_id is not None and ei.from_id in initialized and ei.to_id in initialized:
                new_batch.append(ei); inserted[i] = True; changed = True
        if not changed: break
        ordered.extend(new_batch)
    dropped = sum(1 for v in inserted if not v)
    if dropped: print(f"Warning: could not safely reorder {dropped} edge(s); dropping them as disconnected/unreachable from id0.")
    return Graph(x=g.x.copy(), edges=[Edge(e.type, e.from_id, e.to_id, e.from_idx, e.to_idx, e.measurement.copy(), e.information.copy()) for e in ordered], id_lookup=dict(g.id_lookup), var2node=None if g.var2node is None else g.var2node.copy())
