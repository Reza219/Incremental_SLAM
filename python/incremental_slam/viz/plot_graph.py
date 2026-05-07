from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
from incremental_slam.types import Graph

def plot_graph(g: Graph, ax=None, show_edges: bool = True):
    if ax is None: _, ax = plt.subplots()
    poses = []; lms = []
    for _, ref in sorted(g.id_lookup.items(), key=lambda kv: kv[1].offset):
        if ref.dimension == 3: poses.append(g.x[ref.offset: ref.offset + 2])
        elif ref.dimension == 2: lms.append(g.x[ref.offset: ref.offset + 2])
    if poses:
        P = np.vstack(poses); ax.plot(P[:,0], P[:,1], '-o', markersize=2)
    if lms:
        L = np.vstack(lms); ax.scatter(L[:,0], L[:,1], s=10, marker='x')
    if show_edges:
        for e in g.edges:
            if e.type == 'G' or e.to_id is None: continue
            rf = g.id_lookup[e.from_id]; rt = g.id_lookup[e.to_id]; pf = g.x[rf.offset: rf.offset + 2]; pt = g.x[rt.offset: rt.offset + 2]; ax.plot([pf[0], pt[0]], [pf[1], pt[1]], '-', linewidth=0.4, alpha=0.4)
    ax.set_aspect('equal', adjustable='box'); return ax
