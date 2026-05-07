#include "islam/reorder_edges.hpp"

#include <iostream>
#include <unordered_set>
#include <stdexcept>
#include <vector>

namespace islam {

Graph reorder_edges(const Graph& g) {
    if (g.id_lookup.find(0) == g.id_lookup.end()) {
        throw std::runtime_error("reorder_edges requires seed node id 0");
    }

    Graph out = g;
    out.edges.clear();
    out.node_edges.clear();

    std::unordered_set<int> initialized{0};
    std::vector<bool> inserted(g.edges.size(), false);

    while (true) {
        bool changed = false;
        std::vector<Edge> batch;

        for (size_t i = 0; i < g.edges.size(); ++i) {
            if (inserted[i]) continue;
            const Edge& e = g.edges[i];
            if (e.is_unary()) continue;
            const bool from_ok = initialized.count(e.from_id) > 0;
            const bool to_ok = initialized.count(e.to_id) > 0;
            if (from_ok != to_ok) {
                batch.push_back(e);
                inserted[i] = true;
                changed = true;
                initialized.insert(e.from_id);
                initialized.insert(e.to_id);
                break;
            }
        }

        for (size_t i = 0; i < g.edges.size(); ++i) {
            if (inserted[i]) continue;
            const Edge& e = g.edges[i];
            if (e.is_unary()) {
                if (initialized.count(e.from_id) > 0) {
                    batch.push_back(e);
                    inserted[i] = true;
                    changed = true;
                }
                continue;
            }
            if (e.to_id >= 0 && initialized.count(e.from_id) > 0 && initialized.count(e.to_id) > 0) {
                batch.push_back(e);
                inserted[i] = true;
                changed = true;
            }
        }

        if (!changed) break;
        out.edges.insert(out.edges.end(), batch.begin(), batch.end());
    }

    int dropped = 0;
    for (bool ok : inserted) if (!ok) ++dropped;
    if (dropped > 0) {
        std::cerr << "Warning: could not safely reorder " << dropped
                  << " edge(s); dropping disconnected/unreachable edges from id0.\n";
    }

    return out;
}

} // namespace islam
