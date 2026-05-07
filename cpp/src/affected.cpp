#include "islam/affected.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace islam {
namespace {

AffectedSet expand_nodes_to_vars_and_edges(const Graph& g,
                                           const std::vector<int>& nodes) {
    AffectedSet out;
    out.nodes = nodes;

    std::unordered_set<int> edge_set;
    for (int nid : out.nodes) {
        const auto it_node = g.id_lookup.find(nid);
        if (it_node != g.id_lookup.end()) {
            const int off = it_node->second.offset;
            const int dim = it_node->second.dimension;
            for (int k = 0; k < dim; ++k) out.vars.push_back(off + k);
        }
        const auto it_edges = g.node_edges.find(nid);
        if (it_edges != g.node_edges.end()) {
            for (int eid : it_edges->second) edge_set.insert(eid);
        }
    }

    std::sort(out.vars.begin(), out.vars.end());
    out.vars.erase(std::unique(out.vars.begin(), out.vars.end()), out.vars.end());

    out.edges.assign(edge_set.begin(), edge_set.end());
    std::sort(out.edges.begin(), out.edges.end());
    return out;
}

} // namespace

AffectedSet find_affected(const Graph& g,
                          const Eigen::VectorXd& dx,
                          double threshold) {
    if (dx.size() != g.state_size()) {
        throw std::runtime_error("find_affected: dx size mismatch");
    }
    std::vector<int> all_vars(static_cast<size_t>(g.state_size()));
    for (int i = 0; i < g.state_size(); ++i) all_vars[static_cast<size_t>(i)] = i;
    return find_affected(g, dx, all_vars, threshold);
}

AffectedSet find_affected(const Graph& g,
                          const Eigen::VectorXd& dxa,
                          const std::vector<int>& candidate_vars,
                          double threshold) {
    if (dxa.size() != static_cast<int>(candidate_vars.size())) {
        throw std::runtime_error("find_affected: size mismatch between dxa and candidate_vars");
    }

    // Step 1: prune by threshold at node-block level.
    std::unordered_set<int> node_set;
    for (int k = 0; k < static_cast<int>(candidate_vars.size()); ++k) {
        if (std::abs(dxa[k]) <= threshold) continue;
        const int var = candidate_vars[static_cast<size_t>(k)];
        if (var < 0 || var >= static_cast<int>(g.var2node.size())) {
            throw std::runtime_error("find_affected: var index out of range");
        }
        node_set.insert(g.var2node[static_cast<size_t>(var)]);
    }

    // Step 2: expand once through incident factors. This implements the
    // paper prune-expand SPO rule: if a retained node participates in an
    // edge, all endpoint node blocks of that edge are kept for the next active
    // set so the shared residual can be made consistent.
    std::vector<int> frontier(node_set.begin(), node_set.end());
    for (int nid : frontier) {
        const auto it_edges = g.node_edges.find(nid);
        if (it_edges == g.node_edges.end()) continue;
        for (int eid : it_edges->second) {
            if (eid < 0 || eid >= static_cast<int>(g.edges.size())) {
                throw std::runtime_error("find_affected: edge index out of range");
            }
            const auto& e = g.edges[static_cast<size_t>(eid)];
            node_set.insert(g.var2node.at(static_cast<size_t>(e.from_idx)));
            if (!e.is_unary()) {
                node_set.insert(g.var2node.at(static_cast<size_t>(e.to_idx)));
            }
        }
    }

    std::vector<int> nodes(node_set.begin(), node_set.end());
    std::sort(nodes.begin(), nodes.end());
    return expand_nodes_to_vars_and_edges(g, nodes);
}

} // namespace islam
