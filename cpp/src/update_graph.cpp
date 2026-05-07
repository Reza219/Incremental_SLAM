#include "islam/update_graph.hpp"

#include <set>
#include <stdexcept>
#include <unordered_set>

namespace islam {
namespace {

Graph make_empty_graph() {
    Graph g;
    g.x.resize(0);
    return g;
}

Edge copy_edge_with_current_indices(const Edge& e, const Graph& current) {
    const auto it_from = current.id_lookup.find(e.from_id);
    if (it_from == current.id_lookup.end()) {
        throw std::runtime_error("update_graph missing from_id in current graph");
    }
    Edge out = e;
    out.from_idx = it_from->second.offset;
    if (!e.is_unary()) {
        const auto it_to = current.id_lookup.find(e.to_id);
        if (it_to == current.id_lookup.end()) {
            throw std::runtime_error("update_graph missing to_id in current graph");
        }
        out.to_idx = it_to->second.offset;
    } else {
        out.to_idx = -1;
    }
    return out;
}

void append_node_edge(Graph& g, int node_id, int local_edge_idx) {
    g.node_edges[node_id].push_back(local_edge_idx);
}

} // namespace

UpdateGraphResult update_graph(const std::vector<int>& edge_ids,
                               const Graph& full,
                               const Graph* current,
                               LoopClosureState* lc_state,
                               int lc_gap) {
    Graph gcur = current ? *current : make_empty_graph();

    std::vector<Edge> edges_global;
    edges_global.reserve(edge_ids.size());
    for (int eid : edge_ids) {
        if (eid < 0 || eid >= static_cast<int>(full.edges.size())) {
            throw std::out_of_range("update_graph edge id out of range");
        }
        edges_global.push_back(full.edges[static_cast<size_t>(eid)]);
    }

    std::set<int> touched_set;
    for (const auto& e : edges_global) {
        touched_set.insert(e.from_id);
        if (e.to_id >= 0) touched_set.insert(e.to_id);
    }

    UpdateGraphResult res;
    res.current = std::move(gcur);
    res.edge_nodes.assign(touched_set.begin(), touched_set.end());

    std::unordered_set<int> existing_nodes;
    for (const auto& kv : res.current.id_lookup) existing_nodes.insert(kv.first);

    for (int nid : res.edge_nodes) {
        if (existing_nodes.count(nid) == 0) res.new_nodes.push_back(nid);
    }

    if (lc_state != nullptr) {
        res.loop_closure = detect_loop_closure_unordered(edges_global, *lc_state, lc_gap);
    } else {
        for (const auto& e : edges_global) {
            if (!e.is_unary() && existing_nodes.count(e.from_id) > 0 && existing_nodes.count(e.to_id) > 0) {
                res.loop_closure = true;
                break;
            }
        }
    }

    if (!res.new_nodes.empty()) {
        const int old_size = res.current.state_size();
        int total_new_dim = 0;
        for (int nid : res.new_nodes) {
            const auto it = full.id_lookup.find(nid);
            if (it == full.id_lookup.end()) throw std::runtime_error("update_graph new node missing in full graph");
            total_new_dim += it->second.dimension;
        }

        Eigen::VectorXd new_x(old_size + total_new_dim);
        if (old_size > 0) new_x.head(old_size) = res.current.x;

        int cur_off = old_size;
        for (int nid : res.new_nodes) {
            const auto it = full.id_lookup.find(nid);
            const int src_off = it->second.offset;
            const int dim = it->second.dimension;
            new_x.segment(cur_off, dim) = full.x.segment(src_off, dim);
            res.current.id_lookup[nid] = NodeRef{cur_off, dim};
            for (int k = 0; k < dim; ++k) {
                res.current.var2node.push_back(nid);
                res.new_vars.push_back(cur_off + k);
            }
            cur_off += dim;
        }
        res.current.x = std::move(new_x);
    }

    const int first_new_local_idx = static_cast<int>(res.current.edges.size());
    for (size_t k = 0; k < edges_global.size(); ++k) {
        Edge local_edge = copy_edge_with_current_indices(edges_global[k], res.current);
        const int local_edge_idx = first_new_local_idx + static_cast<int>(k);
        append_node_edge(res.current, local_edge.from_id, local_edge_idx);
        if (!local_edge.is_unary()) append_node_edge(res.current, local_edge.to_id, local_edge_idx);
        res.current.edges.push_back(std::move(local_edge));
    }

    return res;
}

} // namespace islam
