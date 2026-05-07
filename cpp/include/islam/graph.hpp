#pragma once

#include "islam/types.hpp"

#include <Eigen/Dense>

#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace islam {

struct Graph {
    Eigen::VectorXd x;
    std::vector<Edge> edges;
    std::unordered_map<int, NodeRef> id_lookup;
    std::vector<int> var2node;
    std::unordered_map<int, std::vector<int>> node_edges;

    [[nodiscard]] int state_size() const noexcept {
        return static_cast<int>(x.size());
    }

    [[nodiscard]] Eigen::VectorXd block(int start, int dim) const {
        if (start < 0 || dim < 0 || start + dim > x.size()) {
            throw std::out_of_range("Graph::block out of range");
        }
        return x.segment(start, dim);
    }

    [[nodiscard]] Eigen::VectorXd node_block(int node_id) const {
        const auto it = id_lookup.find(node_id);
        if (it == id_lookup.end()) {
            throw std::out_of_range("Graph::node_block missing node_id");
        }
        return block(it->second.offset, it->second.dimension);
    }
};

} // namespace islam
