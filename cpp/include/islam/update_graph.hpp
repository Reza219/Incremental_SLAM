#pragma once

#include "islam/graph.hpp"
#include "islam/loop_closure.hpp"

#include <vector>

namespace islam {

struct UpdateGraphResult {
    Graph current;
    std::vector<int> edge_nodes;
    std::vector<int> new_nodes;
    std::vector<int> new_vars;
    bool loop_closure = false;
};

UpdateGraphResult update_graph(const std::vector<int>& edge_ids,
                               const Graph& full,
                               const Graph* current,
                               LoopClosureState* lc_state = nullptr,
                               int lc_gap = 5);

} // namespace islam
