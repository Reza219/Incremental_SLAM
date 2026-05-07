#pragma once

#include "islam/types.hpp"

#include <unordered_map>
#include <vector>

namespace islam {

struct LoopClosureState {
    std::unordered_map<int, int> pose_order;
    std::unordered_map<int, std::vector<int>> landmark_pose_orders;
    int next_order = 1;
};

bool detect_loop_closure_unordered(const std::vector<Edge>& edges,
                                   LoopClosureState& state,
                                   int order_gap = 5);

} // namespace islam
