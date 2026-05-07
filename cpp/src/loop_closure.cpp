#include "islam/loop_closure.hpp"

#include <algorithm>
#include <cstdlib>

namespace islam {

bool detect_loop_closure_unordered(const std::vector<Edge>& edges,
                                   LoopClosureState& state,
                                   int order_gap) {
    if (order_gap < 0) order_gap = 0;

    for (const auto& e : edges) {
        switch (e.type) {
        case EdgeType::PosePose:
        case EdgeType::PosePose3D: {
            const int i = e.from_id;
            const int j = e.to_id;
            const auto it_i = state.pose_order.find(i);
            const auto it_j = state.pose_order.find(j);
            const bool i_seen = it_i != state.pose_order.end();
            const bool j_seen = it_j != state.pose_order.end();
            if (i_seen && j_seen && std::abs(it_i->second - it_j->second) > order_gap) {
                return true;
            }
            if (!i_seen) state.pose_order[i] = state.next_order++;
            if (!j_seen) state.pose_order[j] = state.next_order++;
            break;
        }
        case EdgeType::PoseLandmark: {
            const int pose_id = e.from_id;
            const int lm_id = e.to_id;
            auto it_pose = state.pose_order.find(pose_id);
            if (it_pose == state.pose_order.end()) {
                it_pose = state.pose_order.emplace(pose_id, state.next_order++).first;
            }
            const int pose_ord = it_pose->second;
            auto& orders = state.landmark_pose_orders[lm_id];
            for (int prev_ord : orders) {
                if (std::abs(prev_ord - pose_ord) > order_gap) {
                    return true;
                }
            }
            orders.push_back(pose_ord);
            break;
        }
        case EdgeType::GpsPrior:
        case EdgeType::GpsPrior3D:
            break;
        }
    }
    return false;
}

} // namespace islam
