#pragma once

#include "islam/graph.hpp"

namespace islam {

Graph make_pose_chain_graph(int num_poses,
                            double dx = 1.0,
                            double information_scale = 1.0,
                            bool add_loop_closure = false);

Graph make_pose_chain_with_periodic_loops(int num_poses,
                                          int loop_period,
                                          double dx = 1.0,
                                          double information_scale = 1.0);

} // namespace islam
