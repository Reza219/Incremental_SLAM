#include "islam/synthetic.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <stdexcept>

namespace islam {
namespace {

Edge make_pose_pose_edge(int from_id, int to_id, int from_idx, int to_idx, double dx, double info_scale) {
    Edge e;
    e.type = EdgeType::PosePose;
    e.from_id = from_id;
    e.to_id = to_id;
    e.from_idx = from_idx;
    e.to_idx = to_idx;
    e.measurement = Eigen::Vector3d(dx, 0.0, 0.0);
    e.information = info_scale * Eigen::Matrix3d::Identity();
    return e;
}

} // namespace

Graph make_pose_chain_graph(int num_poses,
                            double dx,
                            double information_scale,
                            bool add_loop_closure) {
    if (num_poses < 2) {
        throw std::runtime_error("make_pose_chain_graph requires at least 2 poses");
    }

    Graph g;
    g.x = Eigen::VectorXd::Zero(3 * num_poses);
    g.var2node.resize(static_cast<size_t>(3 * num_poses));

    for (int i = 0; i < num_poses; ++i) {
        const int off = 3 * i;
        g.id_lookup[i] = NodeRef{off, 3};
        g.x.segment<3>(off) << dx * i, 0.0, 0.0;
        for (int k = 0; k < 3; ++k) g.var2node[static_cast<size_t>(off + k)] = i;
    }

    for (int i = 0; i < num_poses - 1; ++i) {
        const int from_idx = 3 * i;
        const int to_idx = 3 * (i + 1);
        g.edges.push_back(make_pose_pose_edge(i, i + 1, from_idx, to_idx, dx, information_scale));
    }

    if (add_loop_closure && num_poses >= 4) {
        g.edges.push_back(make_pose_pose_edge(num_poses - 1, 0, 3 * (num_poses - 1), 0,
                                              -dx * (num_poses - 1), information_scale));
    }

    g.node_edges.clear();
    for (int eid = 0; eid < static_cast<int>(g.edges.size()); ++eid) {
        const auto& e = g.edges[static_cast<size_t>(eid)];
        g.node_edges[e.from_id].push_back(eid);
        if (!e.is_unary()) g.node_edges[e.to_id].push_back(eid);
    }

    return g;
}

Graph make_pose_chain_with_periodic_loops(int num_poses,
                                          int loop_period,
                                          double dx,
                                          double information_scale) {
    if (loop_period < 2) {
        throw std::runtime_error("loop_period must be >= 2");
    }

    Graph g = make_pose_chain_graph(num_poses, dx, information_scale, false);

    for (int i = loop_period; i < num_poses; i += loop_period) {
        const int j = std::max(0, i - loop_period);
        if (i != j) {
            g.edges.push_back(make_pose_pose_edge(i, j, 3 * i, 3 * j,
                                                  -dx * (i - j), information_scale));
        }
    }

    g.node_edges.clear();
    for (int eid = 0; eid < static_cast<int>(g.edges.size()); ++eid) {
        const auto& e = g.edges[static_cast<size_t>(eid)];
        g.node_edges[e.from_id].push_back(eid);
        if (!e.is_unary()) g.node_edges[e.to_id].push_back(eid);
    }

    return g;
}

} // namespace islam
