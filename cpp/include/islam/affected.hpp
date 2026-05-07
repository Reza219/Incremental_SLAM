#pragma once

#include "islam/graph.hpp"

#include <Eigen/Dense>

#include <vector>

namespace islam {

struct AffectedSet {
    std::vector<int> vars;
    std::vector<int> nodes;
    std::vector<int> edges;
};

AffectedSet find_affected(const Graph& g,
                          const Eigen::VectorXd& dx,
                          double threshold);

AffectedSet find_affected(const Graph& g,
                          const Eigen::VectorXd& dxa,
                          const std::vector<int>& candidate_vars,
                          double threshold);

} // namespace islam
