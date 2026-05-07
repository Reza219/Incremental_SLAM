#pragma once

#include "islam/factor_manager.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace islam {

struct SelectiveSolveResult {
    Eigen::VectorXd dx;
    std::vector<int> active_vars;
    std::vector<int> active_perm;
    bool fell_back_to_full = false;
    bool used_sparse_partial_solve = false;
    bool used_factor_block_solve = false;
    bool enlarged_to_factor_suffix = false;
};

class SelectiveSolver {
public:
    static SelectiveSolveResult solve_reduced(const FactorManager& fm,
                                              const std::vector<int>& affected_vars,
                                              double alpha = 0.3);

    static std::vector<int> elimination_tree(const Eigen::SparseMatrix<double>& Hperm);
    static std::vector<int> etree_reach(const std::vector<int>& parent,
                                        const std::vector<int>& seeds,
                                        int n);
};

} // namespace islam
