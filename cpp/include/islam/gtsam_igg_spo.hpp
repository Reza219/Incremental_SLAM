#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstddef>
#include <map>
#include <set>
#include <vector>

#ifdef ISLAM_HAS_GTSAM
#include <gtsam/inference/Key.h>
#include <gtsam/inference/Ordering.h>
#include <gtsam/linear/VectorValues.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#endif

namespace islam {

#ifdef ISLAM_HAS_GTSAM

struct GtsamIggSpoParams {
    int max_gn_iterations = 10;
    double eta_threshold = 1.0;
    double dx_threshold = 1e-3;
    double selective_alpha = 0.30;
    double diagonal_damping = 1e-9;
    double sparse_drop_tolerance = 0.0;
    bool force_full_update = false;
    bool always_full_solve = false;
};

struct GtsamIggSpoIterationStats {
    int iteration = 0;
    std::size_t candidate_keys = 0;
    std::size_t active_keys = 0;
    std::size_t active_scalars = 0;
    bool fell_back_to_full = false;
    double dx_inf = 0.0;
};

struct GtsamIggSpoUpdateStats {
    int increment = 0;
    std::size_t num_factors = 0;
    std::size_t num_keys = 0;
    std::size_t tangent_dim = 0;
    bool global_gate = false;
    double eta = 0.0;
    double delta_eta = 0.0;
    int gn_iterations = 0;
    double final_error = 0.0;
    std::vector<GtsamIggSpoIterationStats> iterations;
};

class GtsamIggSpoSolver {
public:
    explicit GtsamIggSpoSolver(GtsamIggSpoParams params = {});

    void reset(const gtsam::NonlinearFactorGraph& graph, const gtsam::Values& values);

    GtsamIggSpoUpdateStats update(const gtsam::NonlinearFactorGraph& new_factors,
                                  const gtsam::Values& new_values = gtsam::Values());

    [[nodiscard]] const gtsam::Values& values() const noexcept { return values_; }
    [[nodiscard]] const gtsam::NonlinearFactorGraph& graph() const noexcept { return graph_; }
    [[nodiscard]] const GtsamIggSpoParams& params() const noexcept { return params_; }
    [[nodiscard]] int increments() const noexcept { return increment_; }

    struct KeyBlock {
        gtsam::Key key = 0;
        int offset = 0;
        int dim = 0;
    };

private:
    struct LinearizedSystem {
        gtsam::Ordering ordering;
        std::vector<KeyBlock> blocks;
        std::map<gtsam::Key, KeyBlock> by_key;
        Eigen::SparseMatrix<double> H;
        Eigen::VectorXd b;
        double eta = 0.0;
    };

    struct PartialSolveResult {
        Eigen::VectorXd dx;
        std::set<gtsam::Key> solved_keys;
        bool fell_back_to_full = false;
    };

    [[nodiscard]] LinearizedSystem linearize_current() const;
    [[nodiscard]] PartialSolveResult solve_active(const LinearizedSystem& sys,
                                                  const std::set<gtsam::Key>& active_keys) const;
    [[nodiscard]] std::set<gtsam::Key> all_keys() const;
    [[nodiscard]] std::set<gtsam::Key> keys_from_factors(const gtsam::NonlinearFactorGraph& factors) const;
    [[nodiscard]] std::set<gtsam::Key> prune_by_increment(const LinearizedSystem& sys,
                                                          const Eigen::VectorXd& dx,
                                                          const std::set<gtsam::Key>& keys) const;
    [[nodiscard]] std::set<gtsam::Key> expand_through_graph(const std::set<gtsam::Key>& keys) const;
    [[nodiscard]] gtsam::VectorValues vector_values_from_dense(const LinearizedSystem& sys,
                                                               const Eigen::VectorXd& dx,
                                                               const std::set<gtsam::Key>& keys) const;
    void insert_new_values(const gtsam::Values& new_values);

    GtsamIggSpoParams params_;
    gtsam::NonlinearFactorGraph graph_;
    gtsam::Values values_;
    double previous_eta_ = 0.0;
    std::size_t previous_dim_ = 0;
    bool have_previous_eta_ = false;
    int increment_ = 0;
};

#endif // ISLAM_HAS_GTSAM

} // namespace islam
