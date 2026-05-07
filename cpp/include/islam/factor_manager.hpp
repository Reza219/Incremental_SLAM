#pragma once

#include "islam/graph.hpp"
#include "islam/dynamic_ccolamd_engine.hpp"
#include "islam/sparse_expanding_cholesky.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace islam {

class SymbolicEngine;
struct SymbolicEngineStats;

struct EdgeContribution {
    int edge_id = -1;
    Eigen::SparseMatrix<double> C;   // Je^T
    Eigen::VectorXd r;
    std::vector<int> touched_vars;
    std::vector<int> touched_nodes;
};

struct NormalEquations {
    Eigen::SparseMatrix<double> H;
    Eigen::VectorXd g;
};

struct PermutedSystem {
    Eigen::SparseMatrix<double> H_perm;
    Eigen::VectorXd g_perm;
    std::vector<int> perm;
    std::vector<int> pinv;
};

struct SymbolicSystem {
    PermutedSystem permuted;
    std::vector<int> etree_parent;
};

struct FactorizationStats {
    std::uint64_t full_refactorizations = 0;
    std::uint64_t same_size_rank_updates = 0;
    std::uint64_t same_size_rank_downdates = 0;
    std::uint64_t dense_growth_extensions = 0;
    std::uint64_t cholmod_growth_refactorizations = 0;
    std::uint64_t custom_sparse_full_factorizations = 0;
    std::uint64_t custom_sparse_suffix_refactorizations = 0;
    std::uint64_t custom_sparse_growth_updates = 0;
    std::uint64_t custom_sparse_dynamic_reorders = 0;
    std::uint64_t custom_sparse_prefix_reuses = 0;
    std::uint64_t custom_sparse_prefix_columns_reused = 0;
    std::uint64_t custom_sparse_l21_recomputes = 0;
    std::uint64_t custom_sparse_left_looking_factorizations = 0;
    std::uint64_t custom_sparse_schur_complements = 0;
    std::uint64_t custom_sparse_triangular_solves = 0;
    std::uint64_t custom_sparse_jitter_regularizations = 0;
    std::uint64_t custom_sparse_affected_closure_refactorizations = 0;
    std::uint64_t custom_sparse_affected_closure_fallbacks = 0;
    std::uint64_t custom_sparse_affected_columns_refactored = 0;
    std::uint64_t custom_sparse_etree_closure_computations = 0;
    std::uint64_t custom_sparse_structural_pattern_changes = 0;
    std::uint64_t custom_sparse_numeric_only_updates = 0;
    std::uint64_t custom_sparse_symbolic_pattern_classifications = 0;
    std::uint64_t custom_sparse_structural_factor_pattern_changes = 0;
    std::uint64_t custom_sparse_structural_factor_columns_changed = 0;
    std::uint64_t custom_sparse_pattern_stable_structural_updates = 0;
    std::uint64_t custom_sparse_structural_closure_attempts = 0;
    std::uint64_t custom_sparse_structural_closure_refactorizations = 0;
    std::uint64_t custom_sparse_structural_closure_fallbacks = 0;
    std::uint64_t custom_sparse_affected_closure_certifications = 0;
    std::uint64_t custom_sparse_affected_closure_certification_failures = 0;
    std::uint64_t custom_sparse_factorization_residual_checks = 0;
    std::uint64_t custom_sparse_column_local_certifications = 0;
    std::uint64_t custom_sparse_column_local_certification_failures = 0;
    std::uint64_t custom_sparse_column_local_certification_columns = 0;
    std::uint64_t custom_sparse_full_certification_fallbacks = 0;
    std::uint64_t custom_sparse_dependency_cache_rebuilds = 0;
    std::uint64_t custom_sparse_dependency_cache_column_refreshes = 0;
    std::uint64_t custom_sparse_dependency_cache_invalidations = 0;
    std::uint64_t custom_sparse_dependency_cache_hits = 0;
    std::uint64_t custom_sparse_dependency_closure_computations = 0;
    std::uint64_t custom_sparse_dependency_closure_columns = 0;
    std::uint64_t custom_sparse_etree_closure_cache_bypasses = 0;
    double custom_sparse_last_certification_residual = 0.0;
    double custom_sparse_last_column_local_certification_residual = 0.0;
    double custom_sparse_max_column_local_certification_residual = 0.0;
    double custom_sparse_max_certification_residual = 0.0;
    std::uint64_t state_expansions = 0;
    std::uint64_t reserve_full_state_calls = 0;
    std::uint64_t incremental_edge_adds = 0;
    std::uint64_t edge_removes = 0;
    std::uint64_t edge_replaces = 0;
};

struct FactorBlockSystem {
    // Lower triangular maintained factor L in the factorization ordering.
    // CHOLMOD stores a factor of the permuted system A(perm,perm) = L L^T.
    // The fallback dense backend exposes its maintained dense L in natural order.
    Eigen::SparseMatrix<double> L_factor;
    Eigen::VectorXd g_factor;
    std::vector<int> perm;   // factor position -> original variable index
    std::vector<int> pinv;   // original variable index -> factor position
    std::vector<int> etree_parent;
    bool available = false;
    bool from_cholmod = false;
};

class FactorManager {
public:
    FactorManager();
    ~FactorManager();

    FactorManager(const FactorManager&) = delete;
    FactorManager& operator=(const FactorManager&) = delete;
    FactorManager(FactorManager&&) noexcept;
    FactorManager& operator=(FactorManager&&) noexcept;

    [[nodiscard]] bool using_cholmod() const noexcept { return using_cholmod_; }
    [[nodiscard]] bool using_sparse_expanding_cholesky() const noexcept { return use_sparse_expanding_backend_; }
    [[nodiscard]] int state_size() const noexcept { return state_size_; }
    [[nodiscard]] bool has_cache() const noexcept { return !edge_cache_.empty(); }
    [[nodiscard]] bool supports_incremental_updates() const noexcept { return supports_incremental_updates_; }
    [[nodiscard]] bool factor_covers_state() const noexcept { return factor_size_ == state_size_ && factorized_; }
    [[nodiscard]] double latent_prior_strength() const noexcept { return latent_prior_strength_; }

    [[nodiscard]] NormalEquations build_normal_equations(
        const Graph& g,
        const std::vector<int>* edge_ids = nullptr,
        double anchor_strength = 1.0) const;

    [[nodiscard]] EdgeContribution build_edge_contribution(const Graph& g, int edge_id) const;

    void clear();
    void enable_sparse_expanding_cholesky(bool enable = true);
    void force_dense_backend(bool enable = true);
    void configure_incremental(double anchor_strength = 1.0,
                               double latent_prior_strength = 0.0,
                               int anchor_dim = 3);
    void ensure_state_size(int n);
    void reserve_full_state(int n,
                            double anchor_strength = 1.0,
                            double latent_prior_strength = 1e-9,
                            int anchor_dim = 3);
    void activate_vars(const std::vector<int>& vars);

    void rebuild_from_graph(
        const Graph& g,
        const std::vector<int>* edge_ids = nullptr,
        double anchor_strength = 1.0);

    void add_edge_contribution(int edge_id, const EdgeContribution& contrib);
    void remove_edge_contribution(int edge_id);
    void replace_edge_contribution(int edge_id, const EdgeContribution& contrib);

    void factorize(const Eigen::SparseMatrix<double>& H);
    [[nodiscard]] Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const;
    [[nodiscard]] Eigen::VectorXd solve_cached() const;
    [[nodiscard]] std::vector<int> active_var_indices() const;
    [[nodiscard]] int active_state_size() const noexcept;
    [[nodiscard]] double maintained_factor_eta_full() const;
    [[nodiscard]] double information_eta(bool active_only = true) const;

    [[nodiscard]] Eigen::VectorXd solve_graph(
        const Graph& g,
        const std::vector<int>* edge_ids = nullptr,
        double anchor_strength = 1.0);

    [[nodiscard]] std::vector<int> get_permutation() const;
    [[nodiscard]] std::vector<int> get_elimination_tree() const;
    [[nodiscard]] PermutedSystem build_permuted_system() const;
    [[nodiscard]] SymbolicSystem build_symbolic_system() const;
    [[nodiscard]] FactorBlockSystem build_factor_block_system() const;

    [[nodiscard]] const Eigen::SparseMatrix<double>& last_H() const noexcept { return last_H_; }
    [[nodiscard]] const Eigen::VectorXd& last_g() const noexcept { return last_g_; }
    [[nodiscard]] const std::unordered_map<int, EdgeContribution>& edge_cache() const noexcept { return edge_cache_; }
    [[nodiscard]] const DynamicCcolamdEngine::Stats& dynamic_ordering_stats() const noexcept;
    [[nodiscard]] const SymbolicEngineStats& symbolic_engine_stats() const noexcept;
    [[nodiscard]] const FactorizationStats& factorization_stats() const noexcept { return factorization_stats_; }

private:
    struct CholmodState;

    void refactorize_current_system();
    void add_anchor_to_current_system();
    void update_cached_system_views();
    void invalidate_numeric_cache() const;
    void invalidate_structure_cache() const;
    void note_structural_change(const std::vector<int>& vars) const;
    [[nodiscard]] bool same_symbolic_pattern(const EdgeContribution& a,
                                             const EdgeContribution& b) const;
    void refresh_symbolic_cache_if_needed() const;
    void refresh_numeric_cache_if_needed() const;
    void apply_contribution_to_system(const EdgeContribution& contrib, bool add);
    void apply_full_size_update_to_factor(const EdgeContribution& contrib, bool add);
    void extend_factor_with_contribution(const EdgeContribution& contrib);
    [[nodiscard]] bool refresh_sparse_expanding_ordering_after_system_update(const std::vector<int>& dirty_vars);
    void apply_sparse_expanding_update_to_factor(const EdgeContribution& contrib, bool add);
    void sync_sparse_expanding_stats();
    void apply_diagonal_rank_update(const std::vector<int>& vars, double weight, bool add);
    [[nodiscard]] static Eigen::SparseMatrix<double> contribution_hessian(const EdgeContribution& contrib);

    int state_size_ = 0;
    int factor_size_ = 0;
    double anchor_strength_ = 1.0;
    double latent_prior_strength_ = 0.0;
    int anchor_dim_ = 3;
    bool factorized_ = false;
    bool using_cholmod_ = false;
    bool use_sparse_expanding_backend_ = false;
    bool supports_incremental_updates_ = true;
    bool cholmod_growth_refactor_pending_ = false;
    FactorizationStats factorization_stats_;

    Eigen::SparseMatrix<double> last_H_;
    Eigen::VectorXd last_g_;
    Eigen::SparseMatrix<double> H_current_;
    Eigen::VectorXd g_current_;
    std::unordered_map<int, EdgeContribution> edge_cache_;
    Eigen::MatrixXd L_dense_;
    std::vector<unsigned char> active_vars_;
    std::unique_ptr<CholmodState> cholmod_state_;
    std::unique_ptr<SparseExpandingCholesky> sparse_expanding_state_;
    std::vector<int> sparse_factor_perm_;
    std::vector<int> sparse_factor_pinv_;

    mutable std::uint64_t numeric_revision_ = 1;
    mutable std::uint64_t cached_numeric_revision_ = 0;
    mutable Eigen::SparseMatrix<double> cached_H_perm_;
    mutable Eigen::VectorXd cached_g_perm_;

    std::unique_ptr<SymbolicEngine> symbolic_engine_;

    void rebuild_symbolic_pattern_from_matrix(const Eigen::SparseMatrix<double>& H);
    void apply_contribution_to_symbolic_pattern(const EdgeContribution& contrib, int delta);
};

} // namespace islam
