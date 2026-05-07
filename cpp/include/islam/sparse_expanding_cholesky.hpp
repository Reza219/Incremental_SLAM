#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstdint>
#include <map>
#include <vector>

namespace islam {

struct SparseExpandingCholeskyOptions {
    double jitter = 1e-12;
    int max_jitter_tries = 6;
    double drop_tolerance = 0.0;
    bool certify_affected_closure_updates = true;
    // M75: first certify only the structurally affected residual columns.
    // This avoids forming L*L^T globally on every accepted local update.
    bool use_column_local_certification = true;
    // If the cheaper local check fails, run the older full residual check before
    // rejecting the update. This keeps the path conservative.
    bool full_certification_fallback = true;
    double certification_tolerance = 1e-8;
    // M76: reuse a maintained factor-dependency graph for numeric-only and
    // pattern-stable local updates. Structural factor-pattern changes still use
    // the fresh symbolic/etree path so newly introduced dependencies are not missed.
    bool use_factor_dependency_cache = true;
};

struct SparseExpandingCholeskyStats {
    std::uint64_t full_factorizations = 0;
    std::uint64_t suffix_refactorizations = 0;
    std::uint64_t same_size_suffix_updates = 0;
    std::uint64_t expansion_suffix_updates = 0;
    std::uint64_t grow_calls = 0;
    std::uint64_t diagonal_updates = 0;
    std::uint64_t dynamic_reorder_refactorizations = 0;
    std::uint64_t reordered_prefix_reuses = 0;
    std::uint64_t prefix_columns_reused = 0;
    std::uint64_t l21_recomputes = 0;
    std::uint64_t sparse_left_looking_factorizations = 0;
    std::uint64_t sparse_schur_complements = 0;
    std::uint64_t sparse_triangular_solves = 0;
    std::uint64_t jitter_regularizations = 0;
    std::uint64_t affected_closure_refactorizations = 0;
    std::uint64_t affected_closure_fallbacks = 0;
    std::uint64_t affected_columns_refactored = 0;
    std::uint64_t etree_closure_computations = 0;
    std::uint64_t structural_pattern_changes = 0;
    std::uint64_t numeric_only_updates = 0;
    std::uint64_t symbolic_pattern_classifications = 0;
    std::uint64_t structural_factor_pattern_changes = 0;
    std::uint64_t structural_factor_columns_changed = 0;
    std::uint64_t pattern_stable_structural_updates = 0;
    std::uint64_t structural_closure_attempts = 0;
    std::uint64_t structural_closure_refactorizations = 0;
    std::uint64_t structural_closure_fallbacks = 0;
    std::uint64_t affected_closure_certifications = 0;
    std::uint64_t affected_closure_certification_failures = 0;
    std::uint64_t factorization_residual_checks = 0;
    std::uint64_t column_local_certifications = 0;
    std::uint64_t column_local_certification_failures = 0;
    std::uint64_t column_local_certification_columns = 0;
    std::uint64_t full_certification_fallbacks = 0;
    std::uint64_t dependency_cache_rebuilds = 0;
    std::uint64_t dependency_cache_column_refreshes = 0;
    std::uint64_t dependency_cache_invalidations = 0;
    std::uint64_t dependency_cache_hits = 0;
    std::uint64_t dependency_closure_computations = 0;
    std::uint64_t dependency_closure_columns = 0;
    std::uint64_t etree_closure_cache_bypasses = 0;
    double last_certification_residual = 0.0;
    double last_column_local_certification_residual = 0.0;
    double max_column_local_certification_residual = 0.0;
    double max_certification_residual = 0.0;
    int last_suffix_start = -1;
    int last_suffix_size = 0;
    int last_reused_prefix = 0;
    int last_affected_closure_size = 0;
    int last_affected_closure_min = -1;
    int last_structural_factor_pattern_changed_columns = 0;
    int last_dependency_closure_size = 0;
};

class SparseExpandingCholesky {
public:
    explicit SparseExpandingCholesky(SparseExpandingCholeskyOptions options = {});

    void clear();
    void set_options(SparseExpandingCholeskyOptions options) { options_ = options; }
    [[nodiscard]] const SparseExpandingCholeskyOptions& options() const noexcept { return options_; }

    void factorize(const Eigen::SparseMatrix<double>& H);
    void refactorize_with_reused_prefix(const Eigen::SparseMatrix<double>& H, int reusable_prefix);
    void grow_to(int n);

    void apply_contribution(const Eigen::SparseMatrix<double>& C,
                            const std::vector<int>& touched_vars,
                            bool add);
    void apply_hessian_delta(const Eigen::SparseMatrix<double>& H_delta,
                             const std::vector<int>& touched_vars,
                             bool add);
    void apply_diagonal_update(const std::vector<int>& vars, double weight, bool add);

    [[nodiscard]] Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const;
    [[nodiscard]] double eta_logdet_half() const;

    [[nodiscard]] bool factorized() const noexcept { return factorized_; }
    [[nodiscard]] int state_size() const noexcept { return n_; }
    [[nodiscard]] int factor_size() const noexcept { return factor_size_; }
    [[nodiscard]] bool covers_state() const noexcept { return factorized_ && factor_size_ == n_; }

    [[nodiscard]] const Eigen::SparseMatrix<double>& normal_matrix() const noexcept { return H_; }
    [[nodiscard]] const Eigen::SparseMatrix<double>& lower_factor() const noexcept { return L_; }
    [[nodiscard]] const SparseExpandingCholeskyStats& stats() const noexcept { return stats_; }
    [[nodiscard]] double scaled_factorization_residual() const;
    [[nodiscard]] bool passes_factorization_residual_check(double* scaled_residual = nullptr) const;

private:
    [[nodiscard]] int suffix_start_from_vars(const std::vector<int>& vars) const;
    [[nodiscard]] int suffix_start_from_matrix(const Eigen::SparseMatrix<double>& A) const;
    [[nodiscard]] std::vector<int> touched_positions_from_vars_and_matrix(
        const std::vector<int>& vars,
        const Eigen::SparseMatrix<double>& A) const;
    [[nodiscard]] std::vector<int> elimination_tree_from_hessian() const;
    [[nodiscard]] std::vector<int> ancestor_closure_from_seeds(
        const std::vector<int>& seeds,
        const std::vector<int>& parent) const;
    void invalidate_factor_dependency_cache();
    void rebuild_factor_dependency_cache();
    void refresh_factor_dependency_cache_columns(const std::vector<int>& columns);
    void refresh_factor_dependency_cache_from(int first_column);
    [[nodiscard]] std::vector<int> factor_dependency_closure_from_seeds(
        const std::vector<int>& seeds);
    [[nodiscard]] std::vector<int> affected_closure_from_seeds(
        const std::vector<int>& seeds,
        bool allow_dependency_cache);
    [[nodiscard]] bool should_use_affected_closure(const std::vector<int>& closure, int suffix_start) const;
    [[nodiscard]] bool delta_changes_normal_pattern(const Eigen::SparseMatrix<double>& H_delta, bool add) const;
    [[nodiscard]] std::vector<std::vector<int>> lower_symbolic_pattern_from_hessian(
        const Eigen::SparseMatrix<double>& H) const;
    [[nodiscard]] std::vector<int> factor_pattern_changed_columns(
        const std::vector<std::vector<int>>& before,
        const std::vector<std::vector<int>>& after) const;
    [[nodiscard]] bool factorization_residual_within_tolerance(double* scaled_residual = nullptr) const;
    [[nodiscard]] std::vector<int> certification_columns_for_factor_columns(
        const std::vector<int>& factor_columns) const;
    [[nodiscard]] bool factorization_residual_columns_within_tolerance(
        const std::vector<int>& columns,
        double* scaled_residual = nullptr) const;
    [[nodiscard]] bool certify_factorization_for_columns(
        const std::vector<int>& factor_columns,
        double* scaled_residual = nullptr);
    [[nodiscard]] bool try_refactor_affected_closure(
        const std::vector<int>& seeds,
        bool allow_dependency_cache);
    [[nodiscard]] bool try_refactor_structural_closure(
        const std::vector<int>& seeds,
        const std::vector<int>& changed_factor_columns);
    void refactor_suffix_from(int suffix_start, bool count_as_full);

    [[nodiscard]] std::vector<std::map<int, double>> lower_columns_from_factor() const;
    [[nodiscard]] std::map<int, double> hessian_lower_column(int col) const;
    void rebuild_lower_from_columns(const std::vector<std::map<int, double>>& lower_cols);

    [[nodiscard]] Eigen::SparseMatrix<double> sparse_suffix_block(int suffix_start) const;
    [[nodiscard]] Eigen::SparseMatrix<double> sparse_a21_block(int suffix_start) const;
    [[nodiscard]] Eigen::SparseMatrix<double> sparse_l11_block(int suffix_start) const;
    [[nodiscard]] Eigen::SparseMatrix<double> solve_lower_prefix_sparse(
        int prefix_size,
        const Eigen::SparseMatrix<double>& rhs) const;
    [[nodiscard]] Eigen::SparseMatrix<double> factorize_sparse_left_looking(
        const Eigen::SparseMatrix<double>& A,
        std::uint64_t* jitter_regularizations) const;

    void rebuild_lower_with_sparse_suffix(int suffix_start,
                                          const Eigen::SparseMatrix<double>& L21,
                                          const Eigen::SparseMatrix<double>& L22);
    static Eigen::SparseMatrix<double> resize_square_preserve(const Eigen::SparseMatrix<double>& A, int n);

    SparseExpandingCholeskyOptions options_;
    int n_ = 0;
    int factor_size_ = 0;
    bool factorized_ = false;
    bool factor_dependency_cache_valid_ = false;
    Eigen::SparseMatrix<double> H_;
    Eigen::SparseMatrix<double> L_;
    std::vector<std::vector<int>> factor_dependency_children_;
    SparseExpandingCholeskyStats stats_;
};

} // namespace islam
