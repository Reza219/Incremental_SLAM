#pragma once

#include <Eigen/Sparse>
#include "dynamic_ccolamd_engine.hpp"
#include "dynamic_exact_ccolamd.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

namespace islam {


struct SymbolicSnapshot {
    std::vector<int> perm;
    std::vector<int> pinv;
    std::vector<int> etree_parent;
    Eigen::SparseMatrix<double> pattern_perm;
};

struct SymbolicEngineStats {
    std::uint64_t etree_exact_recomputes = 0;
    std::uint64_t etree_exact_cache_imported_entries = 0;
    std::uint64_t etree_exact_cache_exported_entries = 0;
    std::uint64_t etree_exact_cache_hits = 0;
    std::uint64_t etree_exact_cache_misses = 0;
    std::uint64_t etree_no_exact_cache_hits = 0;
    std::uint64_t etree_no_exact_cache_misses = 0;
    std::uint64_t etree_local_update_attempts = 0;
    std::uint64_t etree_local_update_accepts = 0;
    std::uint64_t etree_local_update_fallbacks = 0;
};

class SymbolicEngine {
public:
    using OrderingOracle = std::function<std::vector<int>(const Eigen::SparseMatrix<double>&)>;

    SymbolicEngine() = default;

    void clear();
    void reserve_state(int n);
    void rebuild_from_numeric_matrix(const Eigen::SparseMatrix<double>& H);
    void apply_contribution_pattern(const Eigen::SparseMatrix<double>& H_pattern,
                                    const std::vector<int>& touched_vars,
                                    int delta);
    void note_full_refresh();

    void refresh_if_needed(const OrderingOracle& ordering_oracle);

    [[nodiscard]] bool empty() const noexcept { return state_size_ == 0; }
    [[nodiscard]] int state_size() const noexcept { return state_size_; }
    [[nodiscard]] std::uint64_t structure_revision() const noexcept { return structure_revision_; }
    [[nodiscard]] const SymbolicSnapshot& snapshot() const noexcept { return snapshot_; }
    [[nodiscard]] const std::unordered_map<std::uint64_t, int>& pattern_counts() const noexcept { return pattern_counts_; }
    [[nodiscard]] const DynamicCcolamdEngine::Stats& dynamic_ordering_stats() const noexcept;
    [[nodiscard]] const DynamicExactCcolamdPrototype::Stats& dynamic_exact_ordering_stats() const noexcept;
    [[nodiscard]] const SymbolicEngineStats& stats() const noexcept { return stats_; }

private:
    static std::uint64_t pattern_key(int i, int j);
    static std::pair<int, int> decode_pattern_key(std::uint64_t key);
    static std::vector<int> inverse_permutation(const std::vector<int>& perm);
    static std::vector<int> unique_sorted(std::vector<int> vals);
    static Eigen::SparseMatrix<double> principal_submatrix(const Eigen::SparseMatrix<double>& A,
                                                           const std::vector<int>& idx);
    static Eigen::SparseMatrix<double> symmetric_permute_by_order(const Eigen::SparseMatrix<double>& A,
                                                                  const std::vector<int>& perm);
    static std::vector<int> elimination_tree_from_upper(const Eigen::SparseMatrix<double>& Hperm);
    static std::vector<int> expand_pattern_neighborhood(const Eigen::SparseMatrix<double>& H,
                                                        const std::vector<int>& seeds,
                                                        int hops);
    static std::vector<int> incremental_local_amd_permutation(const Eigen::SparseMatrix<double>& H,
                                                              const std::vector<int>& base_perm,
                                                              const std::vector<int>& dirty_vars);
    static std::vector<int> dirty_positions_from_vars(const std::vector<int>& perm,
                                                      const std::vector<int>& dirty_vars);
    static std::vector<int> elimination_tree_from_upper_suffix(const Eigen::SparseMatrix<double>& Hperm,
                                                               int start,
                                                               const std::vector<int>& previous_parent);
    static std::vector<int> incremental_local_etree_update(const Eigen::SparseMatrix<double>& Hperm,
                                                           const std::vector<int>& previous_parent,
                                                           const std::vector<int>& dirty_positions);

    [[nodiscard]] Eigen::SparseMatrix<double> build_pattern_matrix() const;
    void note_structural_change(const std::vector<int>& vars);

    int state_size_ = 0;
    std::uint64_t structure_revision_ = 1;
    std::uint64_t cached_structure_revision_ = 0;
    std::vector<unsigned char> dirty_mask_;
    std::vector<int> dirty_vars_;
    bool force_full_refresh_ = true;
    std::unordered_map<std::uint64_t, int> pattern_counts_;
    mutable bool etree_exact_cache_env_loaded_ = false;
    mutable std::unordered_map<std::uint64_t, std::vector<int>> etree_exact_cache_;
    SymbolicSnapshot snapshot_;
    std::unique_ptr<DynamicCcolamdEngine> dynamic_ordering_;
    std::unique_ptr<DynamicExactCcolamdPrototype> dynamic_exact_ordering_;
    SymbolicEngineStats stats_;
};

} // namespace islam
