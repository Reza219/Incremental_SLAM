#pragma once

#include <Eigen/Sparse>

#include <cstdint>
#include <functional>
#include <vector>
#include <cstddef>

namespace islam {

class DynamicCcolamdEngine {
public:
    using OrderingOracle = std::function<std::vector<int>(const Eigen::SparseMatrix<double>&)>;

    struct Stats {
        std::uint64_t oracle_refreshes = 0;
        std::uint64_t oracle_cache_hits = 0;
        std::uint64_t oracle_cache_misses = 0;
        std::uint64_t oracle_cache_entries = 0;
        std::uint64_t local_pattern_cache_hits = 0;
        std::uint64_t local_pattern_cache_misses = 0;
        std::uint64_t local_pattern_cache_entries = 0;
        std::uint64_t motif_pattern_cache_hits = 0;
        std::uint64_t motif_pattern_cache_misses = 0;
        std::uint64_t motif_pattern_cache_entries = 0;
        std::uint64_t local_attempts = 0;
        std::uint64_t local_accepts = 0;
        std::uint64_t local_rejects = 0;
        std::uint64_t candidate_windows_generated = 0;
        std::uint64_t candidate_windows_tried = 0;
        std::uint64_t one_hop_attempts = 0;
        std::uint64_t two_hop_attempts = 0;
        std::uint64_t interval_attempts = 0;
        std::uint64_t union_attempts = 0;
        std::uint64_t band_attempts = 0;
        std::uint64_t replay_attempts = 0;
        std::uint64_t replay_windows_cached = 0;
        std::uint64_t regime_switches = 0;
        std::uint64_t regime_creations = 0;
        std::uint64_t regime_merges = 0;
        int current_regime_id = -1;
        int num_regimes_discovered = 0;
        std::uint64_t one_hop_accepts = 0;
        std::uint64_t two_hop_accepts = 0;
        std::uint64_t interval_accepts = 0;
        std::uint64_t union_accepts = 0;
        std::uint64_t band_accepts = 0;
        std::uint64_t replay_accepts = 0;
        std::uint64_t adaptive_reorders = 0;
        std::uint64_t overlap_assembly_attempts = 0;
        std::uint64_t overlap_assembly_accepts = 0;
        std::uint64_t overlap_assembly_piece_hits = 0;
        std::uint64_t overlap_assembly_proposals = 0;
        std::uint64_t hierarchical_block_cache_hits = 0;
        std::uint64_t hierarchical_block_cache_misses = 0;
        std::uint64_t hierarchical_block_cache_entries = 0;
        std::uint64_t hierarchical_block_promotions = 0;
        std::uint64_t precedence_cache_hits = 0;
        std::uint64_t precedence_cache_misses = 0;
        std::uint64_t precedence_cache_entries = 0;
        std::uint64_t precedence_cache_promotions = 0;
        std::uint64_t precedence_guided_attempts = 0;
        std::uint64_t precedence_guided_accepts = 0;
        std::uint64_t precedence_consensus_attempts = 0;
        std::uint64_t precedence_consensus_accepts = 0;
        std::uint64_t precedence_consensus_scc_collapses = 0;
        std::uint64_t exact_output_certifications = 0;
        std::uint64_t certified_local_order_accepts = 0;
        std::uint64_t certified_exact_cache_order_accepts = 0;
        std::uint64_t certified_oracle_order_fallbacks = 0;
        std::uint64_t exact_cache_imported_entries = 0;
        std::uint64_t exact_cache_exported_entries = 0;
        std::uint64_t no_oracle_cache_hits = 0;
        std::uint64_t no_oracle_cache_misses = 0;
    };

    void clear();
    void reserve_state(int n);
    void sync_to_exact_reference(const Eigen::SparseMatrix<double>& pattern,
                                 const std::vector<int>& exact_perm);

    [[nodiscard]] std::vector<int> refresh_exact(const Eigen::SparseMatrix<double>& pattern,
                                                 const std::vector<int>& dirty_vars,
                                                 const OrderingOracle& oracle);

    [[nodiscard]] const std::vector<int>& permutation() const noexcept { return perm_; }
    [[nodiscard]] const std::vector<int>& inverse_permutation() const noexcept { return pinv_; }
    [[nodiscard]] const Stats& stats() const noexcept { return stats_; }

private:
    enum class CandidateKind {
        Replay,
        OneHop,
        TwoHop,
        Interval,
        Union,
        Band,
    };


    enum class LocalOrderMethod {
        None,
        ExactCache,
        MotifCache,
        OverlapAssembly,
        PrecedenceCache,
        PrecedenceConsensus,
        GreedyApprox,
    };

    struct ColumnState {
        int degree = 0;
        int external_degree = 0;
        bool dense = false;
        bool alive = true;
    };


    struct LocalOrderResult {
        std::vector<int> ordered_vars;
        LocalOrderMethod method = LocalOrderMethod::None;
        std::uint64_t overlap_piece_hits = 0;
        std::uint64_t overlap_proposals = 0;
    };

    struct CandidateWindow {
        CandidateKind kind = CandidateKind::OneHop;
        std::vector<int> vars;
    };

    struct ReplayWindow {
        std::vector<int> vars;
        std::uint64_t hits = 0;
        std::uint64_t accepts = 0;
        std::uint64_t uses = 0;
    };

    struct KindStats {
        std::uint64_t attempts = 0;
        std::uint64_t accepts = 0;
        std::uint64_t last_success_tick = 0;
    };

    struct RegimeFeatures {
        double log_n = 0.0;
        double avg_degree = 0.0;
        double span_ratio = 0.0;
    };

    struct ExactCacheEntry {
        std::uint64_t key = 0;
        std::vector<int> perm;
        std::uint64_t last_use_tick = 0;
    };

    struct LocalPatternCacheEntry {
        std::uint64_t key = 0;
        std::vector<int> local_perm_idx;
        std::uint64_t hits = 0;
        std::uint64_t last_use_tick = 0;
    };

    struct MotifPatternCacheEntry {
        std::uint64_t key = 0;
        int size = 0;
        std::vector<int> local_perm_idx;
        std::uint64_t hits = 0;
        std::uint64_t last_use_tick = 0;
    };

    struct HierarchicalBlockCacheEntry {
        std::uint64_t key = 0;
        int size = 0;
        std::vector<int> local_perm_idx;
        std::uint64_t hits = 0;
        std::uint64_t last_use_tick = 0;
        double merit = 0.0;
    };

    struct PrecedenceCacheEntry {
        std::uint64_t key = 0;
        int size = 0;
        std::vector<double> pair_scores;
        std::uint64_t hits = 0;
        std::uint64_t last_use_tick = 0;
        double merit = 0.0;
    };

    struct RegimeState {
        int stable_id = -1;
        RegimeFeatures centroid{};
        RegimeFeatures m2{};
        std::vector<ReplayWindow> replay_windows;
        std::vector<KindStats> kind_stats;
        std::uint64_t visits = 0;
    };

    [[nodiscard]] static std::vector<int> inverse_permutation_of(const std::vector<int>& perm);
    [[nodiscard]] static std::vector<int> unique_sorted(std::vector<int> vals);
    [[nodiscard]] static Eigen::SparseMatrix<double> principal_submatrix(const Eigen::SparseMatrix<double>& A,
                                                                         const std::vector<int>& idx);
    [[nodiscard]] static std::vector<int> one_hop_neighborhood(const Eigen::SparseMatrix<double>& pattern,
                                                               const std::vector<int>& seeds);
    [[nodiscard]] static std::vector<std::vector<int>> adjacency_lists(const Eigen::SparseMatrix<double>& A);
    [[nodiscard]] static int dense_degree_threshold(int n);
    [[nodiscard]] static int score_column(const ColumnState& col);
    [[nodiscard]] static std::vector<int> positions_from_vars(const std::vector<int>& perm,
                                                              const std::vector<int>& vars);
    [[nodiscard]] static std::vector<int> interval_window_from_positions(const std::vector<int>& perm,
                                                                         const std::vector<int>& dirty_positions,
                                                                         int pad);
    [[nodiscard]] static std::vector<int> union_sorted_vectors(std::vector<int> a,
                                                               const std::vector<int>& b);
    [[nodiscard]] std::vector<ReplayWindow> replay_windows_from_dirty(const std::vector<int>& dirty_vars) const;
    [[nodiscard]] std::vector<int> band_window_from_dirty(const std::vector<int>& dirty_vars) const;
    [[nodiscard]] static std::vector<int> changed_vars_between_permutations(const std::vector<int>& before,
                                                                            const std::vector<int>& after);
    [[nodiscard]] static std::uint64_t pattern_signature(const Eigen::SparseMatrix<double>& pattern);
    [[nodiscard]] const ExactCacheEntry* lookup_exact_cache(std::uint64_t key) const;
    void store_exact_cache(std::uint64_t key, const std::vector<int>& perm);
    void maybe_load_exact_cache_from_env();
    void maybe_export_exact_cache_entry_to_env(std::uint64_t key, const std::vector<int>& perm);
    [[nodiscard]] const LocalPatternCacheEntry* lookup_local_pattern_cache(std::uint64_t key) const;
    void store_local_pattern_cache(std::uint64_t key, const std::vector<int>& local_perm_idx);
    [[nodiscard]] static std::uint64_t motif_signature(const Eigen::SparseMatrix<double>& pattern);
    [[nodiscard]] static std::uint64_t precedence_signature(const Eigen::SparseMatrix<double>& pattern);
    [[nodiscard]] const MotifPatternCacheEntry* lookup_motif_pattern_cache(std::uint64_t key, int size) const;
    void store_motif_pattern_cache(std::uint64_t key, int size, const std::vector<int>& local_perm_idx);
    [[nodiscard]] const HierarchicalBlockCacheEntry* lookup_hierarchical_block_cache(std::uint64_t key, int size) const;
    void store_hierarchical_block_cache(std::uint64_t key, int size, const std::vector<int>& local_perm_idx, double merit);
    [[nodiscard]] const PrecedenceCacheEntry* lookup_precedence_cache(std::uint64_t key, int size) const;
    void store_precedence_cache(std::uint64_t key, int size, const std::vector<int>& local_perm_idx, double merit);
    void promote_hierarchical_blocks_from_window(const Eigen::SparseMatrix<double>& pattern,
                                                 const std::vector<int>& local_vars,
                                                 const std::vector<int>& local_order_vars);
    [[nodiscard]] static RegimeFeatures extract_regime_features(const Eigen::SparseMatrix<double>& pattern);
    [[nodiscard]] static double regime_distance(const RegimeFeatures& a, const RegimeFeatures& b);
    [[nodiscard]] size_t select_or_create_regime(const RegimeFeatures& features);
    void update_regime_model(size_t regime_index, const RegimeFeatures& features);
    void maybe_merge_regimes();
    [[nodiscard]] std::vector<ReplayWindow>& active_replay_windows();
    [[nodiscard]] const std::vector<ReplayWindow>& active_replay_windows() const;
    [[nodiscard]] std::vector<KindStats>& active_kind_stats();
    [[nodiscard]] const std::vector<KindStats>& active_kind_stats() const;
    [[nodiscard]] double candidate_priority(const CandidateWindow& candidate,
                                            const std::vector<int>& dirty_vars) const;
    [[nodiscard]] std::vector<CandidateWindow> candidate_windows(const Eigen::SparseMatrix<double>& pattern,
                                                                 const std::vector<int>& dirty_vars) const;
    [[nodiscard]] LocalOrderResult order_candidate_window(const Eigen::SparseMatrix<double>& pattern,
                                                           const std::vector<int>& local_vars);
    [[nodiscard]] LocalOrderResult compose_overlapping_window_order(const Eigen::SparseMatrix<double>& local_pattern,
                                                                    const std::vector<int>& local_vars);
    [[nodiscard]] LocalOrderResult precedence_guided_order(const Eigen::SparseMatrix<double>& local_pattern,
                                                         const std::vector<int>& local_vars);
    [[nodiscard]] LocalOrderResult precedence_consensus_order(const Eigen::SparseMatrix<double>& local_pattern,
                                                              const std::vector<int>& local_vars);
    [[nodiscard]] static std::vector<int> stitch_window_order(const std::vector<int>& base_perm,
                                                              const std::vector<int>& local_vars,
                                                              const std::vector<int>& local_order_vars);
    void remember_replay_window(const std::vector<int>& vars);
    void rebuild_state_from_pattern(const Eigen::SparseMatrix<double>& pattern);

    int state_size_ = 0;
    static constexpr int kCandidateKindCount = 6;
    static constexpr size_t kMaxAdaptiveRegimes = 8;
    static constexpr size_t kMaxExactCacheEntries = 64;
    static constexpr size_t kMaxPrecedenceCacheEntries = 96;

    std::vector<int> perm_;
    std::vector<int> pinv_;
    std::vector<std::vector<int>> adjacency_;
    std::vector<RegimeState> regime_states_;
    std::vector<ExactCacheEntry> exact_cache_;
    std::vector<LocalPatternCacheEntry> local_pattern_cache_;
    std::vector<MotifPatternCacheEntry> motif_pattern_cache_;
    std::vector<HierarchicalBlockCacheEntry> hierarchical_block_cache_;
    std::vector<PrecedenceCacheEntry> precedence_cache_;
    size_t current_regime_index_ = 0;
    int next_regime_stable_id_ = 0;
    std::uint64_t tick_ = 0;
    bool exact_cache_env_loaded_ = false;
    Stats stats_;
};

} // namespace islam
