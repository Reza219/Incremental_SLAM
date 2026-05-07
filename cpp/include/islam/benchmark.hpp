#pragma once

#include "graph.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace islam {

enum class GatingMode {
    IGG,
    LCG,
    None,
};

enum class LinearSolverBackend {
    Auto,
    DenseEigen,
    Cholmod,
    CustomSparseExpanding,
};

std::string to_string(GatingMode mode);
GatingMode gating_mode_from_string(const std::string& text);
std::string to_string(LinearSolverBackend backend);
LinearSolverBackend linear_solver_backend_from_string(const std::string& text);

struct BatchStat {
    std::string backend_requested;
    std::string backend_actual;
    int batch = 0;
    int edge_id = -1;
    int state_size = 0;
    int num_edges = 0;
    int gn_iters = 0;
    int candidate_vars = 0;
    int active_perm = 0;
    int affected_vars = 0;
    bool fell_back_to_full = false;
    bool loop_closure = false;
    bool global_gate = false;
    double eta = 0.0;
    double delta_eta = 0.0;
    double dx_inf = 0.0;
    int scalar_measurements = 0;
    double normalized_chi2 = 0.0;
    double ate = 0.0;
    long long solve_flops = 0;
    long long full_solve_flops = 0;
    long long update_flops = 0;
    long long cumulative_solve_flops = 0;
    long long cumulative_full_solve_flops = 0;
    long long cumulative_update_flops = 0;
    long long factor_full_refactorizations = 0;
    long long factor_same_size_rank_updates = 0;
    long long factor_same_size_rank_downdates = 0;
    long long factor_cholmod_growth_refactorizations = 0;
    long long factor_custom_sparse_full_factorizations = 0;
    long long factor_custom_sparse_suffix_refactorizations = 0;
    long long factor_custom_sparse_growth_updates = 0;
    long long factor_custom_sparse_affected_closure_refactorizations = 0;
    long long factor_custom_sparse_structural_closure_refactorizations = 0;
    long long factor_custom_sparse_dependency_cache_hits = 0;
    int order_oracle_refreshes = 0;
    int order_oracle_cache_hits = 0;
    int order_oracle_cache_misses = 0;
    int order_oracle_cache_entries = 0;
    int order_local_pattern_cache_hits = 0;
    int order_local_pattern_cache_misses = 0;
    int order_local_pattern_cache_entries = 0;
    int order_motif_pattern_cache_hits = 0;
    int order_motif_pattern_cache_misses = 0;
    int order_motif_pattern_cache_entries = 0;
    int order_local_attempts = 0;
    int order_local_accepts = 0;
    int order_local_rejects = 0;
    int order_candidate_windows_generated = 0;
    int order_candidate_windows_tried = 0;
    int order_band_attempts = 0;
    int order_replay_attempts = 0;
    int order_replay_windows_cached = 0;
    int order_regime_switches = 0;
    int order_regime_creations = 0;
    int order_regime_merges = 0;
    int order_num_regimes_discovered = 0;
    int order_current_regime_id = -1;
    int order_adaptive_reorders = 0;
    int order_one_hop_accepts = 0;
    int order_two_hop_accepts = 0;
    int order_interval_accepts = 0;
    int order_union_accepts = 0;
    int order_band_accepts = 0;
    int order_replay_accepts = 0;
    int order_overlap_assembly_attempts = 0;
    int order_overlap_assembly_accepts = 0;
    int order_overlap_assembly_piece_hits = 0;
    int order_overlap_assembly_proposals = 0;
    int order_hierarchical_block_cache_hits = 0;
    int order_hierarchical_block_cache_misses = 0;
    int order_hierarchical_block_cache_entries = 0;
    int order_hierarchical_block_promotions = 0;
    int order_precedence_cache_hits = 0;
    int order_precedence_cache_misses = 0;
    int order_precedence_cache_entries = 0;
    int order_precedence_cache_promotions = 0;
    int order_precedence_guided_attempts = 0;
    int order_precedence_guided_accepts = 0;
    int order_precedence_consensus_attempts = 0;
    int order_precedence_consensus_accepts = 0;
    int order_precedence_consensus_scc_collapses = 0;
    int order_exact_output_certifications = 0;
    int order_certified_local_order_accepts = 0;
    int order_certified_exact_cache_order_accepts = 0;
    int order_certified_oracle_order_fallbacks = 0;
    int order_exact_cache_imported_entries = 0;
    int order_exact_cache_exported_entries = 0;
    int order_no_oracle_cache_hits = 0;
    int order_no_oracle_cache_misses = 0;
    double order_local_certification_rate = 0.0;
    double order_oracle_fallback_rate = 0.0;
    int etree_exact_recomputes = 0;
    int etree_exact_cache_imported_entries = 0;
    int etree_exact_cache_exported_entries = 0;
    int etree_exact_cache_hits = 0;
    int etree_exact_cache_misses = 0;
    int etree_no_exact_cache_hits = 0;
    int etree_no_exact_cache_misses = 0;
    int etree_local_update_attempts = 0;
    int etree_local_update_accepts = 0;
    int etree_local_update_fallbacks = 0;
    double etree_local_accept_rate = 0.0;
    int dyn_exact_refreshes = 0;
    int dyn_exact_cold_starts = 0;
    int dyn_exact_prefix_reuse_attempts = 0;
    int dyn_exact_prefix_reuse_successes = 0;
    int dyn_exact_common_prefix = 0;
    int dyn_exact_dirty_boundary = 0;
    int dyn_exact_reusable_checkpoint = 0;
    int dyn_exact_checkpoint_live_vars = 0;
    int dyn_exact_checkpoint_live_edges = 0;
    int dyn_exact_materialized_checkpoint_reuse_attempts = 0;
    int dyn_exact_materialized_checkpoint_reuse_successes = 0;
    int dyn_exact_materialized_checkpoint_reuse_failures = 0;
    int dyn_exact_checkpoint_bank_entries = 0;
    int dyn_exact_checkpoint_bank_insertions = 0;
    int dyn_exact_checkpoint_bank_probes = 0;
    int dyn_exact_checkpoint_bank_hits = 0;
    int dyn_exact_checkpoint_bank_misses = 0;
    int dyn_exact_checkpoint_bank_imported_entries = 0;
    int dyn_exact_checkpoint_bank_import_duplicate_entries = 0;
    int dyn_exact_checkpoint_bank_import_invalid_entries = 0;
    int dyn_exact_checkpoint_bank_import_header_checks = 0;
    int dyn_exact_checkpoint_bank_import_header_failures = 0;
    int dyn_exact_checkpoint_bank_import_digest_mismatches = 0;
    int dyn_exact_checkpoint_bank_import_entry_count_mismatches = 0;
    int dyn_exact_checkpoint_bank_exported_entries = 0;
    int dyn_exact_checkpoint_bank_export_writes = 0;
    int dyn_exact_current_pattern_checkpoint_certification_attempts = 0;
    int dyn_exact_current_pattern_checkpoint_certification_successes = 0;
    int dyn_exact_current_pattern_checkpoint_certification_failures = 0;
    int dyn_exact_structural_checkpoint_certification_attempts = 0;
    int dyn_exact_structural_checkpoint_certification_successes = 0;
    int dyn_exact_structural_checkpoint_certification_failures = 0;
    int dyn_exact_checkpoint_resume_attempts = 0;
    int dyn_exact_checkpoint_resume_successes = 0;
    int dyn_exact_checkpoint_resume_failures = 0;
    int dyn_exact_reference_free_resume_attempts = 0;
    int dyn_exact_reference_free_resume_candidates = 0;
    int dyn_exact_reference_free_resume_failures = 0;
    int dyn_exact_reference_free_resume_reference_matches = 0;
    int dyn_exact_reference_free_resume_reference_mismatches = 0;
    int dyn_exact_dirty_boundary_safety_checks = 0;
    int dyn_exact_dirty_boundary_safe_reuses = 0;
    int dyn_exact_dirty_boundary_unsafe_overestimates = 0;
    int dyn_exact_dirty_boundary_underestimates = 0;
    int dyn_exact_dirty_boundary_exact_matches = 0;
    int dyn_exact_dirty_boundary_candidate_pivots = 0;
    int dyn_exact_dirty_boundary_safe_pivots = 0;
    int dyn_exact_dirty_boundary_unsafe_pivots = 0;
    int dyn_exact_last_dirty_boundary_overestimate = 0;
    int dyn_exact_last_dirty_boundary_underestimate = 0;
    int dyn_exact_dirty_boundary_checkpoint_probe_attempts = 0;
    int dyn_exact_dirty_boundary_checkpoint_probe_structural_successes = 0;
    int dyn_exact_dirty_boundary_checkpoint_probe_structural_failures = 0;
    int dyn_exact_dirty_boundary_checkpoint_probe_resume_matches = 0;
    int dyn_exact_dirty_boundary_checkpoint_probe_resume_mismatches = 0;
    int dyn_exact_dirty_boundary_checkpoint_probe_resume_failures = 0;
    int dyn_exact_last_dirty_boundary_checkpoint_probe_pivots = 0;
    int dyn_exact_rollback_suffix = 0;
    int dyn_exact_suffix_replays = 0;
    int dyn_exact_suffix_replay_failures = 0;
    int dyn_exact_suffix_replay_pivots_reused = 0;
    int dyn_exact_suffix_replay_pivots_recomputed = 0;
    int dyn_exact_reference_checks = 0;
    int dyn_exact_reference_failures = 0;
    int dyn_exact_suitesparse_compat_checks = 0;
    int dyn_exact_suitesparse_compat_failures = 0;
    int dyn_exact_suitesparse_compat_unavailable = 0;
    int dyn_exact_suitesparse_compat_mismatch_positions = 0;
    int dyn_exact_last_suitesparse_compat_mismatches = 0;
};

std::vector<BatchStat> run_selective_benchmark(const Graph& full,
                                               int max_batches,
                                               int max_gn_iter,
                                               double dx_th = 1e-3,
                                               double anchor_strength = 1.0,
                                               double selective_alpha = 0.3,
                                               int lc_gap = 5,
                                               GatingMode gating_mode = GatingMode::IGG,
                                               bool use_spo = true,
                                               double eta_threshold = 1.0,
                                               LinearSolverBackend backend = LinearSolverBackend::Auto);

void write_benchmark_csv(const std::filesystem::path& path,
                         const std::vector<BatchStat>& stats);

} // namespace islam
