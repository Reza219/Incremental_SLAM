#pragma once

#include "islam/deterministic_batch_ccolamd.hpp"

#include <Eigen/Sparse>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace islam {

// Standalone exact dynamic ordering engine for the repo-owned deterministic
// CCOLAMD rule. It never needs SuiteSparse as an oracle. On each structural
// refresh it certifies a reusable eliminated-prefix checkpoint directly on the
// current pattern, resumes from the materialized live fill graph, and falls back
// to a complete in-repo deterministic batch ordering only when no checkpoint can
// be certified.
//
// Exactness is defined relative to DeterministicBatchCcolamd. In this milestone
// that owned rule was moved toward SuiteSparse CCOLAMD semantics: constrained
// member ordering, SuiteSparse-style dense-column thresholds, deterministic
// minimum-degree pivoting, supercolumn handling, and aggressive absorption.
// Optional SuiteSparse compatibility audits remain the guard for any claim of
// byte-for-byte agreement on a particular SuiteSparse build/profile.
class DynamicExactCcolamdPrototype {
public:
    struct Stats {
        std::uint64_t refreshes = 0;
        std::uint64_t cold_starts = 0;
        std::uint64_t prefix_reuse_attempts = 0;
        std::uint64_t prefix_reuse_successes = 0;
        std::uint64_t total_common_prefix = 0;
        std::uint64_t total_rollback_suffix = 0;
        std::uint64_t exact_suffix_recomputes = 0;

        // Reference-free dynamic candidate telemetry. These counters measure
        // whether a checkpoint-resume candidate, computed before using the
        // current batch reference for certification, would have matched the
        // owned deterministic batch ordering. They are a development metric
        // for turning the prototype into a genuinely standalone dynamic engine.
        std::uint64_t reference_free_resume_attempts = 0;
        std::uint64_t reference_free_resume_candidates = 0;
        std::uint64_t reference_free_resume_failures = 0;
        std::uint64_t reference_free_resume_reference_matches = 0;
        std::uint64_t reference_free_resume_reference_mismatches = 0;

        // Dirty-boundary safety telemetry. This checks whether the conservative
        // dirty-var rollback boundary alone would have kept a prefix that is
        // still identical to the owned deterministic reference on the new
        // pattern. Any unsafe overestimate is a blocker for a future truly
        // reference-free dynamic ordering path.
        std::uint64_t dirty_boundary_safety_checks = 0;
        std::uint64_t dirty_boundary_safe_reuses = 0;
        std::uint64_t dirty_boundary_unsafe_overestimates = 0;
        std::uint64_t dirty_boundary_underestimates = 0;
        std::uint64_t dirty_boundary_exact_matches = 0;
        std::uint64_t dirty_boundary_candidate_pivots = 0;
        std::uint64_t dirty_boundary_safe_pivots = 0;
        std::uint64_t dirty_boundary_unsafe_pivots = 0;
        int last_dirty_boundary_overestimate = 0;
        int last_dirty_boundary_underestimate = 0;

        // Dirty-boundary checkpoint probe. Unlike the existing reusable
        // checkpoint, which is clamped by the true common prefix after the
        // current deterministic reference trace is available, this probe uses
        // only the dirty-boundary checkpoint from the previous trace and a
        // structural checkpoint certification on the current pattern. It is a
        // stricter development metric for reference-free dynamic ordering.
        std::uint64_t dirty_boundary_checkpoint_probe_attempts = 0;
        std::uint64_t dirty_boundary_checkpoint_probe_structural_successes = 0;
        std::uint64_t dirty_boundary_checkpoint_probe_structural_failures = 0;
        std::uint64_t dirty_boundary_checkpoint_probe_resume_matches = 0;
        std::uint64_t dirty_boundary_checkpoint_probe_resume_mismatches = 0;
        std::uint64_t dirty_boundary_checkpoint_probe_resume_failures = 0;
        int last_dirty_boundary_checkpoint_probe_pivots = 0;

        std::uint64_t checkpoint_replay_candidates = 0;
        std::uint64_t materialized_checkpoint_reuse_attempts = 0;
        std::uint64_t materialized_checkpoint_reuse_successes = 0;
        std::uint64_t materialized_checkpoint_reuse_failures = 0;
        std::uint64_t checkpoint_bank_entries = 0;
        std::uint64_t checkpoint_bank_insertions = 0;
        std::uint64_t checkpoint_bank_probes = 0;
        std::uint64_t checkpoint_bank_hits = 0;
        std::uint64_t checkpoint_bank_misses = 0;
        std::uint64_t checkpoint_bank_imported_entries = 0;
        std::uint64_t checkpoint_bank_import_duplicate_entries = 0;
        std::uint64_t checkpoint_bank_import_invalid_entries = 0;
        std::uint64_t checkpoint_bank_import_header_checks = 0;
        std::uint64_t checkpoint_bank_import_header_failures = 0;
        std::uint64_t checkpoint_bank_import_digest_mismatches = 0;
        std::uint64_t checkpoint_bank_import_entry_count_mismatches = 0;
        std::uint64_t checkpoint_bank_exported_entries = 0;
        std::uint64_t checkpoint_bank_export_writes = 0;
        std::uint64_t current_pattern_checkpoint_certification_attempts = 0;
        std::uint64_t current_pattern_checkpoint_certification_successes = 0;
        std::uint64_t current_pattern_checkpoint_certification_failures = 0;
        std::uint64_t structural_checkpoint_certification_attempts = 0;
        std::uint64_t structural_checkpoint_certification_successes = 0;
        std::uint64_t structural_checkpoint_certification_failures = 0;
        std::uint64_t checkpoint_resume_attempts = 0;
        std::uint64_t checkpoint_resume_successes = 0;
        std::uint64_t checkpoint_resume_failures = 0;
        std::uint64_t checkpoint_replay_full_prefixes = 0;
        std::uint64_t certified_suffix_replays = 0;
        std::uint64_t certified_suffix_replay_failures = 0;
        std::uint64_t suffix_replay_pivots_reused = 0;
        std::uint64_t suffix_replay_pivots_recomputed = 0;
        std::uint64_t exact_reference_checks = 0;
        std::uint64_t exact_reference_failures = 0;

        // Optional compatibility audit against SuiteSparse/CHOLMOD ordering
        // when SuiteSparse is available and ISLAM_DYNAMIC_EXACT_CHECK_SUITESPARSE=1.
        // This is deliberately not the dynamic engine. correctness oracle; the
        // dynamic backend is exact relative to the owned deterministic reference.
        // These counters measure how far that owned reference is from the
        // production SuiteSparse ordering used in conventional builds.
        std::uint64_t suitesparse_compat_checks = 0;
        std::uint64_t suitesparse_compat_failures = 0;
        std::uint64_t suitesparse_compat_unavailable = 0;
        std::uint64_t suitesparse_compat_mismatch_positions = 0;
        int last_suitesparse_compat_mismatches = 0;

        int last_common_prefix = 0;
        int last_dirty_boundary = 0;
        int last_reusable_checkpoint = 0;
        int last_checkpoint_live_vars = 0;
        int last_checkpoint_live_edges = 0;
        int last_rollback_suffix = 0;
    };

    using CheckpointBank = std::unordered_map<std::uint64_t, std::vector<DeterministicBatchCcolamd::LiveStateCheckpoint>>;

    explicit DynamicExactCcolamdPrototype(DeterministicBatchCcolamd::Options options = {})
        : batch_(options) {}

    [[nodiscard]] std::vector<int> refresh(const Eigen::SparseMatrix<double>& pattern,
                                           const std::vector<int>& dirty_vars = {});

    [[nodiscard]] const DeterministicBatchCcolamd::Trace& current_trace() const noexcept { return current_trace_; }
    [[nodiscard]] const std::vector<int>& last_reused_prefix() const noexcept { return last_reused_prefix_; }
    [[nodiscard]] const Stats& stats() const noexcept { return stats_; }
    [[nodiscard]] bool has_state() const noexcept { return has_state_; }

    void clear();

private:
    DeterministicBatchCcolamd batch_;
    DeterministicBatchCcolamd::Trace current_trace_{};
    std::vector<int> last_reused_prefix_{};
    Stats stats_{};
    CheckpointBank checkpoint_bank_{};
    bool has_state_ = false;
    bool checkpoint_bank_loaded_ = false;
};

using StandaloneExactDynamicCcolamd = DynamicExactCcolamdPrototype;

} // namespace islam
