#pragma once

#include <Eigen/Sparse>

#include <cstdint>
#include <vector>

namespace islam {

// Standalone deterministic batch ordering engine that provides the repo-owned
// reference rule for the dynamic-exact symbolic backend.
//
// The current rule intentionally follows the SuiteSparse COLAMD/CCOLAMD design
// choices that matter for reproducible symbolic ordering: dense column
// postponement using SuiteSparse-style knobs, constrained member ordering,
// deterministic minimum-degree pivoting, supercolumn detection, and aggressive
// absorption.  The dynamic engine is mathematically exact with respect to this
// owned rule because every reused checkpoint is certified by replaying the same
// pivot/absorption rule on the current pattern before suffix replay.
//
// Scope note: this is a clean-room, auditable implementation of a
// SuiteSparse-compatible CCOLAMD-style rule for symmetric normal-equation
// patterns.  It is not claimed to be a vendored copy of SuiteSparse or a proven
// byte-for-byte replacement for every SuiteSparse CCOLAMD input/configuration
// unless the optional SuiteSparse compatibility audit is run and passes on the
// relevant matrices.
class DeterministicBatchCcolamd {
public:
    struct Options {
        // SuiteSparse COLAMD/CCOLAMD-style dense-column threshold.  When
        // use_suitesparse_dense_thresholds is true, columns with live degree
        // greater than max(dense_degree_absolute,
        // suitesparse_dense_col_knob * sqrt(n)) are postponed within their
        // constraint group.  SuiteSparse defaults use knob 10 and absolute floor
        // 16; these are the defaults below.  The legacy ratio threshold is kept
        // only for reproducibility experiments.
        int dense_degree_absolute = 16;
        double suitesparse_dense_col_knob = 10.0;
        double suitesparse_dense_row_knob = 10.0;
        bool use_suitesparse_dense_thresholds = true;
        double dense_degree_ratio = 0.10;

        // CCOLAMD-style behavior switches.
        bool aggressive_absorption = true;
        bool detect_supercolumns = true;
        bool respect_constraint_groups = true;
        bool defer_dense_columns = true;
    };

    struct LiveStateCheckpoint {
        // Materialized symbolic state after an elimination event. This is the
        // checkpoint substrate needed for standalone dynamic rollback/replay.
        int step_after = -1;
        std::uint64_t live_state_signature = 0;
        std::uint64_t prefix_signature = 0;
        std::vector<int> eliminated_prefix;
        std::vector<int> live_variables;
        // Upper-triangular live fill graph as flattened pairs [u0,v0,u1,v1,...].
        std::vector<int> live_edges_upper_flat;
    };

    struct StepRecord {
        int step = 0;
        int pivot = -1;
        int live_degree = 0;
        int constraint_group = 0;
        bool dense_postponed = false;
        std::uint64_t fill_edges_before = 0;
        std::uint64_t fill_edges_after = 0;

        // Checkpoint signatures make rollback decisions auditable. They are
        // structural hashes of the live symbolic state before/after the pivot
        // and of the eliminated prefix after the step. They are not numeric
        // factor hashes and are independent of matrix values.
        std::uint64_t live_state_signature_before = 0;
        std::uint64_t live_state_signature_after = 0;
        std::uint64_t prefix_signature_after = 0;

        std::vector<int> live_neighbors;
        std::vector<int> absorbed_columns;
        std::vector<int> eliminated_prefix_after;
    };

    struct Trace {
        int n = 0;
        std::uint64_t pattern_signature = 0;
        std::uint64_t rule_signature = 0;
        std::vector<StepRecord> steps;
        std::vector<int> permutation;

        // One prefix checkpoint per recorded elimination step. checkpoint[k]
        // corresponds to the state after steps[k]. This is intentionally
        // redundant with StepRecord::eliminated_prefix_after so that dynamic
        // replay code can access checkpoints without scanning records.
        std::vector<std::vector<int>> prefix_checkpoints;

        // Materialized live-state checkpoints aligned with prefix_checkpoints.
        std::vector<LiveStateCheckpoint> live_state_checkpoints;
    };

    struct Stats {
        int n = 0;
        std::uint64_t structural_nnz_upper = 0;
        std::uint64_t elimination_steps = 0;
        std::uint64_t fill_edges_inserted = 0;
        std::uint64_t absorbed_columns = 0;
        std::uint64_t supercolumns_absorbed = 0;
        std::uint64_t dense_columns_postponed = 0;
        std::uint64_t trace_steps_recorded = 0;
        std::uint64_t live_state_checkpoints_recorded = 0;
    };

    DeterministicBatchCcolamd() = default;
    explicit DeterministicBatchCcolamd(Options options) : options_(options) {}

    [[nodiscard]] std::vector<int> order(const Eigen::SparseMatrix<double>& pattern);
    [[nodiscard]] std::vector<int> order_with_constraints(const Eigen::SparseMatrix<double>& pattern,
                                                          const std::vector<int>& cmember);
    [[nodiscard]] Trace order_with_trace(const Eigen::SparseMatrix<double>& pattern);
    [[nodiscard]] Trace order_with_trace_and_constraints(const Eigen::SparseMatrix<double>& pattern,
                                                         const std::vector<int>& cmember);
    [[nodiscard]] Trace order_with_forced_prefix(const Eigen::SparseMatrix<double>& pattern,
                                                 const std::vector<int>& forced_prefix);
    [[nodiscard]] Trace order_with_forced_prefix_and_constraints(const Eigen::SparseMatrix<double>& pattern,
                                                                 const std::vector<int>& forced_prefix,
                                                                 const std::vector<int>& cmember);

    // Resume elimination directly from a materialized live-state checkpoint.
    [[nodiscard]] Trace order_from_checkpoint(const LiveStateCheckpoint& checkpoint);

    // Certify, without consulting a full batch ordering for the current pattern,
    // that the deterministic CCOLAMD-style elimination on `pattern` reaches
    // exactly `checkpoint`. When `prefix_trace` is non-null, it is filled with
    // the certified prefix trace only. This is the mathematical guard that lets
    // the dynamic engine reuse a previous live-state checkpoint safely.
    [[nodiscard]] bool certify_checkpoint_prefix(const Eigen::SparseMatrix<double>& pattern,
                                                 const LiveStateCheckpoint& checkpoint,
                                                 Trace* prefix_trace = nullptr) const;
    [[nodiscard]] const Trace& last_trace() const noexcept { return last_trace_; }
    [[nodiscard]] const Stats& stats() const noexcept { return stats_; }
    [[nodiscard]] const Options& options() const noexcept { return options_; }
    [[nodiscard]] std::uint64_t rule_signature() const noexcept;

    [[nodiscard]] static std::vector<int> order_default(const Eigen::SparseMatrix<double>& pattern) {
        DeterministicBatchCcolamd engine;
        return engine.order(pattern);
    }

    [[nodiscard]] static Trace trace_default(const Eigen::SparseMatrix<double>& pattern) {
        DeterministicBatchCcolamd engine;
        return engine.order_with_trace(pattern);
    }

private:
    Options options_;
    Stats stats_{};
    Trace last_trace_{};
};

// Lightweight deterministic structural hash used for trace/cache keys. It is
// intentionally independent of numeric values and considers only the symmetric
// sparsity pattern.
[[nodiscard]] std::uint64_t deterministic_ccolamd_pattern_signature(const Eigen::SparseMatrix<double>& pattern);

// Compare two traces and return the length of the common elimination prefix.
// This is the primitive needed by rollback/replay dynamic ordering.
[[nodiscard]] int common_trace_prefix_length(const DeterministicBatchCcolamd::Trace& a,
                                             const DeterministicBatchCcolamd::Trace& b);

// Verify that a materialized checkpoint is internally consistent with its
// corresponding trace prefix.
[[nodiscard]] bool checkpoint_matches_trace_prefix(const DeterministicBatchCcolamd::Trace& trace,
                                                   int checkpoint_prefix_length);

// Return the first step in a trace whose pivot, absorbed columns, or live
// neighbor set intersects dirty_vars. If no recorded step is touched, returns
// trace.steps.size(). This provides a conservative rollback boundary for the
// standalone exact dynamic CCOLAMD effort.
[[nodiscard]] int dirty_trace_boundary(const DeterministicBatchCcolamd::Trace& trace,
                                       const std::vector<int>& dirty_vars);

} // namespace islam
