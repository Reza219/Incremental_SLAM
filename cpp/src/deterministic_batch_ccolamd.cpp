#include "islam/deterministic_batch_ccolamd.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace islam {
namespace {

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30u)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27u)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31u);
}

std::uint64_t hash_combine(std::uint64_t h, std::uint64_t x) {
    return mix64(h ^ (x + 0x9e3779b97f4a7c15ull + (h << 6u) + (h >> 2u)));
}

std::vector<std::set<int>> symmetric_adjacency(const Eigen::SparseMatrix<double>& A,
                                               std::uint64_t& structural_nnz_upper) {
    if (A.rows() != A.cols()) {
        throw std::runtime_error("DeterministicBatchCcolamd requires a square pattern");
    }
    const int n = A.rows();
    std::vector<std::set<int>> adj(static_cast<size_t>(n));
    structural_nnz_upper = 0;
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int row = it.row();
            if (row < 0 || row >= n || col < 0 || col >= n) continue;
            if (row == col) {
                ++structural_nnz_upper;
                continue;
            }
            const int i = std::min(row, col);
            const int j = std::max(row, col);
            if (adj[static_cast<size_t>(i)].insert(j).second) {
                adj[static_cast<size_t>(j)].insert(i);
                ++structural_nnz_upper;
            }
        }
    }
    return adj;
}

int live_degree(const std::set<int>& neigh, const std::vector<unsigned char>& alive) {
    int d = 0;
    for (int v : neigh) {
        if (v >= 0 && v < static_cast<int>(alive.size()) && alive[static_cast<size_t>(v)]) ++d;
    }
    return d;
}

int dense_threshold_for_n(const DeterministicBatchCcolamd::Options& options, int n) {
    if (options.use_suitesparse_dense_thresholds) {
        const double knob = std::max(0.0, options.suitesparse_dense_col_knob);
        const int threshold = static_cast<int>(std::ceil(knob * std::sqrt(static_cast<double>(std::max(1, n)))));
        return std::max(options.dense_degree_absolute, threshold);
    }
    return std::max(options.dense_degree_absolute,
                    static_cast<int>(std::ceil(options.dense_degree_ratio * static_cast<double>(std::max(1, n)))));
}

int group_of(const std::vector<int>& cmember, int c, const DeterministicBatchCcolamd::Options& options) {
    if (!options.respect_constraint_groups) return 0;
    if (c < 0 || c >= static_cast<int>(cmember.size())) return 0;
    return std::max(0, cmember[static_cast<size_t>(c)]);
}

std::vector<int> canonical_cmember(int n, const std::vector<int>& cmember) {
    if (cmember.empty()) return std::vector<int>(static_cast<size_t>(n), 0);
    if (static_cast<int>(cmember.size()) != n) {
        throw std::runtime_error("DeterministicBatchCcolamd: cmember size must be zero or equal to n");
    }
    std::vector<int> out = cmember;
    for (int& v : out) v = std::max(0, v);
    return out;
}

std::uint64_t rule_signature_from_options(const DeterministicBatchCcolamd::Options& o) {
    std::uint64_t h = 0xCC01A4D5017E5005ull;
    h = hash_combine(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(o.dense_degree_absolute)));
    h = hash_combine(h, static_cast<std::uint64_t>(std::llround(o.suitesparse_dense_col_knob * 1000000.0)));
    h = hash_combine(h, static_cast<std::uint64_t>(std::llround(o.suitesparse_dense_row_knob * 1000000.0)));
    h = hash_combine(h, o.use_suitesparse_dense_thresholds ? 1ull : 0ull);
    h = hash_combine(h, static_cast<std::uint64_t>(std::llround(o.dense_degree_ratio * 1000000.0)));
    h = hash_combine(h, o.aggressive_absorption ? 1ull : 0ull);
    h = hash_combine(h, o.detect_supercolumns ? 1ull : 0ull);
    h = hash_combine(h, o.respect_constraint_groups ? 1ull : 0ull);
    h = hash_combine(h, o.defer_dense_columns ? 1ull : 0ull);
    return h;
}

std::uint64_t alive_state_signature(const std::vector<std::set<int>>& adj,
                                    const std::vector<unsigned char>& alive) {
    std::uint64_t h = mix64(static_cast<std::uint64_t>(alive.size()));
    for (int i = 0; i < static_cast<int>(alive.size()); ++i) {
        if (!alive[static_cast<size_t>(i)]) continue;
        h ^= mix64(0xA11CE00000000000ull ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)));
        for (int j : adj[static_cast<size_t>(i)]) {
            if (j <= i || j >= static_cast<int>(alive.size())) continue;
            if (!alive[static_cast<size_t>(j)]) continue;
            h ^= mix64((static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 32u)
                     | static_cast<std::uint32_t>(j));
        }
    }
    return h;
}

std::uint64_t prefix_signature(const std::vector<int>& prefix) {
    std::uint64_t h = mix64(static_cast<std::uint64_t>(prefix.size()));
    for (int v : prefix) {
        h = mix64(h ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
    }
    return h;
}

std::vector<int> live_variables_from_alive(const std::vector<unsigned char>& alive) {
    std::vector<int> vars;
    for (int i = 0; i < static_cast<int>(alive.size()); ++i) {
        if (alive[static_cast<size_t>(i)]) vars.push_back(i);
    }
    return vars;
}

std::vector<int> live_edges_flat_from_state(const std::vector<std::set<int>>& adj,
                                            const std::vector<unsigned char>& alive) {
    std::vector<int> flat;
    for (int i = 0; i < static_cast<int>(alive.size()); ++i) {
        if (!alive[static_cast<size_t>(i)]) continue;
        for (int j : adj[static_cast<size_t>(i)]) {
            if (j > i && j < static_cast<int>(alive.size()) && alive[static_cast<size_t>(j)]) {
                flat.push_back(i);
                flat.push_back(j);
            }
        }
    }
    return flat;
}

bool subset_of_live_neighbors(const std::set<int>& a,
                              const std::set<int>& b,
                              const std::vector<unsigned char>& alive,
                              int self_a,
                              int self_b) {
    for (int v : a) {
        if (v == self_b) continue;
        if (!alive[static_cast<size_t>(v)]) continue;
        if (v == self_a) continue;
        if (b.find(v) == b.end() && v != self_b) return false;
    }
    return true;
}

bool equal_live_neighbor_pattern(const std::set<int>& a,
                                 const std::set<int>& b,
                                 const std::vector<unsigned char>& alive,
                                 int self_a,
                                 int self_b) {
    return subset_of_live_neighbors(a, b, alive, self_a, self_b) &&
           subset_of_live_neighbors(b, a, alive, self_b, self_a);
}

struct PivotChoice {
    int column = -1;
    int degree = 0;
    int group = 0;
    bool dense = false;
};

bool pivot_less(const PivotChoice& a, const PivotChoice& b) {
    if (b.column < 0) return true;
    if (a.group != b.group) return a.group < b.group;
    if (a.dense != b.dense) return !a.dense && b.dense;
    if (a.degree != b.degree) return a.degree < b.degree;
    return a.column < b.column;
}

PivotChoice choose_pivot(const std::vector<std::set<int>>& adj,
                         const std::vector<unsigned char>& alive,
                         const std::vector<int>& cmember,
                         const DeterministicBatchCcolamd::Options& options,
                         int dense_threshold,
                         std::vector<unsigned char>& postponed_dense) {
    PivotChoice best;
    const int n = static_cast<int>(alive.size());
    for (int c = 0; c < n; ++c) {
        if (!alive[static_cast<size_t>(c)]) continue;
        const int deg = live_degree(adj[static_cast<size_t>(c)], alive);
        const bool dense = options.defer_dense_columns && deg > dense_threshold;
        if (dense) postponed_dense[static_cast<size_t>(c)] = 1;
        PivotChoice cand{c, deg, group_of(cmember, c, options), dense};
        if (pivot_less(cand, best)) best = cand;
    }
    return best;
}

std::vector<int> live_neighbors_of(const std::vector<std::set<int>>& adj,
                                   const std::vector<unsigned char>& alive,
                                   int c) {
    std::vector<int> neigh;
    for (int v : adj[static_cast<size_t>(c)]) {
        if (alive[static_cast<size_t>(v)]) neigh.push_back(v);
    }
    std::sort(neigh.begin(), neigh.end());
    return neigh;
}

std::vector<int> absorbed_columns_for_pivot(const std::vector<std::set<int>>& adj,
                                            const std::vector<unsigned char>& alive,
                                            const std::vector<int>& neigh,
                                            const std::vector<int>& cmember,
                                            const DeterministicBatchCcolamd::Options& options,
                                            int pivot) {
    std::vector<int> absorbed;
    for (int v : neigh) {
        if (!alive[static_cast<size_t>(v)]) continue;
        if (group_of(cmember, v, options) != group_of(cmember, pivot, options)) continue;
        bool absorb = false;
        if (options.detect_supercolumns) {
            absorb = equal_live_neighbor_pattern(adj[static_cast<size_t>(v)],
                                                 adj[static_cast<size_t>(pivot)],
                                                 alive,
                                                 v,
                                                 pivot);
        }
        if (!absorb && options.aggressive_absorption) {
            absorb = subset_of_live_neighbors(adj[static_cast<size_t>(v)],
                                              adj[static_cast<size_t>(pivot)],
                                              alive,
                                              v,
                                              pivot);
        }
        if (absorb) absorbed.push_back(v);
    }
    std::sort(absorbed.begin(), absorbed.end());
    absorbed.erase(std::unique(absorbed.begin(), absorbed.end()), absorbed.end());
    return absorbed;
}

void append_checkpoint(DeterministicBatchCcolamd::Trace& trace,
                       DeterministicBatchCcolamd::LiveStateCheckpoint ck,
                       const std::vector<std::set<int>>& adj,
                       const std::vector<unsigned char>& alive,
                       const std::vector<int>& perm,
                       int step_after) {
    ck.step_after = step_after;
    ck.live_state_signature = alive_state_signature(adj, alive);
    ck.prefix_signature = prefix_signature(perm);
    ck.eliminated_prefix = perm;
    ck.live_variables = live_variables_from_alive(alive);
    ck.live_edges_upper_flat = live_edges_flat_from_state(adj, alive);
    trace.prefix_checkpoints.push_back(perm);
    trace.live_state_checkpoints.push_back(std::move(ck));
}

} // namespace

std::uint64_t DeterministicBatchCcolamd::rule_signature() const noexcept {
    return rule_signature_from_options(options_);
}

std::uint64_t deterministic_ccolamd_pattern_signature(const Eigen::SparseMatrix<double>& pattern) {
    if (pattern.rows() != pattern.cols()) {
        throw std::runtime_error("deterministic_ccolamd_pattern_signature requires a square matrix");
    }
    std::uint64_t h = mix64(static_cast<std::uint64_t>(pattern.rows()));
    for (int col = 0; col < pattern.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) {
            const int row = it.row();
            if (row < 0 || row >= pattern.rows()) continue;
            const int i = std::min(row, col);
            const int j = std::max(row, col);
            h ^= mix64((static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 32u)
                     | static_cast<std::uint32_t>(j));
        }
    }
    return h;
}

int common_trace_prefix_length(const DeterministicBatchCcolamd::Trace& a,
                               const DeterministicBatchCcolamd::Trace& b) {
    if (a.rule_signature != 0 && b.rule_signature != 0 && a.rule_signature != b.rule_signature) return 0;
    const int n = std::min(static_cast<int>(a.steps.size()), static_cast<int>(b.steps.size()));
    int k = 0;
    for (; k < n; ++k) {
        const auto& sa = a.steps[static_cast<size_t>(k)];
        const auto& sb = b.steps[static_cast<size_t>(k)];
        if (sa.pivot != sb.pivot) break;
        if (sa.live_degree != sb.live_degree) break;
        if (sa.constraint_group != sb.constraint_group) break;
        if (sa.dense_postponed != sb.dense_postponed) break;
        if (sa.live_state_signature_before != sb.live_state_signature_before) break;
        if (sa.live_state_signature_after != sb.live_state_signature_after) break;
        if (sa.prefix_signature_after != sb.prefix_signature_after) break;
        if (sa.live_neighbors != sb.live_neighbors) break;
        if (sa.absorbed_columns != sb.absorbed_columns) break;
        if (sa.eliminated_prefix_after != sb.eliminated_prefix_after) break;
    }
    return k;
}

int dirty_trace_boundary(const DeterministicBatchCcolamd::Trace& trace,
                         const std::vector<int>& dirty_vars) {
    if (dirty_vars.empty()) return static_cast<int>(trace.steps.size());
    std::unordered_set<int> dirty(dirty_vars.begin(), dirty_vars.end());
    for (int k = 0; k < static_cast<int>(trace.steps.size()); ++k) {
        const auto& rec = trace.steps[static_cast<size_t>(k)];
        if (dirty.find(rec.pivot) != dirty.end()) return k;
        for (int v : rec.live_neighbors) {
            if (dirty.find(v) != dirty.end()) return k;
        }
        for (int v : rec.absorbed_columns) {
            if (dirty.find(v) != dirty.end()) return k;
        }
    }
    return static_cast<int>(trace.steps.size());
}

bool checkpoint_matches_trace_prefix(const DeterministicBatchCcolamd::Trace& trace,
                                     int checkpoint_prefix_length) {
    if (checkpoint_prefix_length <= 0) return true;
    const int k = checkpoint_prefix_length - 1;
    if (k < 0 || k >= static_cast<int>(trace.steps.size())) return false;
    if (k >= static_cast<int>(trace.live_state_checkpoints.size())) return false;
    const auto& step = trace.steps[static_cast<size_t>(k)];
    const auto& ck = trace.live_state_checkpoints[static_cast<size_t>(k)];
    return ck.step_after == k &&
           ck.live_state_signature == step.live_state_signature_after &&
           ck.prefix_signature == step.prefix_signature_after &&
           ck.eliminated_prefix == step.eliminated_prefix_after;
}

std::vector<int> DeterministicBatchCcolamd::order(const Eigen::SparseMatrix<double>& pattern) {
    return order_with_trace(pattern).permutation;
}

std::vector<int> DeterministicBatchCcolamd::order_with_constraints(const Eigen::SparseMatrix<double>& pattern,
                                                                   const std::vector<int>& cmember) {
    return order_with_trace_and_constraints(pattern, cmember).permutation;
}

DeterministicBatchCcolamd::Trace DeterministicBatchCcolamd::order_with_forced_prefix(
    const Eigen::SparseMatrix<double>& pattern,
    const std::vector<int>& forced_prefix) {
    return order_with_forced_prefix_and_constraints(pattern, forced_prefix, {});
}

DeterministicBatchCcolamd::Trace DeterministicBatchCcolamd::order_with_forced_prefix_and_constraints(
    const Eigen::SparseMatrix<double>& pattern,
    const std::vector<int>& forced_prefix,
    const std::vector<int>& cmember_in) {
    stats_ = Stats{};
    last_trace_ = Trace{};
    if (pattern.rows() != pattern.cols()) {
        throw std::runtime_error("DeterministicBatchCcolamd::order requires a square matrix");
    }
    const int n = pattern.rows();
    const auto cmember = canonical_cmember(n, cmember_in);
    stats_.n = n;
    last_trace_.n = n;
    last_trace_.pattern_signature = deterministic_ccolamd_pattern_signature(pattern);
    last_trace_.rule_signature = rule_signature();
    std::vector<int> perm;
    perm.reserve(static_cast<size_t>(n));
    if (n <= 0) {
        last_trace_.permutation = perm;
        return last_trace_;
    }

    auto adj = symmetric_adjacency(pattern, stats_.structural_nnz_upper);
    std::vector<unsigned char> alive(static_cast<size_t>(n), 1);
    std::vector<unsigned char> postponed_dense(static_cast<size_t>(n), 0);
    std::vector<unsigned char> forced_seen(static_cast<size_t>(n), 0);
    for (int v : forced_prefix) {
        if (v < 0 || v >= n) {
            throw std::runtime_error("DeterministicBatchCcolamd::order_with_forced_prefix: prefix index out of range");
        }
        if (forced_seen[static_cast<size_t>(v)]) {
            throw std::runtime_error("DeterministicBatchCcolamd::order_with_forced_prefix: duplicate prefix index");
        }
        forced_seen[static_cast<size_t>(v)] = 1;
    }

    const int dense_threshold = dense_threshold_for_n(options_, n);
    int alive_count = n;
    int forced_step = 0;
    while (alive_count > 0) {
        PivotChoice choice;
        bool using_forced = false;
        while (forced_step < static_cast<int>(forced_prefix.size()) &&
               !alive[static_cast<size_t>(forced_prefix[static_cast<size_t>(forced_step)])]) {
            ++forced_step;
        }
        if (forced_step < static_cast<int>(forced_prefix.size())) {
            const int best = forced_prefix[static_cast<size_t>(forced_step++)];
            const int deg = live_degree(adj[static_cast<size_t>(best)], alive);
            const bool dense = options_.defer_dense_columns && deg > dense_threshold;
            if (dense) postponed_dense[static_cast<size_t>(best)] = 1;
            choice = PivotChoice{best, deg, group_of(cmember, best, options_), dense};
            using_forced = true;
        } else {
            choice = choose_pivot(adj, alive, cmember, options_, dense_threshold, postponed_dense);
        }
        if (choice.column < 0) break;
        if (choice.dense) ++stats_.dense_columns_postponed;

        const int best = choice.column;
        const auto neigh = live_neighbors_of(adj, alive, best);
        const auto absorbed = absorbed_columns_for_pivot(adj, alive, neigh, cmember, options_, best);

        StepRecord rec;
        rec.step = static_cast<int>(last_trace_.steps.size());
        rec.pivot = best;
        rec.live_degree = choice.degree;
        rec.constraint_group = choice.group;
        rec.dense_postponed = choice.dense;
        rec.fill_edges_before = stats_.fill_edges_inserted;
        rec.live_state_signature_before = alive_state_signature(adj, alive);
        rec.live_neighbors = neigh;
        rec.absorbed_columns = absorbed;

        for (size_t i = 0; i < neigh.size(); ++i) {
            for (size_t j = i + 1; j < neigh.size(); ++j) {
                const int a = neigh[i];
                const int b = neigh[j];
                if (adj[static_cast<size_t>(a)].insert(b).second) {
                    adj[static_cast<size_t>(b)].insert(a);
                    ++stats_.fill_edges_inserted;
                }
            }
        }
        rec.fill_edges_after = stats_.fill_edges_inserted;
        last_trace_.steps.push_back(rec);

        auto eliminate = [&](int c, bool absorbed_col) {
            if (!alive[static_cast<size_t>(c)]) return;
            alive[static_cast<size_t>(c)] = 0;
            --alive_count;
            perm.push_back(c);
            ++stats_.elimination_steps;
            if (absorbed_col) ++stats_.absorbed_columns;
        };

        eliminate(best, false);
        for (int v : absorbed) {
            if (alive[static_cast<size_t>(v)]) {
                const bool is_super = options_.detect_supercolumns &&
                    equal_live_neighbor_pattern(adj[static_cast<size_t>(v)], adj[static_cast<size_t>(best)], alive, v, best);
                eliminate(v, true);
                if (is_super) ++stats_.supercolumns_absorbed;
            }
        }

        rec.live_state_signature_after = alive_state_signature(adj, alive);
        rec.prefix_signature_after = prefix_signature(perm);
        rec.eliminated_prefix_after = perm;
        last_trace_.steps.back() = rec;
        append_checkpoint(last_trace_, LiveStateCheckpoint{}, adj, alive, perm, rec.step);

        (void)using_forced;
    }

    if (static_cast<int>(perm.size()) != n) {
        for (int i = 0; i < n; ++i) {
            if (std::find(perm.begin(), perm.end(), i) == perm.end()) perm.push_back(i);
        }
    }
    last_trace_.permutation = perm;
    stats_.trace_steps_recorded = static_cast<std::uint64_t>(last_trace_.steps.size());
    stats_.live_state_checkpoints_recorded = static_cast<std::uint64_t>(last_trace_.live_state_checkpoints.size());
    return last_trace_;
}

bool DeterministicBatchCcolamd::certify_checkpoint_prefix(const Eigen::SparseMatrix<double>& pattern,
                                                          const DeterministicBatchCcolamd::LiveStateCheckpoint& checkpoint,
                                                          DeterministicBatchCcolamd::Trace* prefix_trace) const {
    if (prefix_trace != nullptr) *prefix_trace = Trace{};
    if (pattern.rows() != pattern.cols()) return false;
    const int n = pattern.rows();
    const int target_prefix_len = static_cast<int>(checkpoint.eliminated_prefix.size());
    if (target_prefix_len <= 0 || target_prefix_len > n) return false;
    if (checkpoint.step_after < 0) return false;

    std::vector<unsigned char> seen(static_cast<size_t>(n), 0);
    for (int v : checkpoint.eliminated_prefix) {
        if (v < 0 || v >= n || seen[static_cast<size_t>(v)]) return false;
        seen[static_cast<size_t>(v)] = 1;
    }
    std::vector<unsigned char> live_seen(static_cast<size_t>(n), 0);
    for (int v : checkpoint.live_variables) {
        if (v < 0 || v >= n || seen[static_cast<size_t>(v)] || live_seen[static_cast<size_t>(v)]) return false;
        live_seen[static_cast<size_t>(v)] = 1;
    }
    if (static_cast<int>(checkpoint.eliminated_prefix.size() + checkpoint.live_variables.size()) != n) return false;
    if (checkpoint.live_edges_upper_flat.size() % 2 != 0) return false;

    std::uint64_t structural_nnz_upper = 0;
    auto adj = symmetric_adjacency(pattern, structural_nnz_upper);
    std::vector<unsigned char> alive(static_cast<size_t>(n), 1);
    std::vector<unsigned char> postponed_dense(static_cast<size_t>(n), 0);
    std::vector<int> perm;
    perm.reserve(static_cast<size_t>(n));

    Trace trace;
    trace.n = n;
    trace.pattern_signature = deterministic_ccolamd_pattern_signature(pattern);
    trace.rule_signature = rule_signature();

    const auto cmember = canonical_cmember(n, {});
    const int dense_threshold = dense_threshold_for_n(options_, n);
    int alive_count = n;

    while (alive_count > 0 && static_cast<int>(perm.size()) < target_prefix_len) {
        const auto choice = choose_pivot(adj, alive, cmember, options_, dense_threshold, postponed_dense);
        if (choice.column < 0) return false;
        const int best = choice.column;
        const auto neigh = live_neighbors_of(adj, alive, best);
        const auto absorbed = absorbed_columns_for_pivot(adj, alive, neigh, cmember, options_, best);

        StepRecord rec;
        rec.step = static_cast<int>(trace.steps.size());
        rec.pivot = best;
        rec.live_degree = choice.degree;
        rec.constraint_group = choice.group;
        rec.dense_postponed = choice.dense;
        rec.live_state_signature_before = alive_state_signature(adj, alive);
        rec.live_neighbors = neigh;
        rec.absorbed_columns = absorbed;

        std::uint64_t fill_edges_inserted = trace.steps.empty() ? 0 : trace.steps.back().fill_edges_after;
        rec.fill_edges_before = fill_edges_inserted;
        for (size_t i = 0; i < neigh.size(); ++i) {
            for (size_t j = i + 1; j < neigh.size(); ++j) {
                const int a = neigh[i];
                const int b = neigh[j];
                if (adj[static_cast<size_t>(a)].insert(b).second) {
                    adj[static_cast<size_t>(b)].insert(a);
                    ++fill_edges_inserted;
                }
            }
        }
        rec.fill_edges_after = fill_edges_inserted;
        trace.steps.push_back(rec);

        auto eliminate = [&](int c) {
            if (!alive[static_cast<size_t>(c)]) return;
            alive[static_cast<size_t>(c)] = 0;
            --alive_count;
            perm.push_back(c);
        };
        eliminate(best);
        for (int v : absorbed) eliminate(v);

        // A reusable checkpoint must lie on an elimination-event boundary. If an
        // absorption event overshoots the requested prefix, the old checkpoint is
        // not a valid prefix for the current pattern.
        if (static_cast<int>(perm.size()) > target_prefix_len) return false;

        rec.live_state_signature_after = alive_state_signature(adj, alive);
        rec.prefix_signature_after = prefix_signature(perm);
        rec.eliminated_prefix_after = perm;
        trace.steps.back() = rec;
        append_checkpoint(trace, LiveStateCheckpoint{}, adj, alive, perm, rec.step);
    }

    if (static_cast<int>(perm.size()) != target_prefix_len) return false;
    if (perm != checkpoint.eliminated_prefix) return false;
    if (trace.live_state_checkpoints.empty()) return false;
    const auto& reached = trace.live_state_checkpoints.back();
    if (reached.step_after != checkpoint.step_after ||
        reached.live_state_signature != checkpoint.live_state_signature ||
        reached.prefix_signature != checkpoint.prefix_signature ||
        reached.eliminated_prefix != checkpoint.eliminated_prefix ||
        reached.live_variables != checkpoint.live_variables ||
        reached.live_edges_upper_flat != checkpoint.live_edges_upper_flat) {
        return false;
    }
    trace.permutation = perm;
    if (prefix_trace != nullptr) *prefix_trace = std::move(trace);
    return true;
}

DeterministicBatchCcolamd::Trace DeterministicBatchCcolamd::order_from_checkpoint(
    const DeterministicBatchCcolamd::LiveStateCheckpoint& checkpoint) {
    stats_ = Stats{};
    last_trace_ = Trace{};

    int n = 0;
    for (int v : checkpoint.eliminated_prefix) n = std::max(n, v + 1);
    for (int v : checkpoint.live_variables) n = std::max(n, v + 1);
    for (int v : checkpoint.live_edges_upper_flat) n = std::max(n, v + 1);
    if (static_cast<int>(checkpoint.eliminated_prefix.size() + checkpoint.live_variables.size()) > n) {
        n = static_cast<int>(checkpoint.eliminated_prefix.size() + checkpoint.live_variables.size());
    }

    stats_.n = n;
    last_trace_.n = n;
    last_trace_.pattern_signature = checkpoint.live_state_signature;
    last_trace_.rule_signature = rule_signature();

    std::vector<int> perm = checkpoint.eliminated_prefix;
    std::vector<unsigned char> alive(static_cast<size_t>(n), 0);
    for (int v : checkpoint.live_variables) {
        if (v < 0 || v >= n) throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: live variable out of range");
        alive[static_cast<size_t>(v)] = 1;
    }

    std::vector<unsigned char> seen_prefix(static_cast<size_t>(n), 0);
    for (int v : checkpoint.eliminated_prefix) {
        if (v < 0 || v >= n) throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: prefix variable out of range");
        if (seen_prefix[static_cast<size_t>(v)]) throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: duplicate prefix variable");
        if (alive[static_cast<size_t>(v)]) throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: prefix/live overlap");
        seen_prefix[static_cast<size_t>(v)] = 1;
    }

    std::vector<std::set<int>> adj(static_cast<size_t>(n));
    if (checkpoint.live_edges_upper_flat.size() % 2 != 0) {
        throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: malformed live-edge list");
    }
    for (size_t k = 0; k < checkpoint.live_edges_upper_flat.size(); k += 2) {
        const int a = checkpoint.live_edges_upper_flat[k];
        const int b = checkpoint.live_edges_upper_flat[k + 1];
        if (a < 0 || b < 0 || a >= n || b >= n || a == b) {
            throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: live edge out of range");
        }
        if (!alive[static_cast<size_t>(a)] || !alive[static_cast<size_t>(b)]) {
            throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: live edge touches non-live variable");
        }
        adj[static_cast<size_t>(a)].insert(b);
        adj[static_cast<size_t>(b)].insert(a);
        ++stats_.structural_nnz_upper;
    }

    if (alive_state_signature(adj, alive) != checkpoint.live_state_signature) {
        throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: live-state signature mismatch");
    }
    if (prefix_signature(perm) != checkpoint.prefix_signature) {
        throw std::runtime_error("DeterministicBatchCcolamd::order_from_checkpoint: prefix signature mismatch");
    }

    const auto cmember = canonical_cmember(n, {});
    std::vector<unsigned char> postponed_dense(static_cast<size_t>(n), 0);
    const int dense_threshold = dense_threshold_for_n(options_, n);
    int alive_count = static_cast<int>(checkpoint.live_variables.size());

    while (alive_count > 0) {
        const auto choice = choose_pivot(adj, alive, cmember, options_, dense_threshold, postponed_dense);
        if (choice.column < 0) break;
        if (choice.dense) ++stats_.dense_columns_postponed;
        const int best = choice.column;
        const auto neigh = live_neighbors_of(adj, alive, best);
        const auto absorbed = absorbed_columns_for_pivot(adj, alive, neigh, cmember, options_, best);

        StepRecord rec;
        rec.step = static_cast<int>(last_trace_.steps.size());
        rec.pivot = best;
        rec.live_degree = choice.degree;
        rec.constraint_group = choice.group;
        rec.dense_postponed = choice.dense;
        rec.fill_edges_before = stats_.fill_edges_inserted;
        rec.live_state_signature_before = alive_state_signature(adj, alive);
        rec.live_neighbors = neigh;
        rec.absorbed_columns = absorbed;

        for (size_t i = 0; i < neigh.size(); ++i) {
            for (size_t j = i + 1; j < neigh.size(); ++j) {
                const int a = neigh[i];
                const int b = neigh[j];
                if (adj[static_cast<size_t>(a)].insert(b).second) {
                    adj[static_cast<size_t>(b)].insert(a);
                    ++stats_.fill_edges_inserted;
                }
            }
        }
        rec.fill_edges_after = stats_.fill_edges_inserted;
        last_trace_.steps.push_back(rec);

        auto eliminate = [&](int c, bool absorbed_col) {
            if (!alive[static_cast<size_t>(c)]) return;
            alive[static_cast<size_t>(c)] = 0;
            --alive_count;
            perm.push_back(c);
            ++stats_.elimination_steps;
            if (absorbed_col) ++stats_.absorbed_columns;
        };
        eliminate(best, false);
        for (int v : absorbed) {
            if (alive[static_cast<size_t>(v)]) {
                const bool is_super = options_.detect_supercolumns &&
                    equal_live_neighbor_pattern(adj[static_cast<size_t>(v)], adj[static_cast<size_t>(best)], alive, v, best);
                eliminate(v, true);
                if (is_super) ++stats_.supercolumns_absorbed;
            }
        }

        rec.live_state_signature_after = alive_state_signature(adj, alive);
        rec.prefix_signature_after = prefix_signature(perm);
        rec.eliminated_prefix_after = perm;
        last_trace_.steps.back() = rec;
        LiveStateCheckpoint ck;
        append_checkpoint(last_trace_, ck, adj, alive, perm, checkpoint.step_after + 1 + rec.step);
    }

    if (static_cast<int>(perm.size()) != n) {
        for (int i = 0; i < n; ++i) {
            if (std::find(perm.begin(), perm.end(), i) == perm.end()) perm.push_back(i);
        }
    }
    last_trace_.permutation = perm;
    stats_.trace_steps_recorded = static_cast<std::uint64_t>(last_trace_.steps.size());
    stats_.live_state_checkpoints_recorded = static_cast<std::uint64_t>(last_trace_.live_state_checkpoints.size());
    return last_trace_;
}

DeterministicBatchCcolamd::Trace DeterministicBatchCcolamd::order_with_trace(const Eigen::SparseMatrix<double>& pattern) {
    return order_with_trace_and_constraints(pattern, {});
}

DeterministicBatchCcolamd::Trace DeterministicBatchCcolamd::order_with_trace_and_constraints(
    const Eigen::SparseMatrix<double>& pattern,
    const std::vector<int>& cmember) {
    static const std::vector<int> empty_prefix;
    return order_with_forced_prefix_and_constraints(pattern, empty_prefix, cmember);
}

} // namespace islam
