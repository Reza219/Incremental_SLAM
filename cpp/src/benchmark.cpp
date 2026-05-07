#include "islam/benchmark.hpp"

#include "islam/dynamic_exact_ccolamd.hpp"

#include "islam/affected.hpp"
#include "islam/factor_manager.hpp"
#include "islam/linearize.hpp"
#include "islam/selective_solver.hpp"
#include "islam/symbolic_engine.hpp"
#include "islam/update_graph.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace islam {
namespace {

void append_node_vars(const Graph& g, int node_id, std::vector<int>& vars) {
    const auto it = g.id_lookup.find(node_id);
    if (it == g.id_lookup.end()) {
        throw std::runtime_error("append_node_vars: missing node id");
    }
    for (int k = 0; k < it->second.dimension; ++k) {
        vars.push_back(it->second.offset + k);
    }
}

std::vector<int> touched_vars_from_new_edges(const Graph& g, int first_new_local_edge_idx) {
    std::vector<int> vars;
    for (int eid = first_new_local_edge_idx; eid < static_cast<int>(g.edges.size()); ++eid) {
        const auto& e = g.edges[static_cast<size_t>(eid)];
        append_node_vars(g, e.from_id, vars);
        if (!e.is_unary()) append_node_vars(g, e.to_id, vars);
    }
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
    return vars;
}

std::vector<int> all_vars(const Graph& g) {
    std::vector<int> out(static_cast<size_t>(g.state_size()));
    for (int i = 0; i < g.state_size(); ++i) out[static_cast<size_t>(i)] = i;
    return out;
}

int anchor_dim_for_graph(const Graph& g) {
    int best_offset = std::numeric_limits<int>::max();
    int best_dim = 3;
    for (const auto& kv : g.id_lookup) {
        if (kv.second.offset < best_offset) {
            best_offset = kv.second.offset;
            best_dim = kv.second.dimension;
        }
    }
    return std::clamp(best_dim, 1, std::max(1, g.state_size()));
}

std::vector<int> local_to_global_vars(const Graph& current,
                                      const Graph& full,
                                      const std::vector<int>& local_vars) {
    std::vector<int> global;
    global.reserve(local_vars.size());
    for (int lv : local_vars) {
        if (lv < 0 || lv >= current.state_size()) continue;
        const int nid = current.var2node.at(static_cast<size_t>(lv));
        const auto it_cur = current.id_lookup.find(nid);
        const auto it_full = full.id_lookup.find(nid);
        if (it_cur == current.id_lookup.end() || it_full == full.id_lookup.end()) {
            throw std::runtime_error("Failed to map local var to global var");
        }
        const int intra = lv - it_cur->second.offset;
        global.push_back(it_full->second.offset + intra);
    }
    std::sort(global.begin(), global.end());
    global.erase(std::unique(global.begin(), global.end()), global.end());
    return global;
}

std::vector<int> global_to_local_vars(const Graph& current,
                                      const Graph& full,
                                      const std::vector<int>& global_vars) {
    std::vector<int> local;
    local.reserve(global_vars.size());
    for (int gv : global_vars) {
        if (gv < 0 || gv >= full.state_size()) continue;
        const int nid = full.var2node.at(static_cast<size_t>(gv));
        const auto it_full = full.id_lookup.find(nid);
        const auto it_cur = current.id_lookup.find(nid);
        if (it_full == full.id_lookup.end()) {
            throw std::runtime_error("Failed to map global var to local var");
        }
        if (it_cur == current.id_lookup.end()) {
            continue;
        }
        const int intra = gv - it_full->second.offset;
        if (intra >= 0 && intra < it_cur->second.dimension) {
            local.push_back(it_cur->second.offset + intra);
        }
    }
    std::sort(local.begin(), local.end());
    local.erase(std::unique(local.begin(), local.end()), local.end());
    return local;
}

std::vector<int> global_vars_from_edge_ids(const Graph& full,
                                           const std::vector<int>& edge_ids) {
    std::vector<int> vars;
    for (int eid : edge_ids) {
        const auto& e = full.edges.at(static_cast<size_t>(eid));
        append_node_vars(full, e.from_id, vars);
        if (!e.is_unary()) append_node_vars(full, e.to_id, vars);
    }
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
    return vars;
}

Eigen::VectorXd full_dx_to_local(const Graph& current,
                                 const Graph& full,
                                 const Eigen::VectorXd& dx_full) {
    Eigen::VectorXd out = Eigen::VectorXd::Zero(current.state_size());
    for (const auto& [nid, ref_cur] : current.id_lookup) {
        const auto it_full = full.id_lookup.find(nid);
        if (it_full == full.id_lookup.end()) {
            throw std::runtime_error("Missing node in full graph during dx mapping");
        }
        out.segment(ref_cur.offset, ref_cur.dimension) = dx_full.segment(it_full->second.offset, ref_cur.dimension);
    }
    return out;
}

Eigen::VectorXd gather_by_index(const Eigen::VectorXd& x, const std::vector<int>& idx) {
    Eigen::VectorXd out(static_cast<int>(idx.size()));
    for (int i = 0; i < static_cast<int>(idx.size()); ++i) out[i] = x[idx[static_cast<size_t>(i)]];
    return out;
}

void subtract_on_indices(Eigen::VectorXd& x, const Eigen::VectorXd& dx, const std::vector<int>& idx) {
    for (int i : idx) {
        if (i < 0 || i >= x.size() || i >= dx.size()) {
            throw std::runtime_error("subtract_on_indices: index out of range");
        }
        x[i] -= dx[i];
    }
}

std::vector<int> merge_unique(std::vector<int> a, const std::vector<int>& b) {
    a.insert(a.end(), b.begin(), b.end());
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    return a;
}

std::vector<int> unique_sorted_copy(std::vector<int> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

void configure_factor_manager_backend(FactorManager& fm, LinearSolverBackend backend) {
    switch (backend) {
    case LinearSolverBackend::CustomSparseExpanding:
        fm.enable_sparse_expanding_cholesky(true);
        break;
    case LinearSolverBackend::DenseEigen:
        fm.force_dense_backend(true);
        break;
    case LinearSolverBackend::Cholmod:
        fm.enable_sparse_expanding_cholesky(false);
        break;
    case LinearSolverBackend::Auto:
    default:
        break;
    }
}

std::vector<int> factor_column_nnz_by_natural_var(const FactorManager& fm) {
    const FactorBlockSystem fbs = fm.build_factor_block_system();
    const int n = static_cast<int>(fbs.perm.size());
    std::vector<int> kappa(static_cast<size_t>(std::max(0, fm.state_size())), 1);
    if (!fbs.available || n <= 0 || fbs.L_factor.cols() != n || static_cast<int>(fbs.perm.size()) != n) {
        const auto& H = fm.last_H();
        if (H.cols() == static_cast<int>(kappa.size())) {
            for (int c = 0; c < H.outerSize(); ++c) {
                int nnz = 0;
                for (Eigen::SparseMatrix<double>::InnerIterator it(H, c); it; ++it) {
                    if (it.value() != 0.0) ++nnz;
                }
                if (c >= 0 && c < static_cast<int>(kappa.size())) kappa[static_cast<size_t>(c)] = std::max(1, nnz);
            }
        }
        return kappa;
    }

    std::fill(kappa.begin(), kappa.end(), 1);
    for (int pcol = 0; pcol < fbs.L_factor.outerSize(); ++pcol) {
        int nnz = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(fbs.L_factor, pcol); it; ++it) {
            if (it.value() != 0.0) ++nnz;
        }
        const int natural = fbs.perm[static_cast<size_t>(pcol)];
        if (natural >= 0 && natural < static_cast<int>(kappa.size())) {
            kappa[static_cast<size_t>(natural)] = std::max(1, nnz);
        }
    }
    return kappa;
}

long long saturating_ll_from_long_double(long double v) {
    if (!(v > 0.0L)) return 0;
    const long double lim = static_cast<long double>(std::numeric_limits<long long>::max());
    if (v >= lim) return std::numeric_limits<long long>::max();
    return static_cast<long long>(std::llround(v));
}

long long estimate_solve_flops_from_factor(const FactorManager& fm, const std::vector<int>& vars) {
    const auto kappa = factor_column_nnz_by_natural_var(fm);
    if (kappa.empty()) return 0;
    const std::vector<int> use = unique_sorted_copy(vars);
    long double acc = 0.0L;
    for (int v : use) {
        if (v >= 0 && v < static_cast<int>(kappa.size())) acc += 2.0L * static_cast<long double>(kappa[static_cast<size_t>(v)]);
    }
    return saturating_ll_from_long_double(acc);
}

long long estimate_full_solve_flops_from_factor(const FactorManager& fm) {
    const auto kappa = factor_column_nnz_by_natural_var(fm);
    long double acc = 0.0L;
    for (int k : kappa) acc += 2.0L * static_cast<long double>(k);
    return saturating_ll_from_long_double(acc);
}

long long estimate_update_flops_from_factor(const FactorManager& fm,
                                            const std::vector<int>& vars,
                                            bool add_only) {
    const auto kappa = factor_column_nnz_by_natural_var(fm);
    if (kappa.empty()) return 0;
    long double full = 0.0L;
    for (int k : kappa) full += static_cast<long double>(k) * static_cast<long double>(k);
    const std::vector<int> use = unique_sorted_copy(vars);
    long double part = 0.0L;
    for (int v : use) {
        if (v >= 0 && v < static_cast<int>(kappa.size())) {
            const long double kv = static_cast<long double>(kappa[static_cast<size_t>(v)]);
            part += kv * kv;
        }
    }
    const long double scaled_part = add_only ? part : 2.0L * part;
    const long double cap = add_only ? 0.5L * full : full;
    return saturating_ll_from_long_double(std::min(scaled_part, cap));
}

int scalar_measurement_count(const Graph& g) {
    int m = 0;
    for (const auto& edge : g.edges) m += static_cast<int>(edge.measurement.size());
    return m;
}

double normalized_chi_squared(const Graph& g) {
    const int m = scalar_measurement_count(g);
    if (m <= 0) return std::numeric_limits<double>::quiet_NaN();
    long double sum_sq = 0.0L;
    for (int eid = 0; eid < static_cast<int>(g.edges.size()); ++eid) {
        const auto lin = jacobian_edge_jr(g.edges[static_cast<size_t>(eid)], g);
        sum_sq += static_cast<long double>(lin.r.squaredNorm());
    }
    return static_cast<double>(sum_sq / static_cast<long double>(m));
}

int trajectory_position_dim(const Graph& g) {
    for (const auto& kv : g.id_lookup) {
        if (kv.second.dimension == 3) return 2;
        if (kv.second.dimension == 6) return 3;
    }
    return 0;
}

Eigen::VectorXd node_position(const Graph& g, int node_id, int pos_dim) {
    const auto it = g.id_lookup.find(node_id);
    if (it == g.id_lookup.end() || it->second.dimension < pos_dim) return Eigen::VectorXd();
    return g.x.segment(it->second.offset, pos_dim);
}

double aligned_trajectory_rmse(const Graph& estimate, const Graph& reference) {
    const int dim = trajectory_position_dim(estimate);
    if (dim != 2 && dim != 3) return std::numeric_limits<double>::quiet_NaN();
    std::vector<Eigen::VectorXd> x_pts;
    std::vector<Eigen::VectorXd> y_pts;
    for (const auto& kv : estimate.id_lookup) {
        const int nid = kv.first;
        const int node_dim = kv.second.dimension;
        if (!((dim == 2 && node_dim == 3) || (dim == 3 && node_dim == 6))) continue;
        const auto it_ref = reference.id_lookup.find(nid);
        if (it_ref == reference.id_lookup.end() || it_ref->second.dimension != node_dim) continue;
        x_pts.push_back(node_position(estimate, nid, dim));
        y_pts.push_back(node_position(reference, nid, dim));
    }
    const int n = static_cast<int>(x_pts.size());
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    if (n == 1) return 0.0;
    Eigen::MatrixXd X(dim, n);
    Eigen::MatrixXd Y(dim, n);
    for (int i = 0; i < n; ++i) { X.col(i) = x_pts[static_cast<size_t>(i)]; Y.col(i) = y_pts[static_cast<size_t>(i)]; }
    const Eigen::VectorXd mx = X.rowwise().mean();
    const Eigen::VectorXd my = Y.rowwise().mean();
    X.colwise() -= mx;
    Y.colwise() -= my;
    const Eigen::MatrixXd H = X * Y.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::MatrixXd R = V * U.transpose();
    if (R.determinant() < 0.0) { V.col(dim - 1) *= -1.0; R = V * U.transpose(); }
    long double err = 0.0L;
    for (int i = 0; i < n; ++i) {
        const Eigen::VectorXd aligned = R * x_pts[static_cast<size_t>(i)] + my - R * mx;
        err += static_cast<long double>((y_pts[static_cast<size_t>(i)] - aligned).squaredNorm());
    }
    return std::sqrt(static_cast<double>(err / static_cast<long double>(n)));
}

Graph optimize_batch_reference(const Graph& full,
                               int max_gn_iter,
                               double dx_th,
                               double anchor_strength,
                               LinearSolverBackend backend) {
    Graph ref = full;
    if (ref.edges.empty() || ref.state_size() == 0) return ref;
    FactorManager fm;
    configure_factor_manager_backend(fm, backend);
    fm.configure_incremental(anchor_strength, 0.0, anchor_dim_for_graph(ref));
    for (int iter = 0; iter < std::max(1, max_gn_iter); ++iter) {
        fm.rebuild_from_graph(ref, nullptr, anchor_strength);
        const Eigen::VectorXd dx = fm.solve_cached();
        const double dx_inf = dx.size() > 0 ? dx.cwiseAbs().maxCoeff() : 0.0;
        if (dx_inf <= dx_th) break;
        ref.x -= dx;
    }
    return ref;
}

struct GatingDecision {
    std::vector<int> candidate_local_vars;
    bool global_gate = false;
    double eta = 0.0;
    double delta_eta = 0.0;
    int eta_state_size = 0;
};

GatingDecision decide_gating(const FactorManager& fm,
                             const Graph& current,
                             int first_new_local_edge_idx,
                             const UpdateGraphResult& upd,
                             GatingMode mode,
                             double eta_threshold,
                             double prev_eta,
                             int prev_state_size,
                             bool has_prev_eta) {
    GatingDecision out;
    out.eta_state_size = fm.active_state_size();
    out.eta = fm.information_eta(true);
    if (has_prev_eta && prev_state_size > 0 && out.eta_state_size > 0) {
        out.delta_eta = out.eta - (static_cast<double>(prev_state_size) / static_cast<double>(out.eta_state_size)) * prev_eta;
    } else {
        out.delta_eta = 0.0;
    }

    const std::vector<int> local_touched = merge_unique(
        touched_vars_from_new_edges(current, first_new_local_edge_idx), upd.new_vars);

    switch (mode) {
    case GatingMode::None:
        out.global_gate = true;
        out.candidate_local_vars = all_vars(current);
        break;
    case GatingMode::LCG:
        out.global_gate = upd.loop_closure;
        out.candidate_local_vars = out.global_gate ? all_vars(current) : local_touched;
        break;
    case GatingMode::IGG:
    default:
        out.global_gate = out.delta_eta >= eta_threshold;
        out.candidate_local_vars = out.global_gate ? all_vars(current) : local_touched;
        break;
    }

    if (out.candidate_local_vars.empty()) {
        out.candidate_local_vars = all_vars(current);
        out.global_gate = true;
    }
    return out;
}

} // namespace

std::string to_string(GatingMode mode) {
    switch (mode) {
    case GatingMode::IGG:
        return "IGG";
    case GatingMode::LCG:
        return "LCG";
    case GatingMode::None:
    default:
        return "none";
    }
}

GatingMode gating_mode_from_string(const std::string& text) {
    std::string s;
    s.reserve(text.size());
    for (unsigned char ch : text) s.push_back(static_cast<char>(std::tolower(ch)));
    if (s == "igg") return GatingMode::IGG;
    if (s == "lcg") return GatingMode::LCG;
    if (s == "none" || s == "off") return GatingMode::None;
    throw std::runtime_error("Unknown gating mode: " + text + " (expected IGG, LCG, or none)");
}

std::string to_string(LinearSolverBackend backend) {
    switch (backend) {
    case LinearSolverBackend::Auto: return "auto";
    case LinearSolverBackend::DenseEigen: return "dense-eigen";
    case LinearSolverBackend::Cholmod: return "cholmod";
    case LinearSolverBackend::CustomSparseExpanding: return "custom-sparse-expanding";
    default: return "auto";
    }
}

LinearSolverBackend linear_solver_backend_from_string(const std::string& text) {
    std::string s;
    s.reserve(text.size());
    for (unsigned char ch : text) {
        if (ch == '_' || ch == '-') s.push_back('-');
        else s.push_back(static_cast<char>(std::tolower(ch)));
    }
    if (s == "auto" || s == "default") return LinearSolverBackend::Auto;
    if (s == "dense" || s == "eigen" || s == "dense-eigen") return LinearSolverBackend::DenseEigen;
    if (s == "cholmod" || s == "suitesparse") return LinearSolverBackend::Cholmod;
    if (s == "custom" || s == "custom-sparse" || s == "sparse-expanding" || s == "custom-sparse-expanding") {
        return LinearSolverBackend::CustomSparseExpanding;
    }
    throw std::runtime_error("Unknown linear solver backend: " + text +
                             " (expected auto, dense-eigen, cholmod, or custom-sparse-expanding)");
}

std::vector<BatchStat> run_selective_benchmark(const Graph& full,
                                               int max_batches,
                                               int max_gn_iter,
                                               double dx_th,
                                               double anchor_strength,
                                               double selective_alpha,
                                               int lc_gap,
                                               GatingMode gating_mode,
                                               bool use_spo,
                                               double eta_threshold,
                                               LinearSolverBackend backend) {
    std::vector<BatchStat> stats;
    if (full.edges.empty()) return stats;

    FactorManager fm;
    configure_factor_manager_backend(fm, backend);
    fm.configure_incremental(anchor_strength, 0.0, anchor_dim_for_graph(full));
    const std::string backend_requested = to_string(backend);
    const auto actual_backend_name = [&fm]() {
        if (fm.using_sparse_expanding_cholesky()) return std::string("custom-sparse-expanding");
        if (fm.using_cholmod()) return std::string("cholmod");
        return std::string("dense-eigen");
    };

    Graph current;
    bool has_current = false;
    LoopClosureState lc_state;
    double prev_eta = std::numeric_limits<double>::quiet_NaN();
    int prev_state_size = 0;
    bool has_prev_eta = false;
    DynamicCcolamdEngine::Stats prev_order_stats{};
    SymbolicEngineStats prev_symbolic_stats{};
    DynamicExactCcolamdPrototype dyn_exact;
    DynamicExactCcolamdPrototype::Stats prev_dyn_exact_stats{};
    FactorizationStats prev_factor_stats{};
    const Graph reference = optimize_batch_reference(
        full, std::max(max_gn_iter, 50), std::min(dx_th, 1e-9), anchor_strength, backend);
    long long cumulative_update_flops = 0;
    long long cumulative_solve_flops = 0;
    long long cumulative_full_solve_flops = 0;

    const int batches = std::min(max_batches, static_cast<int>(full.edges.size()));
    for (int b = 0; b < batches; ++b) {
        const auto upd = update_graph(std::vector<int>{b}, full, has_current ? &current : nullptr, &lc_state, lc_gap);
        current = upd.current;
        has_current = true;

        const int num_new_edges = 1;
        const int first_new_local_edge_idx = static_cast<int>(current.edges.size()) - num_new_edges;

        fm.ensure_state_size(current.state_size());
        const auto contrib = fm.build_edge_contribution(current, first_new_local_edge_idx);
        const std::vector<int> new_edge_touched_vars = contrib.touched_vars;
        fm.add_edge_contribution(b, contrib);

        const auto gate = decide_gating(
            fm, current, first_new_local_edge_idx, upd, gating_mode, eta_threshold, prev_eta, prev_state_size, has_prev_eta);
        prev_eta = gate.eta;
        prev_state_size = gate.eta_state_size;
        has_prev_eta = true;

        BatchStat st;
        st.backend_requested = backend_requested;
        st.backend_actual = actual_backend_name();
        st.batch = b;
        st.edge_id = b;
        st.state_size = current.state_size();
        st.num_edges = static_cast<int>(current.edges.size());
        st.loop_closure = upd.loop_closure;
        st.global_gate = gate.global_gate;
        st.candidate_vars = static_cast<int>(gate.candidate_local_vars.size());
        st.eta = gate.eta;
        st.delta_eta = gate.delta_eta;
        st.scalar_measurements = scalar_measurement_count(current);
        st.update_flops += estimate_update_flops_from_factor(fm, new_edge_touched_vars, true);

        // Standalone exact dynamic CCOLAMD development telemetry. This does not
        // drive the production ordering path yet; it exercises the owned
        // deterministic-reference rollback/suffix-replay engine on the same
        // evolving symbolic pattern used by the benchmark.
        const std::vector<int> dyn_dirty_vars = touched_vars_from_new_edges(current, first_new_local_edge_idx);
        (void)dyn_exact.refresh(fm.last_H(), dyn_dirty_vars);
        const auto dyn_stats = dyn_exact.stats();
        st.dyn_exact_refreshes = static_cast<int>(dyn_stats.refreshes - prev_dyn_exact_stats.refreshes);
        st.dyn_exact_cold_starts = static_cast<int>(dyn_stats.cold_starts - prev_dyn_exact_stats.cold_starts);
        st.dyn_exact_prefix_reuse_attempts = static_cast<int>(dyn_stats.prefix_reuse_attempts - prev_dyn_exact_stats.prefix_reuse_attempts);
        st.dyn_exact_prefix_reuse_successes = static_cast<int>(dyn_stats.prefix_reuse_successes - prev_dyn_exact_stats.prefix_reuse_successes);
        st.dyn_exact_common_prefix = dyn_stats.last_common_prefix;
        st.dyn_exact_dirty_boundary = dyn_stats.last_dirty_boundary;
        st.dyn_exact_reusable_checkpoint = dyn_stats.last_reusable_checkpoint;
        st.dyn_exact_checkpoint_live_vars = dyn_stats.last_checkpoint_live_vars;
        st.dyn_exact_checkpoint_live_edges = dyn_stats.last_checkpoint_live_edges;
        st.dyn_exact_materialized_checkpoint_reuse_attempts = static_cast<int>(dyn_stats.materialized_checkpoint_reuse_attempts - prev_dyn_exact_stats.materialized_checkpoint_reuse_attempts);
        st.dyn_exact_materialized_checkpoint_reuse_successes = static_cast<int>(dyn_stats.materialized_checkpoint_reuse_successes - prev_dyn_exact_stats.materialized_checkpoint_reuse_successes);
        st.dyn_exact_materialized_checkpoint_reuse_failures = static_cast<int>(dyn_stats.materialized_checkpoint_reuse_failures - prev_dyn_exact_stats.materialized_checkpoint_reuse_failures);
        st.dyn_exact_checkpoint_bank_entries = static_cast<int>(dyn_stats.checkpoint_bank_entries);
        st.dyn_exact_checkpoint_bank_insertions = static_cast<int>(dyn_stats.checkpoint_bank_insertions - prev_dyn_exact_stats.checkpoint_bank_insertions);
        st.dyn_exact_checkpoint_bank_probes = static_cast<int>(dyn_stats.checkpoint_bank_probes - prev_dyn_exact_stats.checkpoint_bank_probes);
        st.dyn_exact_checkpoint_bank_hits = static_cast<int>(dyn_stats.checkpoint_bank_hits - prev_dyn_exact_stats.checkpoint_bank_hits);
        st.dyn_exact_checkpoint_bank_misses = static_cast<int>(dyn_stats.checkpoint_bank_misses - prev_dyn_exact_stats.checkpoint_bank_misses);
        st.dyn_exact_checkpoint_bank_imported_entries = static_cast<int>(dyn_stats.checkpoint_bank_imported_entries);
        st.dyn_exact_checkpoint_bank_import_duplicate_entries = static_cast<int>(dyn_stats.checkpoint_bank_import_duplicate_entries);
        st.dyn_exact_checkpoint_bank_import_invalid_entries = static_cast<int>(dyn_stats.checkpoint_bank_import_invalid_entries);
        st.dyn_exact_checkpoint_bank_import_header_checks = static_cast<int>(dyn_stats.checkpoint_bank_import_header_checks);
        st.dyn_exact_checkpoint_bank_import_header_failures = static_cast<int>(dyn_stats.checkpoint_bank_import_header_failures);
        st.dyn_exact_checkpoint_bank_import_digest_mismatches = static_cast<int>(dyn_stats.checkpoint_bank_import_digest_mismatches);
        st.dyn_exact_checkpoint_bank_import_entry_count_mismatches = static_cast<int>(dyn_stats.checkpoint_bank_import_entry_count_mismatches);
        st.dyn_exact_checkpoint_bank_exported_entries = static_cast<int>(dyn_stats.checkpoint_bank_exported_entries);
        st.dyn_exact_checkpoint_bank_export_writes = static_cast<int>(dyn_stats.checkpoint_bank_export_writes - prev_dyn_exact_stats.checkpoint_bank_export_writes);
        st.dyn_exact_current_pattern_checkpoint_certification_attempts = static_cast<int>(dyn_stats.current_pattern_checkpoint_certification_attempts - prev_dyn_exact_stats.current_pattern_checkpoint_certification_attempts);
        st.dyn_exact_current_pattern_checkpoint_certification_successes = static_cast<int>(dyn_stats.current_pattern_checkpoint_certification_successes - prev_dyn_exact_stats.current_pattern_checkpoint_certification_successes);
        st.dyn_exact_current_pattern_checkpoint_certification_failures = static_cast<int>(dyn_stats.current_pattern_checkpoint_certification_failures - prev_dyn_exact_stats.current_pattern_checkpoint_certification_failures);
        st.dyn_exact_structural_checkpoint_certification_attempts = static_cast<int>(dyn_stats.structural_checkpoint_certification_attempts - prev_dyn_exact_stats.structural_checkpoint_certification_attempts);
        st.dyn_exact_structural_checkpoint_certification_successes = static_cast<int>(dyn_stats.structural_checkpoint_certification_successes - prev_dyn_exact_stats.structural_checkpoint_certification_successes);
        st.dyn_exact_structural_checkpoint_certification_failures = static_cast<int>(dyn_stats.structural_checkpoint_certification_failures - prev_dyn_exact_stats.structural_checkpoint_certification_failures);
        st.dyn_exact_checkpoint_resume_attempts = static_cast<int>(dyn_stats.checkpoint_resume_attempts - prev_dyn_exact_stats.checkpoint_resume_attempts);
        st.dyn_exact_checkpoint_resume_successes = static_cast<int>(dyn_stats.checkpoint_resume_successes - prev_dyn_exact_stats.checkpoint_resume_successes);
        st.dyn_exact_checkpoint_resume_failures = static_cast<int>(dyn_stats.checkpoint_resume_failures - prev_dyn_exact_stats.checkpoint_resume_failures);
        st.dyn_exact_reference_free_resume_attempts = static_cast<int>(dyn_stats.reference_free_resume_attempts - prev_dyn_exact_stats.reference_free_resume_attempts);
        st.dyn_exact_reference_free_resume_candidates = static_cast<int>(dyn_stats.reference_free_resume_candidates - prev_dyn_exact_stats.reference_free_resume_candidates);
        st.dyn_exact_reference_free_resume_failures = static_cast<int>(dyn_stats.reference_free_resume_failures - prev_dyn_exact_stats.reference_free_resume_failures);
        st.dyn_exact_reference_free_resume_reference_matches = static_cast<int>(dyn_stats.reference_free_resume_reference_matches - prev_dyn_exact_stats.reference_free_resume_reference_matches);
        st.dyn_exact_reference_free_resume_reference_mismatches = static_cast<int>(dyn_stats.reference_free_resume_reference_mismatches - prev_dyn_exact_stats.reference_free_resume_reference_mismatches);
        st.dyn_exact_dirty_boundary_safety_checks = static_cast<int>(dyn_stats.dirty_boundary_safety_checks - prev_dyn_exact_stats.dirty_boundary_safety_checks);
        st.dyn_exact_dirty_boundary_safe_reuses = static_cast<int>(dyn_stats.dirty_boundary_safe_reuses - prev_dyn_exact_stats.dirty_boundary_safe_reuses);
        st.dyn_exact_dirty_boundary_unsafe_overestimates = static_cast<int>(dyn_stats.dirty_boundary_unsafe_overestimates - prev_dyn_exact_stats.dirty_boundary_unsafe_overestimates);
        st.dyn_exact_dirty_boundary_underestimates = static_cast<int>(dyn_stats.dirty_boundary_underestimates - prev_dyn_exact_stats.dirty_boundary_underestimates);
        st.dyn_exact_dirty_boundary_exact_matches = static_cast<int>(dyn_stats.dirty_boundary_exact_matches - prev_dyn_exact_stats.dirty_boundary_exact_matches);
        st.dyn_exact_dirty_boundary_candidate_pivots = static_cast<int>(dyn_stats.dirty_boundary_candidate_pivots - prev_dyn_exact_stats.dirty_boundary_candidate_pivots);
        st.dyn_exact_dirty_boundary_safe_pivots = static_cast<int>(dyn_stats.dirty_boundary_safe_pivots - prev_dyn_exact_stats.dirty_boundary_safe_pivots);
        st.dyn_exact_dirty_boundary_unsafe_pivots = static_cast<int>(dyn_stats.dirty_boundary_unsafe_pivots - prev_dyn_exact_stats.dirty_boundary_unsafe_pivots);
        st.dyn_exact_last_dirty_boundary_overestimate = dyn_stats.last_dirty_boundary_overestimate;
        st.dyn_exact_last_dirty_boundary_underestimate = dyn_stats.last_dirty_boundary_underestimate;
        st.dyn_exact_dirty_boundary_checkpoint_probe_attempts = static_cast<int>(dyn_stats.dirty_boundary_checkpoint_probe_attempts - prev_dyn_exact_stats.dirty_boundary_checkpoint_probe_attempts);
        st.dyn_exact_dirty_boundary_checkpoint_probe_structural_successes = static_cast<int>(dyn_stats.dirty_boundary_checkpoint_probe_structural_successes - prev_dyn_exact_stats.dirty_boundary_checkpoint_probe_structural_successes);
        st.dyn_exact_dirty_boundary_checkpoint_probe_structural_failures = static_cast<int>(dyn_stats.dirty_boundary_checkpoint_probe_structural_failures - prev_dyn_exact_stats.dirty_boundary_checkpoint_probe_structural_failures);
        st.dyn_exact_dirty_boundary_checkpoint_probe_resume_matches = static_cast<int>(dyn_stats.dirty_boundary_checkpoint_probe_resume_matches - prev_dyn_exact_stats.dirty_boundary_checkpoint_probe_resume_matches);
        st.dyn_exact_dirty_boundary_checkpoint_probe_resume_mismatches = static_cast<int>(dyn_stats.dirty_boundary_checkpoint_probe_resume_mismatches - prev_dyn_exact_stats.dirty_boundary_checkpoint_probe_resume_mismatches);
        st.dyn_exact_dirty_boundary_checkpoint_probe_resume_failures = static_cast<int>(dyn_stats.dirty_boundary_checkpoint_probe_resume_failures - prev_dyn_exact_stats.dirty_boundary_checkpoint_probe_resume_failures);
        st.dyn_exact_last_dirty_boundary_checkpoint_probe_pivots = dyn_stats.last_dirty_boundary_checkpoint_probe_pivots;
        st.dyn_exact_rollback_suffix = dyn_stats.last_rollback_suffix;
        st.dyn_exact_suffix_replays = static_cast<int>(dyn_stats.certified_suffix_replays - prev_dyn_exact_stats.certified_suffix_replays);
        st.dyn_exact_suffix_replay_failures = static_cast<int>(dyn_stats.certified_suffix_replay_failures - prev_dyn_exact_stats.certified_suffix_replay_failures);
        st.dyn_exact_suffix_replay_pivots_reused = static_cast<int>(dyn_stats.suffix_replay_pivots_reused - prev_dyn_exact_stats.suffix_replay_pivots_reused);
        st.dyn_exact_suffix_replay_pivots_recomputed = static_cast<int>(dyn_stats.suffix_replay_pivots_recomputed - prev_dyn_exact_stats.suffix_replay_pivots_recomputed);
        st.dyn_exact_reference_checks = static_cast<int>(dyn_stats.exact_reference_checks - prev_dyn_exact_stats.exact_reference_checks);
        st.dyn_exact_reference_failures = static_cast<int>(dyn_stats.exact_reference_failures - prev_dyn_exact_stats.exact_reference_failures);
        st.dyn_exact_suitesparse_compat_checks = static_cast<int>(dyn_stats.suitesparse_compat_checks - prev_dyn_exact_stats.suitesparse_compat_checks);
        st.dyn_exact_suitesparse_compat_failures = static_cast<int>(dyn_stats.suitesparse_compat_failures - prev_dyn_exact_stats.suitesparse_compat_failures);
        st.dyn_exact_suitesparse_compat_unavailable = static_cast<int>(dyn_stats.suitesparse_compat_unavailable - prev_dyn_exact_stats.suitesparse_compat_unavailable);
        st.dyn_exact_suitesparse_compat_mismatch_positions = static_cast<int>(dyn_stats.suitesparse_compat_mismatch_positions - prev_dyn_exact_stats.suitesparse_compat_mismatch_positions);
        st.dyn_exact_last_suitesparse_compat_mismatches = dyn_stats.last_suitesparse_compat_mismatches;
        prev_dyn_exact_stats = dyn_stats;

        const auto order_stats = fm.dynamic_ordering_stats();
        st.order_oracle_refreshes = static_cast<int>(order_stats.oracle_refreshes - prev_order_stats.oracle_refreshes);
        st.order_oracle_cache_hits = static_cast<int>(order_stats.oracle_cache_hits - prev_order_stats.oracle_cache_hits);
        st.order_oracle_cache_misses = static_cast<int>(order_stats.oracle_cache_misses - prev_order_stats.oracle_cache_misses);
        st.order_oracle_cache_entries = static_cast<int>(order_stats.oracle_cache_entries);
        st.order_local_pattern_cache_hits = static_cast<int>(order_stats.local_pattern_cache_hits - prev_order_stats.local_pattern_cache_hits);
        st.order_local_pattern_cache_misses = static_cast<int>(order_stats.local_pattern_cache_misses - prev_order_stats.local_pattern_cache_misses);
        st.order_local_pattern_cache_entries = static_cast<int>(order_stats.local_pattern_cache_entries);
        st.order_motif_pattern_cache_hits = static_cast<int>(order_stats.motif_pattern_cache_hits - prev_order_stats.motif_pattern_cache_hits);
        st.order_motif_pattern_cache_misses = static_cast<int>(order_stats.motif_pattern_cache_misses - prev_order_stats.motif_pattern_cache_misses);
        st.order_motif_pattern_cache_entries = static_cast<int>(order_stats.motif_pattern_cache_entries);
        st.order_local_attempts = static_cast<int>(order_stats.local_attempts - prev_order_stats.local_attempts);
        st.order_local_accepts = static_cast<int>(order_stats.local_accepts - prev_order_stats.local_accepts);
        st.order_local_rejects = static_cast<int>(order_stats.local_rejects - prev_order_stats.local_rejects);
        st.order_candidate_windows_generated = static_cast<int>(order_stats.candidate_windows_generated - prev_order_stats.candidate_windows_generated);
        st.order_candidate_windows_tried = static_cast<int>(order_stats.candidate_windows_tried - prev_order_stats.candidate_windows_tried);
        st.order_band_attempts = static_cast<int>(order_stats.band_attempts - prev_order_stats.band_attempts);
        st.order_replay_attempts = static_cast<int>(order_stats.replay_attempts - prev_order_stats.replay_attempts);
        st.order_replay_windows_cached = static_cast<int>(order_stats.replay_windows_cached);
        st.order_regime_switches = static_cast<int>(order_stats.regime_switches - prev_order_stats.regime_switches);
        st.order_regime_creations = static_cast<int>(order_stats.regime_creations - prev_order_stats.regime_creations);
        st.order_regime_merges = static_cast<int>(order_stats.regime_merges - prev_order_stats.regime_merges);
        st.order_num_regimes_discovered = static_cast<int>(order_stats.num_regimes_discovered);
        st.order_current_regime_id = order_stats.current_regime_id;
        st.order_adaptive_reorders = static_cast<int>(order_stats.adaptive_reorders - prev_order_stats.adaptive_reorders);
        st.order_one_hop_accepts = static_cast<int>(order_stats.one_hop_accepts - prev_order_stats.one_hop_accepts);
        st.order_two_hop_accepts = static_cast<int>(order_stats.two_hop_accepts - prev_order_stats.two_hop_accepts);
        st.order_interval_accepts = static_cast<int>(order_stats.interval_accepts - prev_order_stats.interval_accepts);
        st.order_union_accepts = static_cast<int>(order_stats.union_accepts - prev_order_stats.union_accepts);
        st.order_band_accepts = static_cast<int>(order_stats.band_accepts - prev_order_stats.band_accepts);
        st.order_replay_accepts = static_cast<int>(order_stats.replay_accepts - prev_order_stats.replay_accepts);
        st.order_overlap_assembly_attempts = static_cast<int>(order_stats.overlap_assembly_attempts - prev_order_stats.overlap_assembly_attempts);
        st.order_overlap_assembly_accepts = static_cast<int>(order_stats.overlap_assembly_accepts - prev_order_stats.overlap_assembly_accepts);
        st.order_overlap_assembly_piece_hits = static_cast<int>(order_stats.overlap_assembly_piece_hits - prev_order_stats.overlap_assembly_piece_hits);
        st.order_overlap_assembly_proposals = static_cast<int>(order_stats.overlap_assembly_proposals - prev_order_stats.overlap_assembly_proposals);
        st.order_hierarchical_block_cache_hits = static_cast<int>(order_stats.hierarchical_block_cache_hits - prev_order_stats.hierarchical_block_cache_hits);
        st.order_hierarchical_block_cache_misses = static_cast<int>(order_stats.hierarchical_block_cache_misses - prev_order_stats.hierarchical_block_cache_misses);
        st.order_hierarchical_block_cache_entries = static_cast<int>(order_stats.hierarchical_block_cache_entries);
        st.order_hierarchical_block_promotions = static_cast<int>(order_stats.hierarchical_block_promotions - prev_order_stats.hierarchical_block_promotions);
        st.order_precedence_cache_hits = static_cast<int>(order_stats.precedence_cache_hits - prev_order_stats.precedence_cache_hits);
        st.order_precedence_cache_misses = static_cast<int>(order_stats.precedence_cache_misses - prev_order_stats.precedence_cache_misses);
        st.order_precedence_cache_entries = static_cast<int>(order_stats.precedence_cache_entries);
        st.order_precedence_cache_promotions = static_cast<int>(order_stats.precedence_cache_promotions - prev_order_stats.precedence_cache_promotions);
        st.order_precedence_guided_attempts = static_cast<int>(order_stats.precedence_guided_attempts - prev_order_stats.precedence_guided_attempts);
        st.order_precedence_guided_accepts = static_cast<int>(order_stats.precedence_guided_accepts - prev_order_stats.precedence_guided_accepts);
        st.order_precedence_consensus_attempts = static_cast<int>(order_stats.precedence_consensus_attempts - prev_order_stats.precedence_consensus_attempts);
        st.order_precedence_consensus_accepts = static_cast<int>(order_stats.precedence_consensus_accepts - prev_order_stats.precedence_consensus_accepts);
        st.order_precedence_consensus_scc_collapses = static_cast<int>(order_stats.precedence_consensus_scc_collapses - prev_order_stats.precedence_consensus_scc_collapses);
        st.order_exact_output_certifications = static_cast<int>(order_stats.exact_output_certifications - prev_order_stats.exact_output_certifications);
        st.order_certified_local_order_accepts = static_cast<int>(order_stats.certified_local_order_accepts - prev_order_stats.certified_local_order_accepts);
        st.order_certified_exact_cache_order_accepts = static_cast<int>(order_stats.certified_exact_cache_order_accepts - prev_order_stats.certified_exact_cache_order_accepts);
        st.order_certified_oracle_order_fallbacks = static_cast<int>(order_stats.certified_oracle_order_fallbacks - prev_order_stats.certified_oracle_order_fallbacks);
        st.order_exact_cache_imported_entries = static_cast<int>(order_stats.exact_cache_imported_entries - prev_order_stats.exact_cache_imported_entries);
        st.order_exact_cache_exported_entries = static_cast<int>(order_stats.exact_cache_exported_entries - prev_order_stats.exact_cache_exported_entries);
        st.order_no_oracle_cache_hits = static_cast<int>(order_stats.no_oracle_cache_hits - prev_order_stats.no_oracle_cache_hits);
        st.order_no_oracle_cache_misses = static_cast<int>(order_stats.no_oracle_cache_misses - prev_order_stats.no_oracle_cache_misses);
        if (st.order_exact_output_certifications > 0) {
            st.order_local_certification_rate =
                static_cast<double>(st.order_certified_local_order_accepts + st.order_certified_exact_cache_order_accepts) / static_cast<double>(st.order_exact_output_certifications);
            st.order_oracle_fallback_rate =
                static_cast<double>(st.order_certified_oracle_order_fallbacks) / static_cast<double>(st.order_exact_output_certifications);
        }
        const auto symbolic_stats = fm.symbolic_engine_stats();
        st.etree_exact_recomputes = static_cast<int>(symbolic_stats.etree_exact_recomputes - prev_symbolic_stats.etree_exact_recomputes);
        st.etree_exact_cache_imported_entries = static_cast<int>(symbolic_stats.etree_exact_cache_imported_entries - prev_symbolic_stats.etree_exact_cache_imported_entries);
        st.etree_exact_cache_exported_entries = static_cast<int>(symbolic_stats.etree_exact_cache_exported_entries - prev_symbolic_stats.etree_exact_cache_exported_entries);
        st.etree_exact_cache_hits = static_cast<int>(symbolic_stats.etree_exact_cache_hits - prev_symbolic_stats.etree_exact_cache_hits);
        st.etree_exact_cache_misses = static_cast<int>(symbolic_stats.etree_exact_cache_misses - prev_symbolic_stats.etree_exact_cache_misses);
        st.etree_no_exact_cache_hits = static_cast<int>(symbolic_stats.etree_no_exact_cache_hits - prev_symbolic_stats.etree_no_exact_cache_hits);
        st.etree_no_exact_cache_misses = static_cast<int>(symbolic_stats.etree_no_exact_cache_misses - prev_symbolic_stats.etree_no_exact_cache_misses);
        st.etree_local_update_attempts = static_cast<int>(symbolic_stats.etree_local_update_attempts - prev_symbolic_stats.etree_local_update_attempts);
        st.etree_local_update_accepts = static_cast<int>(symbolic_stats.etree_local_update_accepts - prev_symbolic_stats.etree_local_update_accepts);
        st.etree_local_update_fallbacks = static_cast<int>(symbolic_stats.etree_local_update_fallbacks - prev_symbolic_stats.etree_local_update_fallbacks);
        if (st.etree_local_update_attempts > 0) {
            st.etree_local_accept_rate = static_cast<double>(st.etree_local_update_accepts) / static_cast<double>(st.etree_local_update_attempts);
        }
        prev_symbolic_stats = symbolic_stats;
        prev_order_stats = order_stats;

        if (use_spo) {
            AffectedSet affected_local;
            affected_local.vars = gate.candidate_local_vars;

            for (int iter = 0; iter < max_gn_iter; ++iter) {
                const auto sel = SelectiveSolver::solve_reduced(fm, affected_local.vars, selective_alpha);
                const Eigen::VectorXd& dx_local = sel.dx;
                const std::vector<int>& active_local_vars = sel.active_vars;

                st.active_perm = static_cast<int>(sel.active_perm.size());
                st.fell_back_to_full = st.fell_back_to_full || sel.fell_back_to_full;
                st.affected_vars = static_cast<int>(affected_local.vars.size());
                st.solve_flops += estimate_solve_flops_from_factor(fm, active_local_vars);
                st.full_solve_flops += estimate_full_solve_flops_from_factor(fm);
                const Eigen::VectorXd active_dx = gather_by_index(dx_local, active_local_vars);
                st.dx_inf = active_dx.size() > 0 ? active_dx.cwiseAbs().maxCoeff() : 0.0;
                st.gn_iters = iter + 1;

                if (active_local_vars.empty() || st.dx_inf <= dx_th) break;

                affected_local = find_affected(current, active_dx, active_local_vars, dx_th);
                if (affected_local.vars.empty()) break;

                st.update_flops += estimate_update_flops_from_factor(fm, affected_local.vars, false);
                subtract_on_indices(current.x, dx_local, active_local_vars);

                for (int eid : affected_local.edges) {
                    const auto refreshed = fm.build_edge_contribution(current, eid);
                    fm.replace_edge_contribution(eid, refreshed);
                }
            }
        } else {
            const bool do_full_update = (gating_mode == GatingMode::None) || gate.global_gate;
            if (do_full_update) {
                std::vector<int> all_edge_ids(static_cast<size_t>(current.edges.size()));
                std::iota(all_edge_ids.begin(), all_edge_ids.end(), 0);
                for (int iter = 0; iter < max_gn_iter; ++iter) {
                    const Eigen::VectorXd dx_local = fm.solve_cached();
                    st.solve_flops += estimate_full_solve_flops_from_factor(fm);
                    st.full_solve_flops += estimate_full_solve_flops_from_factor(fm);
                    st.active_perm = current.state_size();
                    st.affected_vars = current.state_size();
                    st.dx_inf = dx_local.size() > 0 ? dx_local.cwiseAbs().maxCoeff() : 0.0;
                    st.gn_iters = iter + 1;

                    if (st.dx_inf <= dx_th) break;

                    st.update_flops += estimate_update_flops_from_factor(fm, all_vars(current), false);
                    current.x -= dx_local;
                    for (int eid : all_edge_ids) {
                        const auto refreshed = fm.build_edge_contribution(current, eid);
                        fm.replace_edge_contribution(eid, refreshed);
                    }
                }
            }
        }

        const auto factor_stats = fm.factorization_stats();
        st.factor_full_refactorizations = static_cast<long long>(factor_stats.full_refactorizations - prev_factor_stats.full_refactorizations);
        st.factor_same_size_rank_updates = static_cast<long long>(factor_stats.same_size_rank_updates - prev_factor_stats.same_size_rank_updates);
        st.factor_same_size_rank_downdates = static_cast<long long>(factor_stats.same_size_rank_downdates - prev_factor_stats.same_size_rank_downdates);
        st.factor_cholmod_growth_refactorizations = static_cast<long long>(factor_stats.cholmod_growth_refactorizations - prev_factor_stats.cholmod_growth_refactorizations);
        st.factor_custom_sparse_full_factorizations = static_cast<long long>(factor_stats.custom_sparse_full_factorizations - prev_factor_stats.custom_sparse_full_factorizations);
        st.factor_custom_sparse_suffix_refactorizations = static_cast<long long>(factor_stats.custom_sparse_suffix_refactorizations - prev_factor_stats.custom_sparse_suffix_refactorizations);
        st.factor_custom_sparse_growth_updates = static_cast<long long>(factor_stats.custom_sparse_growth_updates - prev_factor_stats.custom_sparse_growth_updates);
        st.factor_custom_sparse_affected_closure_refactorizations = static_cast<long long>(factor_stats.custom_sparse_affected_closure_refactorizations - prev_factor_stats.custom_sparse_affected_closure_refactorizations);
        st.factor_custom_sparse_structural_closure_refactorizations = static_cast<long long>(factor_stats.custom_sparse_structural_closure_refactorizations - prev_factor_stats.custom_sparse_structural_closure_refactorizations);
        st.factor_custom_sparse_dependency_cache_hits = static_cast<long long>(factor_stats.custom_sparse_dependency_cache_hits - prev_factor_stats.custom_sparse_dependency_cache_hits);
        prev_factor_stats = factor_stats;

        st.normalized_chi2 = normalized_chi_squared(current);
        st.ate = aligned_trajectory_rmse(current, reference);
        cumulative_update_flops += st.update_flops;
        cumulative_solve_flops += st.solve_flops;
        cumulative_full_solve_flops += st.full_solve_flops;
        st.cumulative_update_flops = cumulative_update_flops;
        st.cumulative_solve_flops = cumulative_solve_flops;
        st.cumulative_full_solve_flops = cumulative_full_solve_flops;

        stats.push_back(st);
    }

    return stats;
}

void write_benchmark_csv(const std::filesystem::path& path,
                         const std::vector<BatchStat>& stats) {
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Failed to open benchmark CSV for writing");

    ofs << "backend_requested,backend_actual,batch,edge_id,state_size,num_edges,gn_iters,candidate_vars,active_perm,affected_vars,fell_back_to_full,loop_closure,global_gate,eta,delta_eta,dx_inf,scalar_measurements,normalized_chi2,ate,solve_flops,full_solve_flops,update_flops,cumulative_solve_flops,cumulative_full_solve_flops,cumulative_update_flops,factor_full_refactorizations,factor_same_size_rank_updates,factor_same_size_rank_downdates,factor_cholmod_growth_refactorizations,factor_custom_sparse_full_factorizations,factor_custom_sparse_suffix_refactorizations,factor_custom_sparse_growth_updates,factor_custom_sparse_affected_closure_refactorizations,factor_custom_sparse_structural_closure_refactorizations,factor_custom_sparse_dependency_cache_hits,"
        << "order_oracle_refreshes,order_oracle_cache_hits,order_oracle_cache_misses,order_oracle_cache_entries,"
        << "order_local_pattern_cache_hits,order_local_pattern_cache_misses,order_local_pattern_cache_entries,"
        << "order_motif_pattern_cache_hits,order_motif_pattern_cache_misses,order_motif_pattern_cache_entries,"
        << "order_local_attempts,order_local_accepts,order_local_rejects,order_candidate_windows_generated,order_candidate_windows_tried,"
        << "order_band_attempts,order_replay_attempts,order_replay_windows_cached,order_regime_switches,order_regime_creations,order_regime_merges,order_num_regimes_discovered,order_current_regime_id,order_adaptive_reorders,"
        << "order_one_hop_accepts,order_two_hop_accepts,order_interval_accepts,order_union_accepts,order_band_accepts,order_replay_accepts,"
        << "order_overlap_assembly_attempts,order_overlap_assembly_accepts,order_overlap_assembly_piece_hits,order_overlap_assembly_proposals,"
        << "order_hierarchical_block_cache_hits,order_hierarchical_block_cache_misses,order_hierarchical_block_cache_entries,order_hierarchical_block_promotions,"
        << "order_precedence_cache_hits,order_precedence_cache_misses,order_precedence_cache_entries,order_precedence_cache_promotions,order_precedence_guided_attempts,order_precedence_guided_accepts,"
        << "order_precedence_consensus_attempts,order_precedence_consensus_accepts,order_precedence_consensus_scc_collapses,"
        << "order_exact_output_certifications,order_certified_local_order_accepts,order_certified_exact_cache_order_accepts,order_certified_oracle_order_fallbacks,order_exact_cache_imported_entries,order_exact_cache_exported_entries,order_no_oracle_cache_hits,order_no_oracle_cache_misses,"
        << "order_local_certification_rate,order_oracle_fallback_rate,"
        << "etree_exact_recomputes,etree_exact_cache_imported_entries,etree_exact_cache_exported_entries,etree_exact_cache_hits,etree_exact_cache_misses,etree_no_exact_cache_hits,etree_no_exact_cache_misses,etree_local_update_attempts,etree_local_update_accepts,etree_local_update_fallbacks,etree_local_accept_rate,"
        << "dyn_exact_refreshes,dyn_exact_cold_starts,dyn_exact_prefix_reuse_attempts,dyn_exact_prefix_reuse_successes,dyn_exact_common_prefix,dyn_exact_dirty_boundary,dyn_exact_reusable_checkpoint,dyn_exact_checkpoint_live_vars,dyn_exact_checkpoint_live_edges,dyn_exact_materialized_checkpoint_reuse_attempts,dyn_exact_materialized_checkpoint_reuse_successes,dyn_exact_materialized_checkpoint_reuse_failures,dyn_exact_checkpoint_bank_entries,dyn_exact_checkpoint_bank_insertions,dyn_exact_checkpoint_bank_probes,dyn_exact_checkpoint_bank_hits,dyn_exact_checkpoint_bank_misses,dyn_exact_checkpoint_bank_imported_entries,dyn_exact_checkpoint_bank_import_duplicate_entries,dyn_exact_checkpoint_bank_import_invalid_entries,dyn_exact_checkpoint_bank_import_header_checks,dyn_exact_checkpoint_bank_import_header_failures,dyn_exact_checkpoint_bank_import_digest_mismatches,dyn_exact_checkpoint_bank_import_entry_count_mismatches,dyn_exact_checkpoint_bank_exported_entries,dyn_exact_checkpoint_bank_export_writes,dyn_exact_current_pattern_checkpoint_certification_attempts,dyn_exact_current_pattern_checkpoint_certification_successes,dyn_exact_current_pattern_checkpoint_certification_failures,dyn_exact_structural_checkpoint_certification_attempts,dyn_exact_structural_checkpoint_certification_successes,dyn_exact_structural_checkpoint_certification_failures,dyn_exact_checkpoint_resume_attempts,dyn_exact_checkpoint_resume_successes,dyn_exact_checkpoint_resume_failures,dyn_exact_reference_free_resume_attempts,dyn_exact_reference_free_resume_candidates,dyn_exact_reference_free_resume_failures,dyn_exact_reference_free_resume_reference_matches,dyn_exact_reference_free_resume_reference_mismatches,dyn_exact_dirty_boundary_safety_checks,dyn_exact_dirty_boundary_safe_reuses,dyn_exact_dirty_boundary_unsafe_overestimates,dyn_exact_dirty_boundary_underestimates,dyn_exact_dirty_boundary_exact_matches,dyn_exact_dirty_boundary_candidate_pivots,dyn_exact_dirty_boundary_safe_pivots,dyn_exact_dirty_boundary_unsafe_pivots,dyn_exact_last_dirty_boundary_overestimate,dyn_exact_last_dirty_boundary_underestimate,dyn_exact_dirty_boundary_checkpoint_probe_attempts,dyn_exact_dirty_boundary_checkpoint_probe_structural_successes,dyn_exact_dirty_boundary_checkpoint_probe_structural_failures,dyn_exact_dirty_boundary_checkpoint_probe_resume_matches,dyn_exact_dirty_boundary_checkpoint_probe_resume_mismatches,dyn_exact_dirty_boundary_checkpoint_probe_resume_failures,dyn_exact_last_dirty_boundary_checkpoint_probe_pivots,dyn_exact_rollback_suffix,"
        << "dyn_exact_suffix_replays,dyn_exact_suffix_replay_failures,dyn_exact_suffix_replay_pivots_reused,dyn_exact_suffix_replay_pivots_recomputed,dyn_exact_reference_checks,dyn_exact_reference_failures,"
        << "dyn_exact_suitesparse_compat_checks,dyn_exact_suitesparse_compat_failures,dyn_exact_suitesparse_compat_unavailable,dyn_exact_suitesparse_compat_mismatch_positions,dyn_exact_last_suitesparse_compat_mismatches\n";
    for (const auto& s : stats) {
        ofs << s.backend_requested << ','
            << s.backend_actual << ','
            << s.batch << ','
            << s.edge_id << ','
            << s.state_size << ','
            << s.num_edges << ','
            << s.gn_iters << ','
            << s.candidate_vars << ','
            << s.active_perm << ','
            << s.affected_vars << ','
            << (s.fell_back_to_full ? 1 : 0) << ','
            << (s.loop_closure ? 1 : 0) << ','
            << (s.global_gate ? 1 : 0) << ','
            << s.eta << ','
            << s.delta_eta << ','
            << s.dx_inf << ','
            << s.scalar_measurements << ','
            << s.normalized_chi2 << ','
            << s.ate << ','
            << s.solve_flops << ','
            << s.full_solve_flops << ','
            << s.update_flops << ','
            << s.cumulative_solve_flops << ','
            << s.cumulative_full_solve_flops << ','
            << s.cumulative_update_flops << ','
            << s.factor_full_refactorizations << ','
            << s.factor_same_size_rank_updates << ','
            << s.factor_same_size_rank_downdates << ','
            << s.factor_cholmod_growth_refactorizations << ','
            << s.factor_custom_sparse_full_factorizations << ','
            << s.factor_custom_sparse_suffix_refactorizations << ','
            << s.factor_custom_sparse_growth_updates << ','
            << s.factor_custom_sparse_affected_closure_refactorizations << ','
            << s.factor_custom_sparse_structural_closure_refactorizations << ','
            << s.factor_custom_sparse_dependency_cache_hits << ','
            << s.order_oracle_refreshes << ','
            << s.order_oracle_cache_hits << ','
            << s.order_oracle_cache_misses << ','
            << s.order_oracle_cache_entries << ','
            << s.order_local_pattern_cache_hits << ','
            << s.order_local_pattern_cache_misses << ','
            << s.order_local_pattern_cache_entries << ','
            << s.order_motif_pattern_cache_hits << ','
            << s.order_motif_pattern_cache_misses << ','
            << s.order_motif_pattern_cache_entries << ','
            << s.order_local_attempts << ','
            << s.order_local_accepts << ','
            << s.order_local_rejects << ','
            << s.order_candidate_windows_generated << ','
            << s.order_candidate_windows_tried << ','
            << s.order_band_attempts << ','
            << s.order_replay_attempts << ','
            << s.order_replay_windows_cached << ','
            << s.order_regime_switches << ','
            << s.order_regime_creations << ','
            << s.order_regime_merges << ','
            << s.order_num_regimes_discovered << ','
            << s.order_current_regime_id << ','
            << s.order_adaptive_reorders << ','
            << s.order_one_hop_accepts << ','
            << s.order_two_hop_accepts << ','
            << s.order_interval_accepts << ','
            << s.order_union_accepts << ','
            << s.order_band_accepts << ','
            << s.order_replay_accepts << ','
            << s.order_overlap_assembly_attempts << ','
            << s.order_overlap_assembly_accepts << ','
            << s.order_overlap_assembly_piece_hits << ','
            << s.order_overlap_assembly_proposals << ','
            << s.order_hierarchical_block_cache_hits << ','
            << s.order_hierarchical_block_cache_misses << ','
            << s.order_hierarchical_block_cache_entries << ','
            << s.order_hierarchical_block_promotions << ','
            << s.order_precedence_cache_hits << ','
            << s.order_precedence_cache_misses << ','
            << s.order_precedence_cache_entries << ','
            << s.order_precedence_cache_promotions << ','
            << s.order_precedence_guided_attempts << ','
            << s.order_precedence_guided_accepts << ','
            << s.order_precedence_consensus_attempts << ','
            << s.order_precedence_consensus_accepts << ','
            << s.order_precedence_consensus_scc_collapses << ','
            << s.order_exact_output_certifications << ','
            << s.order_certified_local_order_accepts << ','
            << s.order_certified_exact_cache_order_accepts << ','
            << s.order_certified_oracle_order_fallbacks << ','
            << s.order_exact_cache_imported_entries << ','
            << s.order_exact_cache_exported_entries << ','
            << s.order_no_oracle_cache_hits << ','
            << s.order_no_oracle_cache_misses << ','
            << s.order_local_certification_rate << ','
            << s.order_oracle_fallback_rate << ','
            << s.etree_exact_recomputes << ','
            << s.etree_exact_cache_imported_entries << ','
            << s.etree_exact_cache_exported_entries << ','
            << s.etree_exact_cache_hits << ','
            << s.etree_exact_cache_misses << ','
            << s.etree_no_exact_cache_hits << ','
            << s.etree_no_exact_cache_misses << ','
            << s.etree_local_update_attempts << ','
            << s.etree_local_update_accepts << ','
            << s.etree_local_update_fallbacks << ','
            << s.etree_local_accept_rate << ','
            << s.dyn_exact_refreshes << ','
            << s.dyn_exact_cold_starts << ','
            << s.dyn_exact_prefix_reuse_attempts << ','
            << s.dyn_exact_prefix_reuse_successes << ','
            << s.dyn_exact_common_prefix << ','
            << s.dyn_exact_dirty_boundary << ','
            << s.dyn_exact_reusable_checkpoint << ','
            << s.dyn_exact_checkpoint_live_vars << ','
            << s.dyn_exact_checkpoint_live_edges << ','
            << s.dyn_exact_materialized_checkpoint_reuse_attempts << ','
            << s.dyn_exact_materialized_checkpoint_reuse_successes << ','
            << s.dyn_exact_materialized_checkpoint_reuse_failures << ','
            << s.dyn_exact_checkpoint_bank_entries << ','
            << s.dyn_exact_checkpoint_bank_insertions << ','
            << s.dyn_exact_checkpoint_bank_probes << ','
            << s.dyn_exact_checkpoint_bank_hits << ','
            << s.dyn_exact_checkpoint_bank_misses << ','
            << s.dyn_exact_checkpoint_bank_imported_entries << ','
            << s.dyn_exact_checkpoint_bank_import_duplicate_entries << ','
            << s.dyn_exact_checkpoint_bank_import_invalid_entries << ','
            << s.dyn_exact_checkpoint_bank_import_header_checks << ','
            << s.dyn_exact_checkpoint_bank_import_header_failures << ','
            << s.dyn_exact_checkpoint_bank_import_digest_mismatches << ','
            << s.dyn_exact_checkpoint_bank_import_entry_count_mismatches << ','
            << s.dyn_exact_checkpoint_bank_exported_entries << ','
            << s.dyn_exact_checkpoint_bank_export_writes << ','
            << s.dyn_exact_current_pattern_checkpoint_certification_attempts << ','
            << s.dyn_exact_current_pattern_checkpoint_certification_successes << ','
            << s.dyn_exact_current_pattern_checkpoint_certification_failures << ','
            << s.dyn_exact_structural_checkpoint_certification_attempts << ','
            << s.dyn_exact_structural_checkpoint_certification_successes << ','
            << s.dyn_exact_structural_checkpoint_certification_failures << ','
            << s.dyn_exact_checkpoint_resume_attempts << ','
            << s.dyn_exact_checkpoint_resume_successes << ','
            << s.dyn_exact_checkpoint_resume_failures << ','
            << s.dyn_exact_reference_free_resume_attempts << ','
            << s.dyn_exact_reference_free_resume_candidates << ','
            << s.dyn_exact_reference_free_resume_failures << ','
            << s.dyn_exact_reference_free_resume_reference_matches << ','
            << s.dyn_exact_reference_free_resume_reference_mismatches << ','
            << s.dyn_exact_dirty_boundary_safety_checks << ','
            << s.dyn_exact_dirty_boundary_safe_reuses << ','
            << s.dyn_exact_dirty_boundary_unsafe_overestimates << ','
            << s.dyn_exact_dirty_boundary_underestimates << ','
            << s.dyn_exact_dirty_boundary_exact_matches << ','
            << s.dyn_exact_dirty_boundary_candidate_pivots << ','
            << s.dyn_exact_dirty_boundary_safe_pivots << ','
            << s.dyn_exact_dirty_boundary_unsafe_pivots << ','
            << s.dyn_exact_last_dirty_boundary_overestimate << ','
            << s.dyn_exact_last_dirty_boundary_underestimate << ','
            << s.dyn_exact_dirty_boundary_checkpoint_probe_attempts << ','
            << s.dyn_exact_dirty_boundary_checkpoint_probe_structural_successes << ','
            << s.dyn_exact_dirty_boundary_checkpoint_probe_structural_failures << ','
            << s.dyn_exact_dirty_boundary_checkpoint_probe_resume_matches << ','
            << s.dyn_exact_dirty_boundary_checkpoint_probe_resume_mismatches << ','
            << s.dyn_exact_dirty_boundary_checkpoint_probe_resume_failures << ','
            << s.dyn_exact_last_dirty_boundary_checkpoint_probe_pivots << ','
            << s.dyn_exact_rollback_suffix << ','
            << s.dyn_exact_suffix_replays << ','
            << s.dyn_exact_suffix_replay_failures << ','
            << s.dyn_exact_suffix_replay_pivots_reused << ','
            << s.dyn_exact_suffix_replay_pivots_recomputed << ','
            << s.dyn_exact_reference_checks << ','
            << s.dyn_exact_reference_failures << ','
            << s.dyn_exact_suitesparse_compat_checks << ','
            << s.dyn_exact_suitesparse_compat_failures << ','
            << s.dyn_exact_suitesparse_compat_unavailable << ','
            << s.dyn_exact_suitesparse_compat_mismatch_positions << ','
            << s.dyn_exact_last_suitesparse_compat_mismatches << '\n';
    }
}

} // namespace islam
