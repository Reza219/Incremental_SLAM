#include "islam/affected.hpp"
#include "islam/benchmark.hpp"
#include "islam/factor_manager.hpp"
#include "islam/io_g2o.hpp"
#include "islam/io_toro.hpp"
#include "islam/reorder_edges.hpp"
#include "islam/selective_solver.hpp"
#include "islam/se2.hpp"
#include "islam/update_graph.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

void append_node_vars(const islam::Graph& g, int node_id, std::vector<int>& vars) {
    const auto it = g.id_lookup.find(node_id);
    if (it == g.id_lookup.end()) {
        throw std::runtime_error("append_node_vars: missing node id");
    }
    for (int k = 0; k < it->second.dimension; ++k) vars.push_back(it->second.offset + k);
}

std::vector<int> touched_vars_from_new_edges(const islam::Graph& g, int first_new_local_edge_idx) {
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

std::vector<int> all_vars(const islam::Graph& g) {
    std::vector<int> out(static_cast<size_t>(g.state_size()));
    for (int i = 0; i < g.state_size(); ++i) out[static_cast<size_t>(i)] = i;
    return out;
}

int anchor_dim_for_graph(const islam::Graph& g) {
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

std::vector<int> global_vars_from_edge_ids(const islam::Graph& full,
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

std::vector<int> local_to_global_vars(const islam::Graph& current,
                                      const islam::Graph& full,
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

Eigen::VectorXd full_dx_to_local(const islam::Graph& current,
                                 const islam::Graph& full,
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

struct GatingDecision {
    std::vector<int> candidate_local_vars;
    bool global_gate = false;
    double eta = 0.0;
    double delta_eta = 0.0;
    int eta_state_size = 0;
};

GatingDecision decide_gating(const islam::FactorManager& fm,
                             const islam::Graph& current,
                             int first_new_local_edge_idx,
                             const islam::UpdateGraphResult& upd,
                             islam::GatingMode mode,
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
    case islam::GatingMode::None:
        out.global_gate = true;
        out.candidate_local_vars = all_vars(current);
        break;
    case islam::GatingMode::LCG:
        out.global_gate = upd.loop_closure;
        out.candidate_local_vars = out.global_gate ? all_vars(current) : local_touched;
        break;
    case islam::GatingMode::IGG:
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

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cout << "Usage: islam_demo <graph.{g2o|graph}> [max_batches] [max_gn_iter] [selective_alpha] [lc_gap] [eta_threshold] [gating] [use_spo]\n";
            return 0;
        }

        const fs::path path(argv[1]);
        const int max_batches = (argc >= 3) ? std::max(1, std::stoi(argv[2])) : 10;
        const int max_gn_iter = (argc >= 4) ? std::max(1, std::stoi(argv[3])) : 5;
        const double selective_alpha = (argc >= 5) ? std::stod(argv[4]) : 0.3;
        const int lc_gap = (argc >= 6) ? std::max(0, std::stoi(argv[5])) : 5;
        const double eta_threshold = (argc >= 7) ? std::stod(argv[6]) : 1.0;
        const islam::GatingMode gating_mode = (argc >= 8) ? islam::gating_mode_from_string(argv[7]) : islam::GatingMode::IGG;
        const bool use_spo = (argc >= 9) ? (std::stoi(argv[8]) != 0) : true;
        const double dx_th = 1e-3;
        const double anchor_strength = 1.0;
        const double latent_prior_strength = 1e-9;

        islam::Graph full;
        const auto ext = path.extension().string();
        if (ext == ".g2o") {
            full = islam::read_graph_g2o(path);
        } else if (ext == ".graph") {
            full = islam::read_graph_toro(path);
        } else {
            throw std::runtime_error("Unsupported file extension: " + ext);
        }

        full = islam::reorder_edges(full);

        std::cout << "Loaded graph: " << path << "\n";
        std::cout << "  state size: " << full.state_size() << "\n";
        std::cout << "  nodes:      " << full.id_lookup.size() << "\n";
        std::cout << "  edges:      " << full.edges.size() << "\n";
        std::cout << "  gating:     " << islam::to_string(gating_mode) << "\n";
        std::cout << "  use_spo:    " << (use_spo ? "yes" : "no") << "\n";
        std::cout << "  eta_th:     " << eta_threshold << "\n";

        islam::FactorManager fm;
        fm.configure_incremental(anchor_strength, 0.0, anchor_dim_for_graph(full));

        islam::Graph current;
        bool has_current = false;
        islam::LoopClosureState lc_state;
        double prev_eta = 0.0;
        int prev_state_size = 0;
        bool has_prev_eta = false;

        const int batches = std::min(max_batches, static_cast<int>(full.edges.size()));
        for (int b = 0; b < batches; ++b) {
            const auto upd = islam::update_graph(std::vector<int>{b}, full, has_current ? &current : nullptr, &lc_state, lc_gap);
            current = upd.current;
            has_current = true;

            const int num_new_edges = 1;
            const int first_new_local_edge_idx = static_cast<int>(current.edges.size()) - num_new_edges;

            fm.ensure_state_size(current.state_size());
            const auto contrib = fm.build_edge_contribution(current, first_new_local_edge_idx);
            fm.add_edge_contribution(b, contrib);

            const auto gate = decide_gating(
                fm, current, first_new_local_edge_idx, upd, gating_mode, eta_threshold, prev_eta, prev_state_size, has_prev_eta);
            prev_eta = gate.eta;
            prev_state_size = gate.eta_state_size;
            has_prev_eta = true;

            int gn_used = 0;
            double last_dx_inf = 0.0;
            bool used_full_fallback = false;
            int last_active_perm = 0;
            int last_affected_vars = static_cast<int>(gate.candidate_local_vars.size());

            if (use_spo) {
                islam::AffectedSet affected_local;
                affected_local.vars = gate.candidate_local_vars;

                for (int it = 0; it < max_gn_iter; ++it) {
                    const auto sel = islam::SelectiveSolver::solve_reduced(fm, affected_local.vars, selective_alpha);
                    const Eigen::VectorXd& dx_local = sel.dx;
                    const std::vector<int>& active_local_vars = sel.active_vars;
                    const Eigen::VectorXd active_dx = gather_by_index(dx_local, active_local_vars);
                    last_dx_inf = active_dx.size() > 0 ? active_dx.cwiseAbs().maxCoeff() : 0.0;
                    used_full_fallback = used_full_fallback || sel.fell_back_to_full;
                    last_active_perm = static_cast<int>(sel.active_perm.size());
                    last_affected_vars = static_cast<int>(affected_local.vars.size());
                    gn_used = it + 1;
                    if (active_local_vars.empty() || last_dx_inf <= dx_th) {
                        break;
                    }

                    subtract_on_indices(current.x, dx_local, active_local_vars);

                    const auto next = islam::find_affected(current, active_dx, active_local_vars, dx_th);
                    if (next.edges.empty()) {
                        break;
                    }
                    for (int local_eid : next.edges) {
                        const auto refreshed = fm.build_edge_contribution(current, local_eid);
                        fm.replace_edge_contribution(local_eid, refreshed);
                    }
                    affected_local = next;
                    if (affected_local.vars.empty()) {
                        break;
                    }
                }
            } else {
                const bool do_full_update = (gating_mode == islam::GatingMode::None) || gate.global_gate;
                if (do_full_update) {
                    std::vector<int> all_edge_ids(static_cast<size_t>(current.edges.size()));
                    std::iota(all_edge_ids.begin(), all_edge_ids.end(), 0);
                    for (int it = 0; it < max_gn_iter; ++it) {
                        const Eigen::VectorXd dx_local = fm.solve_cached();
                        last_dx_inf = dx_local.size() > 0 ? dx_local.cwiseAbs().maxCoeff() : 0.0;
                        last_active_perm = current.state_size();
                        last_affected_vars = current.state_size();
                        gn_used = it + 1;
                        if (last_dx_inf <= dx_th) {
                            break;
                        }
                        current.x -= dx_local;
                        for (int eid : all_edge_ids) {
                            const auto refreshed = fm.build_edge_contribution(current, eid);
                            fm.replace_edge_contribution(eid, refreshed);
                        }
                    }
                }
            }

            std::cout << "batch " << (b + 1)
                      << ": state=" << current.state_size()
                      << " edges=" << current.edges.size()
                      << " candidate_vars=" << gate.candidate_local_vars.size()
                      << " loop_closure=" << (upd.loop_closure ? "yes" : "no")
                      << " global_gate=" << (gate.global_gate ? "yes" : "no")
                      << " eta=" << gate.eta
                      << " delta_eta=" << gate.delta_eta
                      << " backend=" << (fm.using_cholmod() ? "cholmod" : "eigen")
                      << " incremental_modify=" << (fm.supports_incremental_updates() ? "yes" : "no")
                      << " active_perm=" << last_active_perm
                      << " affected_vars=" << last_affected_vars
                      << " full_fallback=" << (used_full_fallback ? "yes" : "no")
                      << " max|dx|=" << last_dx_inf
                      << " gn=" << gn_used
                      << "\n";
        }

        const Eigen::Vector3d pose(1.0, 2.0, 0.3);
        const Eigen::Matrix3d T = islam::v2t(pose);
        const Eigen::Vector3d recovered = islam::t2v(T);
        std::cout << "SE2 roundtrip example: " << recovered.transpose() << "\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
}
