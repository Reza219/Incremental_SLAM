#include "islam/affected.hpp"
#include "islam/factor_manager.hpp"
#include "islam/graph.hpp"
#include "islam/linearize.hpp"
#include "islam/selective_solver.hpp"
#include "islam/sparse_expanding_cholesky.hpp"
#include "islam/symbolic_engine.hpp"
#include "islam/dynamic_ccolamd_engine.hpp"
#include "islam/deterministic_batch_ccolamd.hpp"
#include "islam/dynamic_exact_ccolamd.hpp"
#include "islam/se2.hpp"
#include "islam/se3.hpp"
#include "islam/types.hpp"
#include "islam/update_graph.hpp"

#include <Eigen/Dense>
#include <Eigen/OrderingMethods>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool cond, const char* msg) {
    if (!cond) throw std::runtime_error(msg);
}

islam::Graph make_two_pose_graph() {
    islam::Graph g;
    g.x.resize(6);
    g.x << 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0;

    g.id_lookup[0] = islam::NodeRef{0, 3};
    g.id_lookup[1] = islam::NodeRef{3, 3};
    g.var2node = {0, 0, 0, 1, 1, 1};

    islam::Edge e;
    e.type = islam::EdgeType::PosePose;
    e.from_id = 0;
    e.to_id = 1;
    e.from_idx = 0;
    e.to_idx = 3;
    e.measurement = Eigen::Vector3d(1.0, 0.0, 0.0);
    e.information = Eigen::Matrix3d::Identity();
    g.edges.push_back(e);
    g.node_edges[0] = {0};
    g.node_edges[1] = {0};
    return g;
}

islam::Graph make_self_loop_pose_graph() {
    islam::Graph g;
    g.x.resize(3);
    g.x << 0.0, 0.0, 0.0;

    g.id_lookup[0] = islam::NodeRef{0, 3};
    g.var2node = {0, 0, 0};

    islam::Edge e;
    e.type = islam::EdgeType::PosePose;
    e.from_id = 0;
    e.to_id = 0;
    e.from_idx = 0;
    e.to_idx = 0;
    e.measurement = Eigen::Vector3d(0.0, 0.0, 0.0);
    e.information = Eigen::Matrix3d::Identity();
    g.edges.push_back(e);
    g.node_edges[0] = {0};
    return g;
}

islam::Graph make_three_pose_chain() {
    islam::Graph g;
    g.x.resize(9);
    g.x << 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0,
           0.0, 0.0, 0.0;

    g.id_lookup[0] = islam::NodeRef{0, 3};
    g.id_lookup[1] = islam::NodeRef{3, 3};
    g.id_lookup[2] = islam::NodeRef{6, 3};
    g.var2node = {0, 0, 0, 1, 1, 1, 2, 2, 2};

    islam::Edge e01;
    e01.type = islam::EdgeType::PosePose;
    e01.from_id = 0;
    e01.to_id = 1;
    e01.from_idx = 0;
    e01.to_idx = 3;
    e01.measurement = Eigen::Vector3d(1.0, 0.0, 0.0);
    e01.information = Eigen::Matrix3d::Identity();
    g.edges.push_back(e01);

    islam::Edge e12 = e01;
    e12.from_id = 1;
    e12.to_id = 2;
    e12.from_idx = 3;
    e12.to_idx = 6;
    g.edges.push_back(e12);
    g.node_edges[0] = {0};
    g.node_edges[1] = {0,1};
    g.node_edges[2] = {1};
    return g;
}

islam::Graph make_four_pose_chain() {
    islam::Graph g;
    g.x.resize(12);
    g.x << 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0,
           0.0, 0.0, 0.0,
           0.0, 0.0, 0.0;

    g.id_lookup[0] = islam::NodeRef{0, 3};
    g.id_lookup[1] = islam::NodeRef{3, 3};
    g.id_lookup[2] = islam::NodeRef{6, 3};
    g.id_lookup[3] = islam::NodeRef{9, 3};
    g.var2node = {0,0,0,1,1,1,2,2,2,3,3,3};

    islam::Edge e01;
    e01.type = islam::EdgeType::PosePose;
    e01.from_id = 0;
    e01.to_id = 1;
    e01.from_idx = 0;
    e01.to_idx = 3;
    e01.measurement = Eigen::Vector3d(1.0, 0.0, 0.0);
    e01.information = Eigen::Matrix3d::Identity();
    g.edges.push_back(e01);

    islam::Edge e12 = e01;
    e12.from_id = 1;
    e12.to_id = 2;
    e12.from_idx = 3;
    e12.to_idx = 6;
    g.edges.push_back(e12);

    islam::Edge e23 = e01;
    e23.from_id = 2;
    e23.to_id = 3;
    e23.from_idx = 6;
    e23.to_idx = 9;
    g.edges.push_back(e23);

    g.node_edges[0] = {0};
    g.node_edges[1] = {0,1};
    g.node_edges[2] = {1,2};
    g.node_edges[3] = {2};
    return g;
}

islam::Graph make_two_pose3d_graph() {
    islam::Graph g;
    g.x.resize(12);
    g.x.setZero();
    g.x.segment<3>(6) << 1.0, 0.0, 0.0;

    g.id_lookup[0] = islam::NodeRef{0, 6};
    g.id_lookup[1] = islam::NodeRef{6, 6};
    g.var2node = {0,0,0,0,0,0,1,1,1,1,1,1};

    islam::Edge e;
    e.type = islam::EdgeType::PosePose3D;
    e.from_id = 0;
    e.to_id = 1;
    e.from_idx = 0;
    e.to_idx = 6;
    Eigen::Matrix<double, 6, 1> z;
    z << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    e.measurement = z;
    e.information = Eigen::Matrix<double, 6, 6>::Identity();
    g.edges.push_back(e);
    g.node_edges[0] = {0};
    g.node_edges[1] = {0};
    return g;
}

std::vector<int> vars_from_edge(const islam::Graph& g, const islam::Edge& e) {
    std::vector<int> vars;
    const auto it_from = g.id_lookup.find(e.from_id);
    if (it_from == g.id_lookup.end()) throw std::runtime_error("vars_from_edge missing from node");
    for (int k = 0; k < it_from->second.dimension; ++k) vars.push_back(e.from_idx + k);
    if (!e.is_unary()) {
        const auto it_to = g.id_lookup.find(e.to_id);
        if (it_to == g.id_lookup.end()) throw std::runtime_error("vars_from_edge missing to node");
        for (int k = 0; k < it_to->second.dimension; ++k) vars.push_back(e.to_idx + k);
    }
    return vars;
}

bool has_duplicates(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    return std::adjacent_find(values.begin(), values.end()) != values.end();
}

void check_se3_pose_pose_jacobian_against_finite_difference(const Eigen::Matrix<double, 6, 1>& x1,
                                                             const Eigen::Matrix<double, 6, 1>& x2,
                                                             const Eigen::Matrix<double, 6, 1>& z) {
    Eigen::Matrix<double, 6, 1> e;
    Eigen::Matrix<double, 6, 6> A;
    Eigen::Matrix<double, 6, 6> B;
    islam::linearize_pose_pose_3d(x1, x2, z, e, A, B);

    constexpr double eps = 1e-7;
    Eigen::Matrix<double, 6, 6> A_fd;
    Eigen::Matrix<double, 6, 6> B_fd;
    A_fd.setZero();
    B_fd.setZero();
    for (int k = 0; k < 6; ++k) {
        Eigen::Matrix<double, 6, 1> xp = x1;
        Eigen::Matrix<double, 6, 1> xm = x1;
        xp[k] += eps;
        xm[k] -= eps;
        A_fd.col(k) = (islam::se3_relative_error(xp, x2, z) - islam::se3_relative_error(xm, x2, z)) / (2.0 * eps);
    }
    for (int k = 0; k < 6; ++k) {
        Eigen::Matrix<double, 6, 1> xp = x2;
        Eigen::Matrix<double, 6, 1> xm = x2;
        xp[k] += eps;
        xm[k] -= eps;
        B_fd.col(k) = (islam::se3_relative_error(x1, xp, z) - islam::se3_relative_error(x1, xm, z)) / (2.0 * eps);
    }

    require((A - A_fd).norm() < 5e-6, "Analytic SE3 A Jacobian does not match finite difference");
    require((B - B_fd).norm() < 5e-6, "Analytic SE3 B Jacobian does not match finite difference");
}

} // namespace

int main() {
    try {
        {
            const Eigen::Vector3d v(1.0, -2.0, 0.4);
            const Eigen::Matrix3d T = islam::v2t(v);
            const Eigen::Vector3d vv = islam::t2v(T);
            require((v - vv).norm() < 1e-12, "SE2 roundtrip failed");
        }

        {
            Eigen::Matrix<double, 6, 1> v;
            v << 1.0, -2.0, 0.5, 0.1, -0.2, 0.3;
            const Eigen::Isometry3d T = islam::se3_vec_to_iso(v);
            const Eigen::Matrix<double, 6, 1> vv = islam::se3_iso_to_vec(T);
            require((v - vv).norm() < 1e-10, "SE3 roundtrip failed");
        }

        {
            islam::FactorManager fm;
            const islam::Graph g = make_two_pose_graph();
            const auto c = fm.build_edge_contribution(g, 0);
            const std::vector<int> expected = {0, 1, 2, 3, 4, 5};
            require(c.touched_vars == expected, "EdgeContribution touched_vars should contain each endpoint variable exactly once");
            require(!has_duplicates(c.touched_vars), "EdgeContribution touched_vars must not contain duplicates");
        }

        {
            islam::FactorManager fm;
            const islam::Graph g = make_self_loop_pose_graph();
            const auto c = fm.build_edge_contribution(g, 0);
            const std::vector<int> expected = {0, 1, 2};
            require(c.touched_vars == expected, "Self-loop EdgeContribution touched_vars should be unique");
            require(!has_duplicates(c.touched_nodes), "Self-loop EdgeContribution touched_nodes should be unique");
        }


        {
            Eigen::Matrix<double, 6, 1> x1;
            Eigen::Matrix<double, 6, 1> x2;
            Eigen::Matrix<double, 6, 1> z;
            x1 << 0.2, -0.1, 0.3, 0.05, -0.08, 0.12;
            x2 << 1.1, 0.4, -0.2, -0.03, 0.10, -0.04;
            z  << 0.8, 0.6, -0.4, 0.02, -0.03, 0.07;
            check_se3_pose_pose_jacobian_against_finite_difference(x1, x2, z);

            x1 << -0.4, 0.3, 0.7, 0.35, -0.20, 0.15;
            x2 << 0.6, -0.5, 0.1, -0.25, 0.18, 0.22;
            z  << 1.0, -0.2, 0.4, 0.12, 0.05, -0.16;
            check_se3_pose_pose_jacobian_against_finite_difference(x1, x2, z);
        }

        {
            const islam::Graph g = make_two_pose3d_graph();
            const auto lin = islam::jacobian_edge_jr(g.edges.front(), g);
            require(lin.J.rows() == 6, "3D J row count mismatch");
            require(lin.J.cols() == 12, "3D J col count mismatch");
            require(lin.r.size() == 6, "3D residual size mismatch");
            require(lin.r.norm() < 1e-12, "3D pose-pose residual should be zero at measurement-consistent state");

            islam::FactorManager fm;
            const auto ne = fm.build_normal_equations(g, nullptr, 1.0);
            require(ne.H.rows() == 12 && ne.H.cols() == 12, "3D normal matrix wrong size");
            const Eigen::VectorXd dx = fm.solve_graph(g, nullptr, 1.0);
            require(dx.size() == 12, "3D solution size mismatch");
            require(dx.allFinite(), "3D solution contains non-finite values");
        }

        {
            const islam::Graph g = make_two_pose_graph();
            const auto lin = islam::jacobian_edge_jr(g.edges.front(), g);
            require(lin.J.rows() == 3, "Unexpected J row count");
            require(lin.J.cols() == 6, "Unexpected J col count");
            require(lin.r.size() == 3, "Unexpected residual size");
        }

        {
            islam::Graph g = make_two_pose_graph();
            islam::FactorManager fm;
            const auto ne = fm.build_normal_equations(g, nullptr, 1.0);
            require(ne.H.rows() == 6 && ne.H.cols() == 6, "Normal matrix wrong size");
            require(ne.g.size() == 6, "Normal RHS wrong size");

            const Eigen::VectorXd dx = fm.solve_graph(g, nullptr, 1.0);
            require(dx.size() == 6, "Solution size mismatch");
            require(dx.allFinite(), "Solution contains non-finite values");
            require(std::abs(dx[3] + 1.0) < 1e-9, "Unexpected x-update on second pose");
        }

        {
            const islam::Graph g2 = make_two_pose_graph();
            const islam::Graph g3 = make_three_pose_chain();

            islam::FactorManager fm_rebuild;
            fm_rebuild.rebuild_from_graph(g3, nullptr, 1.0);

            islam::FactorManager fm_online;
            fm_online.configure_incremental(1.0, 0.0, 3);
            fm_online.ensure_state_size(g2.state_size());
            fm_online.add_edge_contribution(0, fm_online.build_edge_contribution(g2, 0));
            fm_online.ensure_state_size(g3.state_size());
            fm_online.add_edge_contribution(1, fm_online.build_edge_contribution(g3, 1));

            require(fm_online.factor_covers_state(), "True online expansion should leave a factor covering the grown state");
            require(fm_online.last_H().rows() == g3.state_size(), "True online expansion H rows mismatch");
            require(fm_online.last_H().cols() == g3.state_size(), "True online expansion H cols mismatch");
            require((fm_online.last_g() - fm_rebuild.last_g()).norm() < 1e-10, "True online expansion g mismatch");

            const Eigen::VectorXd dx_full = fm_rebuild.solve_cached();
            const Eigen::VectorXd dx_online = fm_online.solve_cached();
            require((dx_full - dx_online).norm() < 1e-10, "True online expansion solve mismatch");

            const auto stats = fm_online.factorization_stats();
            require(stats.state_expansions == 2, "True online expansion should record two state-growth events");
            require(stats.incremental_edge_adds == 2, "True online expansion should record two incremental edge additions");
        }

        {
            const islam::Graph g2 = make_two_pose_graph();
            const islam::Graph g3 = make_three_pose_chain();

            islam::FactorManager fm_rebuild;
            fm_rebuild.rebuild_from_graph(g3, nullptr, 1.0);

            islam::FactorManager fm_sparse;
            fm_sparse.enable_sparse_expanding_cholesky(true);
            fm_sparse.configure_incremental(1.0, 0.0, 3);
            fm_sparse.ensure_state_size(g2.state_size());
            fm_sparse.add_edge_contribution(0, fm_sparse.build_edge_contribution(g2, 0));
            fm_sparse.ensure_state_size(g3.state_size());
            fm_sparse.add_edge_contribution(1, fm_sparse.build_edge_contribution(g3, 1));

            require(fm_sparse.using_sparse_expanding_cholesky(), "Custom sparse expanding backend should be enabled");
            require(fm_sparse.factor_covers_state(), "Custom sparse expanding backend should cover grown state");
            const Eigen::VectorXd dx_full = fm_rebuild.solve_cached();
            const Eigen::VectorXd dx_sparse = fm_sparse.solve_cached();
            require((dx_full - dx_sparse).norm() < 1e-9, "Custom sparse expanding solve mismatch");
            const auto stats = fm_sparse.factorization_stats();
            require(stats.custom_sparse_growth_updates + stats.custom_sparse_dynamic_reorders >= 1,
                    "Custom sparse backend should record growth update or dynamic-order refactor");
            require(stats.custom_sparse_l21_recomputes >= 1,
                    "Custom sparse backend should recompute L21 for exact bordered suffix factorization");
            require(stats.custom_sparse_left_looking_factorizations >= 1,
                    "Custom sparse backend should use sparse left-looking suffix factorization");
            require(stats.custom_sparse_schur_complements >= 1,
                    "Custom sparse backend should form sparse suffix Schur complements");
            require(stats.custom_sparse_triangular_solves >= 1,
                    "Custom sparse backend should solve sparse prefix triangular systems");
        }

        {
            Eigen::SparseMatrix<double> H(8, 8);
            std::vector<Eigen::Triplet<double>> trips;
            for (int i = 0; i < 8; ++i) {
                trips.emplace_back(i, i, 2.0 + 0.1 * i);
            }
            H.setFromTriplets(trips.begin(), trips.end());
            H.makeCompressed();

            islam::SparseExpandingCholesky dyn;
            dyn.factorize(H);
            dyn.apply_diagonal_update({6}, 0.5, true);

            H.coeffRef(6, 6) += 0.5;
            H.makeCompressed();
            islam::SparseExpandingCholesky full;
            full.factorize(H);

            Eigen::VectorXd rhs(8);
            rhs << 1.0, -0.5, 0.25, 2.0, -1.0, 0.75, 0.4, -0.2;
            require((dyn.solve(rhs) - full.solve(rhs)).norm() < 1e-11,
                    "Affected-closure sparse update solve mismatch");
            require(dyn.stats().affected_closure_refactorizations >= 1,
                    "Sparse backend should record affected etree-closure refactorization");
            require(dyn.stats().affected_columns_refactored == 1,
                    "Diagonal leaf update should refactor only one affected column");
            require(dyn.stats().suffix_refactorizations == 0,
                    "Affected-closure update should avoid suffix refactorization when closure is smaller");
            require(dyn.stats().numeric_only_updates >= 1,
                    "Diagonal update should be classified as numeric-only when the normal pattern is unchanged");
            require(dyn.stats().affected_closure_certifications >= 1,
                    "Affected-closure update should be certified before it is accepted");
            require(dyn.stats().affected_closure_certification_failures == 0,
                    "Simple diagonal closure update should pass certification");
            require(dyn.stats().factorization_residual_checks >= 1,
                    "Affected-closure certification should perform a residual check");
            require(dyn.stats().column_local_certifications >= 1,
                    "Affected-closure certification should use the M75 column-local residual path first");
            require(dyn.stats().column_local_certification_columns >= 1,
                    "Column-local certification should record how many residual columns were checked");
            require(dyn.stats().dependency_cache_rebuilds >= 1,
                    "M76 should build the maintained factor-dependency cache after factorization");
            require(dyn.stats().dependency_cache_hits >= 1,
                    "Numeric-only closure updates should use the maintained dependency cache");
            require(dyn.stats().dependency_closure_computations >= 1,
                    "M76 should compute the affected set from cached factor dependencies");
            require(dyn.stats().last_dependency_closure_size == 1,
                    "Diagonal leaf update should have a one-column cached dependency closure");
        }

        {
            Eigen::SparseMatrix<double> H(6, 6);
            std::vector<Eigen::Triplet<double>> trips;
            for (int i = 0; i < 6; ++i) {
                trips.emplace_back(i, i, 3.0 + 0.1 * i);
            }
            H.setFromTriplets(trips.begin(), trips.end());
            H.makeCompressed();

            islam::SparseExpandingCholesky dyn;
            dyn.factorize(H);

            Eigen::SparseMatrix<double> D(6, 6);
            std::vector<Eigen::Triplet<double>> dtrips{{2, 4, 0.15}, {4, 2, 0.15}};
            D.setFromTriplets(dtrips.begin(), dtrips.end());
            D.makeCompressed();
            dyn.apply_hessian_delta(D, {2, 4}, true);

            H += D;
            H.makeCompressed();
            islam::SparseExpandingCholesky full;
            full.factorize(H);

            Eigen::VectorXd rhs(6);
            rhs << 1.0, -0.5, 0.25, 2.0, -1.0, 0.75;
            require((dyn.solve(rhs) - full.solve(rhs)).norm() < 1e-11,
                    "Structural-pattern fallback sparse update solve mismatch");
            require(dyn.stats().structural_pattern_changes >= 1,
                    "Sparse backend should detect structural normal-pattern changes");
            require(dyn.stats().symbolic_pattern_classifications >= 1,
                    "Structural pattern changes should be classified against the symbolic factor pattern");
            require(dyn.stats().structural_closure_attempts >= 1,
                    "Local structural pattern changes should attempt a certified closure update");
            require(dyn.stats().structural_closure_refactorizations >= 1,
                    "Local structural pattern changes should be accepted only through certified closure recomputation");
            require(dyn.stats().suffix_refactorizations == 0,
                    "Local structural closure update should avoid suffix refactorization when certified");
            require(dyn.stats().column_local_certifications >= 1,
                    "Certified structural closure should use column-local residual certification");
            require(dyn.stats().etree_closure_cache_bypasses >= 1,
                    "Structural factor-pattern changes should bypass the cached dependency graph");
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm_rebuild;
            islam::FactorManager fm_inc;

            fm_rebuild.rebuild_from_graph(g, nullptr, 1.0);

            fm_inc.reserve_full_state(g.state_size(), 1.0, 1e-9);
            fm_inc.activate_vars(vars_from_edge(g, g.edges[0]));
            fm_inc.add_edge_contribution(0, fm_inc.build_edge_contribution(g, 0));
            fm_inc.activate_vars(vars_from_edge(g, g.edges[1]));
            fm_inc.add_edge_contribution(1, fm_inc.build_edge_contribution(g, 1));

            require(fm_inc.last_H().rows() == fm_rebuild.last_H().rows(), "Incremental H rows mismatch");
            require(fm_inc.last_H().cols() == fm_rebuild.last_H().cols(), "Incremental H cols mismatch");
            require((fm_inc.last_g() - fm_rebuild.last_g()).norm() < 1e-10, "Incremental g mismatch");

            const Eigen::VectorXd dx_full = fm_rebuild.solve_cached();
            const Eigen::VectorXd dx_inc = fm_inc.solve_cached();
            require((dx_full - dx_inc).norm() < 1e-10, "Incremental solve mismatch");
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.reserve_full_state(g.state_size(), 1.0, 1e-9);
            fm.activate_vars(vars_from_edge(g, g.edges[0]));
            fm.add_edge_contribution(0, fm.build_edge_contribution(g, 0));
            fm.activate_vars(vars_from_edge(g, g.edges[1]));
            fm.add_edge_contribution(1, fm.build_edge_contribution(g, 1));

            g.x[6] = 0.2;
            const auto c_new = fm.build_edge_contribution(g, 1);
            fm.replace_edge_contribution(1, c_new);
            const Eigen::VectorXd dx = fm.solve_cached();
            require(dx.size() == 9, "Replace contribution solve wrong size");
            require(dx.allFinite(), "Replace contribution solve non-finite");
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.rebuild_from_graph(g, nullptr, 1.0);
            const auto p = fm.get_permutation();
            require(static_cast<int>(p.size()) == g.state_size(), "Permutation size mismatch");
            std::vector<int> seen(static_cast<size_t>(g.state_size()), 0);
            for (int idx : p) {
                require(idx >= 0 && idx < g.state_size(), "Permutation entry out of range");
                seen[static_cast<size_t>(idx)] += 1;
            }
            for (int c : seen) require(c == 1, "Permutation must be bijective");

            const auto ps = fm.build_permuted_system();
            require(ps.H_perm.rows() == g.state_size(), "Permuted H row mismatch");
            require(ps.H_perm.cols() == g.state_size(), "Permuted H col mismatch");
            require(ps.g_perm.size() == g.state_size(), "Permuted g size mismatch");
            const auto etree = fm.get_elimination_tree();
            require(static_cast<int>(etree.size()) == g.state_size(), "Elimination tree size mismatch");
        }


        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.rebuild_from_graph(g, nullptr, 1.0);
            const auto etree_before = fm.get_elimination_tree();

            islam::Graph g_pert = g;
            g_pert.x[6] = 0.15;
            const auto c_new = fm.build_edge_contribution(g_pert, 1);
            fm.replace_edge_contribution(1, c_new);
            const auto etree_after = fm.get_elimination_tree();
            require(etree_before == etree_after, "Numeric-only replace should preserve cached elimination tree structure");
        }

        {
            islam::Graph g = make_four_pose_chain();
            islam::FactorManager fm;
            fm.reserve_full_state(g.state_size(), 1.0, 1e-9);
            fm.activate_vars(vars_from_edge(g, g.edges[0]));
            fm.add_edge_contribution(0, fm.build_edge_contribution(g, 0));
            fm.activate_vars(vars_from_edge(g, g.edges[1]));
            fm.add_edge_contribution(1, fm.build_edge_contribution(g, 1));
            const auto etree_before = fm.get_elimination_tree();

            fm.activate_vars(vars_from_edge(g, g.edges[2]));
            fm.add_edge_contribution(2, fm.build_edge_contribution(g, 2));
            const auto etree_after = fm.get_elimination_tree();
            require(etree_before.size() == etree_after.size(), "Elimination tree size changed unexpectedly");
            int changed = 0;
            for (int i = 0; i < static_cast<int>(etree_before.size()); ++i) {
                if (etree_before[static_cast<size_t>(i)] != etree_after[static_cast<size_t>(i)]) {
                    ++changed;
                }
            }
            require(changed > 0, "Local symbolic structural update should change some etree entries");
            require(changed < g.state_size(), "Local symbolic structural update should not rewrite the entire etree");
        }

        {
            islam::Graph g = make_four_pose_chain();
            islam::FactorManager fm;
            fm.reserve_full_state(g.state_size(), 1.0, 1e-9);
            fm.activate_vars(vars_from_edge(g, g.edges[0]));
            fm.add_edge_contribution(0, fm.build_edge_contribution(g, 0));
            fm.activate_vars(vars_from_edge(g, g.edges[1]));
            fm.add_edge_contribution(1, fm.build_edge_contribution(g, 1));
            const auto p_before = fm.get_permutation();

            fm.activate_vars(vars_from_edge(g, g.edges[2]));
            fm.add_edge_contribution(2, fm.build_edge_contribution(g, 2));
            const auto p_after = fm.get_permutation();
            require(p_before.size() == p_after.size(), "Permutation size changed unexpectedly");
            int changed_positions = 0;
            for (int i = 0; i < static_cast<int>(p_before.size()); ++i) {
                if (p_before[static_cast<size_t>(i)] != p_after[static_cast<size_t>(i)]) {
                    ++changed_positions;
                }
            }
            require(changed_positions > 0, "Local symbolic structural update should change some ordering positions");
            require(changed_positions < g.state_size(), "Local symbolic structural update should not rewrite the entire ordering");

            islam::FactorManager fm_ref;
            fm_ref.rebuild_from_graph(g, nullptr, 1.0);
            const auto p_ref = fm_ref.get_permutation();
            const auto etree_ref = fm_ref.get_elimination_tree();
            require(p_after == p_ref, "Certified symbolic cache permutation must match exact reference recomputation");
            require(fm.get_elimination_tree() == etree_ref, "Certified symbolic cache etree must match exact reference recomputation");
            const auto order_stats = fm.dynamic_ordering_stats();
            require(order_stats.oracle_refreshes > 0 || order_stats.oracle_cache_hits > 0, "Dynamic ordering stats should record exact-order access");
            require(order_stats.oracle_cache_entries >= 1, "Dynamic ordering cache should retain at least one exact pattern entry");
            require(order_stats.local_pattern_cache_hits + order_stats.local_pattern_cache_misses >= order_stats.candidate_windows_tried,
                    "Local pattern cache telemetry should cover tried candidate windows");
            require(order_stats.local_pattern_cache_entries >= 0,
                    "Local pattern cache entry counter should be present");
            require(order_stats.motif_pattern_cache_entries >= 0,
                    "Motif pattern cache entry counter should be present");
            require(order_stats.local_attempts > 0, "Dynamic ordering stats should record local attempts");
            require(order_stats.candidate_windows_generated >= order_stats.candidate_windows_tried,
                    "Generated ordering windows should dominate attempted windows");
            require(order_stats.band_attempts + order_stats.replay_attempts >= 0,
                    "Extended dynamic-ordering counters should be present");
            require(order_stats.hierarchical_block_cache_entries >= 0,
                    "Hierarchical block cache entry counter should be present");
            require(order_stats.hierarchical_block_cache_hits + order_stats.hierarchical_block_cache_misses >= 0,
                    "Hierarchical block cache telemetry should be present");
            require(order_stats.precedence_cache_entries >= 0,
                    "Precedence cache entry counter should be present");
            require(order_stats.precedence_cache_hits + order_stats.precedence_cache_misses >= 0,
                    "Precedence cache telemetry should be present");
            require(order_stats.precedence_guided_attempts >= order_stats.precedence_guided_accepts,
                    "Precedence-guided acceptance count should not exceed attempts");
            require(order_stats.adaptive_reorders >= order_stats.local_accepts,
                    "Adaptive reorder counter should track accepted local proposals");
            require(order_stats.current_regime_id >= 0,
                    "Dynamic ordering should expose a valid current regime id");
            require(order_stats.num_regimes_discovered >= 1,
                    "Self-tuning regime engine should discover at least one regime");
            require(order_stats.regime_creations >= 1,
                    "Self-tuning regime engine should record regime creation");
            require(order_stats.one_hop_accepts + order_stats.two_hop_accepts + order_stats.interval_accepts +
                        order_stats.union_accepts + order_stats.band_accepts + order_stats.replay_accepts
                        == order_stats.local_accepts,
                    "Per-family acceptance counters should sum to local accepts");
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.rebuild_from_graph(g, nullptr, 1.0);
            const auto full = fm.solve_cached();
            const auto sel = islam::SelectiveSolver::solve_reduced(fm, {0,1,2,3,4,5,6,7,8}, 1.0);
            require(!sel.fell_back_to_full, "All-vars selective solve should not fall back when alpha permits");
            require(sel.used_sparse_partial_solve, "Selective all-vars solve should use sparse triangular path");
            require(sel.used_factor_block_solve, "Selective all-vars solve should reuse maintained factor blocks");
            require((sel.dx - full).norm() < 1e-10, "Selective all-vars solve mismatch");
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.rebuild_from_graph(g, nullptr, 1.0);
            const auto sel = islam::SelectiveSolver::solve_reduced(fm, {6,7,8}, 1.0);
            require(!sel.fell_back_to_full, "Subset selective solve should stay selective when alpha permits");
            require(!sel.active_vars.empty(), "Subset selective solve should expose active variables");
            require(sel.used_sparse_partial_solve, "Subset selective solve should use sparse triangular path");
            require(sel.used_factor_block_solve, "Subset selective solve should reuse maintained factor blocks when alpha permits");
            std::vector<int> is_active(static_cast<size_t>(g.state_size()), 0);
            for (int v : sel.active_vars) {
                require(v >= 0 && v < g.state_size(), "Selective active var out of range");
                is_active[static_cast<size_t>(v)] = 1;
            }
            for (int i = 0; i < sel.dx.size(); ++i) {
                if (!is_active[static_cast<size_t>(i)]) {
                    require(std::abs(sel.dx[i]) < 1e-12, "Paper-faithful SPO must leave inactive variables fixed");
                }
            }
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.rebuild_from_graph(g, nullptr, 1.0);
            const auto full = fm.solve_cached();
            const auto sel = islam::SelectiveSolver::solve_reduced(fm, {6,7,8}, 0.05);
            require(sel.fell_back_to_full, "Selective solver should fall back to full for tiny alpha");
            require(!sel.used_sparse_partial_solve, "Fallback full solve should bypass sparse partial path");
            require(!sel.used_factor_block_solve, "Fallback full solve should bypass factor-block partial path");
            require((sel.dx - full).norm() < 1e-10, "Fallback full solve mismatch");
        }

        {
            islam::Graph g = make_three_pose_chain();
            islam::FactorManager fm;
            fm.reserve_full_state(g.state_size(), 1.0, 1e-9);
            fm.activate_vars(vars_from_edge(g, g.edges[0]));
            fm.add_edge_contribution(0, fm.build_edge_contribution(g, 0));
            fm.activate_vars(vars_from_edge(g, g.edges[1]));
            fm.add_edge_contribution(1, fm.build_edge_contribution(g, 1));

            islam::Graph g_pert = g;
            g_pert.x[6] = 0.15;
            const auto c_new = fm.build_edge_contribution(g_pert, 1);
            fm.replace_edge_contribution(1, c_new);

            islam::FactorManager fm_ref;
            fm_ref.rebuild_from_graph(g_pert, nullptr, 1.0);
            require((fm.solve_cached() - fm_ref.solve_cached()).norm() < 1e-10,
                    "Incremental growth+replace solve mismatch");
        }


        {
            islam::SymbolicEngine se;
            se.reserve_state(6);
            Eigen::SparseMatrix<double> p0(6,6);
            std::vector<Eigen::Triplet<double>> t0{{0,0,1.0},{1,1,1.0},{2,2,1.0},{3,3,1.0},{4,4,1.0},{5,5,1.0}};
            p0.setFromTriplets(t0.begin(), t0.end());
            se.rebuild_from_numeric_matrix(p0);
            se.refresh_if_needed([](const Eigen::SparseMatrix<double>& A){
                Eigen::AMDOrdering<int> ordering;
                Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P;
                ordering(A.selfadjointView<Eigen::Upper>(), P);
                std::vector<int> perm(static_cast<size_t>(A.rows()));
                for (int i = 0; i < A.rows(); ++i) perm[static_cast<size_t>(i)] = P.indices()[i];
                return perm;
            });
            const auto before = se.snapshot().etree_parent;

            Eigen::SparseMatrix<double> pat(6,6);
            std::vector<Eigen::Triplet<double>> trips{{3,3,1.0},{4,4,1.0},{5,5,1.0},{3,4,1.0},{4,3,1.0},{4,5,1.0},{5,4,1.0}};
            pat.setFromTriplets(trips.begin(), trips.end());
            se.apply_contribution_pattern(pat, {3,4,5}, +1);
            se.refresh_if_needed([](const Eigen::SparseMatrix<double>& A){
                Eigen::AMDOrdering<int> ordering;
                Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P;
                ordering(A.selfadjointView<Eigen::Upper>(), P);
                std::vector<int> perm(static_cast<size_t>(A.rows()));
                for (int i = 0; i < A.rows(); ++i) perm[static_cast<size_t>(i)] = P.indices()[i];
                return perm;
            });
            const auto after = se.snapshot().etree_parent;
            require(before.size() == after.size(), "SymbolicEngine etree size mismatch");
            int changed = 0;
            for (int i = 0; i < static_cast<int>(before.size()); ++i) if (before[static_cast<size_t>(i)] != after[static_cast<size_t>(i)]) ++changed;
            require(changed > 0, "SymbolicEngine should update etree after structural change");
        }

        {
            islam::Graph g = make_three_pose_chain();
            Eigen::VectorXd dx = Eigen::VectorXd::Zero(9);
            dx[6] = 0.25;
            const auto affected = islam::find_affected(g, dx, 1e-6);
            require(!affected.vars.empty(), "Affected vars should not be empty");
            require(!affected.edges.empty(), "Affected edges should not be empty");
            require(std::find(affected.nodes.begin(), affected.nodes.end(), 2) != affected.nodes.end(), "Affected set should include touched node");
            require(affected.nodes.size() >= 2, "Affected set should expand to neighboring nodes");
        }

        {
            Eigen::SparseMatrix<double> pat(5,5);
            std::vector<Eigen::Triplet<double>> trips{
                {0,0,1.0},{1,1,1.0},{2,2,1.0},{3,3,1.0},{4,4,1.0},
                {0,1,1.0},{1,0,1.0},{1,2,1.0},{2,1,1.0},{2,3,1.0},{3,2,1.0},{3,4,1.0},{4,3,1.0}
            };
            pat.setFromTriplets(trips.begin(), trips.end());
            islam::DeterministicBatchCcolamd engine;
            const auto p1 = engine.order(pat);
            const auto p2 = islam::DeterministicBatchCcolamd::order_default(pat);
            require(p1 == p2, "DeterministicBatchCcolamd should be deterministic");
            require(static_cast<int>(p1.size()) == pat.rows(), "DeterministicBatchCcolamd permutation size mismatch");
            std::vector<int> seen(pat.rows(), 0);
            for (int v : p1) {
                require(v >= 0 && v < pat.rows(), "DeterministicBatchCcolamd permutation entry out of range");
                seen[static_cast<size_t>(v)] += 1;
            }
            for (int c : seen) require(c == 1, "DeterministicBatchCcolamd permutation must be bijective");
            require(engine.stats().elimination_steps == static_cast<std::uint64_t>(pat.rows()), "DeterministicBatchCcolamd elimination count mismatch");
            const auto tr1 = engine.order_with_trace(pat);
            const auto tr2 = islam::DeterministicBatchCcolamd::trace_default(pat);
            require(tr1.permutation == tr2.permutation, "Trace default permutation mismatch");
            require(tr1.pattern_signature == tr2.pattern_signature, "Trace pattern signature mismatch");
            require(islam::common_trace_prefix_length(tr1, tr2) == static_cast<int>(tr1.steps.size()), "Identical traces should share full prefix");
            require(!tr1.steps.empty(), "Trace should contain elimination checkpoints");
            require(tr1.steps.back().eliminated_prefix_after == tr1.permutation, "Final trace checkpoint should equal permutation");
            require(tr1.prefix_checkpoints.size() == tr1.steps.size(), "Trace should have one prefix checkpoint per step");
            require(islam::dirty_trace_boundary(tr1, {}) == static_cast<int>(tr1.steps.size()), "Empty dirty set should preserve the whole trace prefix");
            require(islam::dirty_trace_boundary(tr1, {tr1.steps.front().pivot}) == 0, "Dirty first pivot should force rollback to step zero");

            islam::DynamicExactCcolamdPrototype dyn;
            const auto dp1 = dyn.refresh(pat);
            const auto dp2 = dyn.refresh(pat);
            require(dp1 == p1 && dp2 == p1, "DynamicExactCcolamdPrototype should reproduce deterministic batch ordering");
            require(dyn.stats().refreshes == 2, "DynamicExactCcolamdPrototype refresh counter mismatch");
            require(dyn.stats().prefix_reuse_attempts == 1, "DynamicExactCcolamdPrototype prefix reuse attempt mismatch");
            require(dyn.stats().last_common_prefix == static_cast<int>(tr1.steps.size()), "DynamicExactCcolamdPrototype should detect full prefix reuse for unchanged pattern");
            require(dyn.stats().last_dirty_boundary == static_cast<int>(tr1.steps.size()), "DynamicExactCcolamdPrototype empty dirty set should preserve full prefix");
            require(dyn.stats().last_reusable_checkpoint == static_cast<int>(tr1.steps.size()), "DynamicExactCcolamdPrototype should expose full reusable checkpoint");
            require(dyn.stats().last_rollback_suffix == 0, "DynamicExactCcolamdPrototype unchanged pattern should require zero rollback suffix");
            require(dyn.last_reused_prefix() == tr1.permutation, "Full checkpoint replay prefix should equal permutation for unchanged pattern");
            require(dyn.stats().certified_suffix_replays >= 1, "Unchanged pattern should be certified through suffix replay");
            require(dyn.stats().certified_suffix_replay_failures == 0, "Certified suffix replay should not fail for unchanged pattern");
            require(dyn.stats().suffix_replay_pivots_reused >= static_cast<std::uint64_t>(tr1.permutation.size()), "Suffix replay should report reused pivots");
            require(dyn.stats().checkpoint_resume_attempts >= 1, "Unchanged pattern should attempt checkpoint resume");
            require(dyn.stats().checkpoint_resume_successes >= 1, "Unchanged pattern should resume directly from checkpoint");
            require(dyn.stats().checkpoint_bank_entries > 0, "DynamicExactCcolamdPrototype should populate the checkpoint bank");
            require(dyn.stats().checkpoint_bank_insertions > 0, "DynamicExactCcolamdPrototype should record checkpoint-bank insertions");
            require(dyn.stats().dirty_boundary_safety_checks >= 1, "DynamicExactCcolamdPrototype should audit dirty-boundary safety");
            require(dyn.stats().dirty_boundary_unsafe_overestimates == 0, "Unchanged pattern should not overestimate dirty-boundary safety");
            require(dyn.stats().dirty_boundary_exact_matches >= 1, "Unchanged pattern should exactly match dirty-boundary/common-prefix lengths");

            require(!tr1.live_state_checkpoints.empty(), "Trace should expose materialized live-state checkpoints");
            const auto resumed_from_first_checkpoint = engine.order_from_checkpoint(tr1.live_state_checkpoints.front());
            require(resumed_from_first_checkpoint.permutation == tr1.permutation, "Checkpoint resume should reproduce deterministic permutation");

            const auto replay_full = engine.order_with_forced_prefix(pat, tr1.permutation);
            require(replay_full.permutation == tr1.permutation, "Forced full-prefix replay should reproduce the deterministic permutation");
            const auto partial_prefix = tr1.prefix_checkpoints.empty() ? std::vector<int>{} : tr1.prefix_checkpoints.front();
            const auto replay_partial = engine.order_with_forced_prefix(pat, partial_prefix);
            require(replay_partial.permutation == tr1.permutation, "Forced partial-prefix replay should reproduce the deterministic permutation when the prefix is valid");

            const auto dp3 = dyn.refresh(pat, {tr1.steps.front().pivot});
            require(dp3 == p1, "DynamicExactCcolamdPrototype dirty replay should still reproduce deterministic ordering");
            require(dyn.stats().last_dirty_boundary == 0, "Dirty first pivot should make reusable checkpoint zero");
            require(dyn.stats().last_reusable_checkpoint == 0, "Dirty first pivot should disable prefix reuse");
            require(dyn.stats().last_dirty_boundary_overestimate == 0, "Dirty first pivot should not overestimate reusable prefix");
        }

        std::cout << "All smoke tests passed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failure: " << e.what() << "\n";
        return 1;
    }
}
