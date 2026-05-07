#include "islam/gtsam_igg_spo.hpp"

#ifdef ISLAM_HAS_GTSAM

#include <gtsam/linear/GaussianFactorGraph.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

namespace islam {
namespace {

Eigen::SparseMatrix<double> dense_to_sparse(const Eigen::MatrixXd& dense, double drop_tolerance) {
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<std::size_t>(dense.rows()) * 4U);
    for (int j = 0; j < dense.cols(); ++j) {
        for (int i = 0; i < dense.rows(); ++i) {
            const double v = dense(i, j);
            if (std::abs(v) > drop_tolerance) {
                trips.emplace_back(i, j, v);
            }
        }
    }
    Eigen::SparseMatrix<double> out(dense.rows(), dense.cols());
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

Eigen::SparseMatrix<double> gather_sparse(const Eigen::SparseMatrix<double>& A,
                                          const std::vector<int>& rows,
                                          const std::vector<int>& cols) {
    std::vector<int> row_inv(static_cast<std::size_t>(A.rows()), -1);
    std::vector<int> col_inv(static_cast<std::size_t>(A.cols()), -1);
    for (int i = 0; i < static_cast<int>(rows.size()); ++i) row_inv[static_cast<std::size_t>(rows[i])] = i;
    for (int j = 0; j < static_cast<int>(cols.size()); ++j) col_inv[static_cast<std::size_t>(cols[j])] = j;

    std::vector<Eigen::Triplet<double>> trips;
    for (int col = 0; col < A.outerSize(); ++col) {
        const int new_col = col_inv[static_cast<std::size_t>(col)];
        if (new_col < 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int new_row = row_inv[static_cast<std::size_t>(it.row())];
            if (new_row >= 0) trips.emplace_back(new_row, new_col, it.value());
        }
    }
    Eigen::SparseMatrix<double> out(static_cast<int>(rows.size()), static_cast<int>(cols.size()));
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

Eigen::VectorXd gather_vector(const Eigen::VectorXd& v, const std::vector<int>& idx) {
    Eigen::VectorXd out(static_cast<int>(idx.size()));
    for (int i = 0; i < static_cast<int>(idx.size()); ++i) out[i] = v[idx[static_cast<std::size_t>(i)]];
    return out;
}

std::vector<int> complement_sorted(const std::vector<int>& subset, int n) {
    std::vector<unsigned char> in_subset(static_cast<std::size_t>(n), 0);
    for (int i : subset) {
        if (i >= 0 && i < n) in_subset[static_cast<std::size_t>(i)] = 1;
    }
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(std::max(0, n - static_cast<int>(subset.size()))));
    for (int i = 0; i < n; ++i) {
        if (!in_subset[static_cast<std::size_t>(i)]) out.push_back(i);
    }
    return out;
}

std::vector<int> scalar_indices_for_keys(const std::map<gtsam::Key, GtsamIggSpoSolver::KeyBlock>& by_key,
                                         const std::set<gtsam::Key>& keys) {
    std::vector<int> out;
    for (gtsam::Key key : keys) {
        const auto it = by_key.find(key);
        if (it == by_key.end()) continue;
        for (int k = 0; k < it->second.dim; ++k) out.push_back(it->second.offset + k);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<int> all_scalar_indices(int n) {
    std::vector<int> out(static_cast<std::size_t>(n));
    std::iota(out.begin(), out.end(), 0);
    return out;
}

} // namespace

GtsamIggSpoSolver::GtsamIggSpoSolver(GtsamIggSpoParams params)
    : params_(params) {}

void GtsamIggSpoSolver::reset(const gtsam::NonlinearFactorGraph& graph, const gtsam::Values& values) {
    graph_ = graph;
    values_ = values;
    previous_eta_ = 0.0;
    previous_dim_ = 0;
    have_previous_eta_ = false;
    increment_ = 0;
}

void GtsamIggSpoSolver::insert_new_values(const gtsam::Values& new_values) {
    for (const auto& kv : new_values) {
        if (!values_.exists(kv.key)) {
            values_.insert(kv.key, kv.value);
        }
    }
}

std::set<gtsam::Key> GtsamIggSpoSolver::all_keys() const {
    std::set<gtsam::Key> out;
    for (gtsam::Key key : values_.keys()) out.insert(key);
    return out;
}

std::set<gtsam::Key> GtsamIggSpoSolver::keys_from_factors(const gtsam::NonlinearFactorGraph& factors) const {
    std::set<gtsam::Key> out;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        const auto factor = factors.at(i);
        if (!factor) continue;
        for (gtsam::Key key : factor->keys()) out.insert(key);
    }
    return out;
}

GtsamIggSpoSolver::LinearizedSystem GtsamIggSpoSolver::linearize_current() const {
    LinearizedSystem sys;
    const auto linear = graph_.linearize(values_);
    if (!linear) {
        throw std::runtime_error("GtsamIggSpoSolver: GTSAM linearization returned null");
    }

    sys.ordering = graph_.orderingCOLAMD();
    const auto hess = linear->hessian(sys.ordering);
    Eigen::MatrixXd H_dense = hess.first;
    Eigen::VectorXd b = hess.second;
    if (H_dense.rows() != H_dense.cols() || H_dense.rows() != b.size()) {
        throw std::runtime_error("GtsamIggSpoSolver: inconsistent Hessian dimensions from GTSAM");
    }
    if (params_.diagonal_damping > 0.0) {
        H_dense.diagonal().array() += params_.diagonal_damping;
    }

    const auto dims = linear->getKeyDimMap();
    int offset = 0;
    sys.blocks.reserve(sys.ordering.size());
    for (std::size_t i = 0; i < sys.ordering.size(); ++i) {
        const gtsam::Key key = sys.ordering[i];
        const auto dit = dims.find(key);
        if (dit == dims.end()) {
            throw std::runtime_error("GtsamIggSpoSolver: missing key dimension in GTSAM linear graph");
        }
        KeyBlock block;
        block.key = key;
        block.offset = offset;
        block.dim = static_cast<int>(dit->second);
        offset += block.dim;
        sys.blocks.push_back(block);
        sys.by_key.emplace(key, block);
    }
    if (offset != H_dense.rows()) {
        throw std::runtime_error("GtsamIggSpoSolver: key dimensions do not sum to Hessian size");
    }

    sys.H = dense_to_sparse(H_dense, params_.sparse_drop_tolerance);
    sys.b = b;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
    llt.compute(sys.H);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("GtsamIggSpoSolver: full sparse LLT failed while computing information eta");
    }
    const Eigen::SparseMatrix<double> L = Eigen::SparseMatrix<double>(llt.matrixL());
    double eta = 0.0;
    for (int i = 0; i < std::min(L.rows(), L.cols()); ++i) {
        const double d = L.coeff(i, i);
        if (std::abs(d) > 0.0) eta += std::log(std::abs(d));
    }
    sys.eta = eta;
    return sys;
}

GtsamIggSpoSolver::PartialSolveResult GtsamIggSpoSolver::solve_active(
    const LinearizedSystem& sys,
    const std::set<gtsam::Key>& active_keys) const {
    PartialSolveResult out;
    const int n = static_cast<int>(sys.H.cols());
    out.dx = Eigen::VectorXd::Zero(n);
    if (n == 0 || active_keys.empty()) return out;

    std::vector<int> active_scalars = scalar_indices_for_keys(sys.by_key, active_keys);
    if (active_scalars.empty()) return out;

    if (params_.always_full_solve ||
        (params_.selective_alpha > 0.0 &&
         static_cast<double>(active_scalars.size()) > params_.selective_alpha * static_cast<double>(n))) {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
        llt.compute(sys.H);
        if (llt.info() != Eigen::Success) {
            throw std::runtime_error("GtsamIggSpoSolver: full sparse LLT solve failed");
        }
        out.dx = llt.solve(sys.b);
        out.solved_keys = all_keys();
        out.fell_back_to_full = true;
        return out;
    }

    const std::vector<int> static_scalars = complement_sorted(active_scalars, n);
    std::vector<int> local_order;
    local_order.reserve(static_scalars.size() + active_scalars.size());
    local_order.insert(local_order.end(), static_scalars.begin(), static_scalars.end());
    local_order.insert(local_order.end(), active_scalars.begin(), active_scalars.end());

    const int nu = static_cast<int>(static_scalars.size());
    const int ns = static_cast<int>(active_scalars.size());
    const Eigen::SparseMatrix<double> H_local = gather_sparse(sys.H, local_order, local_order);
    const Eigen::VectorXd b_local = gather_vector(sys.b, local_order);

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> llt;
    llt.compute(H_local);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("GtsamIggSpoSolver: sparse active-set LLT failed");
    }

    Eigen::SparseMatrix<double> L = Eigen::SparseMatrix<double>(llt.matrixL());
    L.makeCompressed();

    Eigen::VectorXd rhsS = b_local.tail(ns);
    if (nu > 0) {
        const Eigen::SparseMatrix<double> LUU = L.block(0, 0, nu, nu);
        const Eigen::SparseMatrix<double> LSU = L.block(nu, 0, ns, nu);
        const Eigen::VectorXd bU = b_local.head(nu);
        const Eigen::VectorXd yU = LUU.triangularView<Eigen::Lower>().solve(bU);
        rhsS.noalias() -= LSU * yU;
    }

    const Eigen::SparseMatrix<double> LSS = L.block(nu, nu, ns, ns);
    const Eigen::VectorXd zS = LSS.triangularView<Eigen::Lower>().solve(rhsS);
    const Eigen::SparseMatrix<double> USS = LSS.transpose();
    const Eigen::VectorXd xS = USS.triangularView<Eigen::Upper>().solve(zS);

    for (int i = 0; i < ns; ++i) {
        out.dx[active_scalars[static_cast<std::size_t>(i)]] = xS[i];
    }
    out.solved_keys = active_keys;
    return out;
}

std::set<gtsam::Key> GtsamIggSpoSolver::prune_by_increment(
    const LinearizedSystem& sys,
    const Eigen::VectorXd& dx,
    const std::set<gtsam::Key>& keys) const {
    std::set<gtsam::Key> kept;
    for (gtsam::Key key : keys) {
        const auto it = sys.by_key.find(key);
        if (it == sys.by_key.end()) continue;
        const auto& b = it->second;
        if (b.offset < 0 || b.offset + b.dim > dx.size()) continue;
        if (dx.segment(b.offset, b.dim).lpNorm<Eigen::Infinity>() > params_.dx_threshold) {
            kept.insert(key);
        }
    }
    return kept;
}

std::set<gtsam::Key> GtsamIggSpoSolver::expand_through_graph(const std::set<gtsam::Key>& keys) const {
    if (keys.empty()) return {};
    std::set<gtsam::Key> out = keys;
    for (std::size_t i = 0; i < graph_.size(); ++i) {
        const auto factor = graph_.at(i);
        if (!factor) continue;
        bool touches = false;
        for (gtsam::Key key : factor->keys()) {
            if (keys.find(key) != keys.end()) {
                touches = true;
                break;
            }
        }
        if (!touches) continue;
        for (gtsam::Key key : factor->keys()) out.insert(key);
    }
    return out;
}

gtsam::VectorValues GtsamIggSpoSolver::vector_values_from_dense(
    const LinearizedSystem& sys,
    const Eigen::VectorXd& dx,
    const std::set<gtsam::Key>& keys) const {
    gtsam::VectorValues delta;
    for (gtsam::Key key : keys) {
        const auto it = sys.by_key.find(key);
        if (it == sys.by_key.end()) continue;
        const auto& b = it->second;
        if (b.offset < 0 || b.offset + b.dim > dx.size()) continue;
        delta.insert(key, dx.segment(b.offset, b.dim));
    }
    return delta;
}

GtsamIggSpoUpdateStats GtsamIggSpoSolver::update(const gtsam::NonlinearFactorGraph& new_factors,
                                                 const gtsam::Values& new_values) {
    for (std::size_t i = 0; i < new_factors.size(); ++i) {
        const auto factor = new_factors.at(i);
        if (factor) graph_.push_back(factor);
    }
    insert_new_values(new_values);

    GtsamIggSpoUpdateStats stats;
    stats.increment = increment_++;
    stats.num_factors = graph_.size();
    stats.num_keys = values_.size();
    stats.tangent_dim = values_.dim();

    LinearizedSystem sys = linearize_current();
    stats.eta = sys.eta;
    if (have_previous_eta_ && previous_dim_ > 0 && stats.tangent_dim > 0) {
        stats.delta_eta = stats.eta -
            (static_cast<double>(previous_dim_) / static_cast<double>(stats.tangent_dim)) * previous_eta_;
    } else {
        stats.delta_eta = 0.0;
    }
    // Do not update previous_eta_ here. The information state carried to the
    // next increment should correspond to the final post-GN linearization
    // point, not just the first linearization after appending new factors.

    std::set<gtsam::Key> active = keys_from_factors(new_factors);
    if (params_.force_full_update || stats.delta_eta >= params_.eta_threshold || active.empty()) {
        active = all_keys();
        stats.global_gate = true;
    }

    for (int iter = 0; iter < params_.max_gn_iterations && !active.empty(); ++iter) {
        sys = linearize_current();
        PartialSolveResult solve = solve_active(sys, active);
        std::set<gtsam::Key> solved_keys = solve.fell_back_to_full ? all_keys() : solve.solved_keys;
        std::set<gtsam::Key> kept = prune_by_increment(sys, solve.dx, solved_keys);

        GtsamIggSpoIterationStats is;
        is.iteration = iter;
        is.candidate_keys = active.size();
        is.active_keys = kept.size();
        is.active_scalars = scalar_indices_for_keys(sys.by_key, kept).size();
        is.fell_back_to_full = solve.fell_back_to_full;
        is.dx_inf = solve.dx.size() > 0 ? solve.dx.lpNorm<Eigen::Infinity>() : 0.0;
        stats.iterations.push_back(is);
        stats.gn_iterations = iter + 1;

        if (kept.empty()) {
            break;
        }

        const gtsam::VectorValues delta = vector_values_from_dense(sys, solve.dx, kept);
        values_ = values_.retract(delta);
        active = expand_through_graph(kept);
    }

    stats.final_error = graph_.error(values_);
    stats.num_factors = graph_.size();
    stats.num_keys = values_.size();
    stats.tangent_dim = values_.dim();

    // Refresh the stored information surrogate from the factor that corresponds
    // to the final maintained Values after selective retraction/relinearization.
    const LinearizedSystem final_sys = linearize_current();
    previous_eta_ = final_sys.eta;
    previous_dim_ = stats.tangent_dim;
    have_previous_eta_ = true;

    return stats;
}

} // namespace islam

#endif // ISLAM_HAS_GTSAM
