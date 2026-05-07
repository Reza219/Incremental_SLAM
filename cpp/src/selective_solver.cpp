#include "islam/selective_solver.hpp"

#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace islam {
namespace {

std::vector<int> invert_perm_checked(const std::vector<int>& p) {
    std::vector<int> pinv(p.size(), -1);
    for (int k = 0; k < static_cast<int>(p.size()); ++k) {
        const int pk = p[static_cast<size_t>(k)];
        if (pk < 0 || pk >= static_cast<int>(p.size())) {
            throw std::runtime_error("SelectiveSolver: invalid permutation entry");
        }
        pinv[static_cast<size_t>(pk)] = k;
    }
    return pinv;
}

Eigen::VectorXd gather_vector(const Eigen::VectorXd& v, const std::vector<int>& idx) {
    Eigen::VectorXd out(static_cast<int>(idx.size()));
    for (int i = 0; i < static_cast<int>(idx.size()); ++i) out[i] = v[idx[static_cast<size_t>(i)]];
    return out;
}

Eigen::SparseMatrix<double> gather_sparse(const Eigen::SparseMatrix<double>& A,
                                          const std::vector<int>& rows,
                                          const std::vector<int>& cols) {
    std::vector<int> row_inv(static_cast<size_t>(A.rows()), -1);
    std::vector<int> col_inv(static_cast<size_t>(A.cols()), -1);
    for (int i = 0; i < static_cast<int>(rows.size()); ++i) row_inv[rows[static_cast<size_t>(i)]] = i;
    for (int j = 0; j < static_cast<int>(cols.size()); ++j) col_inv[cols[static_cast<size_t>(j)]] = j;

    std::vector<Eigen::Triplet<double>> trips;
    for (int col = 0; col < A.outerSize(); ++col) {
        const int new_col = col_inv[static_cast<size_t>(col)];
        if (new_col < 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int new_row = row_inv[static_cast<size_t>(it.row())];
            if (new_row >= 0) trips.emplace_back(new_row, new_col, it.value());
        }
    }

    Eigen::SparseMatrix<double> out(static_cast<int>(rows.size()), static_cast<int>(cols.size()));
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

std::vector<int> complement_sorted(const std::vector<int>& subset, int n) {
    std::vector<char> in_subset(static_cast<size_t>(n), 0);
    for (int i : subset) {
        if (i >= 0 && i < n) in_subset[static_cast<size_t>(i)] = 1;
    }
    std::vector<int> out;
    out.reserve(static_cast<size_t>(std::max(0, n - static_cast<int>(subset.size()))));
    for (int i = 0; i < n; ++i) {
        if (!in_subset[static_cast<size_t>(i)]) out.push_back(i);
    }
    return out;
}

std::vector<int> all_indices(int n) {
    std::vector<int> out(static_cast<size_t>(n));
    std::iota(out.begin(), out.end(), 0);
    return out;
}

bool try_factor_block_suffix_solve(const FactorManager& fm,
                                   const std::vector<int>& affected_vars,
                                   double alpha,
                                   SelectiveSolveResult& out) {
    const FactorBlockSystem fbs = fm.build_factor_block_system();
    const int n = static_cast<int>(fbs.perm.size());
    if (!fbs.available || n == 0 || fbs.L_factor.rows() != n || fbs.L_factor.cols() != n ||
        fbs.g_factor.size() != n || static_cast<int>(fbs.pinv.size()) != n) {
        return false;
    }

    std::vector<int> affected_perm;
    affected_perm.reserve(affected_vars.size());
    for (int v : affected_vars) {
        if (v < 0 || v >= n) continue;
        const int pv = fbs.pinv[static_cast<size_t>(v)];
        if (pv >= 0) affected_perm.push_back(pv);
    }
    std::sort(affected_perm.begin(), affected_perm.end());
    affected_perm.erase(std::unique(affected_perm.begin(), affected_perm.end()), affected_perm.end());
    if (affected_perm.empty()) {
        out.dx = Eigen::VectorXd::Zero(n);
        return true;
    }

    std::vector<int> reach = SelectiveSolver::etree_reach(fbs.etree_parent, affected_perm, n);
    if (reach.empty()) {
        reach = affected_perm;
    }
    const int boundary = std::max(0, *std::min_element(reach.begin(), reach.end()));
    out.active_perm.clear();
    out.active_perm.reserve(static_cast<size_t>(n - boundary));
    for (int p = boundary; p < n; ++p) {
        out.active_perm.push_back(p);
    }
    out.enlarged_to_factor_suffix = (out.active_perm.size() != reach.size());

    out.active_vars.clear();
    out.active_vars.reserve(out.active_perm.size());
    for (int p : out.active_perm) {
        out.active_vars.push_back(fbs.perm[static_cast<size_t>(p)]);
    }
    std::sort(out.active_vars.begin(), out.active_vars.end());
    out.active_vars.erase(std::unique(out.active_vars.begin(), out.active_vars.end()), out.active_vars.end());

    if (alpha > 0.0 && static_cast<double>(out.active_perm.size()) > alpha * static_cast<double>(n)) {
        out.dx = fm.solve_cached();
        out.active_perm = all_indices(n);
        out.active_vars = all_indices(n);
        out.fell_back_to_full = true;
        return true;
    }

    const int nu = boundary;
    const int ns = n - boundary;
    if (ns <= 0) {
        out.dx = Eigen::VectorXd::Zero(n);
        return true;
    }

    const Eigen::VectorXd bS = fbs.g_factor.segment(nu, ns);
    Eigen::VectorXd rhsS = bS;

    if (nu > 0) {
        const Eigen::SparseMatrix<double> LUU = fbs.L_factor.block(0, 0, nu, nu);
        const Eigen::SparseMatrix<double> LSU = fbs.L_factor.block(nu, 0, ns, nu);
        const Eigen::VectorXd bU = fbs.g_factor.head(nu);
        const Eigen::VectorXd yU = LUU.triangularView<Eigen::Lower>().solve(bU);
        rhsS.noalias() -= LSU * yU;
    }

    const Eigen::SparseMatrix<double> LSS = fbs.L_factor.block(nu, nu, ns, ns);
    const Eigen::VectorXd zS = LSS.triangularView<Eigen::Lower>().solve(rhsS);
    const Eigen::SparseMatrix<double> USS = LSS.transpose();
    const Eigen::VectorXd xS = USS.triangularView<Eigen::Upper>().solve(zS);

    Eigen::VectorXd x_factor = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < ns; ++i) {
        x_factor[nu + i] = xS[i];
    }

    out.dx = Eigen::VectorXd::Zero(n);
    for (int p = 0; p < n; ++p) {
        out.dx[fbs.perm[static_cast<size_t>(p)]] = x_factor[p];
    }
    out.used_sparse_partial_solve = true;
    out.used_factor_block_solve = true;
    return true;
}

} // namespace

std::vector<int> SelectiveSolver::elimination_tree(const Eigen::SparseMatrix<double>& Hperm) {
    const int n = Hperm.cols();
    std::vector<int> parent(static_cast<size_t>(n), -1);
    std::vector<int> ancestor(static_cast<size_t>(n), -1);

    for (int k = 0; k < n; ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Hperm, k); it; ++it) {
            int i = it.row();
            if (i >= k) continue;
            while (i != -1 && i < k) {
                const int inext = ancestor[static_cast<size_t>(i)];
                ancestor[static_cast<size_t>(i)] = k;
                if (inext == -1) parent[static_cast<size_t>(i)] = k;
                i = inext;
            }
        }
    }
    return parent;
}

std::vector<int> SelectiveSolver::etree_reach(const std::vector<int>& parent,
                                              const std::vector<int>& seeds,
                                              int n) {
    std::vector<char> marked(static_cast<size_t>(n), 0);
    std::vector<int> reach;
    for (int s : seeds) {
        int j = s;
        while (j >= 0 && j < n && !marked[static_cast<size_t>(j)]) {
            marked[static_cast<size_t>(j)] = 1;
            reach.push_back(j);
            j = parent[static_cast<size_t>(j)];
        }
    }
    std::sort(reach.begin(), reach.end());
    reach.erase(std::unique(reach.begin(), reach.end()), reach.end());
    return reach;
}

SelectiveSolveResult SelectiveSolver::solve_reduced(const FactorManager& fm,
                                                    const std::vector<int>& affected_vars,
                                                    double alpha) {
    // Preferred paper-grade path: reuse the maintained factor itself.  We enlarge
    // the active set to a suffix in the factor ordering so the existing lower
    // triangular factor has the block form [L_UU 0; L_SU L_SS].  This avoids
    // locally refactorizing H_{local}; only sparse triangular subsolves on cached
    // factor blocks are performed.  If the factor is unavailable, the legacy
    // sparse local refactorization path below remains the portable fallback.
    SelectiveSolveResult factor_block_out;
    try {
        if (try_factor_block_suffix_solve(fm, affected_vars, alpha, factor_block_out)) {
            return factor_block_out;
        }
    } catch (const std::exception&) {
        // Fall back to the robust sparse local-refactorization path.  The fallback
        // preserves correctness on platforms where CHOLMOD exposes an unexpected
        // factor representation or when only the portable Eigen backend is used.
    }

    const auto symbolic = fm.build_symbolic_system();
    const auto& ps = symbolic.permuted;
    const int n = static_cast<int>(ps.perm.size());
    SelectiveSolveResult out;
    out.dx = Eigen::VectorXd::Zero(n);

    if (n == 0 || affected_vars.empty()) {
        return out;
    }

    const auto pinv = invert_perm_checked(ps.perm);
    std::vector<int> affected_perm;
    affected_perm.reserve(affected_vars.size());
    for (int v : affected_vars) {
        if (v < 0 || v >= n) continue;
        const int pv = pinv[static_cast<size_t>(v)];
        if (pv >= 0) affected_perm.push_back(pv);
    }
    std::sort(affected_perm.begin(), affected_perm.end());
    affected_perm.erase(std::unique(affected_perm.begin(), affected_perm.end()), affected_perm.end());
    if (affected_perm.empty()) {
        return out;
    }

    out.active_perm = etree_reach(symbolic.etree_parent, affected_perm, n);
    for (int p : out.active_perm) {
        if (p >= 0 && p < n) out.active_vars.push_back(ps.perm[static_cast<size_t>(p)]);
    }
    std::sort(out.active_vars.begin(), out.active_vars.end());
    out.active_vars.erase(std::unique(out.active_vars.begin(), out.active_vars.end()), out.active_vars.end());

    if (alpha > 0.0 && static_cast<double>(out.active_perm.size()) > alpha * static_cast<double>(n)) {
        out.dx = fm.solve_cached();
        out.active_perm = all_indices(n);
        out.active_vars = all_indices(n);
        out.fell_back_to_full = true;
        return out;
    }

    const std::vector<int> static_perm = complement_sorted(out.active_perm, n);

    // Paper-faithful SPO solve.  The paper locally groups the current symbolic
    // ordering into static variables U and active variables S, factors the
    // locally ordered normal matrix as R^T R, and computes only d_S while keeping
    // d_U = 0. The implementation must not reconstruct or return a d_U update.
    // The static block contributes through the triangular correction
    // R_US^T y_U, where R_UU^T y_U = b_U.
    std::vector<int> local_order;
    local_order.reserve(static_perm.size() + out.active_perm.size());
    local_order.insert(local_order.end(), static_perm.begin(), static_perm.end());
    local_order.insert(local_order.end(), out.active_perm.begin(), out.active_perm.end());

    const int nu = static_cast<int>(static_perm.size());
    const int ns = static_cast<int>(out.active_perm.size());

    const Eigen::SparseMatrix<double> H_local_sparse = gather_sparse(ps.H_perm, local_order, local_order);
    const Eigen::VectorXd b_local = gather_vector(ps.g_perm, local_order);

    // Sparse paper-faithful SPO solve. Earlier milestones converted the locally
    // ordered matrix to a dense MatrixXd and used dense LLT. That was useful for
    // semantic validation but defeated the sparse complexity model of the paper.
    // Here we keep the local normal matrix sparse, factor it with a natural-order
    // sparse LLT, and perform the two triangular subsolves directly on sparse
    // blocks of L = R^T. Static variables receive no update.
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> llt;
    llt.compute(H_local_sparse);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("SelectiveSolver::solve_reduced: sparse local LLT factorization failed");
    }

    Eigen::SparseMatrix<double> L = Eigen::SparseMatrix<double>(llt.matrixL());
    L.makeCompressed();

    const Eigen::VectorXd bS = b_local.tail(ns);
    Eigen::VectorXd rhsS = bS;

    if (nu > 0) {
        const Eigen::SparseMatrix<double> LUU = L.block(0, 0, nu, nu);
        const Eigen::SparseMatrix<double> LSU = L.block(nu, 0, ns, nu);
        const Eigen::VectorXd bU = b_local.head(nu);
        const Eigen::VectorXd yU = LUU.triangularView<Eigen::Lower>().solve(bU);
        rhsS.noalias() -= LSU * yU;
    }

    const Eigen::SparseMatrix<double> LSS = L.block(nu, nu, ns, ns);
    const Eigen::VectorXd z = LSS.triangularView<Eigen::Lower>().solve(rhsS);
    const Eigen::SparseMatrix<double> RSS = LSS.transpose();
    const Eigen::VectorXd xS = RSS.triangularView<Eigen::Upper>().solve(z);
    out.used_sparse_partial_solve = true;

    Eigen::VectorXd x_perm = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < ns; ++i) {
        x_perm[out.active_perm[static_cast<size_t>(i)]] = xS[i];
    }

    for (int k = 0; k < n; ++k) {
        out.dx[ps.perm[static_cast<size_t>(k)]] = x_perm[k];
    }

    return out;
}

} // namespace islam
