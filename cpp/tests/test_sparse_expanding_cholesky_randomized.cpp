
#include "islam/sparse_expanding_cholesky.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Sparse = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

[[noreturn]] void fail(const std::string& msg) { throw std::runtime_error(msg); }

Sparse make_initial_hessian(int n) {
    std::vector<Triplet> trips;
    for (int i = 0; i < n; ++i) trips.emplace_back(i, i, 12.0);
    for (int i = 0; i + 1 < n; ++i) {
        trips.emplace_back(i, i, 0.7);
        trips.emplace_back(i + 1, i + 1, 0.7);
        trips.emplace_back(i, i + 1, -0.7);
        trips.emplace_back(i + 1, i, -0.7);
    }
    for (int i = 0; i + 2 < n; i += 2) {
        trips.emplace_back(i, i, 0.2);
        trips.emplace_back(i + 2, i + 2, 0.2);
        trips.emplace_back(i, i + 2, -0.2);
        trips.emplace_back(i + 2, i, -0.2);
    }
    Sparse H(n, n);
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();
    return H;
}


Sparse make_diagonal_hessian(int n, double diag) {
    std::vector<Triplet> trips;
    for (int i = 0; i < n; ++i) trips.emplace_back(i, i, diag);
    Sparse H(n, n);
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();
    return H;
}

Sparse make_column(int n, const std::vector<std::pair<int, double>>& entries) {
    std::vector<Triplet> trips;
    for (const auto& e : entries) trips.emplace_back(e.first, 0, e.second);
    Sparse c(n, 1);
    c.setFromTriplets(trips.begin(), trips.end());
    c.makeCompressed();
    return c;
}

Sparse resize_column_preserve(const Sparse& column, int rows) {
    if (column.rows() == rows) return column;
    std::vector<Triplet> trips;
    for (int col = 0; col < column.outerSize(); ++col) {
        for (Sparse::InnerIterator it(column, col); it; ++it) {
            trips.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
        }
    }
    Sparse out(rows, column.cols());
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

void check_backend(const islam::SparseExpandingCholesky& dyn,
                   const std::string& label,
                   double tol = 2e-8) {
    if (!dyn.covers_state()) fail(label + ": factor does not cover state");

    double residual = 0.0;
    if (!dyn.passes_factorization_residual_check(&residual) || residual > tol) {
        std::ostringstream oss;
        oss << label << ": residual too large: " << residual;
        fail(oss.str());
    }

    const int n = dyn.state_size();
    Eigen::VectorXd rhs(n);
    for (int i = 0; i < n; ++i) rhs[i] = 0.5 + 0.17 * static_cast<double>((i % 5) - 2);

    islam::SparseExpandingCholesky fresh(dyn.options());
    fresh.factorize(dyn.normal_matrix());
    const Eigen::VectorXd x_dyn = dyn.solve(rhs);
    const Eigen::VectorXd x_full = fresh.solve(rhs);
    const double rel = (x_dyn - x_full).norm() / (1.0 + x_full.norm());
    if (!(std::isfinite(rel) && rel <= 2e-8)) {
        std::ostringstream oss;
        oss << label << ": solve disagreement vs fresh refactorization: " << rel;
        fail(oss.str());
    }
}

struct StoredContribution {
    Sparse column;
    std::vector<int> touched;
};
} // namespace

int main() {
    islam::SparseExpandingCholeskyOptions options;
    options.certification_tolerance = 2e-9;
    options.jitter = 1e-12;
    options.max_jitter_tries = 8;
    options.use_factor_dependency_cache = true;
    options.use_column_local_certification = true;
    options.full_certification_fallback = true;

    islam::SparseExpandingCholesky local(options);
    local.factorize(make_diagonal_hessian(8, 10.0));
    local.apply_diagonal_update({2}, 0.25, true);
    check_backend(local, "local closure diagonal update");
    if (local.stats().affected_closure_refactorizations == 0) {
        fail("local closure smoke did not exercise affected-closure refactorization");
    }

    islam::SparseExpandingCholesky dyn(options);
    dyn.factorize(make_initial_hessian(8));
    check_backend(dyn, "initial");

    std::mt19937 rng(78);
    std::uniform_int_distribution<int> op_dist(0, 99);
    std::vector<StoredContribution> removable;

    for (int step = 0; step < 48; ++step) {
        const int n = dyn.state_size();
        const int op = op_dist(rng);

        if (step == 12 || step == 28) {
            const int old_n = dyn.state_size();
            dyn.grow_to(old_n + 1);
            dyn.apply_contribution(make_column(old_n + 1, {{old_n, std::sqrt(9.0)}}), {old_n}, true);
            const Sparse link = make_column(old_n + 1, {{old_n - 1, std::sqrt(0.35)}, {old_n, -std::sqrt(0.35)}});
            dyn.apply_contribution(link, {old_n - 1, old_n}, true);
            removable.push_back({link, {old_n - 1, old_n}});
        } else if (op < 35) {
            std::uniform_int_distribution<int> vdist(0, n - 1);
            dyn.apply_diagonal_update({vdist(rng)}, 0.05 + 0.01 * static_cast<double>(step % 5), true);
        } else if (op < 75 || removable.empty()) {
            std::uniform_int_distribution<int> vdist(0, n - 1);
            int i = vdist(rng);
            int j = vdist(rng);
            if (i == j) j = (j + 1) % n;
            if (i > j) std::swap(i, j);
            const double w = 0.08 + 0.01 * static_cast<double>(step % 7);
            const Sparse c = make_column(n, {{i, std::sqrt(w)}, {j, -std::sqrt(w)}});
            dyn.apply_contribution(c, {i, j}, true);
            removable.push_back({c, {i, j}});
        } else {
            std::uniform_int_distribution<int> rdist(0, static_cast<int>(removable.size()) - 1);
            const int idx = rdist(rng);
            const auto c = removable[static_cast<size_t>(idx)];
            dyn.apply_contribution(resize_column_preserve(c.column, dyn.state_size()), c.touched, false);
            removable.erase(removable.begin() + idx);
        }

        check_backend(dyn, "step " + std::to_string(step));
    }

    const auto& st = dyn.stats();
    if (st.numeric_only_updates == 0) fail("did not exercise numeric-only updates");
    if (st.structural_pattern_changes == 0) fail("did not exercise structural pattern changes");
    if (st.expansion_suffix_updates == 0) fail("did not exercise expansion suffix updates");
    if (st.dependency_cache_hits == 0) fail("did not exercise dependency-cache closure discovery");

    std::cout << "M81 randomized sparse expanding Cholesky differential regression passed\n";
    std::cout << "  final state size: " << dyn.state_size() << "\n";
    std::cout << "  numeric-only updates: " << st.numeric_only_updates << "\n";
    std::cout << "  structural pattern changes: " << st.structural_pattern_changes << "\n";
    std::cout << "  expansion suffix updates: " << st.expansion_suffix_updates << "\n";
    std::cout << "  affected-closure refactorizations: " << st.affected_closure_refactorizations << "\n";
    std::cout << "  suffix refactorizations: " << st.suffix_refactorizations << "\n";
    std::cout << "  dependency cache hits: " << st.dependency_cache_hits << "\n";
    std::cout << "  residual checks: " << st.factorization_residual_checks << "\n";
    std::cout << "  final scaled residual: " << dyn.scaled_factorization_residual() << "\n";
    return 0;
}
