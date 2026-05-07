#include "islam/sparse_expanding_cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace islam {
namespace {

using Sparse = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

void add_sparse_triplets_local(const Sparse& A, std::vector<Triplet>& trips) {
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Sparse::InnerIterator it(A, col); it; ++it) {
            trips.emplace_back(it.row(), it.col(), it.value());
        }
    }
}

std::vector<int> unique_sorted_local(std::vector<int> vals) {
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    return vals;
}

void prune_and_compress(Sparse& A, double tol) {
    if (tol > 0.0) {
        A.prune(tol);
    }
    A.makeCompressed();
}

struct AveragedValue {
    double sum = 0.0;
    int count = 0;
    void add(double v) {
        sum += v;
        ++count;
    }
    [[nodiscard]] double value() const {
        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    }
};

double max_abs_sparse_local(const Sparse& A) {
    double m = 0.0;
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Sparse::InnerIterator it(A, col); it; ++it) {
            m = std::max(m, std::abs(it.value()));
        }
    }
    return m;
}

} // namespace

SparseExpandingCholesky::SparseExpandingCholesky(SparseExpandingCholeskyOptions options)
    : options_(options) {}

void SparseExpandingCholesky::clear() {
    n_ = 0;
    factor_size_ = 0;
    factorized_ = false;
    H_.resize(0, 0);
    L_.resize(0, 0);
    factor_dependency_cache_valid_ = false;
    factor_dependency_children_.clear();
    stats_ = SparseExpandingCholeskyStats{};
}

Sparse SparseExpandingCholesky::resize_square_preserve(const Sparse& A, int n) {
    if (n < A.rows() || n < A.cols()) {
        throw std::runtime_error("SparseExpandingCholesky does not support shrinking");
    }
    std::vector<Triplet> trips;
    trips.reserve(static_cast<size_t>(A.nonZeros()));
    add_sparse_triplets_local(A, trips);
    Sparse out(n, n);
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

void SparseExpandingCholesky::grow_to(int n) {
    if (n < n_) {
        throw std::runtime_error("SparseExpandingCholesky::grow_to does not support shrinking");
    }
    if (n == n_) {
        return;
    }
    H_ = resize_square_preserve(H_, n);
    L_ = resize_square_preserve(L_, n);
    n_ = n;
    invalidate_factor_dependency_cache();
    ++stats_.grow_calls;
}

void SparseExpandingCholesky::factorize(const Sparse& H) {
    if (H.rows() != H.cols()) {
        throw std::runtime_error("SparseExpandingCholesky::factorize requires a square matrix");
    }
    H_ = H;
    H_.makeCompressed();
    n_ = H.rows();
    L_.resize(n_, n_);
    factorized_ = false;
    factor_size_ = 0;
    invalidate_factor_dependency_cache();
    refactor_suffix_from(0, true);
}

void SparseExpandingCholesky::refactorize_with_reused_prefix(const Sparse& H,
                                                             int reusable_prefix) {
    if (H.rows() != H.cols()) {
        throw std::runtime_error("SparseExpandingCholesky::refactorize_with_reused_prefix requires a square matrix");
    }
    const int new_n = H.rows();
    reusable_prefix = std::clamp(reusable_prefix, 0, new_n);
    if (!factorized_ || factor_size_ < reusable_prefix || L_.rows() < reusable_prefix || L_.cols() < reusable_prefix) {
        reusable_prefix = 0;
    }

    H_ = H;
    H_.makeCompressed();
    n_ = new_n;
    if (L_.rows() != n_ || L_.cols() != n_) {
        L_ = resize_square_preserve(L_, n_);
    }

    if (n_ == 0) {
        factorized_ = false;
        factor_size_ = 0;
        return;
    }

    refactor_suffix_from(reusable_prefix, reusable_prefix == 0);
    ++stats_.dynamic_reorder_refactorizations;
    stats_.last_reused_prefix = reusable_prefix;
    if (reusable_prefix > 0) {
        ++stats_.reordered_prefix_reuses;
        stats_.prefix_columns_reused += static_cast<std::uint64_t>(reusable_prefix);
    }
}

int SparseExpandingCholesky::suffix_start_from_vars(const std::vector<int>& vars) const {
    int s = n_;
    for (int v : vars) {
        if (v < 0 || v >= n_) {
            throw std::out_of_range("SparseExpandingCholesky touched variable out of range");
        }
        s = std::min(s, v);
    }
    return s == n_ ? -1 : s;
}

int SparseExpandingCholesky::suffix_start_from_matrix(const Sparse& A) const {
    int s = n_;
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Sparse::InnerIterator it(A, col); it; ++it) {
            if (it.value() == 0.0) continue;
            if (it.row() < 0 || it.row() >= n_ || it.col() < 0 || it.col() >= n_) {
                throw std::out_of_range("SparseExpandingCholesky delta entry out of range");
            }
            s = std::min(s, static_cast<int>(std::min(it.row(), it.col())));
        }
    }
    return s == n_ ? -1 : s;
}

std::vector<int> SparseExpandingCholesky::touched_positions_from_vars_and_matrix(
    const std::vector<int>& vars,
    const Sparse& A) const {
    std::vector<int> touched;
    touched.reserve(vars.size() + static_cast<size_t>(2 * A.nonZeros()));

    for (int v : vars) {
        if (v < 0 || v >= n_) {
            throw std::out_of_range("SparseExpandingCholesky touched variable out of range");
        }
        touched.push_back(v);
    }

    for (int col = 0; col < A.outerSize(); ++col) {
        for (Sparse::InnerIterator it(A, col); it; ++it) {
            if (it.value() == 0.0) continue;
            if (it.row() < 0 || it.row() >= n_ || it.col() < 0 || it.col() >= n_) {
                throw std::out_of_range("SparseExpandingCholesky delta entry out of range");
            }
            touched.push_back(it.row());
            touched.push_back(it.col());
        }
    }

    return unique_sorted_local(std::move(touched));
}

std::vector<int> SparseExpandingCholesky::elimination_tree_from_hessian() const {
    std::vector<int> parent(static_cast<size_t>(n_), -1);
    std::vector<int> ancestor(static_cast<size_t>(n_), -1);

    for (int k = 0; k < n_; ++k) {
        for (Sparse::InnerIterator it(H_, k); it; ++it) {
            int i = it.row();
            if (i >= k) continue;
            while (i != -1 && i < k) {
                const int inext = ancestor[static_cast<size_t>(i)];
                ancestor[static_cast<size_t>(i)] = k;
                if (inext == -1) {
                    parent[static_cast<size_t>(i)] = k;
                }
                i = inext;
            }
        }
    }

    return parent;
}

std::vector<int> SparseExpandingCholesky::ancestor_closure_from_seeds(
    const std::vector<int>& seeds,
    const std::vector<int>& parent) const {
    if (static_cast<int>(parent.size()) != n_) {
        throw std::runtime_error("SparseExpandingCholesky etree size mismatch");
    }

    std::vector<unsigned char> marked(static_cast<size_t>(n_), 0);
    for (int seed : seeds) {
        if (seed < 0 || seed >= n_) {
            throw std::out_of_range("SparseExpandingCholesky seed out of range");
        }
        int v = seed;
        while (v >= 0 && v < n_ && !marked[static_cast<size_t>(v)]) {
            marked[static_cast<size_t>(v)] = 1;
            v = parent[static_cast<size_t>(v)];
        }
    }

    std::vector<int> closure;
    for (int i = 0; i < n_; ++i) {
        if (marked[static_cast<size_t>(i)]) {
            closure.push_back(i);
        }
    }
    return closure;
}


void SparseExpandingCholesky::invalidate_factor_dependency_cache() {
    if (factor_dependency_cache_valid_) {
        ++stats_.dependency_cache_invalidations;
    }
    factor_dependency_cache_valid_ = false;
}

void SparseExpandingCholesky::rebuild_factor_dependency_cache() {
    factor_dependency_children_.assign(static_cast<size_t>(n_), {});
    if (!factorized_ || factor_size_ != n_ || L_.rows() != n_ || L_.cols() != n_) {
        factor_dependency_cache_valid_ = false;
        return;
    }

    for (int col = 0; col < L_.outerSize(); ++col) {
        if (col < 0 || col >= n_) continue;
        auto& children = factor_dependency_children_[static_cast<size_t>(col)];
        children.clear();
        for (Sparse::InnerIterator it(L_, col); it; ++it) {
            const int row = it.row();
            if (row > col && row < n_ && it.value() != 0.0) {
                children.push_back(row);
            }
        }
        children = unique_sorted_local(std::move(children));
    }

    factor_dependency_cache_valid_ = true;
    ++stats_.dependency_cache_rebuilds;
}

void SparseExpandingCholesky::refresh_factor_dependency_cache_columns(
    const std::vector<int>& columns) {
    if (!options_.use_factor_dependency_cache) {
        invalidate_factor_dependency_cache();
        return;
    }
    if (!factorized_ || factor_size_ != n_ || L_.rows() != n_ || L_.cols() != n_) {
        invalidate_factor_dependency_cache();
        return;
    }
    if (!factor_dependency_cache_valid_ || static_cast<int>(factor_dependency_children_.size()) != n_) {
        rebuild_factor_dependency_cache();
        return;
    }

    for (int col : unique_sorted_local(columns)) {
        if (col < 0 || col >= n_) {
            throw std::out_of_range("SparseExpandingCholesky dependency-cache column out of range");
        }
        std::vector<int> children;
        if (col < L_.outerSize()) {
            for (Sparse::InnerIterator it(L_, col); it; ++it) {
                const int row = it.row();
                if (row > col && row < n_ && it.value() != 0.0) {
                    children.push_back(row);
                }
            }
        }
        factor_dependency_children_[static_cast<size_t>(col)] = unique_sorted_local(std::move(children));
        ++stats_.dependency_cache_column_refreshes;
    }
}

void SparseExpandingCholesky::refresh_factor_dependency_cache_from(int first_column) {
    if (!options_.use_factor_dependency_cache) {
        invalidate_factor_dependency_cache();
        return;
    }
    if (!factorized_ || factor_size_ != n_) {
        invalidate_factor_dependency_cache();
        return;
    }
    first_column = std::clamp(first_column, 0, std::max(0, n_ - 1));
    std::vector<int> columns;
    columns.reserve(static_cast<size_t>(n_ - first_column));
    for (int col = first_column; col < n_; ++col) {
        columns.push_back(col);
    }
    refresh_factor_dependency_cache_columns(columns);
}

std::vector<int> SparseExpandingCholesky::factor_dependency_closure_from_seeds(
    const std::vector<int>& seeds) {
    if (static_cast<int>(factor_dependency_children_.size()) != n_ || !factor_dependency_cache_valid_) {
        rebuild_factor_dependency_cache();
    }
    if (!factor_dependency_cache_valid_) {
        return {};
    }

    ++stats_.dependency_cache_hits;
    ++stats_.dependency_closure_computations;

    std::vector<unsigned char> marked(static_cast<size_t>(n_), 0);
    std::vector<int> stack;
    for (int seed : seeds) {
        if (seed < 0 || seed >= n_) {
            throw std::out_of_range("SparseExpandingCholesky dependency seed out of range");
        }
        if (!marked[static_cast<size_t>(seed)]) {
            marked[static_cast<size_t>(seed)] = 1;
            stack.push_back(seed);
        }
    }

    while (!stack.empty()) {
        const int v = stack.back();
        stack.pop_back();
        for (int child : factor_dependency_children_[static_cast<size_t>(v)]) {
            if (child < 0 || child >= n_) {
                throw std::runtime_error("SparseExpandingCholesky dependency-cache child out of range");
            }
            if (!marked[static_cast<size_t>(child)]) {
                marked[static_cast<size_t>(child)] = 1;
                stack.push_back(child);
            }
        }
    }

    std::vector<int> closure;
    for (int i = 0; i < n_; ++i) {
        if (marked[static_cast<size_t>(i)]) closure.push_back(i);
    }
    stats_.dependency_closure_columns += static_cast<std::uint64_t>(closure.size());
    stats_.last_dependency_closure_size = static_cast<int>(closure.size());
    return closure;
}

std::vector<int> SparseExpandingCholesky::affected_closure_from_seeds(
    const std::vector<int>& seeds,
    bool allow_dependency_cache) {
    if (seeds.empty()) return {};
    if (allow_dependency_cache && options_.use_factor_dependency_cache) {
        auto closure = factor_dependency_closure_from_seeds(seeds);
        if (!closure.empty()) {
            return closure;
        }
    }

    ++stats_.etree_closure_cache_bypasses;
    const auto parent = elimination_tree_from_hessian();
    auto closure = ancestor_closure_from_seeds(seeds, parent);
    ++stats_.etree_closure_computations;
    stats_.last_dependency_closure_size = 0;
    return closure;
}

bool SparseExpandingCholesky::should_use_affected_closure(const std::vector<int>& closure,
                                                          int suffix_start) const {
    if (closure.empty() || suffix_start < 0 || suffix_start >= n_) {
        return false;
    }
    const int suffix_size = n_ - suffix_start;
    if (static_cast<int>(closure.size()) >= suffix_size) {
        return false;
    }

    // For very small factors the suffix path is simpler and essentially as cheap.
    // For larger factors, use closure mode only when it avoids a material amount
    // of work.  This keeps the new path conservative rather than merely different.
    return static_cast<int>(closure.size()) <= std::max(1, suffix_size - 2);
}


bool SparseExpandingCholesky::delta_changes_normal_pattern(const Sparse& H_delta,
                                                           bool add) const {
    if (H_delta.rows() != n_ || H_delta.cols() != n_) {
        throw std::runtime_error("SparseExpandingCholesky pattern check dimension mismatch");
    }

    constexpr double zero_tol = 0.0;
    for (int col = 0; col < H_delta.outerSize(); ++col) {
        for (Sparse::InnerIterator it(H_delta, col); it; ++it) {
            const double dv = it.value();
            if (dv == 0.0) continue;
            const double old_v = H_.coeff(it.row(), it.col());
            const double new_v = old_v + (add ? dv : -dv);
            const bool old_nz = std::abs(old_v) > zero_tol;
            const bool new_nz = std::abs(new_v) > zero_tol;
            if (old_nz != new_nz) {
                return true;
            }
        }
    }
    return false;
}


std::vector<std::vector<int>> SparseExpandingCholesky::lower_symbolic_pattern_from_hessian(
    const Sparse& H) const {
    if (H.rows() != H.cols()) {
        throw std::runtime_error("SparseExpandingCholesky symbolic pattern requires a square matrix");
    }
    const int n = H.rows();
    std::vector<std::vector<int>> h_lower(static_cast<size_t>(n));
    std::vector<std::vector<int>> Lpat(static_cast<size_t>(n));
    std::vector<std::vector<int>> row_to_cols(static_cast<size_t>(n));

    for (int col = 0; col < H.outerSize(); ++col) {
        for (Sparse::InnerIterator it(H, col); it; ++it) {
            if (it.value() == 0.0) continue;
            const int r = it.row();
            const int c = it.col();
            if (r < 0 || r >= n || c < 0 || c >= n) {
                throw std::runtime_error("SparseExpandingCholesky symbolic pattern entry out of range");
            }
            const int lower_row = std::max(r, c);
            const int lower_col = std::min(r, c);
            h_lower[static_cast<size_t>(lower_col)].push_back(lower_row);
        }
    }
    for (int k = 0; k < n; ++k) {
        h_lower[static_cast<size_t>(k)].push_back(k);
        h_lower[static_cast<size_t>(k)] = unique_sorted_local(std::move(h_lower[static_cast<size_t>(k)]));
    }

    for (int k = 0; k < n; ++k) {
        std::vector<unsigned char> mark(static_cast<size_t>(n), 0);
        for (int row : h_lower[static_cast<size_t>(k)]) {
            if (row >= k) mark[static_cast<size_t>(row)] = 1;
        }
        mark[static_cast<size_t>(k)] = 1;

        auto previous_cols = unique_sorted_local(std::move(row_to_cols[static_cast<size_t>(k)]));
        for (int p : previous_cols) {
            if (p < 0 || p >= k) continue;
            const auto& pcol = Lpat[static_cast<size_t>(p)];
            if (!std::binary_search(pcol.begin(), pcol.end(), k)) continue;
            auto it = std::lower_bound(pcol.begin(), pcol.end(), k);
            for (; it != pcol.end(); ++it) {
                mark[static_cast<size_t>(*it)] = 1;
            }
        }

        std::vector<int> colpat;
        for (int row = k; row < n; ++row) {
            if (mark[static_cast<size_t>(row)]) colpat.push_back(row);
        }
        for (int row : colpat) {
            if (row > k) row_to_cols[static_cast<size_t>(row)].push_back(k);
        }
        Lpat[static_cast<size_t>(k)] = std::move(colpat);
    }
    return Lpat;
}

std::vector<int> SparseExpandingCholesky::factor_pattern_changed_columns(
    const std::vector<std::vector<int>>& before,
    const std::vector<std::vector<int>>& after) const {
    if (before.size() != after.size()) {
        throw std::runtime_error("SparseExpandingCholesky symbolic pattern comparison size mismatch");
    }
    std::vector<int> changed;
    for (int k = 0; k < static_cast<int>(before.size()); ++k) {
        if (before[static_cast<size_t>(k)] != after[static_cast<size_t>(k)]) {
            changed.push_back(k);
        }
    }
    return changed;
}

bool SparseExpandingCholesky::factorization_residual_within_tolerance(double* scaled_residual) const {
    if (!factorized_ || factor_size_ != n_ || L_.rows() != n_ || L_.cols() != n_) {
        if (scaled_residual != nullptr) *scaled_residual = std::numeric_limits<double>::infinity();
        return false;
    }

    Sparse LLt = (L_ * L_.transpose()).pruned();
    LLt.makeCompressed();
    Sparse R = H_ - LLt;
    R.makeCompressed();

    const double abs_res = max_abs_sparse_local(R);
    const double scale = 1.0 + max_abs_sparse_local(H_);
    const double rel_res = abs_res / scale;
    if (scaled_residual != nullptr) {
        *scaled_residual = rel_res;
    }
    return std::isfinite(rel_res) && rel_res <= options_.certification_tolerance;
}



double SparseExpandingCholesky::scaled_factorization_residual() const {
    double residual = std::numeric_limits<double>::infinity();
    (void)factorization_residual_within_tolerance(&residual);
    return residual;
}

bool SparseExpandingCholesky::passes_factorization_residual_check(double* scaled_residual) const {
    return factorization_residual_within_tolerance(scaled_residual);
}

std::vector<int> SparseExpandingCholesky::certification_columns_for_factor_columns(
    const std::vector<int>& factor_columns) const {
    std::vector<unsigned char> marked(static_cast<size_t>(n_), 0);
    for (int col : factor_columns) {
        if (col < 0 || col >= n_) {
            throw std::out_of_range("SparseExpandingCholesky certification factor column out of range");
        }
        marked[static_cast<size_t>(col)] = 1;
        if (col < L_.outerSize()) {
            for (Sparse::InnerIterator it(L_, col); it; ++it) {
                const int row = it.row();
                if (row >= 0 && row < n_) {
                    marked[static_cast<size_t>(row)] = 1;
                }
            }
        }
    }

    std::vector<int> columns;
    for (int i = 0; i < n_; ++i) {
        if (marked[static_cast<size_t>(i)]) {
            columns.push_back(i);
        }
    }
    return columns;
}

bool SparseExpandingCholesky::factorization_residual_columns_within_tolerance(
    const std::vector<int>& columns,
    double* scaled_residual) const {
    if (!factorized_ || factor_size_ != n_ || L_.rows() != n_ || L_.cols() != n_) {
        if (scaled_residual != nullptr) *scaled_residual = std::numeric_limits<double>::infinity();
        return false;
    }
    if (columns.empty()) {
        if (scaled_residual != nullptr) *scaled_residual = 0.0;
        return true;
    }

    std::vector<int> cert_cols = unique_sorted_local(columns);
    for (int col : cert_cols) {
        if (col < 0 || col >= n_) {
            throw std::out_of_range("SparseExpandingCholesky certification column out of range");
        }
    }

    double max_res = 0.0;
    double max_scale = 0.0;
    for (int j : cert_cols) {
        std::map<int, AveragedValue> h_vals;
        for (int col = 0; col < H_.outerSize(); ++col) {
            for (Sparse::InnerIterator it(H_, col); it; ++it) {
                if (it.value() == 0.0) continue;
                if (col == j) {
                    h_vals[it.row()].add(it.value());
                } else if (it.row() == j) {
                    h_vals[col].add(it.value());
                }
            }
        }

        std::map<int, double> residual_col;
        for (const auto& kv : h_vals) {
            const double v = kv.second.value();
            residual_col[kv.first] = v;
            max_scale = std::max(max_scale, std::abs(v));
        }

        // Column j of L*L^T is sum_k L(:,k) L(j,k), where k runs over
        // the nonzeros in row j of the lower factor.
        for (int k = 0; k < L_.outerSize(); ++k) {
            if (k > j) break;
            double Ljk = 0.0;
            for (Sparse::InnerIterator it(L_, k); it; ++it) {
                if (it.row() == j) {
                    Ljk = it.value();
                    break;
                }
            }
            if (Ljk == 0.0) continue;
            for (Sparse::InnerIterator it(L_, k); it; ++it) {
                const int row = it.row();
                const double v = it.value() * Ljk;
                residual_col[row] -= v;
                max_scale = std::max(max_scale, std::abs(v));
            }
        }

        for (const auto& kv : residual_col) {
            if (kv.first < 0 || kv.first >= n_) {
                throw std::runtime_error("SparseExpandingCholesky residual column index out of range");
            }
            max_res = std::max(max_res, std::abs(kv.second));
        }
    }

    const double rel_res = max_res / (1.0 + max_scale);
    if (scaled_residual != nullptr) {
        *scaled_residual = rel_res;
    }
    return std::isfinite(rel_res) && rel_res <= options_.certification_tolerance;
}

bool SparseExpandingCholesky::certify_factorization_for_columns(
    const std::vector<int>& factor_columns,
    double* scaled_residual) {
    if (!options_.certify_affected_closure_updates) {
        if (scaled_residual != nullptr) *scaled_residual = 0.0;
        return true;
    }

    double residual = std::numeric_limits<double>::infinity();
    bool ok = false;

    if (options_.use_column_local_certification) {
        const auto cert_cols = certification_columns_for_factor_columns(factor_columns);
        ++stats_.column_local_certifications;
        ++stats_.factorization_residual_checks;
        stats_.column_local_certification_columns += static_cast<std::uint64_t>(cert_cols.size());
        ok = factorization_residual_columns_within_tolerance(cert_cols, &residual);
        stats_.last_column_local_certification_residual = residual;
        if (std::isfinite(residual)) {
            stats_.max_column_local_certification_residual =
                std::max(stats_.max_column_local_certification_residual, residual);
        }
        if (!ok) {
            ++stats_.column_local_certification_failures;
        }
    }

    if (!ok && (!options_.use_column_local_certification || options_.full_certification_fallback)) {
        ++stats_.full_certification_fallbacks;
        ++stats_.factorization_residual_checks;
        ok = factorization_residual_within_tolerance(&residual);
    }

    if (scaled_residual != nullptr) {
        *scaled_residual = residual;
    }
    return ok;
}

std::vector<std::map<int, double>> SparseExpandingCholesky::lower_columns_from_factor() const {
    std::vector<std::map<int, double>> cols(static_cast<size_t>(n_));
    for (int col = 0; col < L_.outerSize(); ++col) {
        if (col < 0 || col >= n_) continue;
        for (Sparse::InnerIterator it(L_, col); it; ++it) {
            if (it.row() >= col && it.row() < n_) {
                cols[static_cast<size_t>(col)][it.row()] = it.value();
            }
        }
    }
    return cols;
}

std::map<int, double> SparseExpandingCholesky::hessian_lower_column(int k) const {
    if (k < 0 || k >= n_) {
        throw std::out_of_range("SparseExpandingCholesky Hessian column out of range");
    }

    std::map<int, AveragedValue> vals;
    for (int col = 0; col < H_.outerSize(); ++col) {
        for (Sparse::InnerIterator it(H_, col); it; ++it) {
            const int r = it.row();
            const int c = it.col();
            if (r == k && c >= k) {
                vals[c].add(it.value());
            } else if (c == k && r >= k) {
                vals[r].add(it.value());
            }
        }
    }

    std::map<int, double> out;
    for (const auto& kv : vals) {
        const double v = kv.second.value();
        if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance || kv.first == k) {
            out[kv.first] = v;
        }
    }
    out.emplace(k, 0.0);
    return out;
}

void SparseExpandingCholesky::rebuild_lower_from_columns(
    const std::vector<std::map<int, double>>& lower_cols) {
    if (static_cast<int>(lower_cols.size()) != n_) {
        throw std::runtime_error("SparseExpandingCholesky lower-column size mismatch");
    }

    std::vector<Triplet> trips;
    size_t nnz_est = 0;
    for (const auto& col : lower_cols) nnz_est += col.size();
    trips.reserve(nnz_est);

    for (int col = 0; col < n_; ++col) {
        for (const auto& entry : lower_cols[static_cast<size_t>(col)]) {
            const int row = entry.first;
            const double v = entry.second;
            if (row < col || row >= n_) {
                throw std::runtime_error("SparseExpandingCholesky invalid lower-column entry");
            }
            if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance || row == col) {
                trips.emplace_back(row, col, v);
            }
        }
    }

    L_.resize(n_, n_);
    L_.setFromTriplets(trips.begin(), trips.end());
    prune_and_compress(L_, options_.drop_tolerance);
}

bool SparseExpandingCholesky::try_refactor_affected_closure(const std::vector<int>& seeds,
                                                            bool allow_dependency_cache) {
    if (!factorized_ || factor_size_ != n_ || seeds.empty()) {
        return false;
    }

    const int suffix_start = *std::min_element(seeds.begin(), seeds.end());
    auto closure = affected_closure_from_seeds(seeds, allow_dependency_cache);
    stats_.last_affected_closure_size = static_cast<int>(closure.size());
    stats_.last_affected_closure_min = closure.empty() ? -1 : closure.front();

    if (!should_use_affected_closure(closure, suffix_start)) {
        ++stats_.affected_closure_fallbacks;
        return false;
    }

    const Sparse L_before = L_;
    auto lower_cols = lower_columns_from_factor();
    std::uint64_t jitter_used = 0;

    try {
        for (int k : closure) {
            std::map<int, double> w = hessian_lower_column(k);

            for (int p = 0; p < k; ++p) {
                const auto& pcol = lower_cols[static_cast<size_t>(p)];
                const auto lkp_it = pcol.find(k);
                if (lkp_it == pcol.end()) continue;
                const double Lkp = lkp_it->second;
                if (Lkp == 0.0) continue;
                for (auto it = pcol.lower_bound(k); it != pcol.end(); ++it) {
                    w[it->first] -= it->second * Lkp;
                }
            }

            const double raw_diag = w[k];
            double diag = raw_diag;
            double jitter = options_.jitter;
            int tries = 0;
            while (!(std::isfinite(diag) && diag > 0.0) &&
                   tries < options_.max_jitter_tries &&
                   jitter > 0.0) {
                diag = raw_diag + jitter;
                jitter *= 10.0;
                ++tries;
            }
            if (!(std::isfinite(diag) && diag > 0.0)) {
                ++stats_.affected_closure_fallbacks;
                return false;
            }
            if (tries > 0) {
                ++jitter_used;
            }

            const double Lkk = std::sqrt(diag);
            std::map<int, double> new_col;
            new_col[k] = Lkk;
            for (auto it = w.upper_bound(k); it != w.end(); ++it) {
                const double v = it->second / Lkk;
                if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
                    new_col[it->first] = v;
                }
            }
            lower_cols[static_cast<size_t>(k)] = std::move(new_col);
        }
    } catch (...) {
        ++stats_.affected_closure_fallbacks;
        return false;
    }

    const auto certification_columns_before = certification_columns_for_factor_columns(closure);
    rebuild_lower_from_columns(lower_cols);
    factorized_ = true;
    factor_size_ = n_;

    if (options_.certify_affected_closure_updates) {
        double residual = 0.0;
        ++stats_.affected_closure_certifications;
        std::vector<int> certification_factor_columns = closure;
        certification_factor_columns.insert(certification_factor_columns.end(),
                                            certification_columns_before.begin(),
                                            certification_columns_before.end());
        certification_factor_columns = unique_sorted_local(std::move(certification_factor_columns));
        const bool ok = certify_factorization_for_columns(certification_factor_columns, &residual);
        stats_.last_certification_residual = residual;
        if (std::isfinite(residual)) {
            stats_.max_certification_residual = std::max(stats_.max_certification_residual, residual);
        }
        if (!ok) {
            L_ = L_before;
            factorized_ = true;
            factor_size_ = n_;
            ++stats_.affected_closure_certification_failures;
            ++stats_.affected_closure_fallbacks;
            return false;
        }
    }

    stats_.jitter_regularizations += jitter_used;
    refresh_factor_dependency_cache_columns(closure);
    ++stats_.affected_closure_refactorizations;
    stats_.affected_columns_refactored += static_cast<std::uint64_t>(closure.size());
    stats_.last_suffix_start = -1;
    stats_.last_suffix_size = 0;
    return true;
}


bool SparseExpandingCholesky::try_refactor_structural_closure(
    const std::vector<int>& seeds,
    const std::vector<int>& changed_factor_columns) {
    if (!factorized_ || factor_size_ != n_) {
        return false;
    }

    std::vector<int> structural_seeds = seeds;
    structural_seeds.insert(structural_seeds.end(),
                            changed_factor_columns.begin(),
                            changed_factor_columns.end());
    structural_seeds = unique_sorted_local(std::move(structural_seeds));
    if (structural_seeds.empty()) {
        return false;
    }

    ++stats_.structural_closure_attempts;
    const bool ok = try_refactor_affected_closure(
        structural_seeds, changed_factor_columns.empty());
    if (ok) {
        ++stats_.structural_closure_refactorizations;
    } else {
        ++stats_.structural_closure_fallbacks;
    }
    return ok;
}

void SparseExpandingCholesky::apply_contribution(const Sparse& C,
                                                 const std::vector<int>& touched_vars,
                                                 bool add) {
    if (C.rows() > n_) {
        grow_to(C.rows());
    }
    if (C.rows() != n_) {
        throw std::runtime_error("SparseExpandingCholesky contribution row dimension mismatch");
    }
    Sparse H_delta = C * C.transpose();
    prune_and_compress(H_delta, options_.drop_tolerance);
    apply_hessian_delta(H_delta, touched_vars, add);
}

void SparseExpandingCholesky::apply_hessian_delta(const Sparse& H_delta,
                                                  const std::vector<int>& touched_vars,
                                                  bool add) {
    if (H_delta.rows() > n_ || H_delta.cols() > n_) {
        if (H_delta.rows() != H_delta.cols()) {
            throw std::runtime_error("SparseExpandingCholesky Hessian delta must be square");
        }
        grow_to(H_delta.rows());
    }
    if (H_delta.rows() != n_ || H_delta.cols() != n_) {
        throw std::runtime_error("SparseExpandingCholesky Hessian delta dimension mismatch");
    }

    const bool was_growth = factorized_ && factor_size_ < n_;
    const bool can_classify_symbolic_change =
        factorized_ && factor_size_ == n_ && !was_growth;
    const bool structural_pattern_change = delta_changes_normal_pattern(H_delta, add);

    Sparse H_candidate = H_;
    if (add) {
        H_candidate += H_delta;
    } else {
        H_candidate -= H_delta;
    }
    prune_and_compress(H_candidate, options_.drop_tolerance);

    std::vector<int> changed_factor_columns;
    if (structural_pattern_change) {
        ++stats_.structural_pattern_changes;
        if (can_classify_symbolic_change) {
            ++stats_.symbolic_pattern_classifications;
            const auto before_pattern = lower_symbolic_pattern_from_hessian(H_);
            const auto after_pattern = lower_symbolic_pattern_from_hessian(H_candidate);
            changed_factor_columns = factor_pattern_changed_columns(before_pattern, after_pattern);
            stats_.last_structural_factor_pattern_changed_columns =
                static_cast<int>(changed_factor_columns.size());
            if (changed_factor_columns.empty()) {
                ++stats_.pattern_stable_structural_updates;
            } else {
                ++stats_.structural_factor_pattern_changes;
                stats_.structural_factor_columns_changed +=
                    static_cast<std::uint64_t>(changed_factor_columns.size());
            }
        }
    } else {
        ++stats_.numeric_only_updates;
        stats_.last_structural_factor_pattern_changed_columns = 0;
    }

    H_ = std::move(H_candidate);

    const std::vector<int> seeds = touched_positions_from_vars_and_matrix(touched_vars, H_delta);
    int s = seeds.empty() ? -1 : seeds.front();
    if (s < 0) {
        s = suffix_start_from_matrix(H_delta);
    }
    if (s < 0) {
        return;
    }

    // M73: same-pattern, same-size updates may use the M72 etree-closure
    // update, but accepted closure updates are certified by a residual check.
    // M74: structural normal-pattern changes are no longer forced directly to
    // the suffix path.  We first classify the induced symbolic factor-pattern
    // change.  If the changed symbolic columns are local, they are included in
    // a certified etree-closure update; otherwise the M71 sparse suffix path is
    // used as the conservative recovery path.
    if (!was_growth && factorized_ && factor_size_ == n_) {
        if (!structural_pattern_change && try_refactor_affected_closure(seeds, true)) {
            ++stats_.same_size_suffix_updates;
            return;
        }
        if (structural_pattern_change &&
            can_classify_symbolic_change &&
            try_refactor_structural_closure(seeds, changed_factor_columns)) {
            ++stats_.same_size_suffix_updates;
            return;
        }
        if (structural_pattern_change) {
            ++stats_.affected_closure_fallbacks;
        }
    }

    if (!factorized_ || factor_size_ <= 0 || s > factor_size_) {
        s = 0;
    }
    refactor_suffix_from(s, false);
    if (was_growth) {
        ++stats_.expansion_suffix_updates;
    } else {
        ++stats_.same_size_suffix_updates;
    }
}

void SparseExpandingCholesky::apply_diagonal_update(const std::vector<int>& vars,
                                                    double weight,
                                                    bool add) {
    if (vars.empty() || weight == 0.0) {
        return;
    }
    const auto clean_vars = unique_sorted_local(vars);
    std::vector<Triplet> trips;
    trips.reserve(clean_vars.size());
    for (int v : clean_vars) {
        if (v < 0 || v >= n_) {
            throw std::out_of_range("SparseExpandingCholesky diagonal variable out of range");
        }
        trips.emplace_back(v, v, weight);
    }
    Sparse D(n_, n_);
    D.setFromTriplets(trips.begin(), trips.end());
    D.makeCompressed();
    ++stats_.diagonal_updates;
    apply_hessian_delta(D, clean_vars, add);
}

Sparse SparseExpandingCholesky::sparse_suffix_block(int suffix_start) const {
    const int m = n_ - suffix_start;
    std::map<std::pair<int, int>, AveragedValue> vals;
    for (int col = 0; col < H_.outerSize(); ++col) {
        for (Sparse::InnerIterator it(H_, col); it; ++it) {
            const int r = it.row();
            const int c = it.col();
            if (r < suffix_start || c < suffix_start) continue;
            const int rr = r - suffix_start;
            const int cc = c - suffix_start;
            const auto key = std::make_pair(std::max(rr, cc), std::min(rr, cc));
            vals[key].add(it.value());
        }
    }

    std::vector<Triplet> trips;
    trips.reserve(vals.size() * 2);
    for (const auto& kv : vals) {
        const int r = kv.first.first;
        const int c = kv.first.second;
        const double v = kv.second.value();
        if (options_.drop_tolerance > 0.0 && std::abs(v) <= options_.drop_tolerance) continue;
        trips.emplace_back(r, c, v);
        if (r != c) trips.emplace_back(c, r, v);
    }
    Sparse out(m, m);
    out.setFromTriplets(trips.begin(), trips.end());
    prune_and_compress(out, options_.drop_tolerance);
    return out;
}

Sparse SparseExpandingCholesky::sparse_a21_block(int suffix_start) const {
    const int m = n_ - suffix_start;
    const int p = suffix_start;
    std::map<std::pair<int, int>, AveragedValue> vals;
    if (p <= 0) {
        return Sparse(m, 0);
    }

    for (int col = 0; col < H_.outerSize(); ++col) {
        for (Sparse::InnerIterator it(H_, col); it; ++it) {
            const int r = it.row();
            const int c = it.col();
            if (r >= suffix_start && c < suffix_start) {
                vals[std::make_pair(r - suffix_start, c)].add(it.value());
            } else if (c >= suffix_start && r < suffix_start) {
                vals[std::make_pair(c - suffix_start, r)].add(it.value());
            }
        }
    }

    std::vector<Triplet> trips;
    trips.reserve(vals.size());
    for (const auto& kv : vals) {
        const double v = kv.second.value();
        if (options_.drop_tolerance > 0.0 && std::abs(v) <= options_.drop_tolerance) continue;
        trips.emplace_back(kv.first.first, kv.first.second, v);
    }
    Sparse out(m, p);
    out.setFromTriplets(trips.begin(), trips.end());
    prune_and_compress(out, options_.drop_tolerance);
    return out;
}

Sparse SparseExpandingCholesky::sparse_l11_block(int suffix_start) const {
    Sparse out(suffix_start, suffix_start);
    if (suffix_start <= 0) {
        return out;
    }
    std::vector<Triplet> trips;
    trips.reserve(static_cast<size_t>(L_.nonZeros()));
    for (int col = 0; col < L_.outerSize(); ++col) {
        if (col >= suffix_start) continue;
        for (Sparse::InnerIterator it(L_, col); it; ++it) {
            if (it.row() < suffix_start) {
                trips.emplace_back(it.row(), col, it.value());
            }
        }
    }
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

Sparse SparseExpandingCholesky::solve_lower_prefix_sparse(int prefix_size,
                                                          const Sparse& rhs) const {
    if (rhs.rows() != prefix_size) {
        throw std::runtime_error("SparseExpandingCholesky prefix solve RHS row mismatch");
    }
    if (prefix_size <= 0) {
        return Sparse(0, rhs.cols());
    }

    std::vector<std::vector<std::pair<int, double>>> lower_rows(static_cast<size_t>(prefix_size));
    std::vector<double> diag(static_cast<size_t>(prefix_size), std::numeric_limits<double>::quiet_NaN());

    for (int col = 0; col < L_.outerSize(); ++col) {
        if (col >= prefix_size) continue;
        for (Sparse::InnerIterator it(L_, col); it; ++it) {
            const int row = it.row();
            if (row >= prefix_size) continue;
            if (row == col) {
                diag[static_cast<size_t>(row)] = it.value();
            } else if (row > col) {
                lower_rows[static_cast<size_t>(row)].emplace_back(col, it.value());
            }
        }
    }
    for (int i = 0; i < prefix_size; ++i) {
        if (!(std::isfinite(diag[static_cast<size_t>(i)]) && std::abs(diag[static_cast<size_t>(i)]) > 0.0)) {
            throw std::runtime_error("SparseExpandingCholesky prefix solve missing diagonal");
        }
        std::sort(lower_rows[static_cast<size_t>(i)].begin(), lower_rows[static_cast<size_t>(i)].end());
    }

    std::vector<Triplet> trips;
    std::vector<double> b(static_cast<size_t>(prefix_size), 0.0);
    std::vector<double> x(static_cast<size_t>(prefix_size), 0.0);

    for (int rhs_col = 0; rhs_col < rhs.outerSize(); ++rhs_col) {
        std::fill(b.begin(), b.end(), 0.0);
        std::fill(x.begin(), x.end(), 0.0);
        for (Sparse::InnerIterator it(rhs, rhs_col); it; ++it) {
            if (it.row() < 0 || it.row() >= prefix_size) {
                throw std::runtime_error("SparseExpandingCholesky prefix solve RHS index out of range");
            }
            b[static_cast<size_t>(it.row())] += it.value();
        }
        for (int row = 0; row < prefix_size; ++row) {
            double acc = b[static_cast<size_t>(row)];
            for (const auto& entry : lower_rows[static_cast<size_t>(row)]) {
                acc -= entry.second * x[static_cast<size_t>(entry.first)];
            }
            x[static_cast<size_t>(row)] = acc / diag[static_cast<size_t>(row)];
        }
        for (int row = 0; row < prefix_size; ++row) {
            const double v = x[static_cast<size_t>(row)];
            if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
                trips.emplace_back(row, rhs_col, v);
            }
        }
    }

    Sparse out(prefix_size, rhs.cols());
    out.setFromTriplets(trips.begin(), trips.end());
    prune_and_compress(out, options_.drop_tolerance);
    return out;
}

Sparse SparseExpandingCholesky::factorize_sparse_left_looking(
    const Sparse& A,
    std::uint64_t* jitter_regularizations) const {
    if (A.rows() != A.cols()) {
        throw std::runtime_error("SparseExpandingCholesky sparse factorization requires a square matrix");
    }
    const int n = A.rows();
    if (n <= 0) {
        return Sparse(0, 0);
    }

    std::map<std::pair<int, int>, AveragedValue> avals;
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Sparse::InnerIterator it(A, col); it; ++it) {
            const int r = it.row();
            const int c = it.col();
            if (r < 0 || r >= n || c < 0 || c >= n) {
                throw std::runtime_error("SparseExpandingCholesky sparse factorization entry out of range");
            }
            const auto key = std::make_pair(std::max(r, c), std::min(r, c));
            avals[key].add(it.value());
        }
    }

    std::vector<std::map<int, double>> lower_cols(static_cast<size_t>(n));
    for (const auto& kv : avals) {
        const double v = kv.second.value();
        if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
            lower_cols[static_cast<size_t>(kv.first.second)][kv.first.first] = v;
        }
    }

    std::vector<std::map<int, double>> Lcols(static_cast<size_t>(n));
    std::vector<std::vector<int>> row_to_cols(static_cast<size_t>(n));

    for (int k = 0; k < n; ++k) {
        std::map<int, double> w = lower_cols[static_cast<size_t>(k)];
        w.emplace(k, 0.0);

        auto& previous_cols = row_to_cols[static_cast<size_t>(k)];
        std::sort(previous_cols.begin(), previous_cols.end());
        previous_cols.erase(std::unique(previous_cols.begin(), previous_cols.end()), previous_cols.end());

        for (int p : previous_cols) {
            if (p < 0 || p >= k) continue;
            const auto lk_it = Lcols[static_cast<size_t>(p)].find(k);
            if (lk_it == Lcols[static_cast<size_t>(p)].end()) continue;
            const double Lkp = lk_it->second;
            for (auto it = Lcols[static_cast<size_t>(p)].lower_bound(k); it != Lcols[static_cast<size_t>(p)].end(); ++it) {
                w[it->first] -= it->second * Lkp;
            }
        }

        const double raw_diag = w[k];
        double diag = raw_diag;
        double jitter = options_.jitter;
        int tries = 0;
        while (!(std::isfinite(diag) && diag > 0.0) && tries < options_.max_jitter_tries && jitter > 0.0) {
            diag = raw_diag + jitter;
            jitter *= 10.0;
            ++tries;
        }
        if (!(std::isfinite(diag) && diag > 0.0)) {
            throw std::runtime_error("SparseExpandingCholesky sparse left-looking Cholesky failed; matrix is not SPD");
        }
        if (tries > 0 && jitter_regularizations != nullptr) {
            ++(*jitter_regularizations);
        }

        const double Lkk = std::sqrt(diag);
        std::map<int, double> colmap;
        colmap[k] = Lkk;
        for (auto it = w.upper_bound(k); it != w.end(); ++it) {
            const int row = it->first;
            const double v = it->second / Lkk;
            if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
                colmap[row] = v;
            }
        }
        for (const auto& entry : colmap) {
            if (entry.first > k) {
                row_to_cols[static_cast<size_t>(entry.first)].push_back(k);
            }
        }
        Lcols[static_cast<size_t>(k)] = std::move(colmap);
    }

    std::vector<Triplet> trips;
    size_t nnz_est = 0;
    for (const auto& col : Lcols) nnz_est += col.size();
    trips.reserve(nnz_est);
    for (int col = 0; col < n; ++col) {
        for (const auto& entry : Lcols[static_cast<size_t>(col)]) {
            trips.emplace_back(entry.first, col, entry.second);
        }
    }
    Sparse L(n, n);
    L.setFromTriplets(trips.begin(), trips.end());
    prune_and_compress(L, options_.drop_tolerance);
    return L;
}

void SparseExpandingCholesky::rebuild_lower_with_sparse_suffix(int suffix_start,
                                                               const Sparse& L21,
                                                               const Sparse& L22) {
    const int m = n_ - suffix_start;
    if (L21.rows() != m || L21.cols() != suffix_start) {
        throw std::runtime_error("SparseExpandingCholesky L21 dimension mismatch");
    }
    if (L22.rows() != m || L22.cols() != m) {
        throw std::runtime_error("SparseExpandingCholesky L22 dimension mismatch");
    }

    std::vector<Triplet> trips;
    trips.reserve(static_cast<size_t>(L_.nonZeros() + L21.nonZeros() + L22.nonZeros()));

    // Preserve only L11. L21 and L22 are recomputed from the current matrix.
    for (int col = 0; col < L_.outerSize(); ++col) {
        if (col >= suffix_start) continue;
        for (Sparse::InnerIterator it(L_, col); it; ++it) {
            if (it.row() >= suffix_start) continue;
            const double v = it.value();
            if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
                trips.emplace_back(it.row(), it.col(), v);
            }
        }
    }

    for (int col = 0; col < L21.outerSize(); ++col) {
        for (Sparse::InnerIterator it(L21, col); it; ++it) {
            const double v = it.value();
            if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
                trips.emplace_back(suffix_start + it.row(), col, v);
            }
        }
    }

    for (int col = 0; col < L22.outerSize(); ++col) {
        for (Sparse::InnerIterator it(L22, col); it; ++it) {
            if (it.row() < col) continue;
            const double v = it.value();
            if (options_.drop_tolerance <= 0.0 || std::abs(v) > options_.drop_tolerance) {
                trips.emplace_back(suffix_start + it.row(), suffix_start + col, v);
            }
        }
    }

    L_.resize(n_, n_);
    L_.setFromTriplets(trips.begin(), trips.end());
    prune_and_compress(L_, options_.drop_tolerance);
}

void SparseExpandingCholesky::refactor_suffix_from(int suffix_start, bool count_as_full) {
    if (n_ <= 0) {
        throw std::runtime_error("SparseExpandingCholesky cannot factorize an empty matrix");
    }
    suffix_start = std::clamp(suffix_start, 0, n_ - 1);
    if (!factorized_ || factor_size_ < suffix_start) {
        suffix_start = 0;
    }

    Sparse S = sparse_suffix_block(suffix_start);
    Sparse L21(n_ - suffix_start, suffix_start);

    if (suffix_start > 0) {
        const Sparse A21 = sparse_a21_block(suffix_start);
        const Sparse rhs = A21.transpose(); // L11 * X = A21^T, then L21 = X^T.
        const Sparse X = solve_lower_prefix_sparse(suffix_start, rhs);
        L21 = X.transpose();
        prune_and_compress(L21, options_.drop_tolerance);

        Sparse schur_term = L21 * L21.transpose();
        prune_and_compress(schur_term, options_.drop_tolerance);
        S = S - schur_term;
        prune_and_compress(S, options_.drop_tolerance);

        ++stats_.l21_recomputes;
        ++stats_.sparse_triangular_solves;
        ++stats_.sparse_schur_complements;
    }

    std::uint64_t jitter_used = 0;
    const Sparse L22 = factorize_sparse_left_looking(S, &jitter_used);
    stats_.jitter_regularizations += jitter_used;
    ++stats_.sparse_left_looking_factorizations;

    rebuild_lower_with_sparse_suffix(suffix_start, L21, L22);
    factorized_ = true;
    factor_size_ = n_;
    refresh_factor_dependency_cache_from(suffix_start);
    if (count_as_full || suffix_start == 0) {
        ++stats_.full_factorizations;
    } else {
        ++stats_.suffix_refactorizations;
    }
    stats_.last_suffix_start = suffix_start;
    stats_.last_suffix_size = n_ - suffix_start;
}

Eigen::VectorXd SparseExpandingCholesky::solve(const Eigen::VectorXd& rhs) const {
    if (!covers_state()) {
        throw std::runtime_error("SparseExpandingCholesky::solve called before factor covers state");
    }
    if (rhs.size() != n_) {
        throw std::runtime_error("SparseExpandingCholesky RHS dimension mismatch");
    }
    Eigen::VectorXd y = L_.triangularView<Eigen::Lower>().solve(rhs);
    Sparse Lt = L_.transpose();
    return Lt.triangularView<Eigen::Upper>().solve(y);
}

double SparseExpandingCholesky::eta_logdet_half() const {
    if (!covers_state()) {
        throw std::runtime_error("SparseExpandingCholesky::eta_logdet_half called before factor covers state");
    }
    double eta = 0.0;
    for (int col = 0; col < L_.outerSize(); ++col) {
        bool found = false;
        double diag = 0.0;
        for (Sparse::InnerIterator it(L_, col); it; ++it) {
            if (it.row() == col) {
                diag = it.value();
                found = true;
                break;
            }
        }
        if (!found || !(std::isfinite(diag) && std::abs(diag) > 0.0)) {
            throw std::runtime_error("SparseExpandingCholesky missing/non-positive diagonal");
        }
        eta += std::log(std::abs(diag));
    }
    return eta;
}

} // namespace islam
