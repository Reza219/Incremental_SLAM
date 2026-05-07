#include "islam/factor_manager.hpp"

#include "islam/linearize.hpp"
#include "islam/symbolic_engine.hpp"
#include "islam/deterministic_batch_ccolamd.hpp"

#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <limits>
#include <utility>
#include <vector>

#ifdef ISLAM_HAS_CHOLMOD
#include <cholmod.h>
#endif

namespace islam {
namespace {

void add_sparse_triplets(const Eigen::SparseMatrix<double>& A,
                         std::vector<Eigen::Triplet<double>>& trips,
                         bool upper_only = false) {
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            if (upper_only && it.row() > it.col()) continue;
            trips.emplace_back(it.row(), it.col(), it.value());
        }
    }
}

Eigen::SparseMatrix<double> resize_sparse_preserve(const Eigen::SparseMatrix<double>& A, int n) {
    if (n < A.rows() || n < A.cols()) {
        throw std::runtime_error("resize_sparse_preserve does not support shrinking");
    }
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(A.nonZeros()));
    add_sparse_triplets(A, trips);
    Eigen::SparseMatrix<double> out(n, n);
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

Eigen::SparseMatrix<double> diagonal_sparse(int n, double w) {
    Eigen::SparseMatrix<double> A(n, n);
    if (n <= 0 || w == 0.0) {
        return A;
    }
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        trips.emplace_back(i, i, w);
    }
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    return A;
}

void append_unique_int(std::vector<int>& values, int value) {
    if (std::find(values.begin(), values.end(), value) == values.end()) {
        values.push_back(value);
    }
}

void append_unique_touched_range(std::vector<int>& touched_vars, int start, int dim) {
    for (int k = 0; k < dim; ++k) {
        append_unique_int(touched_vars, start + k);
    }
}

bool chol_rank1_update_lower(Eigen::MatrixXd& L, Eigen::VectorXd x, bool add) {
    const int n = static_cast<int>(x.size());
    if (L.rows() != n || L.cols() != n) {
        throw std::runtime_error("chol_rank1_update_lower dimension mismatch");
    }

    for (int k = 0; k < n; ++k) {
        const double Lkk = L(k, k);
        const double xk = x[k];
        double r = 0.0;
        double c = 0.0;
        double s = 0.0;

        if (add) {
            r = std::hypot(Lkk, xk);
            if (!(r > 0.0)) return false;
            c = r / Lkk;
            s = xk / Lkk;
            L(k, k) = r;
            for (int i = k + 1; i < n; ++i) {
                const double Lik_old = L(i, k);
                const double x_old = x[i];
                const double Lik_new = (Lik_old + s * x_old) / c;
                x[i] = c * x_old - s * Lik_new;
                L(i, k) = Lik_new;
            }
        } else {
            const double delta = Lkk * Lkk - xk * xk;
            if (!(delta > 1e-14)) return false;
            r = std::sqrt(delta);
            c = r / Lkk;
            s = xk / Lkk;
            L(k, k) = r;
            for (int i = k + 1; i < n; ++i) {
                const double Lik_old = L(i, k);
                const double x_old = x[i];
                const double Lik_new = (Lik_old - s * x_old) / c;
                x[i] = c * x_old - s * Lik_new;
                L(i, k) = Lik_new;
            }
        }
    }
    return true;
}

bool chol_rank_k_update_lower(Eigen::MatrixXd& L, const Eigen::MatrixXd& C, bool add) {
    for (int j = 0; j < C.cols(); ++j) {
        if (!chol_rank1_update_lower(L, C.col(j), add)) {
            return false;
        }
    }
    return true;
}

Eigen::MatrixXd sparse_to_dense(const Eigen::SparseMatrix<double>& A) {
    return Eigen::MatrixXd(A);
}

Eigen::SparseMatrix<double> dense_lower_to_sparse(const Eigen::MatrixXd& L) {
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(L.rows() * 4));
    for (int j = 0; j < L.cols(); ++j) {
        for (int i = j; i < L.rows(); ++i) {
            const double v = L(i, j);
            if (v != 0.0) {
                trips.emplace_back(i, j, v);
            }
        }
    }
    Eigen::SparseMatrix<double> out(L.rows(), L.cols());
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}


std::vector<int> inverse_permutation(const std::vector<int>& perm) {
    std::vector<int> pinv(perm.size(), -1);
    for (int k = 0; k < static_cast<int>(perm.size()); ++k) {
        const int pk = perm[static_cast<size_t>(k)];
        if (pk < 0 || pk >= static_cast<int>(perm.size())) {
            throw std::runtime_error("Invalid permutation entry");
        }
        pinv[static_cast<size_t>(pk)] = k;
    }
    return pinv;
}


std::vector<int> elimination_tree_from_upper(const Eigen::SparseMatrix<double>& Hperm) {
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

Eigen::VectorXd permute_vector_by_order(const Eigen::VectorXd& v, const std::vector<int>& perm) {
    Eigen::VectorXd out(static_cast<int>(perm.size()));
    for (int k = 0; k < static_cast<int>(perm.size()); ++k) {
        out[k] = v[perm[static_cast<size_t>(k)]];
    }
    return out;
}

Eigen::VectorXd unpermute_vector_by_order(const Eigen::VectorXd& v_perm, const std::vector<int>& perm) {
    if (v_perm.size() != static_cast<int>(perm.size())) {
        throw std::runtime_error("unpermute_vector_by_order dimension mismatch");
    }
    Eigen::VectorXd out(static_cast<int>(perm.size()));
    for (int k = 0; k < static_cast<int>(perm.size()); ++k) {
        out[perm[static_cast<size_t>(k)]] = v_perm[k];
    }
    return out;
}

Eigen::SparseMatrix<double> row_permute_by_order(const Eigen::SparseMatrix<double>& A,
                                                 const std::vector<int>& perm) {
    if (A.rows() != static_cast<int>(perm.size())) {
        throw std::runtime_error("row_permute_by_order dimension mismatch");
    }
    const auto pinv = inverse_permutation(perm);
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(A.nonZeros()));
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            trips.emplace_back(pinv[static_cast<size_t>(it.row())], col, it.value());
        }
    }
    Eigen::SparseMatrix<double> out(A.rows(), A.cols());
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

std::vector<int> permute_indices_to_positions(const std::vector<int>& vars,
                                              const std::vector<int>& pinv) {
    std::vector<int> out;
    out.reserve(vars.size());
    for (int v : vars) {
        if (v < 0 || v >= static_cast<int>(pinv.size())) {
            throw std::out_of_range("permute_indices_to_positions index out of range");
        }
        out.push_back(pinv[static_cast<size_t>(v)]);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

int common_prefix_size(const std::vector<int>& a, const std::vector<int>& b) {
    const int n = std::min(static_cast<int>(a.size()), static_cast<int>(b.size()));
    int k = 0;
    while (k < n && a[static_cast<size_t>(k)] == b[static_cast<size_t>(k)]) {
        ++k;
    }
    return k;
}

Eigen::SparseMatrix<double> symmetric_permute_by_order(const Eigen::SparseMatrix<double>& A,
                                                       const std::vector<int>& perm) {
    if (A.rows() != A.cols() || A.rows() != static_cast<int>(perm.size())) {
        throw std::runtime_error("symmetric_permute_by_order dimension mismatch");
    }
    const auto pinv = inverse_permutation(perm);

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(A.nonZeros()));
    for (int col = 0; col < A.outerSize(); ++col) {
        const int new_col = pinv[static_cast<size_t>(col)];
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int new_row = pinv[static_cast<size_t>(it.row())];
            trips.emplace_back(new_row, new_col, it.value());
        }
    }

    Eigen::SparseMatrix<double> out(A.rows(), A.cols());
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}


std::vector<int> unique_sorted(std::vector<int> vals) {
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    return vals;
}

int default_anchor_dim_for_graph(const Graph& g) {
    if (g.id_lookup.empty()) {
        return std::min(3, g.state_size());
    }
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

Eigen::SparseMatrix<double> principal_submatrix(const Eigen::SparseMatrix<double>& A,
                                                const std::vector<int>& idx) {
    std::vector<int> inv(static_cast<size_t>(A.rows()), -1);
    for (int k = 0; k < static_cast<int>(idx.size()); ++k) {
        const int v = idx[static_cast<size_t>(k)];
        if (v < 0 || v >= A.rows()) {
            throw std::runtime_error("principal_submatrix index out of range");
        }
        inv[static_cast<size_t>(v)] = k;
    }

    std::vector<Eigen::Triplet<double>> trips;
    for (int col = 0; col < A.outerSize(); ++col) {
        const int new_col = inv[static_cast<size_t>(col)];
        if (new_col < 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int new_row = inv[static_cast<size_t>(it.row())];
            if (new_row >= 0) {
                trips.emplace_back(new_row, new_col, it.value());
            }
        }
    }

    Eigen::SparseMatrix<double> out(static_cast<int>(idx.size()), static_cast<int>(idx.size()));
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

#ifdef ISLAM_HAS_CHOLMOD
std::vector<int> suitesparse_reference_permutation(const Eigen::SparseMatrix<double>& H);
#endif

std::vector<int> amd_permutation(const Eigen::SparseMatrix<double>& H) {
    const int n = H.rows();
    std::vector<int> perm(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) perm[static_cast<size_t>(i)] = i;
    if (n > 1 && H.cols() == n && H.nonZeros() > 0) {
        Eigen::AMDOrdering<int> ordering;
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> perm_matrix;
        ordering(H.selfadjointView<Eigen::Upper>(), perm_matrix);
        for (int k = 0; k < n; ++k) {
            perm[static_cast<size_t>(k)] = perm_matrix.indices()[k];
        }
    }
    return perm;
}

bool standalone_batch_ccolamd_enabled() {
    const char* v = std::getenv("ISLAM_USE_STANDALONE_BATCH_CCOLAMD");
    if (v == nullptr) return false;
    const std::string s(v);
    return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
}

std::vector<int> exact_symbolic_permutation(const Eigen::SparseMatrix<double>& H) {
    if (standalone_batch_ccolamd_enabled()) {
        return DeterministicBatchCcolamd::order_default(H);
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (auto perm = suitesparse_reference_permutation(H); !perm.empty()) {
        return perm;
    }
#endif
    return amd_permutation(H);
}

std::vector<int> expand_pattern_neighborhood(const Eigen::SparseMatrix<double>& H,
                                             const std::vector<int>& seeds,
                                             int hops) {
    const int n = H.rows();
    if (n == 0 || seeds.empty()) return {};

    std::vector<unsigned char> in_set(static_cast<size_t>(n), 0);
    std::vector<int> frontier;
    for (int v : seeds) {
        if (v >= 0 && v < n && !in_set[static_cast<size_t>(v)]) {
            in_set[static_cast<size_t>(v)] = 1;
            frontier.push_back(v);
        }
    }

    for (int hop = 0; hop < hops; ++hop) {
        std::vector<int> next;
        for (int col = 0; col < H.outerSize(); ++col) {
            bool column_touched = false;
            for (Eigen::SparseMatrix<double>::InnerIterator it(H, col); it; ++it) {
                if (in_set[static_cast<size_t>(it.row())]) {
                    column_touched = true;
                    break;
                }
            }
            if (!column_touched && !in_set[static_cast<size_t>(col)]) {
                continue;
            }
            if (!in_set[static_cast<size_t>(col)]) {
                in_set[static_cast<size_t>(col)] = 1;
                next.push_back(col);
            }
            for (Eigen::SparseMatrix<double>::InnerIterator it(H, col); it; ++it) {
                const int r = it.row();
                if (!in_set[static_cast<size_t>(r)]) {
                    in_set[static_cast<size_t>(r)] = 1;
                    next.push_back(r);
                }
            }
        }
        if (next.empty()) break;
        frontier = std::move(next);
    }

    std::vector<int> out;
    for (int i = 0; i < n; ++i) if (in_set[static_cast<size_t>(i)]) out.push_back(i);
    return out;
}

std::vector<int> incremental_local_amd_permutation(const Eigen::SparseMatrix<double>& H,
                                                   const std::vector<int>& base_perm,
                                                   const std::vector<int>& dirty_vars) {
    const int n = H.rows();
    if (static_cast<int>(base_perm.size()) != n) {
        return amd_permutation(H);
    }

    auto dirty = unique_sorted(dirty_vars);
    if (dirty.empty()) {
        return base_perm;
    }

    const auto expanded = expand_pattern_neighborhood(H, dirty, 1);
    if (expanded.empty()) {
        return base_perm;
    }

    if (static_cast<int>(expanded.size()) >= std::max(8, n / 3)) {
        return amd_permutation(H);
    }

    const auto base_pinv = inverse_permutation(base_perm);
    std::vector<int> local_positions;
    local_positions.reserve(expanded.size());
    for (int v : expanded) {
        if (v >= 0 && v < n) {
            const int p = base_pinv[static_cast<size_t>(v)];
            if (p >= 0) local_positions.push_back(p);
        }
    }
    local_positions = unique_sorted(std::move(local_positions));
    if (local_positions.size() < 2) {
        return base_perm;
    }

    std::vector<int> local_vars;
    local_vars.reserve(local_positions.size());
    for (int pos : local_positions) {
        local_vars.push_back(base_perm[static_cast<size_t>(pos)]);
    }

    const Eigen::SparseMatrix<double> H_local = principal_submatrix(H, local_vars);
    const auto local_perm = amd_permutation(H_local);

    std::vector<int> new_perm = base_perm;
    for (int k = 0; k < static_cast<int>(local_positions.size()); ++k) {
        new_perm[static_cast<size_t>(local_positions[static_cast<size_t>(k)])] =
            local_vars[static_cast<size_t>(local_perm[static_cast<size_t>(k)])];
    }
    return new_perm;
}



std::uint64_t pattern_key(int i, int j) {
    if (i > j) std::swap(i, j);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 32u)
         | static_cast<std::uint32_t>(j);
}

std::pair<int,int> decode_pattern_key(std::uint64_t key) {
    const int i = static_cast<int>(static_cast<std::uint32_t>(key >> 32u));
    const int j = static_cast<int>(static_cast<std::uint32_t>(key & 0xffffffffu));
    return {i, j};
}

std::vector<int> dirty_positions_from_vars(const std::vector<int>& perm,
                                           const std::vector<int>& dirty_vars) {
    if (perm.empty() || dirty_vars.empty()) return {};
    const auto pinv = inverse_permutation(perm);
    std::vector<int> pos;
    pos.reserve(dirty_vars.size());
    for (int v : dirty_vars) {
        if (v < 0 || v >= static_cast<int>(pinv.size())) continue;
        const int p = pinv[static_cast<size_t>(v)];
        if (p >= 0) pos.push_back(p);
    }
    return unique_sorted(std::move(pos));
}

std::vector<int> elimination_tree_from_upper_suffix(const Eigen::SparseMatrix<double>& Hperm,
                                                    int start,
                                                    const std::vector<int>& previous_parent) {
    const int n = Hperm.cols();
    if (start <= 0 || static_cast<int>(previous_parent.size()) != n) {
        return elimination_tree_from_upper(Hperm);
    }

    std::vector<int> parent = previous_parent;
    std::vector<int> ancestor(static_cast<size_t>(n), -1);

    const int prefix = std::min(start, n);
    for (int k = 0; k < prefix; ++k) {
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

    for (int k = start; k < n; ++k) {
        parent[static_cast<size_t>(k)] = -1;
    }
    for (int k = start; k < n; ++k) {
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

std::vector<int> incremental_local_etree_update(const Eigen::SparseMatrix<double>& Hperm,
                                                const std::vector<int>& previous_parent,
                                                const std::vector<int>& dirty_positions) {
    const int n = Hperm.cols();
    if (n == 0) return {};
    if (dirty_positions.empty() || static_cast<int>(previous_parent.size()) != n) {
        return elimination_tree_from_upper(Hperm);
    }
    auto expanded = expand_pattern_neighborhood(Hperm, dirty_positions, 1);
    expanded = unique_sorted(std::move(expanded));
    if (expanded.empty()) {
        return previous_parent;
    }
    if (static_cast<int>(expanded.size()) >= std::max(8, n / 3)) {
        return elimination_tree_from_upper(Hperm);
    }
    const int start = std::max(0, expanded.front() - 1);
    return elimination_tree_from_upper_suffix(Hperm, start, previous_parent);
}
#ifdef ISLAM_HAS_CHOLMOD
using ss_long = SuiteSparse_long;

cholmod_sparse* eigen_to_cholmod_sparse(const Eigen::SparseMatrix<double>& in,
                                        int stype,
                                        cholmod_common* common) {
    Eigen::SparseMatrix<double> A = in;
    A.makeCompressed();

    cholmod_sparse* out = cholmod_l_allocate_sparse(
        static_cast<size_t>(A.rows()),
        static_cast<size_t>(A.cols()),
        static_cast<size_t>(A.nonZeros()),
        1,
        1,
        stype,
        CHOLMOD_REAL,
        common);
    if (out == nullptr) {
        throw std::runtime_error("cholmod allocate_sparse failed");
    }

    auto* p = static_cast<ss_long*>(out->p);
    auto* i = static_cast<ss_long*>(out->i);
    auto* x = static_cast<double*>(out->x);

    for (int col = 0; col <= A.cols(); ++col) {
        p[col] = static_cast<ss_long>(A.outerIndexPtr()[col]);
    }
    for (int nz = 0; nz < A.nonZeros(); ++nz) {
        i[nz] = static_cast<ss_long>(A.innerIndexPtr()[nz]);
        x[nz] = A.valuePtr()[nz];
    }
    out->sorted = 1;
    out->packed = 1;
    return out;
}

cholmod_dense* eigen_to_cholmod_dense(const Eigen::VectorXd& rhs,
                                      cholmod_common* common) {
    cholmod_dense* out = cholmod_l_allocate_dense(
        static_cast<size_t>(rhs.size()),
        1,
        static_cast<size_t>(rhs.size()),
        CHOLMOD_REAL,
        common);
    if (out == nullptr) {
        throw std::runtime_error("cholmod allocate_dense failed");
    }
    auto* x = static_cast<double*>(out->x);
    std::memcpy(x, rhs.data(), static_cast<size_t>(rhs.size()) * sizeof(double));
    return out;
}

Eigen::VectorXd cholmod_dense_to_eigen(const cholmod_dense* x) {
    if (x == nullptr || x->ncol != 1) {
        throw std::runtime_error("Invalid cholmod dense vector");
    }
    Eigen::VectorXd out(static_cast<int>(x->nrow));
    const auto* px = static_cast<const double*>(x->x);
    for (int r = 0; r < static_cast<int>(x->nrow); ++r) {
        out[r] = px[r];
    }
    return out;
}

bool is_valid_permutation_vector(const std::vector<int>& perm, int n) {
    if (static_cast<int>(perm.size()) != n) {
        return false;
    }
    std::vector<unsigned char> seen(static_cast<size_t>(n), 0);
    for (int v : perm) {
        if (v < 0 || v >= n || seen[static_cast<size_t>(v)]) {
            return false;
        }
        seen[static_cast<size_t>(v)] = 1;
    }
    return true;
}

std::vector<int> identity_permutation(int n) {
    std::vector<int> perm(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) perm[static_cast<size_t>(i)] = i;
    return perm;
}

std::vector<int> suitesparse_extract_factor_perm(const cholmod_factor* factor, int n) {
    if (factor == nullptr) {
        return {};
    }
    if (factor->Perm == nullptr) {
        return identity_permutation(n);
    }
    std::vector<int> perm(static_cast<size_t>(n), -1);
    const auto* p = static_cast<const ss_long*>(factor->Perm);
    for (int i = 0; i < n; ++i) {
        perm[static_cast<size_t>(i)] = static_cast<int>(p[i]);
    }
    if (!is_valid_permutation_vector(perm, n)) {
        return {};
    }
    return perm;
}

Eigen::SparseMatrix<double> cholmod_factor_to_eigen_lower(const cholmod_factor* factor, int n) {
    if (factor == nullptr || n <= 0) {
        return Eigen::SparseMatrix<double>();
    }
    if (factor->is_super != 0) {
        throw std::runtime_error("CHOLMOD factor-block solve requires simplicial factor storage");
    }
    if (factor->xtype != CHOLMOD_REAL) {
        throw std::runtime_error("CHOLMOD factor-block solve requires real-valued factor storage");
    }
    if (factor->is_ll == 0) {
        throw std::runtime_error("CHOLMOD factor-block solve requires an LL^T factor; set common.final_ll = 1");
    }
    const auto* p = static_cast<const ss_long*>(factor->p);
    const auto* i = static_cast<const ss_long*>(factor->i);
    const auto* x = static_cast<const double*>(factor->x);
    if (p == nullptr || i == nullptr || x == nullptr) {
        throw std::runtime_error("Invalid CHOLMOD factor storage");
    }

    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(p[n]));
    for (int col = 0; col < n; ++col) {
        for (ss_long kk = p[col]; kk < p[col + 1]; ++kk) {
            const int row = static_cast<int>(i[kk]);
            if (row < 0 || row >= n) {
                throw std::runtime_error("CHOLMOD factor row index out of range");
            }
            if (row >= col) {
                trips.emplace_back(row, col, x[kk]);
            }
        }
    }
    Eigen::SparseMatrix<double> L(n, n);
    L.setFromTriplets(trips.begin(), trips.end());
    L.makeCompressed();
    return L;
}

std::vector<int> suitesparse_reference_permutation_try(const Eigen::SparseMatrix<double>& H,
                                                       int ordering_code) {
    const int n = H.rows();
    if (n <= 0 || H.cols() != n) {
        return {};
    }

    cholmod_common common{};
    cholmod_l_start(&common);
    common.nmethods = 1;
    common.method[0].ordering = ordering_code;
    common.postorder = 0;
    common.supernodal = CHOLMOD_SIMPLICIAL;

    cholmod_sparse* A = nullptr;
    cholmod_factor* F = nullptr;
    std::vector<int> perm;
    try {
        A = eigen_to_cholmod_sparse(H, 1, &common);
        F = cholmod_l_analyze(A, &common);
        perm = suitesparse_extract_factor_perm(F, n);
    } catch (...) {
        perm.clear();
    }

    if (F != nullptr) {
        cholmod_l_free_factor(&F, &common);
    }
    if (A != nullptr) {
        cholmod_l_free_sparse(&A, &common);
    }
    cholmod_l_finish(&common);
    return perm;
}

std::vector<int> suitesparse_reference_permutation(const Eigen::SparseMatrix<double>& H) {
#if defined(CHOLMOD_COLAMD)
    if (auto perm = suitesparse_reference_permutation_try(H, CHOLMOD_COLAMD); !perm.empty()) {
        return perm;
    }
#endif
#if defined(CHOLMOD_AMD)
    if (auto perm = suitesparse_reference_permutation_try(H, CHOLMOD_AMD); !perm.empty()) {
        return perm;
    }
#endif
#if defined(CHOLMOD_NATURAL)
    if (auto perm = suitesparse_reference_permutation_try(H, CHOLMOD_NATURAL); !perm.empty()) {
        return perm;
    }
#endif
    return {};
}

cholmod_sparse* standard_basis_columns(const std::vector<int>& vars,
                                       int n,
                                       double scale,
                                       cholmod_common* common) {
    cholmod_sparse* out = cholmod_l_allocate_sparse(
        static_cast<size_t>(n),
        static_cast<size_t>(vars.size()),
        static_cast<size_t>(vars.size()),
        1,
        1,
        0,
        CHOLMOD_REAL,
        common);
    if (out == nullptr) {
        throw std::runtime_error("cholmod allocate_sparse failed for basis columns");
    }
    auto* p = static_cast<ss_long*>(out->p);
    auto* i = static_cast<ss_long*>(out->i);
    auto* x = static_cast<double*>(out->x);
    for (int col = 0; col < static_cast<int>(vars.size()); ++col) {
        p[col] = static_cast<ss_long>(col);
        i[col] = static_cast<ss_long>(vars[static_cast<size_t>(col)]);
        x[col] = scale;
    }
    p[vars.size()] = static_cast<ss_long>(vars.size());
    out->sorted = 1;
    out->packed = 1;
    return out;
}
#endif

} // namespace

struct FactorManager::CholmodState {
#ifdef ISLAM_HAS_CHOLMOD
    cholmod_common common{};
    cholmod_factor* factor = nullptr;
    bool started = false;
#endif
};

FactorManager::FactorManager()
    : cholmod_state_(std::make_unique<CholmodState>()),
      sparse_expanding_state_(std::make_unique<SparseExpandingCholesky>()),
      symbolic_engine_(std::make_unique<SymbolicEngine>()) {
#ifdef ISLAM_HAS_CHOLMOD
    cholmod_l_start(&cholmod_state_->common);
    cholmod_state_->started = true;
    cholmod_state_->common.nmethods = 1;
    cholmod_state_->common.method[0].ordering = CHOLMOD_NATURAL;
    cholmod_state_->common.postorder = 0;
    cholmod_state_->common.supernodal = CHOLMOD_SIMPLICIAL;
    cholmod_state_->common.final_ll = 1;
    using_cholmod_ = true;
    supports_incremental_updates_ = true;
#else
    supports_incremental_updates_ = true;
#endif
}

FactorManager::~FactorManager() {
#ifdef ISLAM_HAS_CHOLMOD
    if (cholmod_state_ && cholmod_state_->started) {
        if (cholmod_state_->factor != nullptr) {
            cholmod_l_free_factor(&cholmod_state_->factor, &cholmod_state_->common);
        }
        cholmod_l_finish(&cholmod_state_->common);
        cholmod_state_->started = false;
    }
#endif
}

FactorManager::FactorManager(FactorManager&& other) noexcept = default;
FactorManager& FactorManager::operator=(FactorManager&& other) noexcept = default;

void FactorManager::clear() {
    state_size_ = 0;
    factor_size_ = 0;
    anchor_strength_ = 1.0;
    latent_prior_strength_ = 0.0;
    factorized_ = false;
    cholmod_growth_refactor_pending_ = false;
    if (sparse_expanding_state_) {
        sparse_expanding_state_->clear();
    }
    factorization_stats_ = FactorizationStats{};
    last_H_.resize(0, 0);
    last_g_.resize(0);
    H_current_.resize(0, 0);
    g_current_.resize(0);
    edge_cache_.clear();
    L_dense_.resize(0, 0);
    active_vars_.clear();
    sparse_factor_perm_.clear();
    sparse_factor_pinv_.clear();
    cached_H_perm_.resize(0, 0);
    cached_g_perm_.resize(0);
    numeric_revision_ = 1;
    cached_numeric_revision_ = 0;
    if (symbolic_engine_) {
        symbolic_engine_->clear();
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (cholmod_state_ && cholmod_state_->started && cholmod_state_->factor != nullptr) {
        cholmod_l_free_factor(&cholmod_state_->factor, &cholmod_state_->common);
    }
#endif
}


void FactorManager::enable_sparse_expanding_cholesky(bool enable) {
    if (state_size_ != 0 || factorized_ || !edge_cache_.empty()) {
        throw std::runtime_error("enable_sparse_expanding_cholesky must be called before adding variables or factors");
    }
    use_sparse_expanding_backend_ = enable;
    if (enable) {
        using_cholmod_ = false;
        if (!sparse_expanding_state_) {
            sparse_expanding_state_ = std::make_unique<SparseExpandingCholesky>();
        }
    } else {
#ifdef ISLAM_HAS_CHOLMOD
        using_cholmod_ = true;
#else
        using_cholmod_ = false;
#endif
    }
}

void FactorManager::force_dense_backend(bool enable) {
    if (state_size_ != 0 || factorized_ || !edge_cache_.empty()) {
        throw std::runtime_error("force_dense_backend must be called before adding variables or factors");
    }
    if (enable) {
        use_sparse_expanding_backend_ = false;
        using_cholmod_ = false;
    } else {
#ifdef ISLAM_HAS_CHOLMOD
        using_cholmod_ = true;
#else
        using_cholmod_ = false;
#endif
    }
}

void FactorManager::configure_incremental(double anchor_strength,
                                          double latent_prior_strength,
                                          int anchor_dim) {
    if (state_size_ != 0 || factorized_ || !edge_cache_.empty()) {
        throw std::runtime_error("configure_incremental must be called before adding variables or factors");
    }
    anchor_strength_ = anchor_strength;
    latent_prior_strength_ = std::max(0.0, latent_prior_strength);
    anchor_dim_ = std::max(0, anchor_dim);
}

void FactorManager::reserve_full_state(int n,
                                       double anchor_strength,
                                       double latent_prior_strength,
                                       int anchor_dim) {
    if (n < 0) {
        throw std::runtime_error("Negative state size in reserve_full_state");
    }
    clear();
    ++factorization_stats_.reserve_full_state_calls;
    state_size_ = n;
    factor_size_ = n;
    anchor_strength_ = anchor_strength;
    latent_prior_strength_ = std::max(0.0, latent_prior_strength);
    anchor_dim_ = std::clamp(anchor_dim, 0, std::max(0, n));
    active_vars_.assign(static_cast<size_t>(n), 0);
    H_current_ = diagonal_sparse(n, latent_prior_strength_);
    g_current_ = Eigen::VectorXd::Zero(n);
    if (symbolic_engine_) {
        symbolic_engine_->reserve_state(n);
    }
    invalidate_structure_cache();
    update_cached_system_views();
    rebuild_symbolic_pattern_from_matrix(last_H_);
    factorize(last_H_);
}

void FactorManager::ensure_state_size(int n) {
    if (n < 0) {
        throw std::runtime_error("Negative state size in FactorManager");
    }
    if (n < state_size_) {
        throw std::runtime_error("FactorManager::ensure_state_size does not support shrinking");
    }
    if (n == state_size_) {
        return;
    }

    ++factorization_stats_.state_expansions;

    if (state_size_ == 0) {
        state_size_ = n;
        H_current_.resize(n, n);
        H_current_.setZero();
        g_current_ = Eigen::VectorXd::Zero(n);
        last_H_.resize(n, n);
        last_H_.setZero();
        last_g_ = Eigen::VectorXd::Zero(n);
        factor_size_ = 0;
        factorized_ = false;
        active_vars_.assign(static_cast<size_t>(n), 1);
        if (symbolic_engine_) {
            symbolic_engine_->reserve_state(n);
        }
        invalidate_structure_cache();
        return;
    }

    H_current_ = resize_sparse_preserve(H_current_, n);
    last_H_ = resize_sparse_preserve(last_H_, n);

    const int old_n = state_size_;
    g_current_.conservativeResize(n);
    g_current_.segment(old_n, n - old_n).setZero();
    last_g_.conservativeResize(n);
    last_g_.segment(old_n, n - old_n).setZero();

    state_size_ = n;
    active_vars_.resize(static_cast<size_t>(n), 1);
    if (use_sparse_expanding_backend_ && sparse_expanding_state_ && factorized_) {
        sparse_expanding_state_->grow_to(n);
        factor_size_ = old_n;
    }
    anchor_dim_ = std::clamp(anchor_dim_, 0, std::max(0, state_size_));

#ifdef ISLAM_HAS_CHOLMOD
    // CHOLMOD cannot safely apply an up/down-date to a factor whose dimension
    // has changed.  We therefore preserve the numeric/symbolic system, mark the
    // numeric factor stale, and refactorize the current-size system after the
    // next newly arrived contribution has been inserted.  This gives true
    // unknown-size online expansion semantics; the dense fallback can still use
    // the explicit block-extension routine below when possible.
    if (using_cholmod_ && factorized_) {
        if (cholmod_state_ && cholmod_state_->started && cholmod_state_->factor != nullptr) {
            cholmod_l_free_factor(&cholmod_state_->factor, &cholmod_state_->common);
        }
        factorized_ = false;
        factor_size_ = 0;
        cholmod_growth_refactor_pending_ = true;
    }
#endif

    if (symbolic_engine_) {
        symbolic_engine_->reserve_state(n);
        symbolic_engine_->rebuild_from_numeric_matrix(last_H_);
    }
    invalidate_structure_cache();
    update_cached_system_views();
}

NormalEquations FactorManager::build_normal_equations(
    const Graph& g,
    const std::vector<int>* edge_ids,
    double anchor_strength) const {

    if (g.state_size() <= 0) {
        throw std::runtime_error("Cannot build normal equations for empty graph");
    }

    std::vector<int> ids;
    if (edge_ids != nullptr) {
        ids = *edge_ids;
    } else {
        ids.resize(static_cast<int>(g.edges.size()));
        for (int i = 0; i < static_cast<int>(g.edges.size()); ++i) ids[i] = i;
    }

    Eigen::VectorXd gvec = Eigen::VectorXd::Zero(g.state_size());
    std::vector<Eigen::Triplet<double>> trips;

    for (int eid : ids) {
        if (eid < 0 || eid >= static_cast<int>(g.edges.size())) {
            throw std::out_of_range("Edge id out of range in build_normal_equations");
        }
        const auto contrib = build_edge_contribution(g, eid);
        const Eigen::SparseMatrix<double> He = contribution_hessian(contrib);
        const Eigen::VectorXd ge = contrib.C * contrib.r;
        gvec += ge;
        add_sparse_triplets(He, trips);
    }

    const int anchor_dim = default_anchor_dim_for_graph(g);
    if (anchor_strength > 0.0) {
        for (int i = 0; i < anchor_dim; ++i) {
            trips.emplace_back(i, i, anchor_strength);
        }
    }

    Eigen::SparseMatrix<double> H(g.state_size(), g.state_size());
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();

    return NormalEquations{std::move(H), std::move(gvec)};
}

EdgeContribution FactorManager::build_edge_contribution(const Graph& g, int edge_id) const {
    if (edge_id < 0 || edge_id >= static_cast<int>(g.edges.size())) {
        throw std::out_of_range("Edge id out of range in build_edge_contribution");
    }
    const auto& edge = g.edges[static_cast<size_t>(edge_id)];
    const auto lin = jacobian_edge_jr(edge, g);

    EdgeContribution contrib;
    contrib.edge_id = edge_id;
    contrib.C = lin.J.transpose();
    contrib.r = lin.r;

    if (edge.from_idx >= 0) {
        const auto it_from = g.id_lookup.find(edge.from_id);
        if (it_from == g.id_lookup.end()) {
            throw std::runtime_error("build_edge_contribution: missing from node");
        }
        append_unique_touched_range(contrib.touched_vars, edge.from_idx, it_from->second.dimension);
        append_unique_int(contrib.touched_nodes, edge.from_id);
    }
    if (!edge.is_unary() && edge.to_idx >= 0) {
        const auto it_to = g.id_lookup.find(edge.to_id);
        if (it_to == g.id_lookup.end()) {
            throw std::runtime_error("build_edge_contribution: missing to node");
        }
        append_unique_touched_range(contrib.touched_vars, edge.to_idx, it_to->second.dimension);
        append_unique_int(contrib.touched_nodes, edge.to_id);
    }
    return contrib;
}

Eigen::SparseMatrix<double> FactorManager::contribution_hessian(const EdgeContribution& contrib) {
    return (contrib.C * contrib.C.transpose()).eval();
}

void FactorManager::add_anchor_to_current_system() {
    if (state_size_ <= 0) {
        last_H_.resize(0, 0);
        last_g_.resize(0);
        return;
    }

    std::vector<Eigen::Triplet<double>> trips;
    const int anchor_dim = std::clamp(anchor_dim_, 0, state_size_);
    trips.reserve(static_cast<size_t>(H_current_.nonZeros() + anchor_dim));
    add_sparse_triplets(H_current_, trips);
    for (int i = 0; i < anchor_dim; ++i) {
        trips.emplace_back(i, i, anchor_strength_);
    }
    Eigen::SparseMatrix<double> Hanch(state_size_, state_size_);
    Hanch.setFromTriplets(trips.begin(), trips.end());
    Hanch.makeCompressed();
    last_H_ = std::move(Hanch);
    last_g_ = g_current_;
}

void FactorManager::update_cached_system_views() {
    add_anchor_to_current_system();
    invalidate_numeric_cache();
}

void FactorManager::rebuild_symbolic_pattern_from_matrix(const Eigen::SparseMatrix<double>& H) {
    if (symbolic_engine_) {
        symbolic_engine_->rebuild_from_numeric_matrix(H);
    }
}

void FactorManager::apply_contribution_to_symbolic_pattern(const EdgeContribution& contrib, int delta) {
    if (delta == 0 || !symbolic_engine_) return;
    symbolic_engine_->apply_contribution_pattern(contribution_hessian(contrib), contrib.touched_vars, delta);
    invalidate_numeric_cache();
}

void FactorManager::invalidate_numeric_cache() const {
    ++numeric_revision_;
    cached_numeric_revision_ = 0;
}

void FactorManager::invalidate_structure_cache() const {
    if (symbolic_engine_) {
        symbolic_engine_->note_full_refresh();
    }
    invalidate_numeric_cache();
}

void FactorManager::note_structural_change(const std::vector<int>& vars) const {
    (void)vars;
    if (symbolic_engine_) {
        symbolic_engine_->note_full_refresh();
    }
    invalidate_numeric_cache();
}

bool FactorManager::same_symbolic_pattern(const EdgeContribution& a,
                                          const EdgeContribution& b) const {
    return a.C.rows() == b.C.rows() && a.C.cols() == b.C.cols() && a.touched_vars == b.touched_vars;
}

void FactorManager::refresh_symbolic_cache_if_needed() const {
    if (!symbolic_engine_) {
        return;
    }
    symbolic_engine_->refresh_if_needed([](const Eigen::SparseMatrix<double>& pattern) {
        return exact_symbolic_permutation(pattern);
    });
}

void FactorManager::refresh_numeric_cache_if_needed() const {
    refresh_symbolic_cache_if_needed();
    if (cached_numeric_revision_ == numeric_revision_) {
        return;
    }
    const auto& snap = symbolic_engine_->snapshot();
    cached_H_perm_ = symmetric_permute_by_order(last_H_, snap.perm);
    cached_g_perm_ = permute_vector_by_order(last_g_, snap.perm);
    cached_numeric_revision_ = numeric_revision_;
}

void FactorManager::refactorize_current_system() {
    if (state_size_ <= 0) {
        throw std::runtime_error("Cannot factorize an empty system");
    }
    update_cached_system_views();
    factorize(last_H_);
    last_g_ = g_current_;
}

void FactorManager::rebuild_from_graph(
    const Graph& g,
    const std::vector<int>* edge_ids,
    double anchor_strength) {

    clear();
    anchor_strength_ = anchor_strength;
    anchor_dim_ = default_anchor_dim_for_graph(g);
    ensure_state_size(g.state_size());
    H_current_.setZero();
    g_current_.setZero();
    edge_cache_.clear();
    active_vars_.assign(static_cast<size_t>(g.state_size()), 1);

    std::vector<int> ids;
    if (edge_ids != nullptr) {
        ids = *edge_ids;
    } else {
        ids.resize(static_cast<int>(g.edges.size()));
        for (int i = 0; i < static_cast<int>(g.edges.size()); ++i) ids[i] = i;
    }

    for (int eid : ids) {
        auto contrib = build_edge_contribution(g, eid);
        H_current_ += contribution_hessian(contrib);
        g_current_ += contrib.C * contrib.r;
        edge_cache_[eid] = std::move(contrib);
    }

    invalidate_structure_cache();
    refactorize_current_system();
}

void FactorManager::apply_contribution_to_system(const EdgeContribution& contrib, bool add) {
    const Eigen::SparseMatrix<double> He = contribution_hessian(contrib);
    const Eigen::VectorXd ge = contrib.C * contrib.r;
    if (add) {
        H_current_ += He;
        g_current_ += ge;
    } else {
        H_current_ -= He;
        g_current_ -= ge;
    }
    update_cached_system_views();
}

void FactorManager::apply_diagonal_rank_update(const std::vector<int>& vars, double weight, bool add) {
    if (vars.empty() || weight <= 0.0) {
        return;
    }
    if (use_sparse_expanding_backend_) {
        if (!factorized_ || factor_size_ != state_size_) {
            throw std::runtime_error("apply_diagonal_rank_update requires a full-size factor");
        }
        if (refresh_sparse_expanding_ordering_after_system_update(vars)) {
            return;
        }
        if (sparse_factor_pinv_.empty()) {
            sparse_factor_perm_ = identity_permutation(state_size_);
            sparse_factor_pinv_ = inverse_permutation(sparse_factor_perm_);
        }
        const std::vector<int> factor_vars = permute_indices_to_positions(vars, sparse_factor_pinv_);
        sparse_expanding_state_->apply_diagonal_update(factor_vars, weight, add);
        factor_size_ = sparse_expanding_state_->factor_size();
        factorized_ = sparse_expanding_state_->factorized();
        sync_sparse_expanding_stats();
        return;
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (using_cholmod_) {
        auto* C = standard_basis_columns(vars, state_size_, std::sqrt(weight), &cholmod_state_->common);
        const int ok = cholmod_l_updown(add ? 1 : 0, C, cholmod_state_->factor, &cholmod_state_->common);
        cholmod_l_free_sparse(&C, &cholmod_state_->common);
        if (!ok) {
            throw std::runtime_error(add ? "CHOLMOD diagonal update failed" : "CHOLMOD diagonal downdate failed");
        }
        return;
    }
#endif

    if (!factorized_ || factor_size_ != state_size_) {
        throw std::runtime_error("apply_diagonal_rank_update requires a full-size factor");
    }
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(state_size_, static_cast<int>(vars.size()));
    const double scale = std::sqrt(weight);
    for (int col = 0; col < static_cast<int>(vars.size()); ++col) {
        C(vars[static_cast<size_t>(col)], col) = scale;
    }
    if (!chol_rank_k_update_lower(L_dense_, C, add)) {
        throw std::runtime_error(add ? "Dense diagonal update failed" : "Dense diagonal downdate failed");
    }
}

void FactorManager::activate_vars(const std::vector<int>& vars) {
    if (vars.empty()) {
        return;
    }
    if (state_size_ <= 0) {
        throw std::runtime_error("activate_vars called before reserve/factorize");
    }

    std::vector<int> newly_active;
    newly_active.reserve(vars.size());
    for (int v : vars) {
        if (v < 0 || v >= state_size_) {
            throw std::out_of_range("activate_vars index out of range");
        }
        if (!active_vars_[static_cast<size_t>(v)]) {
            active_vars_[static_cast<size_t>(v)] = 1;
            newly_active.push_back(v);
        }
    }
    if (newly_active.empty()) {
        return;
    }

    if (latent_prior_strength_ > 0.0) {
        note_structural_change(newly_active);
        for (int v : newly_active) {
            H_current_.coeffRef(v, v) -= latent_prior_strength_;
        }
        H_current_.makeCompressed();
        update_cached_system_views();
        if (factorized_) {
            apply_diagonal_rank_update(newly_active, latent_prior_strength_, false);
        }
    }
}

void FactorManager::apply_full_size_update_to_factor(const EdgeContribution& contrib, bool add) {
    if (!use_sparse_expanding_backend_) {
        if (add) {
            ++factorization_stats_.same_size_rank_updates;
        } else {
            ++factorization_stats_.same_size_rank_downdates;
        }
    }
    if (!factorized_ || factor_size_ != state_size_) {
        throw std::runtime_error("apply_full_size_update_to_factor requires a full-size factor");
    }
    if (use_sparse_expanding_backend_) {
        apply_sparse_expanding_update_to_factor(contrib, add);
        return;
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (using_cholmod_) {
        auto* C = eigen_to_cholmod_sparse(contrib.C, 0, &cholmod_state_->common);
        const int ok = cholmod_l_updown(add ? 1 : 0, C, cholmod_state_->factor, &cholmod_state_->common);
        cholmod_l_free_sparse(&C, &cholmod_state_->common);
        if (!ok) {
            throw std::runtime_error(add ? "CHOLMOD incremental update failed" : "CHOLMOD incremental downdate failed");
        }
        return;
    }
#endif
    const Eigen::MatrixXd C = sparse_to_dense(contrib.C);
    if (!chol_rank_k_update_lower(L_dense_, C, add)) {
        throw std::runtime_error(add
            ? "Incremental Cholesky update failed"
            : "Incremental Cholesky downdate failed");
    }
}

void FactorManager::extend_factor_with_contribution(const EdgeContribution& contrib) {
    if (!factorized_) {
        throw std::runtime_error("extend_factor_with_contribution requires an existing factor");
    }
    if (factor_size_ >= state_size_) {
        throw std::runtime_error("extend_factor_with_contribution called without pending growth");
    }

    const int old_n = factor_size_;
    const int new_n = state_size_;
    const int k = new_n - old_n;
    const Eigen::MatrixXd Cfull = sparse_to_dense(contrib.C);
    const Eigen::MatrixXd U = Cfull.topRows(old_n);
    const Eigen::MatrixXd V = Cfull.bottomRows(k);

    if (U.cols() > 0 && !chol_rank_k_update_lower(L_dense_, U, true)) {
        throw std::runtime_error("Incremental update of old block failed during growth");
    }

    const Eigen::MatrixXd B = U * V.transpose();
    const Eigen::MatrixXd Cblk = V * V.transpose();
    Eigen::MatrixXd Y(old_n, k);
    if (old_n > 0) {
        Y = L_dense_.triangularView<Eigen::Lower>().solve(B);
    } else {
        Y.resize(0, k);
    }

    Eigen::MatrixXd S = Cblk;
    if (old_n > 0) {
        S -= Y.transpose() * Y;
    }
    S = 0.5 * (S + S.transpose());

    Eigen::LLT<Eigen::MatrixXd> llt(S);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("State growth Schur factorization failed; new vars are not sufficiently constrained");
    }
    const Eigen::MatrixXd L22 = llt.matrixL();

    Eigen::MatrixXd Lnew = Eigen::MatrixXd::Zero(new_n, new_n);
    if (old_n > 0) {
        Lnew.topLeftCorner(old_n, old_n) = L_dense_;
        Lnew.bottomLeftCorner(k, old_n) = Y.transpose();
    }
    Lnew.bottomRightCorner(k, k) = L22;
    L_dense_ = std::move(Lnew);
    factor_size_ = new_n;
    factorized_ = true;
    ++factorization_stats_.dense_growth_extensions;
}

void FactorManager::sync_sparse_expanding_stats() {
    if (!sparse_expanding_state_) {
        return;
    }
    const auto& st = sparse_expanding_state_->stats();
    factorization_stats_.custom_sparse_full_factorizations = st.full_factorizations;
    factorization_stats_.custom_sparse_suffix_refactorizations = st.suffix_refactorizations;
    factorization_stats_.custom_sparse_growth_updates = st.expansion_suffix_updates;
    factorization_stats_.custom_sparse_dynamic_reorders = st.dynamic_reorder_refactorizations;
    factorization_stats_.custom_sparse_prefix_reuses = st.reordered_prefix_reuses;
    factorization_stats_.custom_sparse_prefix_columns_reused = st.prefix_columns_reused;
    factorization_stats_.custom_sparse_l21_recomputes = st.l21_recomputes;
    factorization_stats_.custom_sparse_left_looking_factorizations = st.sparse_left_looking_factorizations;
    factorization_stats_.custom_sparse_schur_complements = st.sparse_schur_complements;
    factorization_stats_.custom_sparse_triangular_solves = st.sparse_triangular_solves;
    factorization_stats_.custom_sparse_jitter_regularizations = st.jitter_regularizations;
    factorization_stats_.custom_sparse_affected_closure_refactorizations = st.affected_closure_refactorizations;
    factorization_stats_.custom_sparse_affected_closure_fallbacks = st.affected_closure_fallbacks;
    factorization_stats_.custom_sparse_affected_columns_refactored = st.affected_columns_refactored;
    factorization_stats_.custom_sparse_etree_closure_computations = st.etree_closure_computations;
    factorization_stats_.custom_sparse_structural_pattern_changes = st.structural_pattern_changes;
    factorization_stats_.custom_sparse_numeric_only_updates = st.numeric_only_updates;
    factorization_stats_.custom_sparse_symbolic_pattern_classifications = st.symbolic_pattern_classifications;
    factorization_stats_.custom_sparse_structural_factor_pattern_changes = st.structural_factor_pattern_changes;
    factorization_stats_.custom_sparse_structural_factor_columns_changed = st.structural_factor_columns_changed;
    factorization_stats_.custom_sparse_pattern_stable_structural_updates = st.pattern_stable_structural_updates;
    factorization_stats_.custom_sparse_structural_closure_attempts = st.structural_closure_attempts;
    factorization_stats_.custom_sparse_structural_closure_refactorizations = st.structural_closure_refactorizations;
    factorization_stats_.custom_sparse_structural_closure_fallbacks = st.structural_closure_fallbacks;
    factorization_stats_.custom_sparse_affected_closure_certifications = st.affected_closure_certifications;
    factorization_stats_.custom_sparse_affected_closure_certification_failures = st.affected_closure_certification_failures;
    factorization_stats_.custom_sparse_factorization_residual_checks = st.factorization_residual_checks;
    factorization_stats_.custom_sparse_column_local_certifications = st.column_local_certifications;
    factorization_stats_.custom_sparse_column_local_certification_failures = st.column_local_certification_failures;
    factorization_stats_.custom_sparse_column_local_certification_columns = st.column_local_certification_columns;
    factorization_stats_.custom_sparse_full_certification_fallbacks = st.full_certification_fallbacks;
    factorization_stats_.custom_sparse_dependency_cache_rebuilds = st.dependency_cache_rebuilds;
    factorization_stats_.custom_sparse_dependency_cache_column_refreshes = st.dependency_cache_column_refreshes;
    factorization_stats_.custom_sparse_dependency_cache_invalidations = st.dependency_cache_invalidations;
    factorization_stats_.custom_sparse_dependency_cache_hits = st.dependency_cache_hits;
    factorization_stats_.custom_sparse_dependency_closure_computations = st.dependency_closure_computations;
    factorization_stats_.custom_sparse_dependency_closure_columns = st.dependency_closure_columns;
    factorization_stats_.custom_sparse_etree_closure_cache_bypasses = st.etree_closure_cache_bypasses;
    factorization_stats_.custom_sparse_last_certification_residual = st.last_certification_residual;
    factorization_stats_.custom_sparse_last_column_local_certification_residual = st.last_column_local_certification_residual;
    factorization_stats_.custom_sparse_max_column_local_certification_residual = st.max_column_local_certification_residual;
    factorization_stats_.custom_sparse_max_certification_residual = st.max_certification_residual;
    factorization_stats_.custom_sparse_last_dependency_closure_size = st.last_dependency_closure_size;
}

bool FactorManager::refresh_sparse_expanding_ordering_after_system_update(const std::vector<int>& dirty_vars) {
    if (!use_sparse_expanding_backend_ || !sparse_expanding_state_ || !factorized_) {
        return false;
    }
    if (last_H_.rows() != state_size_) {
        return false;
    }

    refresh_symbolic_cache_if_needed();
    std::vector<int> new_perm;
    if (symbolic_engine_ && static_cast<int>(symbolic_engine_->snapshot().perm.size()) == state_size_) {
        new_perm = symbolic_engine_->snapshot().perm;
    } else {
        new_perm = identity_permutation(state_size_);
    }

    if (sparse_factor_perm_.empty()) {
        sparse_factor_perm_ = new_perm;
        sparse_factor_pinv_ = inverse_permutation(sparse_factor_perm_);
        return false;
    }

    if (new_perm == sparse_factor_perm_) {
        return false;
    }

    int reusable_prefix = common_prefix_size(sparse_factor_perm_, new_perm);
    if (!dirty_vars.empty()) {
        // A variable touched by the latest structural/numeric change cannot be
        // kept in the supposedly reusable eliminated prefix.  This is a
        // conservative exactness guard for the custom dynamic sparse path.
        for (int v : dirty_vars) {
            if (v >= 0 && v < static_cast<int>(new_perm.size())) {
                const auto it = std::find(new_perm.begin(), new_perm.end(), v);
                if (it != new_perm.end()) {
                    const int pos = static_cast<int>(std::distance(new_perm.begin(), it));
                    reusable_prefix = std::min(reusable_prefix, pos);
                }
            }
        }
    }

    const Eigen::SparseMatrix<double> H_perm = symmetric_permute_by_order(last_H_, new_perm);
    sparse_expanding_state_->refactorize_with_reused_prefix(H_perm, reusable_prefix);
    sparse_factor_perm_ = std::move(new_perm);
    sparse_factor_pinv_ = inverse_permutation(sparse_factor_perm_);
    factor_size_ = sparse_expanding_state_->factor_size();
    factorized_ = sparse_expanding_state_->factorized();
    sync_sparse_expanding_stats();
    return true;
}

void FactorManager::apply_sparse_expanding_update_to_factor(const EdgeContribution& contrib, bool add) {
    if (!use_sparse_expanding_backend_ || !sparse_expanding_state_) {
        throw std::runtime_error("Sparse expanding backend is not enabled");
    }
    if (!factorized_) {
        throw std::runtime_error("Sparse expanding update requires an existing factor");
    }
    if (contrib.C.rows() != state_size_) {
        throw std::runtime_error("Sparse expanding contribution dimension mismatch");
    }

    // If the exact dynamic symbolic engine changed the factor ordering after
    // this structural/numeric update, rebuild the custom sparse factor in the
    // new ordering while reusing the certified common prefix where safe.  The
    // factor is then already current with last_H_, so applying this contribution
    // again would double-count it.
    if (refresh_sparse_expanding_ordering_after_system_update(contrib.touched_vars)) {
        return;
    }

    if (sparse_factor_perm_.empty()) {
        sparse_factor_perm_ = identity_permutation(state_size_);
        sparse_factor_pinv_ = inverse_permutation(sparse_factor_perm_);
    }
    if (static_cast<int>(sparse_factor_perm_.size()) != state_size_) {
        throw std::runtime_error("Sparse expanding factor permutation is stale");
    }

    const Eigen::SparseMatrix<double> C_factor = row_permute_by_order(contrib.C, sparse_factor_perm_);
    const std::vector<int> touched_factor =
        permute_indices_to_positions(contrib.touched_vars, sparse_factor_pinv_);
    sparse_expanding_state_->apply_contribution(C_factor, touched_factor, add);
    factor_size_ = sparse_expanding_state_->factor_size();
    factorized_ = sparse_expanding_state_->factorized();
    sync_sparse_expanding_stats();
}

void FactorManager::add_edge_contribution(int edge_id, const EdgeContribution& contrib) {
    ++factorization_stats_.incremental_edge_adds;
    if (state_size_ <= 0) {
        throw std::runtime_error("Call rebuild_from_graph or reserve_full_state before add_edge_contribution");
    }
    if (contrib.C.rows() != state_size_) {
        throw std::runtime_error("Contribution state dimension mismatch in add_edge_contribution");
    }
    if (edge_cache_.count(edge_id) > 0) {
        throw std::runtime_error("Duplicate edge id in add_edge_contribution");
    }

    apply_contribution_to_symbolic_pattern(contrib, +1);
    apply_contribution_to_system(contrib, true);
    edge_cache_[edge_id] = contrib;

    if (!factorized_) {
        refactorize_current_system();
        return;
    }

    if (factor_size_ < state_size_) {
        if (use_sparse_expanding_backend_) {
            apply_sparse_expanding_update_to_factor(contrib, true);
        } else {
            extend_factor_with_contribution(contrib);
        }
    } else {
        apply_full_size_update_to_factor(contrib, true);
    }
}

void FactorManager::remove_edge_contribution(int edge_id) {
    ++factorization_stats_.edge_removes;
    const auto it = edge_cache_.find(edge_id);
    if (it == edge_cache_.end()) {
        throw std::runtime_error("Missing edge id in remove_edge_contribution");
    }
    const auto old = it->second;
    if (factor_size_ != state_size_) {
        throw std::runtime_error("Cannot remove edge while factor growth is pending");
    }

    apply_contribution_to_symbolic_pattern(old, -1);
    apply_contribution_to_system(old, false);
    edge_cache_.erase(it);
    apply_full_size_update_to_factor(old, false);
}

void FactorManager::replace_edge_contribution(int edge_id, const EdgeContribution& contrib) {
    ++factorization_stats_.edge_replaces;
    const auto it = edge_cache_.find(edge_id);
    if (it == edge_cache_.end()) {
        throw std::runtime_error("Missing edge id in replace_edge_contribution");
    }
    if (contrib.C.rows() != state_size_) {
        throw std::runtime_error("Contribution state dimension mismatch in replace_edge_contribution");
    }
    if (factor_size_ != state_size_) {
        throw std::runtime_error("Cannot replace edge while factor growth is pending");
    }

    const auto old = it->second;
    if (!same_symbolic_pattern(old, contrib)) {
        apply_contribution_to_symbolic_pattern(old, -1);
        apply_contribution_to_symbolic_pattern(contrib, +1);
    }
    apply_contribution_to_system(old, false);
    apply_full_size_update_to_factor(old, false);

    apply_contribution_to_system(contrib, true);
    apply_full_size_update_to_factor(contrib, true);
    it->second = contrib;
}

void FactorManager::factorize(const Eigen::SparseMatrix<double>& H) {
    ++factorization_stats_.full_refactorizations;
    if (cholmod_growth_refactor_pending_) {
        ++factorization_stats_.cholmod_growth_refactorizations;
        cholmod_growth_refactor_pending_ = false;
    }
    if (H.rows() != H.cols()) {
        throw std::runtime_error("Normal matrix must be square");
    }

    last_H_ = H;
    rebuild_symbolic_pattern_from_matrix(last_H_);
    invalidate_structure_cache();
    if (use_sparse_expanding_backend_) {
        if (!sparse_expanding_state_) {
            sparse_expanding_state_ = std::make_unique<SparseExpandingCholesky>();
        }
        refresh_symbolic_cache_if_needed();
        if (symbolic_engine_ && static_cast<int>(symbolic_engine_->snapshot().perm.size()) == H.rows()) {
            sparse_factor_perm_ = symbolic_engine_->snapshot().perm;
        } else {
            sparse_factor_perm_ = identity_permutation(H.rows());
        }
        sparse_factor_pinv_ = inverse_permutation(sparse_factor_perm_);
        const Eigen::SparseMatrix<double> H_factor = symmetric_permute_by_order(H, sparse_factor_perm_);
        sparse_expanding_state_->factorize(H_factor);
        factor_size_ = H.rows();
        factorized_ = true;
        using_cholmod_ = false;
        sync_sparse_expanding_stats();
        return;
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (using_cholmod_) {
        // Analyze with the repo-owned symbolic permutation rather than letting
        // CHOLMOD silently choose an unrelated ordering.  This keeps the
        // maintained numeric factor, the exposed factor-block SPO path, and the
        // dynamic CCOLAMD/etree symbolic state aligned with the paper's \pi_t.
        refresh_symbolic_cache_if_needed();
        const auto& snap = symbolic_engine_->snapshot();
        std::vector<ss_long> user_perm;
        user_perm.reserve(static_cast<size_t>(H.rows()));
        if (static_cast<int>(snap.perm.size()) == H.rows()) {
            for (int v : snap.perm) {
                user_perm.push_back(static_cast<ss_long>(v));
            }
        }

        auto* A = eigen_to_cholmod_sparse(H, 1, &cholmod_state_->common);
        if (cholmod_state_->factor != nullptr) {
            cholmod_l_free_factor(&cholmod_state_->factor, &cholmod_state_->common);
        }
        if (!user_perm.empty()) {
            cholmod_state_->factor = cholmod_l_analyze_p(
                A,
                user_perm.data(),
                nullptr,
                0,
                &cholmod_state_->common);
        } else {
            cholmod_state_->factor = cholmod_l_analyze(A, &cholmod_state_->common);
        }
        if (cholmod_state_->factor == nullptr) {
            cholmod_l_free_sparse(&A, &cholmod_state_->common);
            throw std::runtime_error("CHOLMOD analyze failed");
        }
        const int ok = cholmod_l_factorize(A, cholmod_state_->factor, &cholmod_state_->common);
        cholmod_l_free_sparse(&A, &cholmod_state_->common);
        if (!ok) {
            throw std::runtime_error("CHOLMOD factorize failed");
        }
        factor_size_ = H.rows();
        factorized_ = true;
        return;
    }
#endif

    const Eigen::MatrixXd Hd = sparse_to_dense(H);
    Eigen::LLT<Eigen::MatrixXd> llt(Hd);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("Dense LLT factorization failed");
    }
    L_dense_ = llt.matrixL();
    factor_size_ = H.rows();
    factorized_ = true;
    using_cholmod_ = false;
}

Eigen::VectorXd FactorManager::solve(const Eigen::VectorXd& rhs) const {
    if (!factorized_) {
        throw std::runtime_error("FactorManager::solve called before factorize");
    }
    if (rhs.size() != last_H_.rows()) {
        throw std::runtime_error("RHS size does not match factorized system");
    }
    if (factor_size_ != last_H_.rows()) {
        throw std::runtime_error("FactorManager::solve called while factor growth is pending");
    }
    if (use_sparse_expanding_backend_) {
        if (static_cast<int>(sparse_factor_perm_.size()) != rhs.size()) {
            throw std::runtime_error("Sparse expanding solve permutation is unavailable");
        }
        const Eigen::VectorXd rhs_factor = permute_vector_by_order(rhs, sparse_factor_perm_);
        const Eigen::VectorXd x_factor = sparse_expanding_state_->solve(rhs_factor);
        return unpermute_vector_by_order(x_factor, sparse_factor_perm_);
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (using_cholmod_) {
        auto* b = eigen_to_cholmod_dense(rhs, &cholmod_state_->common);
        cholmod_dense* x = cholmod_l_solve(CHOLMOD_A, cholmod_state_->factor, b, &cholmod_state_->common);
        cholmod_l_free_dense(&b, &cholmod_state_->common);
        if (x == nullptr) {
            throw std::runtime_error("CHOLMOD solve failed");
        }
        Eigen::VectorXd out = cholmod_dense_to_eigen(x);
        cholmod_l_free_dense(&x, &cholmod_state_->common);
        return out;
    }
#endif
    const Eigen::VectorXd y = L_dense_.triangularView<Eigen::Lower>().solve(rhs);
    return L_dense_.transpose().triangularView<Eigen::Upper>().solve(y);
}

Eigen::VectorXd FactorManager::solve_cached() const {
    if (last_g_.size() == 0) {
        throw std::runtime_error("solve_cached called without cached RHS");
    }
    return solve(last_g_);
}

std::vector<int> FactorManager::active_var_indices() const {
    std::vector<int> idx;
    idx.reserve(active_vars_.size());
    for (int i = 0; i < static_cast<int>(active_vars_.size()); ++i) {
        if (active_vars_[static_cast<size_t>(i)] != 0) idx.push_back(i);
    }
    return idx;
}

int FactorManager::active_state_size() const noexcept {
    int count = 0;
    for (unsigned char v : active_vars_) count += (v != 0);
    return count;
}

double FactorManager::maintained_factor_eta_full() const {
    if (!factorized_ || factor_size_ <= 0) {
        throw std::runtime_error("maintained_factor_eta_full called without a factorized system");
    }
    if (use_sparse_expanding_backend_) {
        return sparse_expanding_state_->eta_logdet_half();
    }
#ifdef ISLAM_HAS_CHOLMOD
    if (using_cholmod_ && cholmod_state_ && cholmod_state_->factor != nullptr) {
        const auto* F = cholmod_state_->factor;
        const auto* p = static_cast<const ss_long*>(F->p);
        const auto* i = static_cast<const ss_long*>(F->i);
        const auto* x = static_cast<const double*>(F->x);
        if (p == nullptr || i == nullptr || x == nullptr) {
            throw std::runtime_error("maintained_factor_eta_full: invalid CHOLMOD factor storage");
        }
        double eta = 0.0;
        for (int col = 0; col < factor_size_; ++col) {
            double diag = std::numeric_limits<double>::quiet_NaN();
            for (ss_long kk = p[col]; kk < p[col + 1]; ++kk) {
                if (static_cast<int>(i[kk]) == col) {
                    diag = x[kk];
                    break;
                }
            }
            if (!(std::isfinite(diag) && std::abs(diag) > 0.0)) {
                throw std::runtime_error("maintained_factor_eta_full: non-positive/missing CHOLMOD factor diagonal");
            }
            eta += std::log(std::abs(diag));
        }
        return eta;
    }
#endif
    if (L_dense_.rows() != factor_size_ || L_dense_.cols() != factor_size_) {
        throw std::runtime_error("maintained_factor_eta_full: dense factor is unavailable");
    }
    double eta = 0.0;
    for (int i = 0; i < factor_size_; ++i) {
        const double diag = L_dense_(i, i);
        if (!(std::isfinite(diag) && std::abs(diag) > 0.0)) {
            throw std::runtime_error("maintained_factor_eta_full: non-positive dense factor diagonal");
        }
        eta += std::log(std::abs(diag));
    }
    return eta;
}

double FactorManager::information_eta(bool active_only) const {
    const double full_eta = maintained_factor_eta_full();
    if (!active_only) return full_eta;

    // In the reserved-state incremental path, inactive future variables are
    // represented by independent latent priors. Their factor-diagonal contribution
    // is known exactly and should not influence the paper's current-state eta_t.
    const int active_n = active_state_size();
    const int inactive_n = std::max(0, factor_size_ - active_n);
    if (inactive_n == 0) return full_eta;
    if (!(latent_prior_strength_ > 0.0)) return full_eta;
    return full_eta - 0.5 * static_cast<double>(inactive_n) * std::log(latent_prior_strength_);
}

const DynamicCcolamdEngine::Stats& FactorManager::dynamic_ordering_stats() const noexcept {
    static const DynamicCcolamdEngine::Stats empty_stats{};
    return symbolic_engine_ ? symbolic_engine_->dynamic_ordering_stats() : empty_stats;
}

const SymbolicEngineStats& FactorManager::symbolic_engine_stats() const noexcept {
    static const SymbolicEngineStats empty_stats{};
    return symbolic_engine_ ? symbolic_engine_->stats() : empty_stats;
}

std::vector<int> FactorManager::get_permutation() const {
    if (last_H_.rows() == 0) {
        return {};
    }
    refresh_symbolic_cache_if_needed();
    return symbolic_engine_ ? symbolic_engine_->snapshot().perm : std::vector<int>{};
}

std::vector<int> FactorManager::get_elimination_tree() const {
    if (last_H_.rows() == 0) {
        return {};
    }
    refresh_symbolic_cache_if_needed();
    return symbolic_engine_ ? symbolic_engine_->snapshot().etree_parent : std::vector<int>{};
}

PermutedSystem FactorManager::build_permuted_system() const {
    if (last_H_.rows() == 0 || last_g_.size() == 0) {
        throw std::runtime_error("build_permuted_system called without cached system");
    }
    refresh_numeric_cache_if_needed();
    const auto& snap = symbolic_engine_->snapshot();
    return PermutedSystem{
        cached_H_perm_,
        cached_g_perm_,
        snap.perm,
        snap.pinv
    };
}

SymbolicSystem FactorManager::build_symbolic_system() const {
    if (last_H_.rows() == 0 || last_g_.size() == 0) {
        throw std::runtime_error("build_symbolic_system called without cached system");
    }
    refresh_numeric_cache_if_needed();
    const auto& snap = symbolic_engine_->snapshot();
    return SymbolicSystem{
        PermutedSystem{cached_H_perm_, cached_g_perm_, snap.perm, snap.pinv},
        snap.etree_parent
    };
}

FactorBlockSystem FactorManager::build_factor_block_system() const {
    FactorBlockSystem out;
    if (!factorized_ || factor_size_ <= 0 || last_g_.size() != factor_size_) {
        return out;
    }
    if (factor_size_ != state_size_) {
        return out;
    }

    if (use_sparse_expanding_backend_ && sparse_expanding_state_) {
        if (static_cast<int>(sparse_factor_perm_.size()) == factor_size_) {
            out.perm = sparse_factor_perm_;
        } else {
            out.perm = identity_permutation(factor_size_);
        }
        out.pinv = inverse_permutation(out.perm);
        out.g_factor = permute_vector_by_order(last_g_, out.perm);
        out.L_factor = sparse_expanding_state_->lower_factor();
        Eigen::SparseMatrix<double> H_factor = (out.L_factor * out.L_factor.transpose()).pruned();
        H_factor.makeCompressed();
        out.etree_parent = elimination_tree_from_upper(H_factor);
        out.available = true;
        out.from_cholmod = false;
        return out;
    }

#ifdef ISLAM_HAS_CHOLMOD
    if (using_cholmod_ && cholmod_state_ && cholmod_state_->factor != nullptr) {
        out.perm = suitesparse_extract_factor_perm(cholmod_state_->factor, factor_size_);
        if (out.perm.empty()) {
            out.perm = identity_permutation(factor_size_);
        }
        out.pinv = inverse_permutation(out.perm);
        out.g_factor = permute_vector_by_order(last_g_, out.perm);
        out.L_factor = cholmod_factor_to_eigen_lower(cholmod_state_->factor, factor_size_);
        Eigen::SparseMatrix<double> H_factor = (out.L_factor * out.L_factor.transpose()).pruned();
        H_factor.makeCompressed();
        out.etree_parent = elimination_tree_from_upper(H_factor);
        out.available = true;
        out.from_cholmod = true;
        return out;
    }
#endif

    if (L_dense_.rows() == factor_size_ && L_dense_.cols() == factor_size_) {
        out.perm.resize(static_cast<size_t>(factor_size_));
        for (int i = 0; i < factor_size_; ++i) out.perm[static_cast<size_t>(i)] = i;
        out.pinv = inverse_permutation(out.perm);
        out.g_factor = last_g_;
        out.L_factor = dense_lower_to_sparse(L_dense_);
        Eigen::SparseMatrix<double> H_factor = (out.L_factor * out.L_factor.transpose()).pruned();
        H_factor.makeCompressed();
        out.etree_parent = elimination_tree_from_upper(H_factor);
        out.available = true;
        out.from_cholmod = false;
    }
    return out;
}

Eigen::VectorXd FactorManager::solve_graph(
    const Graph& g,
    const std::vector<int>* edge_ids,
    double anchor_strength) {
    rebuild_from_graph(g, edge_ids, anchor_strength);
    return solve_cached();
}

} // namespace islam
