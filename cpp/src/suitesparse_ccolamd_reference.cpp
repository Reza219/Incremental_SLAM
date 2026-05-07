#include "islam/suitesparse_ccolamd_reference.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <sstream>
#include <stdexcept>

#ifdef ISLAM_HAS_CCOLAMD
#include <ccolamd.h>
#endif

namespace islam {
namespace {

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30u)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27u)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31u);
}

std::vector<int> canonical_cmember(int n, const std::vector<int>& cmember) {
    if (cmember.empty()) return std::vector<int>(static_cast<size_t>(n), 0);
    if (static_cast<int>(cmember.size()) != n) {
        throw std::runtime_error("SuiteSparse CCOLAMD reference: cmember size must be zero or equal to n_col");
    }
    std::vector<int> out = cmember;
    for (int& v : out) v = std::max(0, v);
    return out;
}

Eigen::SparseMatrix<double> compressed_binary_pattern(const Eigen::SparseMatrix<double>& in) {
    Eigen::SparseMatrix<double> A = in;
    A.makeCompressed();
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(static_cast<size_t>(A.nonZeros()));
    for (int col = 0; col < A.outerSize(); ++col) {
        std::set<int> rows;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            if (it.row() >= 0 && it.row() < A.rows()) rows.insert(it.row());
        }
        for (int r : rows) trips.emplace_back(r, col, 1.0);
    }
    Eigen::SparseMatrix<double> B(A.rows(), A.cols());
    B.setFromTriplets(trips.begin(), trips.end());
    B.makeCompressed();
    return B;
}

} // namespace

bool is_valid_permutation(const std::vector<int>& perm, int n) {
    if (static_cast<int>(perm.size()) != n) return false;
    std::vector<unsigned char> seen(static_cast<size_t>(n), 0);
    for (int v : perm) {
        if (v < 0 || v >= n || seen[static_cast<size_t>(v)]) return false;
        seen[static_cast<size_t>(v)] = 1;
    }
    return true;
}

int permutation_mismatch_positions(const std::vector<int>& a, const std::vector<int>& b) {
    const int n = std::min(static_cast<int>(a.size()), static_cast<int>(b.size()));
    int mismatches = std::abs(static_cast<int>(a.size()) - static_cast<int>(b.size()));
    for (int i = 0; i < n; ++i) {
        if (a[static_cast<size_t>(i)] != b[static_cast<size_t>(i)]) ++mismatches;
    }
    return mismatches;
}

std::uint64_t sparse_pattern_audit_signature(const Eigen::SparseMatrix<double>& pattern) {
    Eigen::SparseMatrix<double> A = compressed_binary_pattern(pattern);
    std::uint64_t h = mix64(static_cast<std::uint64_t>(static_cast<std::uint32_t>(A.rows())));
    h ^= mix64(static_cast<std::uint64_t>(static_cast<std::uint32_t>(A.cols())) << 1u);
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const std::uint64_t x = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(it.row())) << 32u)
                                | static_cast<std::uint32_t>(col);
            h ^= mix64(x);
        }
    }
    return h;
}

bool suitesparse_ccolamd_reference_available() {
#ifdef ISLAM_HAS_CCOLAMD
    return true;
#else
    return false;
#endif
}

SuiteSparseCcolamdReferenceResult suitesparse_ccolamd_reference_order(
    const Eigen::SparseMatrix<double>& pattern,
    const std::vector<int>& cmember_in) {
    SuiteSparseCcolamdReferenceResult result;
#ifdef ISLAM_HAS_CCOLAMD
    result.available = true;
    result.backend = "SuiteSparse CCOLAMD direct API";
    if (pattern.rows() < 0 || pattern.cols() < 0) {
        result.message = "invalid matrix dimensions";
        return result;
    }
    const int n_row_i = pattern.rows();
    const int n_col_i = pattern.cols();
    if (n_col_i == 0) {
        result.ok = true;
        return result;
    }

    Eigen::SparseMatrix<double> A_eigen = compressed_binary_pattern(pattern);
    const auto n_row = static_cast<SuiteSparse_long>(n_row_i);
    const auto n_col = static_cast<SuiteSparse_long>(n_col_i);
    const auto nnz = static_cast<SuiteSparse_long>(A_eigen.nonZeros());
    const SuiteSparse_long recommended = ccolamd_l_recommended(nnz, n_row, n_col);
    if (recommended <= nnz || recommended <= 0) {
        result.message = "ccolamd_l_recommended returned an invalid workspace length";
        return result;
    }

    std::vector<SuiteSparse_long> A(static_cast<size_t>(recommended), 0);
    std::vector<SuiteSparse_long> p(static_cast<size_t>(n_col_i + 1), 0);
    for (int col = 0; col <= n_col_i; ++col) {
        p[static_cast<size_t>(col)] = static_cast<SuiteSparse_long>(A_eigen.outerIndexPtr()[col]);
    }
    for (int nz = 0; nz < A_eigen.nonZeros(); ++nz) {
        A[static_cast<size_t>(nz)] = static_cast<SuiteSparse_long>(A_eigen.innerIndexPtr()[nz]);
    }

    double knobs[CCOLAMD_KNOBS];
    ccolamd_l_set_defaults(knobs);
    SuiteSparse_long stats[CCOLAMD_STATS];
    std::fill(stats, stats + CCOLAMD_STATS, static_cast<SuiteSparse_long>(0));

    std::vector<int> cmember_canon = canonical_cmember(n_col_i, cmember_in);
    std::vector<SuiteSparse_long> cmember(static_cast<size_t>(n_col_i), 0);
    for (int c = 0; c < n_col_i; ++c) {
        cmember[static_cast<size_t>(c)] = static_cast<SuiteSparse_long>(cmember_canon[static_cast<size_t>(c)]);
    }

    const int ok = ccolamd_l(
        n_row,
        n_col,
        recommended,
        A.data(),
        p.data(),
        knobs,
        stats,
        cmember.empty() ? nullptr : cmember.data());

    result.stats.reserve(CCOLAMD_STATS);
    for (SuiteSparse_long s : stats) result.stats.push_back(static_cast<long long>(s));
    if (!ok) {
        std::ostringstream oss;
        oss << "ccolamd_l failed";
#ifdef CCOLAMD_STATUS
        oss << " with status " << static_cast<long long>(stats[CCOLAMD_STATUS]);
#endif
        result.message = oss.str();
        return result;
    }

    result.permutation.assign(static_cast<size_t>(n_col_i), -1);
    for (int k = 0; k < n_col_i; ++k) {
        result.permutation[static_cast<size_t>(k)] = static_cast<int>(p[static_cast<size_t>(k)]);
    }
    result.ok = is_valid_permutation(result.permutation, n_col_i);
    result.message = result.ok ? "ok" : "ccolamd_l returned an invalid permutation";
    return result;
#else
    (void)pattern;
    (void)cmember_in;
    result.available = false;
    result.ok = false;
    result.backend = "unavailable";
    result.message = "project was not compiled with ISLAM_HAS_CCOLAMD=1";
    return result;
#endif
}

} // namespace islam
