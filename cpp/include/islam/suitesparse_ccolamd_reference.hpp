#pragma once

#include <Eigen/Sparse>

#include <cstdint>
#include <string>
#include <vector>

namespace islam {

struct SuiteSparseCcolamdReferenceResult {
    bool available = false;
    bool ok = false;
    std::string backend;
    std::string message;
    std::vector<int> permutation;
    std::vector<long long> stats;
};

// True only when the project was compiled with the SuiteSparse CCOLAMD header
// and library.  This is deliberately narrower than ISLAM_HAS_CHOLMOD: CHOLMOD
// may produce a fill-reducing ordering internally, but a byte-for-byte CCOLAMD
// audit must call the CCOLAMD API directly.
[[nodiscard]] bool suitesparse_ccolamd_reference_available();

// Direct SuiteSparse CCOLAMD reference ordering for a sparse column pattern.
// The input is interpreted as the matrix whose columns are to be ordered.  For
// the current SLAM symbolic path this is normally the scalar normal-equation
// pattern H = J^T J used by the in-repo deterministic batch rule.  `cmember`
// may be empty, in which case all columns belong to constraint group zero.
[[nodiscard]] SuiteSparseCcolamdReferenceResult suitesparse_ccolamd_reference_order(
    const Eigen::SparseMatrix<double>& pattern,
    const std::vector<int>& cmember = {});

[[nodiscard]] bool is_valid_permutation(const std::vector<int>& perm, int n);
[[nodiscard]] int permutation_mismatch_positions(const std::vector<int>& a,
                                                 const std::vector<int>& b);

// Compact structural signature for audit reports.  This intentionally mirrors
// deterministic_ccolamd_pattern_signature without forcing callers to include the
// dynamic-exact headers.
[[nodiscard]] std::uint64_t sparse_pattern_audit_signature(const Eigen::SparseMatrix<double>& pattern);

} // namespace islam
