#include "islam/deterministic_batch_ccolamd.hpp"
#include "islam/dynamic_exact_ccolamd.hpp"
#include "islam/suitesparse_ccolamd_reference.hpp"

#include <Eigen/Sparse>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct CaseSpec {
    std::string name;
    int poses = 0;
    int block_dim = 3;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> priors;
};

void add_scalar_edge(std::set<std::pair<int, int>>& entries, int a, int b) {
    entries.insert({a, b});
    entries.insert({b, a});
}

void add_block_self(std::set<std::pair<int, int>>& entries, int node, int dim) {
    const int base = node * dim;
    for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < dim; ++c) add_scalar_edge(entries, base + r, base + c);
    }
}

void add_block_edge(std::set<std::pair<int, int>>& entries, int a, int b, int dim) {
    add_block_self(entries, a, dim);
    add_block_self(entries, b, dim);
    const int ba = a * dim;
    const int bb = b * dim;
    for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < dim; ++c) add_scalar_edge(entries, ba + r, bb + c);
    }
}

Eigen::SparseMatrix<double> make_pattern(const CaseSpec& spec) {
    std::set<std::pair<int, int>> entries;
    for (int i = 0; i < spec.poses; ++i) add_block_self(entries, i, spec.block_dim);
    for (auto [a, b] : spec.edges) add_block_edge(entries, a, b, spec.block_dim);
    for (int p : spec.priors) add_block_self(entries, p, spec.block_dim);
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(entries.size());
    for (auto [r, c] : entries) trips.emplace_back(r, c, 1.0);
    Eigen::SparseMatrix<double> H(spec.poses * spec.block_dim, spec.poses * spec.block_dim);
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();
    return H;
}

std::vector<CaseSpec> representative_cases() {
    std::vector<CaseSpec> cases;
    auto chain = [](int n) {
        std::vector<std::pair<int,int>> e;
        for (int i = 1; i < n; ++i) e.push_back({i - 1, i});
        return e;
    };
    {
        CaseSpec c{"chain_20", 20, 3, chain(20), {0}};
        cases.push_back(c);
    }
    {
        CaseSpec c{"chain_50_loop_closures", 50, 3, chain(50), {0}};
        for (int i = 10; i < 50; i += 10) c.edges.push_back({0, i});
        for (int i = 12; i < 45; i += 8) c.edges.push_back({i, std::min(49, i + 5)});
        cases.push_back(c);
    }
    {
        CaseSpec c{"dense_local_windows_60", 60, 3, chain(60), {0}};
        for (int i = 0; i + 5 < 60; i += 3) c.edges.push_back({i, i + 5});
        for (int i = 0; i + 12 < 60; i += 12) c.edges.push_back({i, i + 12});
        cases.push_back(c);
    }
    {
        CaseSpec c{"intermittent_priors_80", 80, 3, chain(80), {0}};
        for (int i = 10; i < 80; i += 10) c.priors.push_back(i);
        for (int i = 20; i < 80; i += 20) c.edges.push_back({std::max(0, i - 17), i});
        cases.push_back(c);
    }
    {
        CaseSpec c{"multi_loop_120", 120, 3, chain(120), {0}};
        for (int i = 20; i < 120; i += 20) c.edges.push_back({0, i});
        for (int i = 30; i + 30 < 120; i += 15) c.edges.push_back({i, i + 30});
        for (int i = 5; i + 7 < 120; i += 11) c.edges.push_back({i, i + 7});
        cases.push_back(c);
    }
    return cases;
}

std::vector<int> dirty_vars_for_edge(std::pair<int,int> e, int dim) {
    std::vector<int> out;
    for (int n : {e.first, e.second}) {
        for (int k = 0; k < dim; ++k) out.push_back(n * dim + k);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

int first_mismatch_position(const std::vector<int>& a, const std::vector<int>& b) {
    const int n = std::min(static_cast<int>(a.size()), static_cast<int>(b.size()));
    for (int i = 0; i < n; ++i) if (a[static_cast<size_t>(i)] != b[static_cast<size_t>(i)]) return i;
    if (a.size() != b.size()) return n;
    return -1;
}

std::string head_perm(const std::vector<int>& p, int count = 16) {
    std::ostringstream oss;
    oss << "[";
    for (int i = 0; i < std::min(count, static_cast<int>(p.size())); ++i) {
        if (i) oss << ' ';
        oss << p[static_cast<size_t>(i)];
    }
    if (static_cast<int>(p.size()) > count) oss << " ...";
    oss << "]";
    return oss.str();
}

struct Totals {
    int cases = 0;
    int unavailable = 0;
    int batch_failures = 0;
    int dynamic_failures = 0;
    long long batch_mismatch_positions = 0;
    long long dynamic_mismatch_positions = 0;
};

void audit_case(const CaseSpec& spec, bool verbose, Totals& totals) {
    ++totals.cases;
    const auto H = make_pattern(spec);
    islam::DeterministicBatchCcolamd owned;
    const auto owned_perm = owned.order(H);
    const auto ref = islam::suitesparse_ccolamd_reference_order(H);
    std::cout << "case=" << spec.name
              << " n=" << H.rows()
              << " nnz=" << H.nonZeros()
              << " sig=" << islam::sparse_pattern_audit_signature(H);

    if (!ref.available || !ref.ok) {
        ++totals.unavailable;
        std::cout << " suitesparse=unavailable message=\"" << ref.message << "\"\n";
        return;
    }

    const int mism = islam::permutation_mismatch_positions(owned_perm, ref.permutation);
    totals.batch_mismatch_positions += mism;
    if (mism != 0) ++totals.batch_failures;
    std::cout << " batch_mismatches=" << mism;
    if (mism != 0) {
        std::cout << " first=" << first_mismatch_position(owned_perm, ref.permutation);
    }

    islam::DynamicExactCcolamdPrototype dyn;
    std::vector<int> dirty;
    if (!spec.edges.empty()) dirty = dirty_vars_for_edge(spec.edges.back(), spec.block_dim);
    const auto dyn_perm = dyn.refresh(H, dirty);
    const int dyn_mism = islam::permutation_mismatch_positions(dyn_perm, ref.permutation);
    totals.dynamic_mismatch_positions += dyn_mism;
    if (dyn_mism != 0) ++totals.dynamic_failures;
    std::cout << " dynamic_mismatches=" << dyn_mism;
    if (dyn_mism != 0) {
        std::cout << " dynamic_first=" << first_mismatch_position(dyn_perm, ref.permutation);
    }
    std::cout << "\n";

    if (verbose && (mism != 0 || dyn_mism != 0)) {
        std::cout << "  owned_head=" << head_perm(owned_perm) << "\n";
        std::cout << "  dyn_head  =" << head_perm(dyn_perm) << "\n";
        std::cout << "  ss_head   =" << head_perm(ref.permutation) << "\n";
    }
}

} // namespace

int main(int argc, char** argv) {
    bool strict = false;
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--strict") strict = true;
        else if (arg == "--verbose") verbose = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: islam_suitesparse_compat_audit [--strict] [--verbose]\n";
            return 0;
        }
    }

    Totals totals;
    for (const auto& c : representative_cases()) audit_case(c, verbose, totals);

    std::cout << "summary cases=" << totals.cases
              << " unavailable=" << totals.unavailable
              << " batch_failures=" << totals.batch_failures
              << " dynamic_failures=" << totals.dynamic_failures
              << " batch_mismatch_positions=" << totals.batch_mismatch_positions
              << " dynamic_mismatch_positions=" << totals.dynamic_mismatch_positions
              << "\n";

    if (strict) {
        if (totals.unavailable != 0) return 3;
        if (totals.batch_failures != 0 || totals.dynamic_failures != 0) return 2;
    }
    return 0;
}
