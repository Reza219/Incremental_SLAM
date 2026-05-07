#include "islam/symbolic_engine.hpp"

#include <Eigen/OrderingMethods>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>

namespace islam {
namespace {

bool strict_symbolic_runtime_enabled() {
    const char* v = std::getenv("ISLAM_STRICT_SYMBOLIC");
    if (v == nullptr) return false;
    const std::string s(v);
    return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
}


bool use_dynamic_exact_ccolamd_backend() {
    const char* v = std::getenv("ISLAM_USE_DYNAMIC_EXACT_CCOLAMD_BACKEND");
    if (v == nullptr) return false;
    const std::string ss(v);
    return ss == "1" || ss == "true" || ss == "TRUE" || ss == "on" || ss == "ON";
}

bool no_exact_etree_enabled() {
    const char* v = std::getenv("ISLAM_ETREE_NO_EXACT");
    if (v == nullptr) return false;
    const std::string s(v);
    return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
}

std::uint64_t sparse_pattern_signature(const Eigen::SparseMatrix<double>& pattern) {
    const std::uint64_t kOffset = 1469598103934665603ull;
    const std::uint64_t kPrime = 1099511628211ull;
    std::uint64_t h = kOffset;
    auto mix = [&](std::uint64_t x) { h ^= x; h *= kPrime; };
    mix(static_cast<std::uint64_t>(pattern.rows()));
    mix(static_cast<std::uint64_t>(pattern.cols()));
    for (int col = 0; col < pattern.outerSize(); ++col) {
        mix(static_cast<std::uint64_t>(col + 1));
        for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) {
            mix(static_cast<std::uint64_t>(it.row() + 1));
            mix(static_cast<std::uint64_t>(col + 1));
        }
    }
    return h;
}

std::vector<int> parse_int_list_semicolon(const std::string& text) {
    std::vector<int> out;
    std::stringstream ss(text);
    std::string tok;
    while (std::getline(ss, tok, ';')) {
        if (!tok.empty()) out.push_back(std::stoi(tok));
    }
    return out;
}

std::string int_list_to_semicolon(const std::vector<int>& xs) {
    std::ostringstream os;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (i) os << ';';
        os << xs[i];
    }
    return os.str();
}

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

} // namespace

void SymbolicEngine::clear() {
    if (dynamic_ordering_) dynamic_ordering_->clear();
    if (dynamic_exact_ordering_) dynamic_exact_ordering_->clear();
    state_size_ = 0;
    ++structure_revision_;
    cached_structure_revision_ = 0;
    dirty_mask_.clear();
    dirty_vars_.clear();
    force_full_refresh_ = true;
    pattern_counts_.clear();
    etree_exact_cache_env_loaded_ = false;
    etree_exact_cache_.clear();
    snapshot_ = SymbolicSnapshot{};
    stats_ = SymbolicEngineStats{};
}

const DynamicCcolamdEngine::Stats& SymbolicEngine::dynamic_ordering_stats() const noexcept {
    static const DynamicCcolamdEngine::Stats empty_stats{};
    return dynamic_ordering_ ? dynamic_ordering_->stats() : empty_stats;
}

const DynamicExactCcolamdPrototype::Stats& SymbolicEngine::dynamic_exact_ordering_stats() const noexcept {
    static const DynamicExactCcolamdPrototype::Stats empty_stats{};
    return dynamic_exact_ordering_ ? dynamic_exact_ordering_->stats() : empty_stats;
}

void SymbolicEngine::reserve_state(int n) {
    if (!dynamic_ordering_) dynamic_ordering_ = std::make_unique<DynamicCcolamdEngine>();
    if (n < 0) {
        throw std::runtime_error("Negative state size in SymbolicEngine::reserve_state");
    }
    clear();
    state_size_ = n;
    dirty_mask_.assign(static_cast<size_t>(n), 0);
    dynamic_ordering_->reserve_state(n);
    if (dynamic_exact_ordering_) dynamic_exact_ordering_->clear();
    for (int i = 0; i < n; ++i) {
        pattern_counts_[pattern_key(i, i)] = 1;
    }
}

void SymbolicEngine::rebuild_from_numeric_matrix(const Eigen::SparseMatrix<double>& H) {
    reserve_state(H.rows());
    for (int col = 0; col < H.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H, col); it; ++it) {
            if (it.row() > it.col()) continue;
            const auto key = pattern_key(it.row(), it.col());
            pattern_counts_[key] = std::max(pattern_counts_[key], 1);
        }
    }
    note_full_refresh();
}

void SymbolicEngine::apply_contribution_pattern(const Eigen::SparseMatrix<double>& H_pattern,
                                                const std::vector<int>& touched_vars,
                                                int delta) {
    if (delta == 0) return;
    std::vector<int> changed;
    changed.reserve(static_cast<size_t>(H_pattern.nonZeros()) + touched_vars.size());
    for (int col = 0; col < H_pattern.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H_pattern, col); it; ++it) {
            if (it.row() > it.col()) continue;
            const auto key = pattern_key(it.row(), it.col());
            int next = 0;
            const auto found = pattern_counts_.find(key);
            if (found != pattern_counts_.end()) next = found->second;
            next += delta;
            if (next <= 0) {
                if (it.row() != it.col()) {
                    pattern_counts_.erase(key);
                } else {
                    pattern_counts_[key] = 1;
                }
            } else {
                pattern_counts_[key] = next;
            }
            changed.push_back(it.row());
            changed.push_back(it.col());
        }
    }
    changed.insert(changed.end(), touched_vars.begin(), touched_vars.end());
    note_structural_change(unique_sorted(std::move(changed)));
}

void SymbolicEngine::note_full_refresh() {
    ++structure_revision_;
    cached_structure_revision_ = 0;
    force_full_refresh_ = true;
    dirty_mask_.assign(static_cast<size_t>(state_size_), 0);
    dirty_vars_.clear();
}

void SymbolicEngine::refresh_if_needed(const OrderingOracle& ordering_oracle) {
    if (cached_structure_revision_ == structure_revision_) {
        return;
    }

    if (state_size_ <= 0) {
        snapshot_ = SymbolicSnapshot{};
        cached_structure_revision_ = structure_revision_;
        force_full_refresh_ = false;
        return;
    }

    const Eigen::SparseMatrix<double> pattern = build_pattern_matrix();
    const auto dirty_before = dirty_vars_;
    const bool can_do_local = !force_full_refresh_
        && static_cast<int>(snapshot_.perm.size()) == state_size_
        && !dirty_before.empty();

    if (use_dynamic_exact_ccolamd_backend()) {
        if (!dynamic_exact_ordering_) dynamic_exact_ordering_ = std::make_unique<DynamicExactCcolamdPrototype>();
        snapshot_.perm = dynamic_exact_ordering_->refresh(pattern, dirty_before);
    } else {
        if (!dynamic_ordering_) dynamic_ordering_ = std::make_unique<DynamicCcolamdEngine>();
        snapshot_.perm = dynamic_ordering_->refresh_exact(pattern, dirty_before, ordering_oracle ? ordering_oracle : OrderingOracle(amd_permutation));
    }
    snapshot_.pinv = inverse_permutation(snapshot_.perm);
    snapshot_.pattern_perm = symmetric_permute_by_order(pattern, snapshot_.perm);
    const auto etree_cache_key = sparse_pattern_signature(snapshot_.pattern_perm);
    if (!etree_exact_cache_env_loaded_) {
        etree_exact_cache_env_loaded_ = true;
        const char* in_path = std::getenv("ISLAM_ETREE_EXACT_CACHE_IN");
        if (in_path != nullptr && std::string(in_path).size() > 0) {
            std::ifstream ifs(in_path);
            if (!ifs) {
                if (no_exact_etree_enabled()) {
                    throw std::runtime_error("ISLAM_ETREE_NO_EXACT=1: failed to open ISLAM_ETREE_EXACT_CACHE_IN");
                }
            } else {
                std::string line;
                while (std::getline(ifs, line)) {
                    if (line.empty() || line[0] == '#') continue;
                    const auto comma = line.find(',');
                    if (comma == std::string::npos) continue;
                    const auto key = static_cast<std::uint64_t>(std::stoull(line.substr(0, comma)));
                    const auto etree = parse_int_list_semicolon(line.substr(comma + 1));
                    if (!etree.empty() && etree_exact_cache_.find(key) == etree_exact_cache_.end()) {
                        etree_exact_cache_.emplace(key, etree);
                        ++stats_.etree_exact_cache_imported_entries;
                    }
                }
            }
        }
    }

    std::vector<int> exact_etree;
    const auto cached_etree = etree_exact_cache_.find(etree_cache_key);
    if (cached_etree != etree_exact_cache_.end()) {
        exact_etree = cached_etree->second;
        ++stats_.etree_exact_cache_hits;
        if (no_exact_etree_enabled()) ++stats_.etree_no_exact_cache_hits;
    } else {
        ++stats_.etree_exact_cache_misses;
        if (no_exact_etree_enabled()) {
            ++stats_.etree_no_exact_cache_misses;
            throw std::runtime_error("ISLAM_ETREE_NO_EXACT=1: exact etree cache miss; run once with ISLAM_ETREE_EXACT_CACHE_OUT to precompute the cache");
        }
        exact_etree = elimination_tree_from_upper(snapshot_.pattern_perm);
        ++stats_.etree_exact_recomputes;
        etree_exact_cache_[etree_cache_key] = exact_etree;
        const char* out_path = std::getenv("ISLAM_ETREE_EXACT_CACHE_OUT");
        if (out_path != nullptr && std::string(out_path).size() > 0) {
            std::ofstream ofs(out_path, std::ios::app);
            if (ofs) {
                ofs << etree_cache_key << ',' << int_list_to_semicolon(exact_etree) << '\n';
                ++stats_.etree_exact_cache_exported_entries;
            }
        }
    }

    const auto dirty_positions = dirty_positions_from_vars(snapshot_.perm, dirty_before);
    if (can_do_local && static_cast<int>(snapshot_.etree_parent.size()) == state_size_) {
        ++stats_.etree_local_update_attempts;
        auto local_etree = incremental_local_etree_update(snapshot_.pattern_perm,
                                                          snapshot_.etree_parent,
                                                          dirty_positions);
        if (local_etree == exact_etree) {
            ++stats_.etree_local_update_accepts;
            snapshot_.etree_parent = std::move(local_etree);
        } else {
            ++stats_.etree_local_update_fallbacks;
            if (strict_symbolic_runtime_enabled()) {
                throw std::runtime_error("Strict symbolic mode failed: local elimination-tree update did not certify against the exact etree");
            }
            snapshot_.etree_parent = exact_etree;
        }
    } else {
        snapshot_.etree_parent = exact_etree;
    }

    cached_structure_revision_ = structure_revision_;
    dirty_mask_.assign(static_cast<size_t>(state_size_), 0);
    dirty_vars_.clear();
    force_full_refresh_ = false;
}

std::uint64_t SymbolicEngine::pattern_key(int i, int j) {
    if (i > j) std::swap(i, j);
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 32u)
         | static_cast<std::uint32_t>(j);
}

std::pair<int, int> SymbolicEngine::decode_pattern_key(std::uint64_t key) {
    const int i = static_cast<int>(static_cast<std::uint32_t>(key >> 32u));
    const int j = static_cast<int>(static_cast<std::uint32_t>(key & 0xffffffffu));
    return {i, j};
}

std::vector<int> SymbolicEngine::inverse_permutation(const std::vector<int>& perm) {
    std::vector<int> pinv(perm.size(), -1);
    for (int k = 0; k < static_cast<int>(perm.size()); ++k) {
        const int pk = perm[static_cast<size_t>(k)];
        if (pk < 0 || pk >= static_cast<int>(perm.size())) {
            throw std::runtime_error("Invalid permutation entry in SymbolicEngine");
        }
        pinv[static_cast<size_t>(pk)] = k;
    }
    return pinv;
}

std::vector<int> SymbolicEngine::unique_sorted(std::vector<int> vals) {
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    return vals;
}

Eigen::SparseMatrix<double> SymbolicEngine::principal_submatrix(const Eigen::SparseMatrix<double>& A,
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
            if (new_row >= 0) trips.emplace_back(new_row, new_col, it.value());
        }
    }

    Eigen::SparseMatrix<double> out(static_cast<int>(idx.size()), static_cast<int>(idx.size()));
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

Eigen::SparseMatrix<double> SymbolicEngine::symmetric_permute_by_order(const Eigen::SparseMatrix<double>& A,
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

std::vector<int> SymbolicEngine::elimination_tree_from_upper(const Eigen::SparseMatrix<double>& Hperm) {
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

std::vector<int> SymbolicEngine::expand_pattern_neighborhood(const Eigen::SparseMatrix<double>& H,
                                                             const std::vector<int>& seeds,
                                                             int hops) {
    const int n = H.rows();
    if (n == 0 || seeds.empty()) return {};

    std::vector<unsigned char> in_set(static_cast<size_t>(n), 0);
    for (int v : seeds) {
        if (v >= 0 && v < n) in_set[static_cast<size_t>(v)] = 1;
    }

    for (int hop = 0; hop < hops; ++hop) {
        std::vector<int> next;
        for (int col = 0; col < H.outerSize(); ++col) {
            bool column_touched = in_set[static_cast<size_t>(col)] != 0;
            if (!column_touched) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(H, col); it; ++it) {
                    if (in_set[static_cast<size_t>(it.row())]) {
                        column_touched = true;
                        break;
                    }
                }
            }
            if (!column_touched) continue;
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
    }

    std::vector<int> out;
    for (int i = 0; i < n; ++i) if (in_set[static_cast<size_t>(i)]) out.push_back(i);
    return out;
}

std::vector<int> SymbolicEngine::incremental_local_amd_permutation(const Eigen::SparseMatrix<double>& H,
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

std::vector<int> SymbolicEngine::dirty_positions_from_vars(const std::vector<int>& perm,
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

std::vector<int> SymbolicEngine::elimination_tree_from_upper_suffix(const Eigen::SparseMatrix<double>& Hperm,
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

std::vector<int> SymbolicEngine::incremental_local_etree_update(const Eigen::SparseMatrix<double>& Hperm,
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

Eigen::SparseMatrix<double> SymbolicEngine::build_pattern_matrix() const {
    Eigen::SparseMatrix<double> A(state_size_, state_size_);
    if (state_size_ <= 0) return A;
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(pattern_counts_.size() * 2 + static_cast<size_t>(state_size_));
    for (const auto& kv : pattern_counts_) {
        if (kv.second <= 0) continue;
        const auto ij = decode_pattern_key(kv.first);
        trips.emplace_back(ij.first, ij.second, 1.0);
        if (ij.first != ij.second) {
            trips.emplace_back(ij.second, ij.first, 1.0);
        }
    }
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    return A;
}

void SymbolicEngine::note_structural_change(const std::vector<int>& vars) {
    ++structure_revision_;
    cached_structure_revision_ = 0;
    if (force_full_refresh_) {
        return;
    }
    if (static_cast<int>(dirty_mask_.size()) != state_size_) {
        dirty_mask_.assign(static_cast<size_t>(state_size_), 0);
        dirty_vars_.clear();
    }
    for (int v : vars) {
        if (v < 0 || v >= state_size_) continue;
        if (!dirty_mask_[static_cast<size_t>(v)]) {
            dirty_mask_[static_cast<size_t>(v)] = 1;
            dirty_vars_.push_back(v);
        }
    }
}

} // namespace islam
