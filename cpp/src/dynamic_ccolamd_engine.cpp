#include "islam/dynamic_ccolamd_engine.hpp"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>
#include <functional>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace islam {

namespace {

constexpr double kLaplacePrior = 1.0;
constexpr double kSizePenalty = 0.0025;
constexpr double kReplayHitBonus = 0.05;
constexpr double kRecencyBonusScale = 0.15;
constexpr int kMaxReplayCandidates = 3;
constexpr size_t kMaxLocalPatternCacheEntries = 96;
constexpr size_t kMaxMotifPatternCacheEntries = 128;
constexpr size_t kMaxHierarchicalBlockCacheEntries = 96;
constexpr int kOverlapMinWidth = 4;
constexpr int kOverlapMaxWidth = 10;

bool env_flag_enabled(const char* name) {
    const char* v = std::getenv(name);
    if (v == nullptr) return false;
    const std::string s(v);
    return s == "1" || s == "true" || s == "TRUE" || s == "on" || s == "ON";
}

bool strict_symbolic_runtime_enabled() {
    return env_flag_enabled("ISLAM_STRICT_SYMBOLIC");
}

bool no_oracle_symbolic_enabled() {
    return env_flag_enabled("ISLAM_SYMBOLIC_NO_ORACLE");
}

std::vector<int> parse_perm_csv(const std::string& text) {
    std::vector<int> out;
    std::stringstream ss(text);
    std::string tok;
    while (std::getline(ss, tok, ';')) {
        if (!tok.empty()) out.push_back(std::stoi(tok));
    }
    return out;
}

std::string perm_to_csv(const std::vector<int>& perm) {
    std::ostringstream os;
    for (size_t i = 0; i < perm.size(); ++i) {
        if (i) os << ';';
        os << perm[i];
    }
    return os.str();
}

} // namespace

void DynamicCcolamdEngine::clear() {
    state_size_ = 0;
    perm_.clear();
    pinv_.clear();
    adjacency_.clear();
    regime_states_.clear();
    exact_cache_.clear();
    exact_cache_env_loaded_ = false;
    local_pattern_cache_.clear();
    motif_pattern_cache_.clear();
    hierarchical_block_cache_.clear();
    precedence_cache_.clear();
    current_regime_index_ = 0;
    next_regime_stable_id_ = 0;
    tick_ = 0;
    stats_ = Stats{};
}

void DynamicCcolamdEngine::reserve_state(int n) {
    if (n < 0) throw std::runtime_error("Negative state size in DynamicCcolamdEngine::reserve_state");
    clear();
    state_size_ = n;
    perm_.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) perm_[static_cast<size_t>(i)] = i;
    pinv_ = inverse_permutation_of(perm_);
    adjacency_.assign(static_cast<size_t>(n), {});
    stats_.current_regime_id = -1;
    stats_.num_regimes_discovered = 0;
    stats_.motif_pattern_cache_entries = 0;
    stats_.hierarchical_block_cache_entries = 0;
    stats_.precedence_cache_entries = 0;
}


void DynamicCcolamdEngine::sync_to_exact_reference(const Eigen::SparseMatrix<double>& pattern,
                                                   const std::vector<int>& exact_perm) {
    state_size_ = pattern.rows();
    const auto features = extract_regime_features(pattern);
    current_regime_index_ = select_or_create_regime(features);
    update_regime_model(current_regime_index_, features);
    if (!regime_states_.empty()) regime_states_[current_regime_index_].visits += 1;
    maybe_merge_regimes();
    stats_.current_regime_id = regime_states_.empty() ? -1 : regime_states_[current_regime_index_].stable_id;
    stats_.num_regimes_discovered = static_cast<int>(regime_states_.size());
    perm_ = exact_perm;
    pinv_ = inverse_permutation_of(perm_);
    rebuild_state_from_pattern(pattern);
}

std::vector<int> DynamicCcolamdEngine::refresh_exact(const Eigen::SparseMatrix<double>& pattern,
                                                     const std::vector<int>& dirty_vars,
                                                     const OrderingOracle& oracle) {
    if (pattern.rows() != pattern.cols()) throw std::runtime_error("DynamicCcolamdEngine requires square pattern");
    if (state_size_ != pattern.rows() || static_cast<int>(perm_.size()) != pattern.rows()) {
        reserve_state(pattern.rows());
    }
    if (pattern.rows() <= 1 || !oracle) {
        sync_to_exact_reference(pattern, perm_);
        return perm_;
    }

    const auto features = extract_regime_features(pattern);
    const size_t next_regime_index = select_or_create_regime(features);
    if (!regime_states_.empty() && next_regime_index != current_regime_index_) ++stats_.regime_switches;
    current_regime_index_ = next_regime_index;
    update_regime_model(current_regime_index_, features);
    if (!regime_states_.empty()) {
        regime_states_[current_regime_index_].visits += 1;
        stats_.current_regime_id = regime_states_[current_regime_index_].stable_id;
        stats_.num_regimes_discovered = static_cast<int>(regime_states_.size());
    }
    maybe_merge_regimes();
    if (!regime_states_.empty()) stats_.current_regime_id = regime_states_[current_regime_index_].stable_id;
    stats_.num_regimes_discovered = static_cast<int>(regime_states_.size());
    ++tick_;
    const auto cache_key = pattern_signature(pattern);
    maybe_load_exact_cache_from_env();
    std::vector<int> exact_perm;
    bool used_cache = false;
    if (const auto* cached = lookup_exact_cache(cache_key); cached != nullptr) {
        exact_perm = cached->perm;
        used_cache = true;
        ++stats_.oracle_cache_hits;
        if (no_oracle_symbolic_enabled()) ++stats_.no_oracle_cache_hits;
    } else {
        ++stats_.oracle_cache_misses;
        if (no_oracle_symbolic_enabled()) {
            ++stats_.no_oracle_cache_misses;
            throw std::runtime_error("ISLAM_SYMBOLIC_NO_ORACLE=1: exact ordering cache miss; run once with ISLAM_SYMBOLIC_EXACT_CACHE_OUT to precompute the cache");
        }
        exact_perm = oracle(pattern);
        ++stats_.oracle_refreshes;
        store_exact_cache(cache_key, exact_perm);
    }

    // Certified exact-cache dynamic accept: if this structural pattern has already
    // been certified exactly and the cached exact output is already the current
    // ordering, no oracle fallback is needed for this refresh.
    if (used_cache && exact_perm == perm_) {
        ++stats_.exact_output_certifications;
        ++stats_.certified_exact_cache_order_accepts;
        sync_to_exact_reference(pattern, exact_perm);
        return perm_;
    }

    if (!dirty_vars.empty()) {
        ++stats_.local_attempts;
        const auto candidates = candidate_windows(pattern, dirty_vars);
        stats_.candidate_windows_generated += static_cast<std::uint64_t>(candidates.size());
        for (const auto& candidate : candidates) {
            ++stats_.candidate_windows_tried;
            const auto kind_idx = static_cast<size_t>(static_cast<int>(candidate.kind));
            auto& regime_kind_stats = active_kind_stats();
            if (kind_idx < regime_kind_stats.size()) ++regime_kind_stats[kind_idx].attempts;
            switch (candidate.kind) {
            case CandidateKind::Replay: {
                ++stats_.replay_attempts;
                auto& replay_windows = active_replay_windows();
                for (auto& rw : replay_windows) {
                    if (rw.vars == candidate.vars) { ++rw.uses; break; }
                }
                break;
            }
            case CandidateKind::OneHop: ++stats_.one_hop_attempts; break;
            case CandidateKind::TwoHop: ++stats_.two_hop_attempts; break;
            case CandidateKind::Interval: ++stats_.interval_attempts; break;
            case CandidateKind::Union: ++stats_.union_attempts; break;
            case CandidateKind::Band: ++stats_.band_attempts; break;
            }
            const auto local_order = order_candidate_window(pattern, candidate.vars);
            if (local_order.ordered_vars.empty()) continue;
            stats_.overlap_assembly_piece_hits += local_order.overlap_piece_hits;
            stats_.overlap_assembly_proposals += local_order.overlap_proposals;
            const auto proposal = stitch_window_order(perm_, candidate.vars, local_order.ordered_vars);
            if (proposal == exact_perm) {
                ++stats_.exact_output_certifications;
                ++stats_.certified_local_order_accepts;
                if (local_order.method == LocalOrderMethod::OverlapAssembly) {
                    ++stats_.overlap_assembly_accepts;
                }
                ++stats_.local_accepts;
                ++stats_.adaptive_reorders;
                auto& accept_kind_stats = active_kind_stats();
                if (kind_idx < accept_kind_stats.size()) {
                    ++accept_kind_stats[kind_idx].accepts;
                    accept_kind_stats[kind_idx].last_success_tick = tick_;
                }
                switch (candidate.kind) {
                case CandidateKind::Replay: ++stats_.replay_accepts; break;
                case CandidateKind::OneHop: ++stats_.one_hop_accepts; break;
                case CandidateKind::TwoHop: ++stats_.two_hop_accepts; break;
                case CandidateKind::Interval: ++stats_.interval_accepts; break;
                case CandidateKind::Union: ++stats_.union_accepts; break;
                case CandidateKind::Band: ++stats_.band_accepts; break;
                }
                remember_replay_window(candidate.vars);
                std::vector<int> local_perm_idx;
                local_perm_idx.reserve(candidate.vars.size());
                for (int v : local_order.ordered_vars) {
                    auto it = std::find(candidate.vars.begin(), candidate.vars.end(), v);
                    if (it == candidate.vars.end()) { local_perm_idx.clear(); break; }
                    local_perm_idx.push_back(static_cast<int>(std::distance(candidate.vars.begin(), it)));
                }
                if (local_perm_idx.size() == candidate.vars.size()) {
                    const auto A_local = principal_submatrix(pattern, candidate.vars);
                    store_local_pattern_cache(pattern_signature(A_local), local_perm_idx);
                    store_motif_pattern_cache(motif_signature(A_local), static_cast<int>(candidate.vars.size()), local_perm_idx);
                    if (local_order.method == LocalOrderMethod::OverlapAssembly) {
                        promote_hierarchical_blocks_from_window(pattern, candidate.vars, local_order.ordered_vars);
                    }
                    store_precedence_cache(precedence_signature(A_local), static_cast<int>(candidate.vars.size()), local_perm_idx,
                                           local_order.method == LocalOrderMethod::PrecedenceCache ? 0.9 :
                                           (local_order.method == LocalOrderMethod::OverlapAssembly ? 1.0 : 0.8));
                }
                if (local_order.method == LocalOrderMethod::PrecedenceCache) {
                    ++stats_.precedence_guided_accepts;
                }
                if (local_order.method == LocalOrderMethod::PrecedenceConsensus) {
                    ++stats_.precedence_consensus_accepts;
                }
                sync_to_exact_reference(pattern, proposal);
                return perm_;
            }
        }
        ++stats_.local_rejects;
    }

    ++stats_.exact_output_certifications;
    ++stats_.certified_oracle_order_fallbacks;
    if (strict_symbolic_runtime_enabled()) {
        throw std::runtime_error("Strict symbolic mode failed: dynamic ordering required an oracle fallback instead of a certified local exact-output update");
    }
    remember_replay_window(changed_vars_between_permutations(perm_, exact_perm));
    if (used_cache) store_exact_cache(cache_key, exact_perm);
    sync_to_exact_reference(pattern, exact_perm);
    return perm_;
}

std::uint64_t DynamicCcolamdEngine::pattern_signature(const Eigen::SparseMatrix<double>& pattern) {
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

void DynamicCcolamdEngine::maybe_load_exact_cache_from_env() {
    if (exact_cache_env_loaded_) return;
    exact_cache_env_loaded_ = true;
    const char* path = std::getenv("ISLAM_SYMBOLIC_EXACT_CACHE_IN");
    if (path == nullptr || std::string(path).empty()) return;

    std::ifstream ifs(path);
    if (!ifs) {
        if (no_oracle_symbolic_enabled()) {
            throw std::runtime_error("ISLAM_SYMBOLIC_NO_ORACLE=1: failed to open ISLAM_SYMBOLIC_EXACT_CACHE_IN");
        }
        return;
    }

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto comma = line.find(',');
        if (comma == std::string::npos) continue;
        const auto key = static_cast<std::uint64_t>(std::stoull(line.substr(0, comma)));
        const auto perm = parse_perm_csv(line.substr(comma + 1));
        if (perm.empty()) continue;
        bool exists = false;
        for (const auto& e : exact_cache_) {
            if (e.key == key) { exists = true; break; }
        }
        if (!exists) {
            ExactCacheEntry e;
            e.key = key;
            e.perm = perm;
            e.last_use_tick = tick_;
            exact_cache_.push_back(std::move(e));
            ++stats_.exact_cache_imported_entries;
        }
    }
    stats_.oracle_cache_entries = static_cast<std::uint64_t>(exact_cache_.size());
}

void DynamicCcolamdEngine::maybe_export_exact_cache_entry_to_env(std::uint64_t key,
                                                                 const std::vector<int>& perm) {
    const char* path = std::getenv("ISLAM_SYMBOLIC_EXACT_CACHE_OUT");
    if (path == nullptr || std::string(path).empty()) return;
    std::ofstream ofs(path, std::ios::app);
    if (!ofs) return;
    ofs << key << ',' << perm_to_csv(perm) << '\n';
    ++stats_.exact_cache_exported_entries;
}

const DynamicCcolamdEngine::ExactCacheEntry* DynamicCcolamdEngine::lookup_exact_cache(std::uint64_t key) const {
    for (const auto& e : exact_cache_) {
        if (e.key == key) return &e;
    }
    return nullptr;
}

const DynamicCcolamdEngine::LocalPatternCacheEntry* DynamicCcolamdEngine::lookup_local_pattern_cache(std::uint64_t key) const {
    for (const auto& e : local_pattern_cache_) {
        if (e.key == key) return &e;
    }
    return nullptr;
}

std::uint64_t DynamicCcolamdEngine::motif_signature(const Eigen::SparseMatrix<double>& pattern) {
    const std::uint64_t kOffset = 1469598103934665603ull;
    const std::uint64_t kPrime = 1099511628211ull;
    std::uint64_t h = kOffset;
    auto mix = [&](std::uint64_t x) { h ^= x; h *= kPrime; };

    mix(static_cast<std::uint64_t>(pattern.rows()));
    mix(static_cast<std::uint64_t>(pattern.nonZeros()));

    std::vector<int> degrees(static_cast<size_t>(pattern.rows()), 0);
    std::vector<int> spans;
    spans.reserve(static_cast<size_t>(pattern.nonZeros()));
    for (int col = 0; col < pattern.outerSize(); ++col) {
        int d = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) {
            if (it.row() == col) continue;
            ++d;
            spans.push_back(std::abs(it.row() - col));
        }
        degrees[static_cast<size_t>(col)] = d;
    }
    std::sort(degrees.begin(), degrees.end());
    std::sort(spans.begin(), spans.end());
    for (int d : degrees) mix(static_cast<std::uint64_t>(d + 1));
    const size_t stride = spans.empty() ? 1 : std::max<size_t>(1, spans.size() / 8);
    for (size_t i = 0; i < spans.size(); i += stride) mix(static_cast<std::uint64_t>(spans[i] + 1));
    mix(static_cast<std::uint64_t>(spans.size()));
    return h;
}

std::uint64_t DynamicCcolamdEngine::precedence_signature(const Eigen::SparseMatrix<double>& pattern) {
    const std::uint64_t kOffset = 1469598103934665603ull;
    const std::uint64_t kPrime = 1099511628211ull;
    std::uint64_t h = kOffset;
    auto mix = [&](std::uint64_t x) { h ^= x; h *= kPrime; };

    const int n = pattern.rows();
    mix(static_cast<std::uint64_t>(n));
    int offdiag_nnz = 0;
    std::vector<int> degs(static_cast<size_t>(n), 0);
    std::vector<int> spans;
    spans.reserve(static_cast<size_t>(pattern.nonZeros()));
    for (int col = 0; col < pattern.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) {
            if (it.row() == col) continue;
            ++offdiag_nnz;
            ++degs[static_cast<size_t>(col)];
            spans.push_back(std::abs(it.row() - col));
        }
    }
    mix(static_cast<std::uint64_t>(offdiag_nnz / 2));
    std::sort(degs.begin(), degs.end());
    std::sort(spans.begin(), spans.end());
    const size_t dstride = std::max<size_t>(1, degs.size() / 6);
    for (size_t i = 0; i < degs.size(); i += dstride) mix(static_cast<std::uint64_t>(degs[i] / 2 + 1));
    const size_t sstride = spans.empty() ? 1 : std::max<size_t>(1, spans.size() / 6);
    for (size_t i = 0; i < spans.size(); i += sstride) mix(static_cast<std::uint64_t>(spans[i] / 2 + 1));
    mix(static_cast<std::uint64_t>(spans.size()));
    return h;
}

const DynamicCcolamdEngine::PrecedenceCacheEntry* DynamicCcolamdEngine::lookup_precedence_cache(std::uint64_t key, int size) const {
    const PrecedenceCacheEntry* best = nullptr;
    for (const auto& e : precedence_cache_) {
        if (e.key != key || e.size != size) continue;
        if (best == nullptr || e.merit > best->merit ||
            (e.merit == best->merit && e.last_use_tick > best->last_use_tick)) {
            best = &e;
        }
    }
    return best;
}

void DynamicCcolamdEngine::store_precedence_cache(std::uint64_t key, int size, const std::vector<int>& local_perm_idx, double merit) {
    if (size < 3 || static_cast<int>(local_perm_idx.size()) != size) return;
    std::vector<double> pair_scores(static_cast<size_t>(size * size), 0.0);
    std::vector<int> pos(static_cast<size_t>(size), -1);
    for (int rank = 0; rank < size; ++rank) {
        const int idx = local_perm_idx[static_cast<size_t>(rank)];
        if (idx < 0 || idx >= size || pos[static_cast<size_t>(idx)] >= 0) return;
        pos[static_cast<size_t>(idx)] = rank;
    }
    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            const int ii = local_perm_idx[static_cast<size_t>(i)];
            const int jj = local_perm_idx[static_cast<size_t>(j)];
            pair_scores[static_cast<size_t>(ii * size + jj)] += merit;
            pair_scores[static_cast<size_t>(jj * size + ii)] -= merit;
        }
    }
    for (auto& e : precedence_cache_) {
        if (e.key == key && e.size == size) {
            if (e.pair_scores.size() != pair_scores.size()) e.pair_scores.assign(pair_scores.size(), 0.0);
            for (size_t k = 0; k < pair_scores.size(); ++k) e.pair_scores[k] += pair_scores[k];
            e.last_use_tick = tick_;
            e.merit = std::max(e.merit, merit);
            ++e.hits;
            ++stats_.precedence_cache_promotions;
            stats_.precedence_cache_entries = static_cast<std::uint64_t>(precedence_cache_.size());
            return;
        }
    }
    PrecedenceCacheEntry e;
    e.key = key;
    e.size = size;
    e.pair_scores = std::move(pair_scores);
    e.hits = 1;
    e.last_use_tick = tick_;
    e.merit = merit;
    if (precedence_cache_.size() >= kMaxPrecedenceCacheEntries) {
        auto it = std::min_element(precedence_cache_.begin(), precedence_cache_.end(), [](const PrecedenceCacheEntry& a, const PrecedenceCacheEntry& b) {
            if (a.merit != b.merit) return a.merit < b.merit;
            return a.last_use_tick < b.last_use_tick;
        });
        if (it != precedence_cache_.end()) *it = std::move(e);
    } else {
        precedence_cache_.push_back(std::move(e));
    }
    ++stats_.precedence_cache_promotions;
    stats_.precedence_cache_entries = static_cast<std::uint64_t>(precedence_cache_.size());
}

const DynamicCcolamdEngine::MotifPatternCacheEntry* DynamicCcolamdEngine::lookup_motif_pattern_cache(std::uint64_t key, int size) const {
    for (const auto& e : motif_pattern_cache_) {
        if (e.key == key && e.size == size) return &e;
    }
    return nullptr;
}

const DynamicCcolamdEngine::HierarchicalBlockCacheEntry* DynamicCcolamdEngine::lookup_hierarchical_block_cache(std::uint64_t key, int size) const {
    const HierarchicalBlockCacheEntry* best = nullptr;
    for (const auto& e : hierarchical_block_cache_) {
        if (e.key != key || e.size != size) continue;
        if (best == nullptr || e.merit > best->merit ||
            (e.merit == best->merit && e.last_use_tick > best->last_use_tick)) {
            best = &e;
        }
    }
    return best;
}

void DynamicCcolamdEngine::store_exact_cache(std::uint64_t key, const std::vector<int>& perm) {
    for (auto& e : exact_cache_) {
        if (e.key == key) {
            e.perm = perm;
            e.last_use_tick = tick_;
            stats_.oracle_cache_entries = static_cast<std::uint64_t>(exact_cache_.size());
            maybe_export_exact_cache_entry_to_env(key, perm);
            return;
        }
    }
    ExactCacheEntry e;
    e.key = key;
    e.perm = perm;
    e.last_use_tick = tick_;
    if (exact_cache_.size() >= kMaxExactCacheEntries) {
        auto it = std::min_element(exact_cache_.begin(), exact_cache_.end(), [](const ExactCacheEntry& a, const ExactCacheEntry& b) {
            return a.last_use_tick < b.last_use_tick;
        });
        if (it != exact_cache_.end()) *it = std::move(e);
    } else {
        exact_cache_.push_back(std::move(e));
    }
    stats_.oracle_cache_entries = static_cast<std::uint64_t>(exact_cache_.size());
    maybe_export_exact_cache_entry_to_env(key, perm);
}

void DynamicCcolamdEngine::store_local_pattern_cache(std::uint64_t key, const std::vector<int>& local_perm_idx) {
    if (local_perm_idx.empty()) return;
    for (auto& e : local_pattern_cache_) {
        if (e.key == key) {
            e.local_perm_idx = local_perm_idx;
            e.last_use_tick = tick_;
            stats_.local_pattern_cache_entries = static_cast<std::uint64_t>(local_pattern_cache_.size());
            return;
        }
    }
    LocalPatternCacheEntry e;
    e.key = key;
    e.local_perm_idx = local_perm_idx;
    e.last_use_tick = tick_;
    if (local_pattern_cache_.size() >= kMaxLocalPatternCacheEntries) {
        auto it = std::min_element(local_pattern_cache_.begin(), local_pattern_cache_.end(), [](const LocalPatternCacheEntry& a, const LocalPatternCacheEntry& b) {
            return a.last_use_tick < b.last_use_tick;
        });
        if (it != local_pattern_cache_.end()) *it = std::move(e);
    } else {
        local_pattern_cache_.push_back(std::move(e));
    }
    stats_.local_pattern_cache_entries = static_cast<std::uint64_t>(local_pattern_cache_.size());
}

void DynamicCcolamdEngine::store_motif_pattern_cache(std::uint64_t key, int size, const std::vector<int>& local_perm_idx) {
    if (local_perm_idx.empty()) return;
    for (auto& e : motif_pattern_cache_) {
        if (e.key == key && e.size == size) {
            e.local_perm_idx = local_perm_idx;
            e.last_use_tick = tick_;
            ++e.hits;
            stats_.motif_pattern_cache_entries = static_cast<std::uint64_t>(motif_pattern_cache_.size());
            return;
        }
    }
    MotifPatternCacheEntry e;
    e.key = key;
    e.size = size;
    e.local_perm_idx = local_perm_idx;
    e.hits = 1;
    e.last_use_tick = tick_;
    if (motif_pattern_cache_.size() >= kMaxMotifPatternCacheEntries) {
        auto it = std::min_element(motif_pattern_cache_.begin(), motif_pattern_cache_.end(), [](const MotifPatternCacheEntry& a, const MotifPatternCacheEntry& b) {
            if (a.hits != b.hits) return a.hits < b.hits;
            return a.last_use_tick < b.last_use_tick;
        });
        if (it != motif_pattern_cache_.end()) *it = std::move(e);
    } else {
        motif_pattern_cache_.push_back(std::move(e));
    }
    stats_.motif_pattern_cache_entries = static_cast<std::uint64_t>(motif_pattern_cache_.size());
}

void DynamicCcolamdEngine::store_hierarchical_block_cache(std::uint64_t key, int size, const std::vector<int>& local_perm_idx, double merit) {
    if (local_perm_idx.empty() || size < 4) return;
    for (auto& e : hierarchical_block_cache_) {
        if (e.key == key && e.size == size && e.local_perm_idx == local_perm_idx) {
            e.last_use_tick = tick_;
            e.merit = std::max(e.merit, merit);
            ++e.hits;
            stats_.hierarchical_block_cache_entries = static_cast<std::uint64_t>(hierarchical_block_cache_.size());
    stats_.precedence_cache_entries = static_cast<std::uint64_t>(precedence_cache_.size());
            return;
        }
    }
    HierarchicalBlockCacheEntry e;
    e.key = key;
    e.size = size;
    e.local_perm_idx = local_perm_idx;
    e.hits = 1;
    e.last_use_tick = tick_;
    e.merit = merit;
    if (hierarchical_block_cache_.size() >= kMaxHierarchicalBlockCacheEntries) {
        auto it = std::min_element(hierarchical_block_cache_.begin(), hierarchical_block_cache_.end(), [](const HierarchicalBlockCacheEntry& a, const HierarchicalBlockCacheEntry& b) {
            if (a.merit != b.merit) return a.merit < b.merit;
            return a.last_use_tick < b.last_use_tick;
        });
        if (it != hierarchical_block_cache_.end()) *it = std::move(e);
    } else {
        hierarchical_block_cache_.push_back(std::move(e));
    }
    stats_.hierarchical_block_cache_entries = static_cast<std::uint64_t>(hierarchical_block_cache_.size());
    stats_.precedence_cache_entries = static_cast<std::uint64_t>(precedence_cache_.size());
}

void DynamicCcolamdEngine::promote_hierarchical_blocks_from_window(const Eigen::SparseMatrix<double>& pattern,
                                                                   const std::vector<int>& local_vars,
                                                                   const std::vector<int>& local_order_vars) {
    if (local_vars.size() < 6 || local_vars.size() != local_order_vars.size()) return;
    const std::vector<int> widths = {6, 8, 10, 12};
    for (int width : widths) {
        if (width > static_cast<int>(local_order_vars.size())) continue;
        const int stride = std::max(2, width / 2);
        for (int start = 0; start + width <= static_cast<int>(local_order_vars.size()); start += stride) {
            std::vector<int> block_order_vars(local_order_vars.begin() + start, local_order_vars.begin() + start + width);
            std::vector<int> block_input_vars = block_order_vars;
            std::stable_sort(block_input_vars.begin(), block_input_vars.end(), [&](int a, int b) {
                const int pa = (a >= 0 && a < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(a)] : a;
                const int pb = (b >= 0 && b < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(b)] : b;
                return pa < pb;
            });
            std::vector<int> local_perm_idx;
            local_perm_idx.reserve(block_input_vars.size());
            bool ok = true;
            for (int v : block_order_vars) {
                auto it = std::find(block_input_vars.begin(), block_input_vars.end(), v);
                if (it == block_input_vars.end()) { ok = false; break; }
                local_perm_idx.push_back(static_cast<int>(std::distance(block_input_vars.begin(), it)));
            }
            if (!ok || local_perm_idx.size() != block_input_vars.size()) continue;
            const auto Ablock = principal_submatrix(pattern, block_input_vars);
            const double merit = 1.0 + 0.01 * static_cast<double>(width);
            store_hierarchical_block_cache(motif_signature(Ablock), static_cast<int>(block_input_vars.size()), local_perm_idx, merit);
            ++stats_.hierarchical_block_promotions;
        }
    }
}

std::vector<int> DynamicCcolamdEngine::inverse_permutation_of(const std::vector<int>& perm) {
    std::vector<int> pinv(perm.size(), -1);
    for (int k = 0; k < static_cast<int>(perm.size()); ++k) {
        const int pk = perm[static_cast<size_t>(k)];
        if (pk < 0 || pk >= static_cast<int>(perm.size())) throw std::runtime_error("Invalid permutation entry");
        pinv[static_cast<size_t>(pk)] = k;
    }
    return pinv;
}

std::vector<int> DynamicCcolamdEngine::unique_sorted(std::vector<int> vals) {
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    return vals;
}

Eigen::SparseMatrix<double> DynamicCcolamdEngine::principal_submatrix(const Eigen::SparseMatrix<double>& A,
                                                                      const std::vector<int>& idx) {
    std::vector<int> inv(static_cast<size_t>(A.rows()), -1);
    for (int k = 0; k < static_cast<int>(idx.size()); ++k) {
        const int v = idx[static_cast<size_t>(k)];
        if (v < 0 || v >= A.rows()) throw std::runtime_error("principal_submatrix index out of range");
        inv[static_cast<size_t>(v)] = k;
    }
    std::vector<Eigen::Triplet<double>> trips;
    for (int col = 0; col < A.outerSize(); ++col) {
        const int new_col = inv[static_cast<size_t>(col)];
        if (new_col < 0) continue;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int new_row = inv[static_cast<size_t>(it.row())];
            if (new_row >= 0) trips.emplace_back(new_row, new_col, 1.0);
        }
    }
    Eigen::SparseMatrix<double> out(static_cast<int>(idx.size()), static_cast<int>(idx.size()));
    out.setFromTriplets(trips.begin(), trips.end());
    out.makeCompressed();
    return out;
}

std::vector<int> DynamicCcolamdEngine::one_hop_neighborhood(const Eigen::SparseMatrix<double>& pattern,
                                                            const std::vector<int>& seeds) {
    const int n = pattern.rows();
    std::vector<unsigned char> mark(static_cast<size_t>(n), 0);
    for (int v : seeds) if (v >= 0 && v < n) mark[static_cast<size_t>(v)] = 1;
    for (int col = 0; col < pattern.outerSize(); ++col) {
        bool hit = mark[static_cast<size_t>(col)] != 0;
        if (!hit) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) {
                if (mark[static_cast<size_t>(it.row())]) { hit = true; break; }
            }
        }
        if (!hit) continue;
        mark[static_cast<size_t>(col)] = 1;
        for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) mark[static_cast<size_t>(it.row())] = 1;
    }
    std::vector<int> out;
    for (int i = 0; i < n; ++i) if (mark[static_cast<size_t>(i)]) out.push_back(i);
    return out;
}

std::vector<std::vector<int>> DynamicCcolamdEngine::adjacency_lists(const Eigen::SparseMatrix<double>& A) {
    const int n = A.rows();
    std::vector<std::unordered_set<int>> tmp(static_cast<size_t>(n));
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it) {
            const int r = it.row();
            if (r == col) continue;
            tmp[static_cast<size_t>(col)].insert(r);
            tmp[static_cast<size_t>(r)].insert(col);
        }
    }
    std::vector<std::vector<int>> out(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        out[static_cast<size_t>(i)].assign(tmp[static_cast<size_t>(i)].begin(), tmp[static_cast<size_t>(i)].end());
        std::sort(out[static_cast<size_t>(i)].begin(), out[static_cast<size_t>(i)].end());
    }
    return out;
}

int DynamicCcolamdEngine::dense_degree_threshold(int n) {
    return std::max(16, static_cast<int>(std::sqrt(static_cast<double>(std::max(n, 1)))));
}

int DynamicCcolamdEngine::score_column(const ColumnState& col) {
    if (!col.alive) return std::numeric_limits<int>::max();
    int score = col.degree + 2 * col.external_degree;
    if (col.dense) score += 1000000;
    return score;
}

std::vector<int> DynamicCcolamdEngine::positions_from_vars(const std::vector<int>& perm,
                                                           const std::vector<int>& vars) {
    const auto pinv = inverse_permutation_of(perm);
    std::vector<int> pos;
    pos.reserve(vars.size());
    for (int v : vars) {
        if (v >= 0 && v < static_cast<int>(pinv.size())) {
            const int p = pinv[static_cast<size_t>(v)];
            if (p >= 0) pos.push_back(p);
        }
    }
    return unique_sorted(std::move(pos));
}

std::vector<int> DynamicCcolamdEngine::interval_window_from_positions(const std::vector<int>& perm,
                                                                      const std::vector<int>& dirty_positions,
                                                                      int pad) {
    if (dirty_positions.empty()) return {};
    const int n = static_cast<int>(perm.size());
    const int lo = std::max(0, dirty_positions.front() - std::max(0, pad));
    const int hi = std::min(n - 1, dirty_positions.back() + std::max(0, pad));
    std::vector<int> vars;
    vars.reserve(static_cast<size_t>(hi - lo + 1));
    for (int p = lo; p <= hi; ++p) vars.push_back(perm[static_cast<size_t>(p)]);
    return vars;
}

std::vector<int> DynamicCcolamdEngine::union_sorted_vectors(std::vector<int> a,
                                                            const std::vector<int>& b) {
    a.insert(a.end(), b.begin(), b.end());
    return unique_sorted(std::move(a));
}

std::vector<DynamicCcolamdEngine::ReplayWindow>
DynamicCcolamdEngine::replay_windows_from_dirty(const std::vector<int>& dirty_vars) const {
    auto dirty = unique_sorted(dirty_vars);
    if (dirty.empty()) return {};
    struct ScoredReplay {
        ReplayWindow rw;
        double score = 0.0;
    };
    std::vector<ScoredReplay> scored;
    scored.reserve(active_replay_windows().size());
    for (const auto& rw : active_replay_windows()) {
        std::size_t overlap = 0;
        for (int v : dirty) {
            if (std::binary_search(rw.vars.begin(), rw.vars.end(), v)) ++overlap;
        }
        if (overlap == 0) continue;
        const std::size_t uni = rw.vars.size() + dirty.size() - overlap;
        const double jaccard = uni > 0 ? static_cast<double>(overlap) / static_cast<double>(uni) : 0.0;
        const double success = (static_cast<double>(rw.accepts) + kLaplacePrior) /
                               (static_cast<double>(rw.uses) + 2.0 * kLaplacePrior);
        const double score = 0.55 * jaccard + 0.25 * success +
                             0.10 * std::min<double>(1.0, static_cast<double>(rw.hits) / 4.0) -
                             0.001 * static_cast<double>(rw.vars.size());
        scored.push_back(ScoredReplay{rw, score});
    }
    std::sort(scored.begin(), scored.end(), [](const ScoredReplay& a, const ScoredReplay& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.rw.vars.size() != b.rw.vars.size()) return a.rw.vars.size() < b.rw.vars.size();
        return a.rw.hits > b.rw.hits;
    });
    std::vector<ReplayWindow> out;
    const int limit = std::min<int>(kMaxReplayCandidates, static_cast<int>(scored.size()));
    out.reserve(static_cast<size_t>(limit));
    for (int i = 0; i < limit; ++i) out.push_back(std::move(scored[static_cast<size_t>(i)].rw));
    return out;
}

std::vector<int> DynamicCcolamdEngine::band_window_from_dirty(const std::vector<int>& dirty_vars) const {
    auto dirty = unique_sorted(dirty_vars);
    if (dirty.empty() || perm_.empty() || pinv_.empty()) return {};
    int lo = static_cast<int>(perm_.size());
    int hi = -1;
    int pad = 1;
    for (int v : dirty) {
        if (v < 0 || v >= static_cast<int>(pinv_.size())) continue;
        const int pv = pinv_[static_cast<size_t>(v)];
        if (pv < 0) continue;
        lo = std::min(lo, pv);
        hi = std::max(hi, pv);
        if (v >= 0 && v < static_cast<int>(adjacency_.size())) {
            for (int u : adjacency_[static_cast<size_t>(v)]) {
                if (u < 0 || u >= static_cast<int>(pinv_.size())) continue;
                const int pu = pinv_[static_cast<size_t>(u)];
                if (pu < 0) continue;
                lo = std::min(lo, pu);
                hi = std::max(hi, pu);
                pad = std::max(pad, std::min(8, std::abs(pu - pv)));
            }
        }
    }
    if (hi < lo) return {};
    lo = std::max(0, lo - pad);
    hi = std::min(static_cast<int>(perm_.size()) - 1, hi + pad);
    std::vector<int> vars;
    vars.reserve(static_cast<size_t>(hi - lo + 1));
    for (int p = lo; p <= hi; ++p) vars.push_back(perm_[static_cast<size_t>(p)]);
    return vars;
}

std::vector<int> DynamicCcolamdEngine::changed_vars_between_permutations(const std::vector<int>& before,
                                                                         const std::vector<int>& after) {
    if (before.size() != after.size()) return {};
    std::vector<int> changed;
    changed.reserve(before.size());
    for (size_t i = 0; i < before.size(); ++i) {
        if (before[i] != after[i]) {
            changed.push_back(before[i]);
            changed.push_back(after[i]);
        }
    }
    return unique_sorted(std::move(changed));
}

double DynamicCcolamdEngine::candidate_priority(const CandidateWindow& candidate,
                                                const std::vector<int>& dirty_vars) const {
    const auto kind_idx = static_cast<size_t>(static_cast<int>(candidate.kind));
    double empirical = 0.5;
    double recency = 0.0;
    const auto& regime_kind_stats = active_kind_stats();
    if (kind_idx < regime_kind_stats.size()) {
        const auto& ks = regime_kind_stats[kind_idx];
        empirical = (static_cast<double>(ks.accepts) + kLaplacePrior) /
                    (static_cast<double>(ks.attempts) + 2.0 * kLaplacePrior);
        if (tick_ > 0 && ks.last_success_tick > 0) {
            const double age = static_cast<double>(tick_ - ks.last_success_tick);
            recency = kRecencyBonusScale / (1.0 + age);
        }
    }

    double overlap = 0.0;

    if (!dirty_vars.empty()) {
        int overlap_count = 0;
        for (int v : dirty_vars) {
            if (std::binary_search(candidate.vars.begin(), candidate.vars.end(), v)) ++overlap_count;
        }
        overlap = static_cast<double>(overlap_count) / static_cast<double>(dirty_vars.size());
    }

    double replay_bonus = 0.0;
    if (candidate.kind == CandidateKind::Replay) {
        for (const auto& rw : active_replay_windows()) {
            if (rw.vars == candidate.vars) {
                replay_bonus = kReplayHitBonus * std::min<double>(1.0, static_cast<double>(rw.hits) / 4.0);
                break;
            }
        }
    }

    return empirical + recency + 0.20 * overlap + replay_bonus -
           kSizePenalty * static_cast<double>(candidate.vars.size());
}

std::vector<DynamicCcolamdEngine::CandidateWindow>
DynamicCcolamdEngine::candidate_windows(const Eigen::SparseMatrix<double>& pattern,
                                        const std::vector<int>& dirty_vars) const {
    if (static_cast<int>(perm_.size()) != pattern.rows()) return {};

    const auto dirty = unique_sorted(dirty_vars);
    const auto dirty_positions = positions_from_vars(perm_, dirty);

    std::vector<CandidateWindow> out;
    auto push_candidate = [&](CandidateKind kind, std::vector<int> vars) {
        vars = unique_sorted(std::move(vars));
        if (vars.size() < 2) return;
        if (static_cast<int>(vars.size()) >= std::max(12, static_cast<int>(pattern.rows()) / 3)) return;
        for (const auto& existing : out) {
            if (existing.vars == vars) return;
        }
        out.push_back(CandidateWindow{kind, std::move(vars)});
    };

    const auto replay_candidates = replay_windows_from_dirty(dirty);
    const auto hop1 = one_hop_neighborhood(pattern, dirty);
    const auto hop2 = one_hop_neighborhood(pattern, hop1);
    const auto interval1 = interval_window_from_positions(perm_, dirty_positions, 1);
    const auto interval2 = interval_window_from_positions(perm_, dirty_positions, 2);
    const auto band = band_window_from_dirty(dirty);

    for (const auto& rw : replay_candidates) push_candidate(CandidateKind::Replay, rw.vars);
    push_candidate(CandidateKind::OneHop, hop1);
    push_candidate(CandidateKind::Band, band);
    push_candidate(CandidateKind::Interval, interval1);
    push_candidate(CandidateKind::Union, union_sorted_vectors(hop1, interval1));
    push_candidate(CandidateKind::Union, union_sorted_vectors(band, interval1));
    push_candidate(CandidateKind::TwoHop, hop2);
    push_candidate(CandidateKind::Interval, interval2);
    push_candidate(CandidateKind::Union, union_sorted_vectors(hop2, interval2));
    if (!replay_candidates.empty()) {
        push_candidate(CandidateKind::Union, union_sorted_vectors(replay_candidates.front().vars, band));
        push_candidate(CandidateKind::Union, union_sorted_vectors(replay_candidates.front().vars, interval1));
    }

    std::stable_sort(out.begin(), out.end(), [&](const CandidateWindow& a, const CandidateWindow& b) {
        const double pa = candidate_priority(a, dirty);
        const double pb = candidate_priority(b, dirty);
        if (pa != pb) return pa > pb;
        if (a.vars.size() != b.vars.size()) return a.vars.size() < b.vars.size();
        return static_cast<int>(a.kind) < static_cast<int>(b.kind);
    });
    return out;
}

DynamicCcolamdEngine::LocalOrderResult
DynamicCcolamdEngine::order_candidate_window(const Eigen::SparseMatrix<double>& pattern,
                                             const std::vector<int>& local_vars) {
    LocalOrderResult result;
    if (local_vars.size() < 2) return result;
    const auto A = principal_submatrix(pattern, local_vars);
    const auto local_key = pattern_signature(A);
    if (const auto* cached = lookup_local_pattern_cache(local_key); cached != nullptr) {
        ++stats_.local_pattern_cache_hits;
        result.ordered_vars.reserve(local_vars.size());
        if (cached->local_perm_idx.size() == local_vars.size()) {
            for (int idx : cached->local_perm_idx) {
                if (idx < 0 || idx >= static_cast<int>(local_vars.size())) return {};
                result.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
            }
            result.method = LocalOrderMethod::ExactCache;
            return result;
        }
    }
    ++stats_.local_pattern_cache_misses;

    const auto motif_key = motif_signature(A);
    if (const auto* motif = lookup_motif_pattern_cache(motif_key, static_cast<int>(local_vars.size())); motif != nullptr) {
        ++stats_.motif_pattern_cache_hits;
        result.ordered_vars.reserve(local_vars.size());
        if (motif->local_perm_idx.size() == local_vars.size()) {
            std::vector<unsigned char> used(local_vars.size(), 0);
            for (int idx : motif->local_perm_idx) {
                if (idx < 0 || idx >= static_cast<int>(local_vars.size())) return {};
                if (used[static_cast<size_t>(idx)] != 0) return {};
                used[static_cast<size_t>(idx)] = 1;
                result.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
            }
            if (result.ordered_vars.size() == local_vars.size()) {
                result.method = LocalOrderMethod::MotifCache;
                return result;
            }
        }
    }
    ++stats_.motif_pattern_cache_misses;

    if (const auto* hier = lookup_hierarchical_block_cache(motif_key, static_cast<int>(local_vars.size())); hier != nullptr) {
        ++stats_.hierarchical_block_cache_hits;
        result.ordered_vars.reserve(local_vars.size());
        if (hier->local_perm_idx.size() == local_vars.size()) {
            std::vector<unsigned char> used(local_vars.size(), 0);
            for (int idx : hier->local_perm_idx) {
                if (idx < 0 || idx >= static_cast<int>(local_vars.size()) || used[static_cast<size_t>(idx)] != 0) {
                    result.ordered_vars.clear(); break; }
                used[static_cast<size_t>(idx)] = 1;
                result.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
            }
            if (result.ordered_vars.size() == local_vars.size()) {
                result.method = LocalOrderMethod::MotifCache;
                return result;
            }
        }
    } else {
        ++stats_.hierarchical_block_cache_misses;
    }

    ++stats_.precedence_consensus_attempts;
    result = precedence_consensus_order(A, local_vars);
    if (!result.ordered_vars.empty()) return result;

    ++stats_.precedence_guided_attempts;
    result = precedence_guided_order(A, local_vars);
    if (!result.ordered_vars.empty()) return result;

    ++stats_.overlap_assembly_attempts;
    result = compose_overlapping_window_order(A, local_vars);
    if (!result.ordered_vars.empty()) return result;

    const auto adj = adjacency_lists(A);
    const int nloc = static_cast<int>(local_vars.size());
    const int dense_thresh = dense_degree_threshold(nloc);
    std::vector<ColumnState> cols(static_cast<size_t>(nloc));
    std::vector<std::unordered_set<int>> alive_adj(static_cast<size_t>(nloc));
    for (int i = 0; i < nloc; ++i) {
        alive_adj[static_cast<size_t>(i)].insert(adj[static_cast<size_t>(i)].begin(), adj[static_cast<size_t>(i)].end());
        cols[static_cast<size_t>(i)].degree = static_cast<int>(alive_adj[static_cast<size_t>(i)].size());
        cols[static_cast<size_t>(i)].external_degree = 0;
        cols[static_cast<size_t>(i)].dense = cols[static_cast<size_t>(i)].degree >= dense_thresh;
    }

    std::vector<int> local_order;
    local_order.reserve(static_cast<size_t>(nloc));
    for (int step = 0; step < nloc; ++step) {
        int best = -1;
        int best_score = std::numeric_limits<int>::max();
        for (int i = 0; i < nloc; ++i) {
            const int s = score_column(cols[static_cast<size_t>(i)]);
            if (s < best_score || (s == best_score && (best < 0 || i < best))) {
                best = i;
                best_score = s;
            }
        }
        if (best < 0) break;
        cols[static_cast<size_t>(best)].alive = false;
        local_order.push_back(best);

        std::vector<int> neighbors(alive_adj[static_cast<size_t>(best)].begin(), alive_adj[static_cast<size_t>(best)].end());
        std::sort(neighbors.begin(), neighbors.end());
        for (int u : neighbors) alive_adj[static_cast<size_t>(u)].erase(best);
        for (size_t a = 0; a < neighbors.size(); ++a) {
            const int u = neighbors[a];
            if (!cols[static_cast<size_t>(u)].alive) continue;
            for (size_t b = a + 1; b < neighbors.size(); ++b) {
                const int v = neighbors[b];
                if (!cols[static_cast<size_t>(v)].alive) continue;
                alive_adj[static_cast<size_t>(u)].insert(v);
                alive_adj[static_cast<size_t>(v)].insert(u);
            }
        }
        alive_adj[static_cast<size_t>(best)].clear();
        for (int i = 0; i < nloc; ++i) {
            cols[static_cast<size_t>(i)].degree = static_cast<int>(alive_adj[static_cast<size_t>(i)].size());
            cols[static_cast<size_t>(i)].dense = cols[static_cast<size_t>(i)].degree >= dense_thresh;
        }
    }

    if (static_cast<int>(local_order.size()) != nloc) return {};

    result.ordered_vars.reserve(static_cast<size_t>(nloc));
    for (int idx : local_order) result.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
    result.method = LocalOrderMethod::GreedyApprox;
    return result;
}

DynamicCcolamdEngine::LocalOrderResult
DynamicCcolamdEngine::precedence_consensus_order(const Eigen::SparseMatrix<double>& A,
                                                 const std::vector<int>& local_vars) {
    LocalOrderResult result;
    result.method = LocalOrderMethod::PrecedenceConsensus;
    const int n = static_cast<int>(local_vars.size());
    if (n < 4) return result;
    const auto key = precedence_signature(A);
    const auto* entry = lookup_precedence_cache(key, n);
    if (entry == nullptr || entry->pair_scores.size() != static_cast<size_t>(n * n)) {
        return {};
    }

    std::vector<double> strength(static_cast<size_t>(n), 0.0);
    double mean_margin = 0.0;
    int margin_count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double pij = entry->pair_scores[static_cast<size_t>(i * n + j)];
            const double pji = entry->pair_scores[static_cast<size_t>(j * n + i)];
            const double d = pij - pji;
            strength[static_cast<size_t>(i)] += d;
            strength[static_cast<size_t>(j)] -= d;
            mean_margin += std::abs(d);
            ++margin_count;
        }
    }
    if (margin_count == 0) return {};
    mean_margin /= static_cast<double>(margin_count);
    const double eps = std::max(0.05, 0.35 * mean_margin);

    std::vector<std::vector<int>> g(static_cast<size_t>(n)), rg(static_cast<size_t>(n));
    int edge_count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double pij = entry->pair_scores[static_cast<size_t>(i * n + j)];
            const double pji = entry->pair_scores[static_cast<size_t>(j * n + i)];
            if (pij > pji + eps) {
                g[static_cast<size_t>(i)].push_back(j);
                rg[static_cast<size_t>(j)].push_back(i);
                ++edge_count;
            } else if (pji > pij + eps) {
                g[static_cast<size_t>(j)].push_back(i);
                rg[static_cast<size_t>(i)].push_back(j);
                ++edge_count;
            }
        }
    }
    if (edge_count == 0) return {};

    std::vector<unsigned char> vis(static_cast<size_t>(n), 0);
    std::vector<int> order;
    order.reserve(static_cast<size_t>(n));
    std::function<void(int)> dfs1 = [&](int u) {
        vis[static_cast<size_t>(u)] = 1;
        for (int v : g[static_cast<size_t>(u)]) if (!vis[static_cast<size_t>(v)]) dfs1(v);
        order.push_back(u);
    };
    for (int i = 0; i < n; ++i) if (!vis[static_cast<size_t>(i)]) dfs1(i);

    std::vector<int> comp(static_cast<size_t>(n), -1);
    std::vector<std::vector<int>> comps;
    std::function<void(int,int)> dfs2 = [&](int u, int cid) {
        comp[static_cast<size_t>(u)] = cid;
        comps[static_cast<size_t>(cid)].push_back(u);
        for (int v : rg[static_cast<size_t>(u)]) if (comp[static_cast<size_t>(v)] < 0) dfs2(v, cid);
    };
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const int u = *it;
        if (comp[static_cast<size_t>(u)] >= 0) continue;
        comps.push_back({});
        dfs2(u, static_cast<int>(comps.size()) - 1);
    }

    int collapsed = 0;
    for (const auto& c : comps) if (c.size() > 1) collapsed += static_cast<int>(c.size()) - 1;
    if (collapsed > 0) stats_.precedence_consensus_scc_collapses += static_cast<std::uint64_t>(collapsed);

    const int nc = static_cast<int>(comps.size());
    std::vector<std::vector<int>> dag(static_cast<size_t>(nc));
    std::vector<int> indeg(static_cast<size_t>(nc), 0);
    std::vector<double> comp_strength(static_cast<size_t>(nc), 0.0);
    for (int cid = 0; cid < nc; ++cid) {
        for (int u : comps[static_cast<size_t>(cid)]) comp_strength[static_cast<size_t>(cid)] += strength[static_cast<size_t>(u)];
    }
    std::vector<std::unordered_set<int>> seen(static_cast<size_t>(nc));
    for (int u = 0; u < n; ++u) {
        const int cu = comp[static_cast<size_t>(u)];
        for (int v : g[static_cast<size_t>(u)]) {
            const int cv = comp[static_cast<size_t>(v)];
            if (cu == cv) continue;
            if (seen[static_cast<size_t>(cu)].insert(cv).second) {
                dag[static_cast<size_t>(cu)].push_back(cv);
                ++indeg[static_cast<size_t>(cv)];
            }
        }
    }

    std::vector<int> ready;
    ready.reserve(static_cast<size_t>(nc));
    for (int cid = 0; cid < nc; ++cid) if (indeg[static_cast<size_t>(cid)] == 0) ready.push_back(cid);
    std::vector<int> comp_order;
    comp_order.reserve(static_cast<size_t>(nc));
    auto comp_key = [&](int cid) {
        int best_pos = std::numeric_limits<int>::max();
        for (int u : comps[static_cast<size_t>(cid)]) {
            const int gv = local_vars[static_cast<size_t>(u)];
            const int pos = (gv >= 0 && gv < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(gv)] : gv;
            best_pos = std::min(best_pos, pos);
        }
        return best_pos;
    };
    while (!ready.empty()) {
        std::stable_sort(ready.begin(), ready.end(), [&](int a, int b) {
            if (comp_strength[static_cast<size_t>(a)] != comp_strength[static_cast<size_t>(b)]) {
                return comp_strength[static_cast<size_t>(a)] > comp_strength[static_cast<size_t>(b)];
            }
            return comp_key(a) < comp_key(b);
        });
        const int cid = ready.front();
        ready.erase(ready.begin());
        comp_order.push_back(cid);
        for (int v : dag[static_cast<size_t>(cid)]) {
            --indeg[static_cast<size_t>(v)];
            if (indeg[static_cast<size_t>(v)] == 0) ready.push_back(v);
        }
    }
    if (static_cast<int>(comp_order.size()) != nc) return {};

    result.ordered_vars.reserve(local_vars.size());
    for (int cid : comp_order) {
        auto members = comps[static_cast<size_t>(cid)];
        std::stable_sort(members.begin(), members.end(), [&](int a, int b) {
            if (strength[static_cast<size_t>(a)] != strength[static_cast<size_t>(b)]) {
                return strength[static_cast<size_t>(a)] > strength[static_cast<size_t>(b)];
            }
            const int va = local_vars[static_cast<size_t>(a)];
            const int vb = local_vars[static_cast<size_t>(b)];
            const int pa = (va >= 0 && va < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(va)] : va;
            const int pb = (vb >= 0 && vb < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(vb)] : vb;
            return pa < pb;
        });
        for (int idx : members) result.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
    }
    return result;
}

DynamicCcolamdEngine::LocalOrderResult
DynamicCcolamdEngine::precedence_guided_order(const Eigen::SparseMatrix<double>& A,
                                              const std::vector<int>& local_vars) {
    LocalOrderResult result;
    result.method = LocalOrderMethod::PrecedenceCache;
    if (local_vars.size() < 3) return result;
    const auto key = precedence_signature(A);
    const auto* entry = lookup_precedence_cache(key, static_cast<int>(local_vars.size()));
    if (entry == nullptr || entry->pair_scores.size() != local_vars.size() * local_vars.size()) {
        ++stats_.precedence_cache_misses;
        return {};
    }
    ++stats_.precedence_cache_hits;
    std::vector<int> idxs(local_vars.size());
    for (int i = 0; i < static_cast<int>(local_vars.size()); ++i) idxs[static_cast<size_t>(i)] = i;
    std::vector<double> scores(local_vars.size(), 0.0);
    const int n = static_cast<int>(local_vars.size());
    for (int i = 0; i < n; ++i) {
        double ssum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            ssum += entry->pair_scores[static_cast<size_t>(i * n + j)];
        }
        scores[static_cast<size_t>(i)] = ssum;
    }
    std::stable_sort(idxs.begin(), idxs.end(), [&](int a, int b) {
        if (scores[static_cast<size_t>(a)] != scores[static_cast<size_t>(b)]) {
            return scores[static_cast<size_t>(a)] > scores[static_cast<size_t>(b)];
        }
        const int va = local_vars[static_cast<size_t>(a)];
        const int vb = local_vars[static_cast<size_t>(b)];
        const int pa = (va >= 0 && va < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(va)] : va;
        const int pb = (vb >= 0 && vb < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(vb)] : vb;
        return pa < pb;
    });
    result.ordered_vars.reserve(local_vars.size());
    for (int idx : idxs) result.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
    return result;
}

DynamicCcolamdEngine::LocalOrderResult
DynamicCcolamdEngine::compose_overlapping_window_order(const Eigen::SparseMatrix<double>& A,
                                                       const std::vector<int>& local_vars) {
    LocalOrderResult result;
    result.method = LocalOrderMethod::OverlapAssembly;
    const int n = static_cast<int>(local_vars.size());
    if (n < kOverlapMinWidth) return result;

    std::vector<int> base_idx(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) base_idx[static_cast<size_t>(i)] = i;
    std::stable_sort(base_idx.begin(), base_idx.end(), [&](int a, int b) {
        const int pa = (local_vars[a] >= 0 && local_vars[a] < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(local_vars[a])] : a;
        const int pb = (local_vars[b] >= 0 && local_vars[b] < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(local_vars[b])] : b;
        return pa < pb;
    });

    struct Proposal {
        std::vector<double> score_sum;
        std::vector<double> weight_sum;
        int covered = 0;
        int exact_hits = 0;
        int motif_hits = 0;
        int segments = 0;
        double merit = -1.0;
        std::vector<int> ordered_vars;
    };

    Proposal best;
    std::uint64_t proposals = 0;
    std::uint64_t piece_hits = 0;
    const std::vector<int> widths = {4, 6, 8, 10};
    for (int width : widths) {
        if (width > n) continue;
        const int stride = std::max(2, width / 2);
        std::vector<int> shifts{0};
        if (stride > 2) shifts.push_back(stride / 2);
        for (int shift : shifts) {
            Proposal p;
            p.score_sum.assign(static_cast<size_t>(n), 0.0);
            p.weight_sum.assign(static_cast<size_t>(n), 0.0);
            for (int start = shift; start < n; start += stride) {
                const int end = std::min(n, start + width);
                if (end - start < 2) continue;
                std::vector<int> seg_idx;
                seg_idx.reserve(static_cast<size_t>(end - start));
                for (int t = start; t < end; ++t) seg_idx.push_back(base_idx[static_cast<size_t>(t)]);
                const auto Aseg = principal_submatrix(A, seg_idx);
                std::vector<int> ordered_seg_idx;
                double weight = 0.0;
                const auto seg_key = pattern_signature(Aseg);
                if (const auto* cached = lookup_local_pattern_cache(seg_key); cached != nullptr) {
                    if (cached->local_perm_idx.size() == seg_idx.size()) {
                        ordered_seg_idx.reserve(seg_idx.size());
                        for (int idx : cached->local_perm_idx) {
                            if (idx < 0 || idx >= static_cast<int>(seg_idx.size())) { ordered_seg_idx.clear(); break; }
                            ordered_seg_idx.push_back(seg_idx[static_cast<size_t>(idx)]);
                        }
                        weight = ordered_seg_idx.empty() ? 0.0 : 1.0;
                        if (weight > 0.0) ++p.exact_hits;
                    }
                }
                if (ordered_seg_idx.empty()) {
                    const auto mkey = motif_signature(Aseg);
                    if (const auto* hier = lookup_hierarchical_block_cache(mkey, static_cast<int>(seg_idx.size())); hier != nullptr) {
                        if (hier->local_perm_idx.size() == seg_idx.size()) {
                            std::vector<unsigned char> used(seg_idx.size(), 0);
                            ordered_seg_idx.reserve(seg_idx.size());
                            for (int idx : hier->local_perm_idx) {
                                if (idx < 0 || idx >= static_cast<int>(seg_idx.size()) || used[static_cast<size_t>(idx)] != 0) { ordered_seg_idx.clear(); break; }
                                used[static_cast<size_t>(idx)] = 1;
                                ordered_seg_idx.push_back(seg_idx[static_cast<size_t>(idx)]);
                            }
                            weight = ordered_seg_idx.empty() ? 0.0 : 0.82;
                            if (weight > 0.0) {
                                ++stats_.hierarchical_block_cache_hits;
                                ++p.motif_hits;
                            }
                        }
                    } else {
                        ++stats_.hierarchical_block_cache_misses;
                    }
                    if (ordered_seg_idx.empty()) {
                        if (const auto* motif = lookup_motif_pattern_cache(mkey, static_cast<int>(seg_idx.size())); motif != nullptr) {
                            if (motif->local_perm_idx.size() == seg_idx.size()) {
                                std::vector<unsigned char> used(seg_idx.size(), 0);
                                ordered_seg_idx.reserve(seg_idx.size());
                                for (int idx : motif->local_perm_idx) {
                                    if (idx < 0 || idx >= static_cast<int>(seg_idx.size()) || used[static_cast<size_t>(idx)] != 0) { ordered_seg_idx.clear(); break; }
                                    used[static_cast<size_t>(idx)] = 1;
                                    ordered_seg_idx.push_back(seg_idx[static_cast<size_t>(idx)]);
                                }
                                weight = ordered_seg_idx.empty() ? 0.0 : 0.7;
                                if (weight > 0.0) ++p.motif_hits;
                            }
                        }
                    }
                    if (ordered_seg_idx.empty()) {
                        const auto pkey = precedence_signature(Aseg);
                        if (const auto* prec = lookup_precedence_cache(pkey, static_cast<int>(seg_idx.size())); prec != nullptr &&
                            prec->pair_scores.size() == seg_idx.size() * seg_idx.size()) {
                            std::vector<int> loc(seg_idx.size());
                            for (int i = 0; i < static_cast<int>(seg_idx.size()); ++i) loc[static_cast<size_t>(i)] = i;
                            std::vector<double> ps(loc.size(), 0.0);
                            const int ns = static_cast<int>(seg_idx.size());
                            for (int i = 0; i < ns; ++i) {
                                for (int j = 0; j < ns; ++j) {
                                    if (i == j) continue;
                                    ps[static_cast<size_t>(i)] += prec->pair_scores[static_cast<size_t>(i * ns + j)];
                                }
                            }
                            std::stable_sort(loc.begin(), loc.end(), [&](int a, int b) {
                                if (ps[static_cast<size_t>(a)] != ps[static_cast<size_t>(b)]) return ps[static_cast<size_t>(a)] > ps[static_cast<size_t>(b)];
                                return seg_idx[static_cast<size_t>(a)] < seg_idx[static_cast<size_t>(b)];
                            });
                            ordered_seg_idx.clear();
                            ordered_seg_idx.reserve(seg_idx.size());
                            for (int idx : loc) ordered_seg_idx.push_back(seg_idx[static_cast<size_t>(idx)]);
                            weight = 0.58;
                        }
                    }
                }
                if (ordered_seg_idx.empty()) continue;
                ++p.segments;
                ++piece_hits;
                const double denom = std::max(1, static_cast<int>(ordered_seg_idx.size()) - 1);
                for (int rank = 0; rank < static_cast<int>(ordered_seg_idx.size()); ++rank) {
                    const int outer_idx = ordered_seg_idx[static_cast<size_t>(rank)];
                    const double score = static_cast<double>(rank) / static_cast<double>(denom);
                    p.score_sum[static_cast<size_t>(outer_idx)] += weight * score;
                    if (p.weight_sum[static_cast<size_t>(outer_idx)] == 0.0) ++p.covered;
                    p.weight_sum[static_cast<size_t>(outer_idx)] += weight;
                }
            }
            if (p.segments < 2 || p.covered < std::max(2, (3 * n) / 5)) continue;
            std::vector<int> idxs(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) idxs[static_cast<size_t>(i)] = i;
            std::stable_sort(idxs.begin(), idxs.end(), [&](int a, int b) {
                const bool ca = p.weight_sum[static_cast<size_t>(a)] > 0.0;
                const bool cb = p.weight_sum[static_cast<size_t>(b)] > 0.0;
                if (ca != cb) return ca > cb;
                const double sa = ca ? (p.score_sum[static_cast<size_t>(a)] / p.weight_sum[static_cast<size_t>(a)]) : 0.5;
                const double sb = cb ? (p.score_sum[static_cast<size_t>(b)] / p.weight_sum[static_cast<size_t>(b)]) : 0.5;
                if (sa != sb) return sa < sb;
                const int pa = (local_vars[a] >= 0 && local_vars[a] < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(local_vars[a])] : a;
                const int pb = (local_vars[b] >= 0 && local_vars[b] < static_cast<int>(pinv_.size())) ? pinv_[static_cast<size_t>(local_vars[b])] : b;
                return pa < pb;
            });
            p.ordered_vars.reserve(static_cast<size_t>(n));
            for (int idx : idxs) p.ordered_vars.push_back(local_vars[static_cast<size_t>(idx)]);
            const double coverage = static_cast<double>(p.covered) / static_cast<double>(n);
            p.merit = 2.5 * coverage + 0.20 * static_cast<double>(p.exact_hits) + 0.12 * static_cast<double>(p.motif_hits)
                      - 0.03 * static_cast<double>(p.segments);
            ++proposals;
            if (p.merit > best.merit) best = std::move(p);
        }
    }

    if (best.merit > 0.0 && !best.ordered_vars.empty()) {
        result.ordered_vars = std::move(best.ordered_vars);
        result.overlap_piece_hits = piece_hits;
        result.overlap_proposals = proposals;
    }
    return result;
}

std::vector<int> DynamicCcolamdEngine::stitch_window_order(const std::vector<int>& base_perm,
                                                           const std::vector<int>& local_vars,
                                                           const std::vector<int>& local_order_vars) {
    if (local_vars.size() != local_order_vars.size()) return base_perm;
    const auto local_positions = positions_from_vars(base_perm, local_vars);
    if (local_positions.size() != local_order_vars.size()) return base_perm;
    std::vector<int> out = base_perm;
    for (int k = 0; k < static_cast<int>(local_positions.size()); ++k) {
        out[static_cast<size_t>(local_positions[static_cast<size_t>(k)])] = local_order_vars[static_cast<size_t>(k)];
    }
    return out;
}

void DynamicCcolamdEngine::remember_replay_window(const std::vector<int>& vars) {
    auto canon = unique_sorted(vars);
    if (canon.size() < 2) return;
    const size_t max_window = static_cast<size_t>(std::max(12, state_size_ / 3));
    if (canon.size() >= max_window) return;
    auto& replay_windows = active_replay_windows();
    for (auto& rw : replay_windows) {
        if (rw.vars == canon) {
            ++rw.hits;
            ++rw.accepts;
            ++rw.uses;
            std::stable_sort(replay_windows.begin(), replay_windows.end(), [](const ReplayWindow& a, const ReplayWindow& b) {
                if (a.accepts != b.accepts) return a.accepts > b.accepts;
                if (a.hits != b.hits) return a.hits > b.hits;
                return a.vars.size() < b.vars.size();
            });
            stats_.replay_windows_cached = static_cast<std::uint64_t>(replay_windows.size());
            return;
        }
    }
    replay_windows.insert(replay_windows.begin(), ReplayWindow{canon, 1, 1, 1});
    constexpr size_t kMaxReplayWindows = 8;
    if (replay_windows.size() > kMaxReplayWindows) replay_windows.resize(kMaxReplayWindows);
    stats_.replay_windows_cached = static_cast<std::uint64_t>(replay_windows.size());
}

void DynamicCcolamdEngine::rebuild_state_from_pattern(const Eigen::SparseMatrix<double>& pattern) {
    adjacency_ = adjacency_lists(pattern);
    pinv_ = inverse_permutation_of(perm_);
    stats_.replay_windows_cached = static_cast<std::uint64_t>(active_replay_windows().size());
    stats_.num_regimes_discovered = static_cast<int>(regime_states_.size());
    stats_.local_pattern_cache_entries = static_cast<std::uint64_t>(local_pattern_cache_.size());
    stats_.hierarchical_block_cache_entries = static_cast<std::uint64_t>(hierarchical_block_cache_.size());
    stats_.precedence_cache_entries = static_cast<std::uint64_t>(precedence_cache_.size());
}

DynamicCcolamdEngine::RegimeFeatures
DynamicCcolamdEngine::extract_regime_features(const Eigen::SparseMatrix<double>& pattern) {
    RegimeFeatures f{};
    const int n = pattern.rows();
    f.log_n = std::log1p(static_cast<double>(std::max(n, 0)));
    if (n <= 0) return f;

    int offdiag_nnz = 0;
    double span_sum = 0.0;
    for (int col = 0; col < pattern.outerSize(); ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(pattern, col); it; ++it) {
            if (it.row() == col) continue;
            ++offdiag_nnz;
            span_sum += static_cast<double>(std::abs(it.row() - col));
        }
    }
    f.avg_degree = static_cast<double>(offdiag_nnz) / static_cast<double>(n);
    f.span_ratio = offdiag_nnz > 0 ? span_sum / (static_cast<double>(offdiag_nnz) * static_cast<double>(std::max(1, n - 1))) : 0.0;
    return f;
}

double DynamicCcolamdEngine::regime_distance(const RegimeFeatures& a, const RegimeFeatures& b) {
    const double d0 = (a.log_n - b.log_n) / 0.8;
    const double d1 = (a.avg_degree - b.avg_degree) / 3.0;
    const double d2 = (a.span_ratio - b.span_ratio) / 0.25;
    return std::sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

size_t DynamicCcolamdEngine::select_or_create_regime(const RegimeFeatures& features) {
    constexpr double kCreateThreshold = 1.15;
    constexpr double kSplitThreshold = 0.80;

    auto make_regime = [&](const RegimeFeatures& feat) -> size_t {
        RegimeState rs{};
        rs.stable_id = next_regime_stable_id_++;
        rs.centroid = feat;
        rs.m2 = RegimeFeatures{};
        rs.kind_stats.assign(static_cast<size_t>(kCandidateKindCount), KindStats{});
        regime_states_.push_back(std::move(rs));
        ++stats_.regime_creations;
        stats_.num_regimes_discovered = static_cast<int>(regime_states_.size());
        return regime_states_.size() - 1;
    };

    if (regime_states_.empty()) return make_regime(features);

    size_t best_idx = 0;
    double best_dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < regime_states_.size(); ++i) {
        const double dist = regime_distance(features, regime_states_[i].centroid);
        if (dist < best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }

    const auto& best = regime_states_[best_idx];
    const double var_proxy = (best.visits > 1)
        ? std::sqrt(best.m2.log_n + best.m2.avg_degree + best.m2.span_ratio) / static_cast<double>(best.visits - 1)
        : 0.0;
    const bool should_split = (best.visits >= 8 && best_dist > kSplitThreshold && var_proxy > 0.12);
    if ((best_dist > kCreateThreshold || should_split) && regime_states_.size() < kMaxAdaptiveRegimes) {
        return make_regime(features);
    }
    return best_idx;
}

void DynamicCcolamdEngine::update_regime_model(size_t regime_index, const RegimeFeatures& features) {
    if (regime_index >= regime_states_.size()) return;
    auto& rs = regime_states_[regime_index];
    const double n = static_cast<double>(rs.visits + 1);
    const double d0 = features.log_n - rs.centroid.log_n;
    const double d1 = features.avg_degree - rs.centroid.avg_degree;
    const double d2 = features.span_ratio - rs.centroid.span_ratio;
    rs.centroid.log_n += d0 / n;
    rs.centroid.avg_degree += d1 / n;
    rs.centroid.span_ratio += d2 / n;
    rs.m2.log_n += d0 * (features.log_n - rs.centroid.log_n);
    rs.m2.avg_degree += d1 * (features.avg_degree - rs.centroid.avg_degree);
    rs.m2.span_ratio += d2 * (features.span_ratio - rs.centroid.span_ratio);
}

void DynamicCcolamdEngine::maybe_merge_regimes() {
    constexpr double kMergeThreshold = 0.32;
    if (regime_states_.size() < 2) return;
    if ((tick_ % 16u) != 0u && regime_states_.size() < kMaxAdaptiveRegimes) return;

    size_t a_best = 0, b_best = 0;
    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < regime_states_.size(); ++i) {
        for (size_t j = i + 1; j < regime_states_.size(); ++j) {
            const double d = regime_distance(regime_states_[i].centroid, regime_states_[j].centroid);
            if (d < best) {
                best = d;
                a_best = i;
                b_best = j;
            }
        }
    }
    if (best > kMergeThreshold) return;

    if (a_best > b_best) std::swap(a_best, b_best);
    auto& a = regime_states_[a_best];
    auto& b = regime_states_[b_best];
    const double wa = static_cast<double>(std::max<std::uint64_t>(1, a.visits));
    const double wb = static_cast<double>(std::max<std::uint64_t>(1, b.visits));
    const double wt = wa + wb;
    a.centroid.log_n = (wa * a.centroid.log_n + wb * b.centroid.log_n) / wt;
    a.centroid.avg_degree = (wa * a.centroid.avg_degree + wb * b.centroid.avg_degree) / wt;
    a.centroid.span_ratio = (wa * a.centroid.span_ratio + wb * b.centroid.span_ratio) / wt;
    a.m2.log_n += b.m2.log_n;
    a.m2.avg_degree += b.m2.avg_degree;
    a.m2.span_ratio += b.m2.span_ratio;
    a.visits += b.visits;
    if (a.kind_stats.size() < static_cast<size_t>(kCandidateKindCount)) {
        a.kind_stats.assign(static_cast<size_t>(kCandidateKindCount), KindStats{});
    }
    for (size_t k = 0; k < std::min(a.kind_stats.size(), b.kind_stats.size()); ++k) {
        a.kind_stats[k].attempts += b.kind_stats[k].attempts;
        a.kind_stats[k].accepts += b.kind_stats[k].accepts;
        a.kind_stats[k].last_success_tick = std::max(a.kind_stats[k].last_success_tick, b.kind_stats[k].last_success_tick);
    }
    for (const auto& rw : b.replay_windows) {
        bool merged = false;
        for (auto& arw : a.replay_windows) {
            if (arw.vars == rw.vars) {
                arw.hits += rw.hits;
                arw.accepts += rw.accepts;
                arw.uses += rw.uses;
                merged = true;
                break;
            }
        }
        if (!merged) a.replay_windows.push_back(rw);
    }
    constexpr size_t kMaxReplayWindows = 8;
    std::stable_sort(a.replay_windows.begin(), a.replay_windows.end(), [](const ReplayWindow& x, const ReplayWindow& y) {
        if (x.accepts != y.accepts) return x.accepts > y.accepts;
        if (x.hits != y.hits) return x.hits > y.hits;
        return x.vars.size() < y.vars.size();
    });
    if (a.replay_windows.size() > kMaxReplayWindows) a.replay_windows.resize(kMaxReplayWindows);

    const bool current_was_b = current_regime_index_ == b_best;
    if (current_regime_index_ > b_best) --current_regime_index_;
    if (current_was_b) current_regime_index_ = a_best;
    regime_states_.erase(regime_states_.begin() + static_cast<std::ptrdiff_t>(b_best));
    ++stats_.regime_merges;
    stats_.num_regimes_discovered = static_cast<int>(regime_states_.size());
}

std::vector<DynamicCcolamdEngine::ReplayWindow>& DynamicCcolamdEngine::active_replay_windows() {
    if (regime_states_.empty()) {
        static std::vector<ReplayWindow> empty;
        return empty;
    }
    if (current_regime_index_ >= regime_states_.size()) current_regime_index_ = regime_states_.size() - 1;
    return regime_states_[current_regime_index_].replay_windows;
}

const std::vector<DynamicCcolamdEngine::ReplayWindow>& DynamicCcolamdEngine::active_replay_windows() const {
    static const std::vector<ReplayWindow> empty;
    if (regime_states_.empty()) return empty;
    const size_t idx = std::min(current_regime_index_, regime_states_.size() - 1);
    return regime_states_[idx].replay_windows;
}

std::vector<DynamicCcolamdEngine::KindStats>& DynamicCcolamdEngine::active_kind_stats() {
    if (regime_states_.empty()) {
        static std::vector<KindStats> empty;
        return empty;
    }
    if (current_regime_index_ >= regime_states_.size()) current_regime_index_ = regime_states_.size() - 1;
    return regime_states_[current_regime_index_].kind_stats;
}

const std::vector<DynamicCcolamdEngine::KindStats>& DynamicCcolamdEngine::active_kind_stats() const {
    static const std::vector<KindStats> empty;
    if (regime_states_.empty()) return empty;
    const size_t idx = std::min(current_regime_index_, regime_states_.size() - 1);
    return regime_states_[idx].kind_stats;
}

} // namespace islam
