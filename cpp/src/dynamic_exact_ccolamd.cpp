#include "islam/dynamic_exact_ccolamd.hpp"
#include "islam/suitesparse_ccolamd_reference.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <string>

#ifdef ISLAM_HAS_CHOLMOD
#include <cholmod.h>
#endif

namespace islam {
namespace {

bool env_flag_enabled(const char* name) {
    const char* v = std::getenv(name);
    return v && std::string(v) != "0" && std::string(v) != "false" && std::string(v) != "FALSE";
}

bool strict_dynamic_exact_self_verify_enabled() {
    return env_flag_enabled("ISLAM_DYNAMIC_EXACT_SELF_VERIFY");
}

bool suitesparse_compat_check_enabled() {
    return env_flag_enabled("ISLAM_DYNAMIC_EXACT_CHECK_SUITESPARSE");
}

bool strict_suitesparse_compat_enabled() {
    return env_flag_enabled("ISLAM_DYNAMIC_EXACT_STRICT_SUITESPARSE_COMPAT");
}

void certify_against_reference(DynamicExactCcolamdPrototype::Stats& stats,
                               const std::vector<int>& produced,
                               const std::vector<int>& reference) {
    ++stats.exact_reference_checks;
    if (produced == reference) return;
    ++stats.exact_reference_failures;
    if (strict_dynamic_exact_self_verify_enabled()) {
        throw std::runtime_error("DynamicExactCcolamdPrototype: produced ordering differs from deterministic batch reference");
    }
}

int permutation_mismatch_positions(const std::vector<int>& a, const std::vector<int>& b) {
    const int n = std::min(static_cast<int>(a.size()), static_cast<int>(b.size()));
    int mismatches = std::abs(static_cast<int>(a.size()) - static_cast<int>(b.size()));
    for (int i = 0; i < n; ++i) {
        if (a[static_cast<size_t>(i)] != b[static_cast<size_t>(i)]) ++mismatches;
    }
    return mismatches;
}

bool live_checkpoint_equivalent(const DeterministicBatchCcolamd::LiveStateCheckpoint& a,
                                const DeterministicBatchCcolamd::LiveStateCheckpoint& b) {
    return a.step_after == b.step_after &&
           a.live_state_signature == b.live_state_signature &&
           a.prefix_signature == b.prefix_signature &&
           a.eliminated_prefix == b.eliminated_prefix &&
           a.live_variables == b.live_variables &&
           a.live_edges_upper_flat == b.live_edges_upper_flat;
}

std::uint64_t checkpoint_bank_key(const DeterministicBatchCcolamd::LiveStateCheckpoint& ck) {
    std::uint64_t h = 0xC0A11DCC01A4Dull;
    h ^= ck.live_state_signature + 0x9e3779b97f4a7c15ull + (h << 6u) + (h >> 2u);
    h ^= ck.prefix_signature + 0x9e3779b97f4a7c15ull + (h << 6u) + (h >> 2u);
    h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(ck.step_after)) +
         0x9e3779b97f4a7c15ull + (h << 6u) + (h >> 2u);
    h ^= static_cast<std::uint64_t>(ck.live_variables.size()) +
         0x9e3779b97f4a7c15ull + (h << 6u) + (h >> 2u);
    h ^= static_cast<std::uint64_t>(ck.live_edges_upper_flat.size()) +
         0x9e3779b97f4a7c15ull + (h << 6u) + (h >> 2u);
    return h;
}


#ifdef ISLAM_HAS_CHOLMOD
using ss_long = SuiteSparse_long;

bool is_valid_permutation_vector(const std::vector<int>& perm, int n) {
    if (static_cast<int>(perm.size()) != n) return false;
    std::vector<unsigned char> seen(static_cast<size_t>(n), 0);
    for (int v : perm) {
        if (v < 0 || v >= n || seen[static_cast<size_t>(v)]) return false;
        seen[static_cast<size_t>(v)] = 1;
    }
    return true;
}

cholmod_sparse* eigen_pattern_to_cholmod_sparse(const Eigen::SparseMatrix<double>& in,
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
        throw std::runtime_error("DynamicExactCcolamdPrototype: CHOLMOD sparse allocation failed");
    }

    auto* p = static_cast<ss_long*>(out->p);
    auto* i = static_cast<ss_long*>(out->i);
    auto* x = static_cast<double*>(out->x);
    for (int col = 0; col <= A.cols(); ++col) {
        p[col] = static_cast<ss_long>(A.outerIndexPtr()[col]);
    }
    for (int nz = 0; nz < A.nonZeros(); ++nz) {
        i[nz] = static_cast<ss_long>(A.innerIndexPtr()[nz]);
        x[nz] = 1.0;
    }
    out->sorted = 1;
    out->packed = 1;
    return out;
}

std::vector<int> extract_cholmod_perm(const cholmod_factor* factor, int n) {
    if (factor == nullptr) return {};
    if (factor->Perm == nullptr) {
        std::vector<int> id(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) id[static_cast<size_t>(i)] = i;
        return id;
    }
    std::vector<int> perm(static_cast<size_t>(n), -1);
    const auto* p = static_cast<const ss_long*>(factor->Perm);
    for (int i = 0; i < n; ++i) perm[static_cast<size_t>(i)] = static_cast<int>(p[i]);
    if (!is_valid_permutation_vector(perm, n)) return {};
    return perm;
}

std::vector<int> suitesparse_cholmod_reference_permutation(const Eigen::SparseMatrix<double>& pattern) {
    const int n = pattern.rows();
    if (n <= 0 || pattern.cols() != n) return {};

    cholmod_common common{};
    cholmod_l_start(&common);
    common.nmethods = 1;
#if defined(CHOLMOD_COLAMD)
    common.method[0].ordering = CHOLMOD_COLAMD;
#elif defined(CHOLMOD_AMD)
    common.method[0].ordering = CHOLMOD_AMD;
#else
    common.method[0].ordering = CHOLMOD_NATURAL;
#endif
    common.postorder = 0;
    common.supernodal = CHOLMOD_SIMPLICIAL;

    cholmod_sparse* A = nullptr;
    cholmod_factor* F = nullptr;
    std::vector<int> perm;
    try {
        A = eigen_pattern_to_cholmod_sparse(pattern, 1, &common);
        F = cholmod_l_analyze(A, &common);
        perm = extract_cholmod_perm(F, n);
    } catch (...) {
        perm.clear();
    }

    if (F != nullptr) cholmod_l_free_factor(&F, &common);
    if (A != nullptr) cholmod_l_free_sparse(&A, &common);
    cholmod_l_finish(&common);
    return perm;
}
#endif


void append_dynamic_exact_trace_event(const DynamicExactCcolamdPrototype::Stats& stats,
                                      const Eigen::SparseMatrix<double>& pattern,
                                      const std::vector<int>& dirty_vars,
                                      const char* outcome,
                                      int reused_prefix_size) {
    const char* path = std::getenv("ISLAM_DYNAMIC_EXACT_TRACE_CSV");
    if (path == nullptr || std::string(path).empty()) return;

    bool need_header = true;
    {
        std::ifstream probe(path);
        need_header = !probe.good() || probe.peek() == std::ifstream::traits_type::eof();
    }

    std::ofstream os(path, std::ios::app);
    if (!os) {
        throw std::runtime_error("DynamicExactCcolamdPrototype: failed to open ISLAM_DYNAMIC_EXACT_TRACE_CSV for append");
    }
    if (need_header) {
        os << "refresh,outcome,n,pattern_signature,dirty_count,reused_prefix_size,"
              "common_prefix,dirty_boundary,reusable_checkpoint,checkpoint_live_vars,checkpoint_live_edges,rollback_suffix,"
              "current_pattern_checkpoint_certification_attempts,current_pattern_checkpoint_certification_successes,current_pattern_checkpoint_certification_failures,"
              "structural_checkpoint_certification_attempts,structural_checkpoint_certification_successes,structural_checkpoint_certification_failures,"
              "checkpoint_bank_entries,checkpoint_bank_insertions,checkpoint_bank_probes,checkpoint_bank_hits,checkpoint_bank_misses,checkpoint_bank_imported_entries,checkpoint_bank_import_duplicate_entries,checkpoint_bank_import_invalid_entries,checkpoint_bank_import_header_checks,checkpoint_bank_import_header_failures,checkpoint_bank_import_digest_mismatches,checkpoint_bank_import_entry_count_mismatches,checkpoint_bank_exported_entries,checkpoint_bank_export_writes,"
              "checkpoint_resume_attempts,checkpoint_resume_successes,checkpoint_resume_failures,"
              "reference_free_resume_attempts,reference_free_resume_candidates,reference_free_resume_failures,"
              "reference_free_resume_reference_matches,reference_free_resume_reference_mismatches,"
              "dirty_boundary_safety_checks,dirty_boundary_safe_reuses,dirty_boundary_unsafe_overestimates,"
              "dirty_boundary_underestimates,dirty_boundary_exact_matches,dirty_boundary_candidate_pivots,"
              "dirty_boundary_safe_pivots,dirty_boundary_unsafe_pivots,last_dirty_boundary_overestimate,last_dirty_boundary_underestimate,"
              "dirty_boundary_checkpoint_probe_attempts,dirty_boundary_checkpoint_probe_structural_successes,dirty_boundary_checkpoint_probe_structural_failures,"
              "dirty_boundary_checkpoint_probe_resume_matches,dirty_boundary_checkpoint_probe_resume_mismatches,dirty_boundary_checkpoint_probe_resume_failures,last_dirty_boundary_checkpoint_probe_pivots,"
              "certified_suffix_replays,certified_suffix_replay_failures,"
              "suffix_replay_pivots_reused,suffix_replay_pivots_recomputed,"
              "exact_reference_checks,exact_reference_failures,"
              "suitesparse_compat_checks,suitesparse_compat_failures,"
              "suitesparse_compat_unavailable,suitesparse_compat_mismatch_positions,"
              "last_suitesparse_compat_mismatches\n";
    }
    os << stats.refreshes << ','
       << outcome << ','
       << pattern.rows() << ','
       << deterministic_ccolamd_pattern_signature(pattern) << ','
       << dirty_vars.size() << ','
       << reused_prefix_size << ','
       << stats.last_common_prefix << ','
       << stats.last_dirty_boundary << ','
       << stats.last_reusable_checkpoint << ','
       << stats.last_checkpoint_live_vars << ','
       << stats.last_checkpoint_live_edges << ','
       << stats.last_rollback_suffix << ','
       << stats.current_pattern_checkpoint_certification_attempts << ','
       << stats.current_pattern_checkpoint_certification_successes << ','
       << stats.current_pattern_checkpoint_certification_failures << ','
       << stats.structural_checkpoint_certification_attempts << ','
       << stats.structural_checkpoint_certification_successes << ','
       << stats.structural_checkpoint_certification_failures << ','
       << stats.checkpoint_bank_entries << ','
       << stats.checkpoint_bank_insertions << ','
       << stats.checkpoint_bank_probes << ','
       << stats.checkpoint_bank_hits << ','
       << stats.checkpoint_bank_misses << ','
       << stats.checkpoint_bank_imported_entries << ','
       << stats.checkpoint_bank_import_duplicate_entries << ','
       << stats.checkpoint_bank_import_invalid_entries << ','
       << stats.checkpoint_bank_import_header_checks << ','
       << stats.checkpoint_bank_import_header_failures << ','
       << stats.checkpoint_bank_import_digest_mismatches << ','
       << stats.checkpoint_bank_import_entry_count_mismatches << ','
       << stats.checkpoint_bank_exported_entries << ','
       << stats.checkpoint_bank_export_writes << ','
       << stats.checkpoint_resume_attempts << ','
       << stats.checkpoint_resume_successes << ','
       << stats.checkpoint_resume_failures << ','
       << stats.reference_free_resume_attempts << ','
       << stats.reference_free_resume_candidates << ','
       << stats.reference_free_resume_failures << ','
       << stats.reference_free_resume_reference_matches << ','
       << stats.reference_free_resume_reference_mismatches << ','
       << stats.dirty_boundary_safety_checks << ','
       << stats.dirty_boundary_safe_reuses << ','
       << stats.dirty_boundary_unsafe_overestimates << ','
       << stats.dirty_boundary_underestimates << ','
       << stats.dirty_boundary_exact_matches << ','
       << stats.dirty_boundary_candidate_pivots << ','
       << stats.dirty_boundary_safe_pivots << ','
       << stats.dirty_boundary_unsafe_pivots << ','
       << stats.last_dirty_boundary_overestimate << ','
       << stats.last_dirty_boundary_underestimate << ','
       << stats.dirty_boundary_checkpoint_probe_attempts << ','
       << stats.dirty_boundary_checkpoint_probe_structural_successes << ','
       << stats.dirty_boundary_checkpoint_probe_structural_failures << ','
       << stats.dirty_boundary_checkpoint_probe_resume_matches << ','
       << stats.dirty_boundary_checkpoint_probe_resume_mismatches << ','
       << stats.dirty_boundary_checkpoint_probe_resume_failures << ','
       << stats.last_dirty_boundary_checkpoint_probe_pivots << ','
       << stats.certified_suffix_replays << ','
       << stats.certified_suffix_replay_failures << ','
       << stats.suffix_replay_pivots_reused << ','
       << stats.suffix_replay_pivots_recomputed << ','
       << stats.exact_reference_checks << ','
       << stats.exact_reference_failures << ','
       << stats.suitesparse_compat_checks << ','
       << stats.suitesparse_compat_failures << ','
       << stats.suitesparse_compat_unavailable << ','
       << stats.suitesparse_compat_mismatch_positions << ','
       << stats.last_suitesparse_compat_mismatches << '\n';
}

void audit_suitesparse_compatibility(DynamicExactCcolamdPrototype::Stats& stats,
                                     const std::vector<int>& produced,
                                     const Eigen::SparseMatrix<double>& pattern) {
    stats.last_suitesparse_compat_mismatches = 0;
    if (!suitesparse_compat_check_enabled()) return;

    const auto ref = suitesparse_ccolamd_reference_order(pattern);
    if (!ref.available || !ref.ok || ref.permutation.empty()) {
        ++stats.suitesparse_compat_unavailable;
        if (strict_suitesparse_compat_enabled()) {
            std::string why = ref.message.empty() ? "SuiteSparse CCOLAMD direct reference unavailable" : ref.message;
            throw std::runtime_error("DynamicExactCcolamdPrototype: " + why);
        }
        return;
    }

    ++stats.suitesparse_compat_checks;
    const int mismatches = permutation_mismatch_positions(produced, ref.permutation);
    stats.last_suitesparse_compat_mismatches = mismatches;
    stats.suitesparse_compat_mismatch_positions += static_cast<std::uint64_t>(mismatches);
    if (mismatches != 0) {
        ++stats.suitesparse_compat_failures;
        if (strict_suitesparse_compat_enabled()) {
            throw std::runtime_error("DynamicExactCcolamdPrototype: ordering differs from direct SuiteSparse CCOLAMD reference");
        }
    }
}

} // namespace

void DynamicExactCcolamdPrototype::clear() {
    checkpoint_bank_loaded_ = false;
    current_trace_ = DeterministicBatchCcolamd::Trace{};
    last_reused_prefix_.clear();
    stats_ = Stats{};
    checkpoint_bank_.clear();
    has_state_ = false;
}

namespace {

void store_trace_checkpoints(DynamicExactCcolamdPrototype::Stats& stats,
                             DynamicExactCcolamdPrototype::CheckpointBank& bank,
                             const DeterministicBatchCcolamd::Trace& trace) {
    for (const auto& ck : trace.live_state_checkpoints) {
        const auto key = checkpoint_bank_key(ck);
        auto& bucket = bank[key];
        bool exists = false;
        for (const auto& old : bucket) {
            if (live_checkpoint_equivalent(old, ck)) {
                exists = true;
                break;
            }
        }
        if (!exists) {
            bucket.push_back(ck);
            ++stats.checkpoint_bank_insertions;
        }
    }
    std::uint64_t entries = 0;
    for (const auto& kv : bank) entries += static_cast<std::uint64_t>(kv.second.size());
    stats.checkpoint_bank_entries = entries;
}

bool lookup_checkpoint_bank(DynamicExactCcolamdPrototype::Stats& stats,
                            const DynamicExactCcolamdPrototype::CheckpointBank& bank,
                            const DeterministicBatchCcolamd::LiveStateCheckpoint& target,
                            DeterministicBatchCcolamd::LiveStateCheckpoint& out) {
    ++stats.checkpoint_bank_probes;
    const auto key = checkpoint_bank_key(target);
    const auto it = bank.find(key);
    if (it != bank.end()) {
        for (const auto& ck : it->second) {
            if (live_checkpoint_equivalent(ck, target)) {
                out = ck;
                ++stats.checkpoint_bank_hits;
                return true;
            }
        }
    }
    ++stats.checkpoint_bank_misses;
    return false;
}

bool structurally_certify_checkpoint_against_pattern(DeterministicBatchCcolamd& batch,
                                                     const Eigen::SparseMatrix<double>& pattern,
                                                     const DeterministicBatchCcolamd::LiveStateCheckpoint& checkpoint) {
    // Reference-free checkpoint validity test: force only the checkpoint prefix
    // on the current pattern, then compare the materialized live symbolic state
    // at that prefix. This does not compare against the current full batch
    // ordering suffix and is therefore a stepping stone toward a truly
    // reference-free dynamic ordering backend.
    auto forced_trace = batch.order_with_forced_prefix(pattern, checkpoint.eliminated_prefix);
    for (const auto& ck : forced_trace.live_state_checkpoints) {
        if (ck.eliminated_prefix == checkpoint.eliminated_prefix) {
            return live_checkpoint_equivalent(ck, checkpoint);
        }
    }
    return false;
}


std::string checkpoint_bank_import_path() {
    const char* p = std::getenv("ISLAM_DYNAMIC_EXACT_CHECKPOINT_BANK_IN");
    return p ? std::string(p) : std::string();
}

std::string checkpoint_bank_export_path() {
    const char* p = std::getenv("ISLAM_DYNAMIC_EXACT_CHECKPOINT_BANK_OUT");
    return p ? std::string(p) : std::string();
}

bool strict_checkpoint_bank_import_enabled() {
    return env_flag_enabled("ISLAM_DYNAMIC_EXACT_CHECKPOINT_BANK_STRICT_IMPORT");
}

bool require_canonical_checkpoint_bank_import_enabled() {
    return env_flag_enabled("ISLAM_DYNAMIC_EXACT_CHECKPOINT_BANK_REQUIRE_CANONICAL");
}

void write_int_vector(std::ostream& os, const std::vector<int>& values) {
    os << values.size();
    for (int v : values) os << ' ' << v;
}

bool read_int_vector(std::istream& is, std::vector<int>& values) {
    size_t n = 0;
    if (!(is >> n)) return false;
    values.clear();
    values.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        int v = 0;
        if (!(is >> v)) return false;
        values.push_back(v);
    }
    return true;
}

std::uint64_t checkpoint_mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30u)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27u)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31u);
}

std::uint64_t checkpoint_prefix_signature(const std::vector<int>& prefix) {
    std::uint64_t h = checkpoint_mix64(static_cast<std::uint64_t>(prefix.size()));
    for (int v : prefix) {
        h = checkpoint_mix64(h ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
    }
    return h;
}

bool is_strictly_increasing_nonnegative(const std::vector<int>& values) {
    int prev = -1;
    for (int v : values) {
        if (v < 0 || v <= prev) return false;
        prev = v;
    }
    return true;
}

bool checkpoint_vectors_valid(const DeterministicBatchCcolamd::LiveStateCheckpoint& ck) {
    if (ck.step_after < 0) return false;
    // step_after is the zero-based elimination-event index of the checkpoint,
    // not the number of eliminated variables. Absorption can eliminate more
    // than one variable in one event, so the prefix length may exceed
    // step_after + 1.
    if (ck.eliminated_prefix.empty()) return false;
    if (ck.step_after >= static_cast<int>(ck.eliminated_prefix.size())) return false;
    std::vector<int> prefix = ck.eliminated_prefix;
    std::sort(prefix.begin(), prefix.end());
    if (std::adjacent_find(prefix.begin(), prefix.end()) != prefix.end()) return false;
    for (int v : prefix) if (v < 0) return false;
    if (!is_strictly_increasing_nonnegative(ck.live_variables)) return false;
    for (int v : ck.live_variables) {
        if (std::binary_search(prefix.begin(), prefix.end(), v)) return false;
    }
    if (ck.live_edges_upper_flat.size() % 2 != 0) return false;
    int prev_a = -1;
    int prev_b = -1;
    for (size_t k = 0; k < ck.live_edges_upper_flat.size(); k += 2) {
        const int a = ck.live_edges_upper_flat[k];
        const int b = ck.live_edges_upper_flat[k + 1];
        if (a < 0 || b < 0 || a >= b) return false;
        if (!std::binary_search(ck.live_variables.begin(), ck.live_variables.end(), a)) return false;
        if (!std::binary_search(ck.live_variables.begin(), ck.live_variables.end(), b)) return false;
        if (prev_a > a || (prev_a == a && prev_b >= b)) return false;
        prev_a = a;
        prev_b = b;
    }
    if (checkpoint_prefix_signature(ck.eliminated_prefix) != ck.prefix_signature) return false;
    return true;
}

bool checkpoint_bank_digest_less(const DeterministicBatchCcolamd::LiveStateCheckpoint& a,
                                 const DeterministicBatchCcolamd::LiveStateCheckpoint& b) {
    if (a.step_after != b.step_after) return a.step_after < b.step_after;
    if (a.live_state_signature != b.live_state_signature) return a.live_state_signature < b.live_state_signature;
    if (a.prefix_signature != b.prefix_signature) return a.prefix_signature < b.prefix_signature;
    if (a.eliminated_prefix != b.eliminated_prefix) return a.eliminated_prefix < b.eliminated_prefix;
    if (a.live_variables != b.live_variables) return a.live_variables < b.live_variables;
    return a.live_edges_upper_flat < b.live_edges_upper_flat;
}

struct CheckpointBankHeader {
    bool present = false;
    std::uint64_t declared_entries = 0;
    std::uint64_t declared_digest64 = 0;
};

std::uint64_t checkpoint_bank_canonical_digest(const std::vector<DeterministicBatchCcolamd::LiveStateCheckpoint>& entries);

bool parse_checkpoint_bank_canonical_header(const std::string& line, CheckpointBankHeader& header) {
    std::istringstream ss(line);
    std::string hash;
    std::string tag;
    std::string entries_label;
    std::string digest_label;
    std::uint64_t entries = 0;
    std::uint64_t digest = 0;
    if (!(ss >> hash >> tag >> entries_label >> entries >> digest_label >> digest)) return false;
    if (hash != "#" || tag != "canonical_dynamic_exact_checkpoint_bank_v1" ||
        entries_label != "entries" || digest_label != "digest64") {
        return false;
    }
    std::string trailing;
    if (ss >> trailing) return false;
    header.present = true;
    header.declared_entries = entries;
    header.declared_digest64 = digest;
    return true;
}


void load_checkpoint_bank_once(DynamicExactCcolamdPrototype::Stats& stats,
                               DynamicExactCcolamdPrototype::CheckpointBank& bank,
                               bool& loaded_flag) {
    if (loaded_flag) return;
    loaded_flag = true;

    const std::string path = checkpoint_bank_import_path();
    if (path.empty()) return;

    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("DynamicExactCcolamdPrototype: failed to open checkpoint-bank import file");
    }

    std::string line;
    std::uint64_t imported = 0;
    std::uint64_t duplicates = 0;
    std::uint64_t invalid = 0;
    std::uint64_t data_rows = 0;
    CheckpointBankHeader header;
    std::vector<DeterministicBatchCcolamd::LiveStateCheckpoint> valid_input_entries;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            CheckpointBankHeader parsed;
            if (parse_checkpoint_bank_canonical_header(line, parsed)) {
                header = parsed;
            }
            continue;
        }
        ++data_rows;
        std::istringstream ss(line);
        DeterministicBatchCcolamd::LiveStateCheckpoint ck;
        if (!(ss >> ck.step_after >> ck.live_state_signature >> ck.prefix_signature)) {
            throw std::runtime_error("DynamicExactCcolamdPrototype: malformed checkpoint-bank import header");
        }
        if (!read_int_vector(ss, ck.eliminated_prefix) ||
            !read_int_vector(ss, ck.live_variables) ||
            !read_int_vector(ss, ck.live_edges_upper_flat)) {
            throw std::runtime_error("DynamicExactCcolamdPrototype: malformed checkpoint-bank import vectors");
        }

        if (!checkpoint_vectors_valid(ck)) {
            ++invalid;
            if (strict_checkpoint_bank_import_enabled()) {
                throw std::runtime_error("DynamicExactCcolamdPrototype: invalid checkpoint-bank import entry");
            }
            continue;
        }

        valid_input_entries.push_back(ck);
        const auto key = checkpoint_bank_key(ck);
        auto& bucket = bank[key];
        bool exists = false;
        for (const auto& old : bucket) {
            if (live_checkpoint_equivalent(old, ck)) {
                exists = true;
                break;
            }
        }
        if (!exists) {
            bucket.push_back(std::move(ck));
            ++imported;
        } else {
            ++duplicates;
        }
    }

    if (header.present || require_canonical_checkpoint_bank_import_enabled()) {
        ++stats.checkpoint_bank_import_header_checks;
        bool header_failed = false;
        if (!header.present) {
            header_failed = true;
        } else {
            if (header.declared_entries != data_rows) {
                ++stats.checkpoint_bank_import_entry_count_mismatches;
                header_failed = true;
            }
            std::sort(valid_input_entries.begin(), valid_input_entries.end(), checkpoint_bank_digest_less);
            const std::uint64_t digest = checkpoint_bank_canonical_digest(valid_input_entries);
            if (digest != header.declared_digest64) {
                ++stats.checkpoint_bank_import_digest_mismatches;
                header_failed = true;
            }
        }
        if (header_failed) {
            ++stats.checkpoint_bank_import_header_failures;
            if (strict_checkpoint_bank_import_enabled() || require_canonical_checkpoint_bank_import_enabled()) {
                throw std::runtime_error("DynamicExactCcolamdPrototype: checkpoint-bank canonical import header validation failed");
            }
        }
    }

    stats.checkpoint_bank_imported_entries += imported;
    stats.checkpoint_bank_import_duplicate_entries += duplicates;
    stats.checkpoint_bank_import_invalid_entries += invalid;
    std::uint64_t entries = 0;
    for (const auto& kv : bank) entries += static_cast<std::uint64_t>(kv.second.size());
    stats.checkpoint_bank_entries = entries;
}

bool checkpoint_bank_export_less(const DeterministicBatchCcolamd::LiveStateCheckpoint& a,
                                  const DeterministicBatchCcolamd::LiveStateCheckpoint& b) {
    if (a.step_after != b.step_after) return a.step_after < b.step_after;
    if (a.live_state_signature != b.live_state_signature) return a.live_state_signature < b.live_state_signature;
    if (a.prefix_signature != b.prefix_signature) return a.prefix_signature < b.prefix_signature;
    if (a.eliminated_prefix != b.eliminated_prefix) return a.eliminated_prefix < b.eliminated_prefix;
    if (a.live_variables != b.live_variables) return a.live_variables < b.live_variables;
    return a.live_edges_upper_flat < b.live_edges_upper_flat;
}

std::uint64_t checkpoint_bank_canonical_digest(const std::vector<DeterministicBatchCcolamd::LiveStateCheckpoint>& entries) {
    std::uint64_t h = 0xD1A6E57C0FFEE48ull;
    for (const auto& ck : entries) {
        h = checkpoint_mix64(h ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(ck.step_after)));
        h = checkpoint_mix64(h ^ ck.live_state_signature);
        h = checkpoint_mix64(h ^ ck.prefix_signature);
        for (int v : ck.eliminated_prefix) h = checkpoint_mix64(h ^ 0xE1000000ull ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
        for (int v : ck.live_variables) h = checkpoint_mix64(h ^ 0xE2000000ull ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
        for (int v : ck.live_edges_upper_flat) h = checkpoint_mix64(h ^ 0xE3000000ull ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
    }
    return h;
}

void export_checkpoint_bank_if_requested(DynamicExactCcolamdPrototype::Stats& stats,
                                         const DynamicExactCcolamdPrototype::CheckpointBank& bank) {
    const std::string path = checkpoint_bank_export_path();
    if (path.empty()) return;

    std::vector<DeterministicBatchCcolamd::LiveStateCheckpoint> entries;
    for (const auto& kv : bank) {
        entries.insert(entries.end(), kv.second.begin(), kv.second.end());
    }
    std::sort(entries.begin(), entries.end(), checkpoint_bank_export_less);
    entries.erase(std::unique(entries.begin(), entries.end(), live_checkpoint_equivalent), entries.end());

    std::ofstream out(path, std::ios::trunc);
    if (!out) {
        throw std::runtime_error("DynamicExactCcolamdPrototype: failed to open checkpoint-bank export file");
    }
    out << "# canonical_dynamic_exact_checkpoint_bank_v1 entries " << entries.size()
        << " digest64 " << checkpoint_bank_canonical_digest(entries) << "\n";
    out << "# step_after live_state_signature prefix_signature eliminated_prefix live_variables live_edges_upper_flat\n";
    std::uint64_t exported = 0;
    for (const auto& ck : entries) {
        out << ck.step_after << ' ' << ck.live_state_signature << ' ' << ck.prefix_signature << ' ';
        write_int_vector(out, ck.eliminated_prefix);
        out << ' ';
        write_int_vector(out, ck.live_variables);
        out << ' ';
        write_int_vector(out, ck.live_edges_upper_flat);
        out << '\n';
        ++exported;
    }
    stats.checkpoint_bank_exported_entries = exported;
    ++stats.checkpoint_bank_export_writes;
}
} // namespace

std::vector<int> DynamicExactCcolamdPrototype::refresh(const Eigen::SparseMatrix<double>& pattern,
                                                       const std::vector<int>& dirty_vars) {
    load_checkpoint_bank_once(stats_, checkpoint_bank_, checkpoint_bank_loaded_);
    ++stats_.refreshes;

    auto maybe_self_verify = [&](const std::vector<int>& produced) {
        if (!strict_dynamic_exact_self_verify_enabled()) return;
        auto reference_trace = batch_.order_with_trace(pattern);
        certify_against_reference(stats_, produced, reference_trace.permutation);
    };

    auto install_trace = [&](DeterministicBatchCcolamd::Trace trace,
                             const char* outcome,
                             int reused_prefix_size) -> std::vector<int> {
        current_trace_ = std::move(trace);
        has_state_ = true;
        store_trace_checkpoints(stats_, checkpoint_bank_, current_trace_);
        export_checkpoint_bank_if_requested(stats_, checkpoint_bank_);
        maybe_self_verify(current_trace_.permutation);
        audit_suitesparse_compatibility(stats_, current_trace_.permutation, pattern);
        append_dynamic_exact_trace_event(stats_, pattern, dirty_vars, outcome, reused_prefix_size);
        return current_trace_.permutation;
    };

    auto combine_prefix_and_suffix = [&](DeterministicBatchCcolamd::Trace prefix_trace,
                                         DeterministicBatchCcolamd::Trace suffix_trace) {
        DeterministicBatchCcolamd::Trace out;
        out.n = pattern.rows();
        out.pattern_signature = deterministic_ccolamd_pattern_signature(pattern);
        out.rule_signature = suffix_trace.rule_signature ? suffix_trace.rule_signature : prefix_trace.rule_signature;
        out.steps = std::move(prefix_trace.steps);
        out.prefix_checkpoints = std::move(prefix_trace.prefix_checkpoints);
        out.live_state_checkpoints = std::move(prefix_trace.live_state_checkpoints);
        const int prefix_events = static_cast<int>(out.steps.size());
        for (int k = 0; k < static_cast<int>(suffix_trace.steps.size()); ++k) {
            auto rec = suffix_trace.steps[static_cast<size_t>(k)];
            rec.step = prefix_events + k;
            out.steps.push_back(std::move(rec));
        }
        for (int k = 0; k < static_cast<int>(suffix_trace.prefix_checkpoints.size()); ++k) {
            out.prefix_checkpoints.push_back(std::move(suffix_trace.prefix_checkpoints[static_cast<size_t>(k)]));
        }
        for (int k = 0; k < static_cast<int>(suffix_trace.live_state_checkpoints.size()); ++k) {
            auto ck = suffix_trace.live_state_checkpoints[static_cast<size_t>(k)];
            ck.step_after = prefix_events + k;
            out.live_state_checkpoints.push_back(std::move(ck));
        }
        out.permutation = std::move(suffix_trace.permutation);
        return out;
    };

    if (!has_state_) {
        ++stats_.cold_starts;
        stats_.last_common_prefix = 0;
        stats_.last_dirty_boundary = 0;
        stats_.last_reusable_checkpoint = 0;
        stats_.last_checkpoint_live_vars = 0;
        stats_.last_checkpoint_live_edges = 0;
        last_reused_prefix_.clear();
        ++stats_.exact_suffix_recomputes;
        auto trace = batch_.order_with_trace(pattern);
        stats_.last_rollback_suffix = static_cast<int>(trace.steps.size());
        stats_.total_rollback_suffix += static_cast<std::uint64_t>(stats_.last_rollback_suffix);
        return install_trace(std::move(trace), "cold_start", 0);
    }

    ++stats_.prefix_reuse_attempts;
    const int dirty_boundary = dirty_trace_boundary(current_trace_, dirty_vars);
    const int max_reusable = std::min(dirty_boundary, static_cast<int>(current_trace_.live_state_checkpoints.size()));

    int certified_events = 0;
    DeterministicBatchCcolamd::LiveStateCheckpoint certified_checkpoint{};
    DeterministicBatchCcolamd::Trace certified_prefix_trace{};
    bool have_certified_checkpoint = false;

    if (max_reusable > 0) {
        ++stats_.checkpoint_replay_candidates;
    }

    for (int events = max_reusable; events > 0; --events) {
        ++stats_.materialized_checkpoint_reuse_attempts;
        ++stats_.structural_checkpoint_certification_attempts;
        const auto& candidate = current_trace_.live_state_checkpoints[static_cast<size_t>(events - 1)];
        DeterministicBatchCcolamd::Trace prefix_trace;
        bool ok = false;
        try {
            ok = batch_.certify_checkpoint_prefix(pattern, candidate, &prefix_trace);
        } catch (...) {
            ok = false;
        }
        if (ok) {
            ++stats_.structural_checkpoint_certification_successes;
            ++stats_.materialized_checkpoint_reuse_successes;
            certified_events = events;
            certified_checkpoint = candidate;
            certified_prefix_trace = std::move(prefix_trace);
            have_certified_checkpoint = true;
            break;
        }
        ++stats_.structural_checkpoint_certification_failures;
        ++stats_.materialized_checkpoint_reuse_failures;
    }

    stats_.last_common_prefix = certified_events;
    stats_.last_dirty_boundary = dirty_boundary;
    stats_.last_reusable_checkpoint = certified_events;
    stats_.last_dirty_boundary_overestimate = std::max(0, dirty_boundary - certified_events);
    stats_.last_dirty_boundary_underestimate = 0;
    ++stats_.dirty_boundary_safety_checks;
    stats_.dirty_boundary_candidate_pivots += static_cast<std::uint64_t>(std::max(0, dirty_boundary));
    if (certified_events == dirty_boundary) {
        ++stats_.dirty_boundary_safe_reuses;
        ++stats_.dirty_boundary_exact_matches;
        stats_.dirty_boundary_safe_pivots += static_cast<std::uint64_t>(std::max(0, certified_events));
    } else {
        ++stats_.dirty_boundary_unsafe_overestimates;
        stats_.dirty_boundary_safe_pivots += static_cast<std::uint64_t>(std::max(0, certified_events));
        stats_.dirty_boundary_unsafe_pivots += static_cast<std::uint64_t>(std::max(0, dirty_boundary - certified_events));
    }

    if (have_certified_checkpoint) {
        last_reused_prefix_ = certified_checkpoint.eliminated_prefix;
        stats_.last_checkpoint_live_vars = static_cast<int>(certified_checkpoint.live_variables.size());
        stats_.last_checkpoint_live_edges = static_cast<int>(certified_checkpoint.live_edges_upper_flat.size() / 2);
        ++stats_.checkpoint_resume_attempts;
        try {
            auto suffix_trace = batch_.order_from_checkpoint(certified_checkpoint);
            ++stats_.checkpoint_resume_successes;
            ++stats_.certified_suffix_replays;
            stats_.suffix_replay_pivots_reused += static_cast<std::uint64_t>(last_reused_prefix_.size());
            stats_.suffix_replay_pivots_recomputed +=
                static_cast<std::uint64_t>(std::max(0, static_cast<int>(suffix_trace.permutation.size()) -
                                                       static_cast<int>(last_reused_prefix_.size())));
            auto full_trace = combine_prefix_and_suffix(std::move(certified_prefix_trace), std::move(suffix_trace));
            stats_.last_rollback_suffix = std::max(0, static_cast<int>(full_trace.steps.size()) - certified_events);
            stats_.total_common_prefix += static_cast<std::uint64_t>(certified_events);
            stats_.total_rollback_suffix += static_cast<std::uint64_t>(stats_.last_rollback_suffix);
            if (certified_events > 0) ++stats_.prefix_reuse_successes;
            if (stats_.last_rollback_suffix == 0) ++stats_.checkpoint_replay_full_prefixes;
            return install_trace(std::move(full_trace), "certified_checkpoint_resume", static_cast<int>(last_reused_prefix_.size()));
        } catch (...) {
            ++stats_.checkpoint_resume_failures;
            ++stats_.certified_suffix_replay_failures;
            last_reused_prefix_.clear();
            stats_.last_checkpoint_live_vars = 0;
            stats_.last_checkpoint_live_edges = 0;
        }
    } else {
        last_reused_prefix_.clear();
        stats_.last_checkpoint_live_vars = 0;
        stats_.last_checkpoint_live_edges = 0;
    }

    // Exact fallback: if no previous checkpoint can be certified on the current
    // pattern, recompute the owned deterministic CCOLAMD ordering from scratch.
    // This preserves mathematical exactness and is the unavoidable worst case
    // for a tie-sensitive minimum-degree ordering under arbitrary structural edits.
    ++stats_.exact_suffix_recomputes;
    auto trace = batch_.order_with_trace(pattern);
    stats_.last_rollback_suffix = static_cast<int>(trace.steps.size());
    stats_.total_common_prefix += static_cast<std::uint64_t>(certified_events);
    stats_.total_rollback_suffix += static_cast<std::uint64_t>(stats_.last_rollback_suffix);
    return install_trace(std::move(trace), "exact_full_recompute", static_cast<int>(last_reused_prefix_.size()));
}

} // namespace islam
