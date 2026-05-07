#include "islam/benchmark.hpp"
#include "islam/io_g2o.hpp"
#include "islam/io_toro.hpp"
#include "islam/reorder_edges.hpp"
#include "islam/synthetic.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (char ch : line) {
        if (ch == '"') {
            in_quotes = !in_quotes;
        } else if (ch == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(cur);
    return out;
}

long long get_ll(const std::map<std::string, int>& col,
                 const std::vector<std::string>& row,
                 const std::string& name) {
    const auto it = col.find(name);
    if (it == col.end() || it->second < 0 || it->second >= static_cast<int>(row.size())) return 0;
    if (row[static_cast<size_t>(it->second)].empty()) return 0;
    return std::stoll(row[static_cast<size_t>(it->second)]);
}

double get_double(const std::map<std::string, int>& col,
                  const std::vector<std::string>& row,
                  const std::string& name) {
    const auto it = col.find(name);
    if (it == col.end() || it->second < 0 || it->second >= static_cast<int>(row.size())) return 0.0;
    if (row[static_cast<size_t>(it->second)].empty()) return 0.0;
    return std::stod(row[static_cast<size_t>(it->second)]);
}

int run_symbolic_audit(const fs::path& csv_path,
                       double max_order_oracle_fallback_rate,
                       double max_etree_fallback_rate) {
    std::ifstream ifs(csv_path);
    if (!ifs) throw std::runtime_error("Failed to open benchmark CSV for symbolic audit: " + csv_path.string());

    std::string header_line;
    if (!std::getline(ifs, header_line)) throw std::runtime_error("Empty benchmark CSV: " + csv_path.string());
    const auto header = split_csv_line(header_line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[static_cast<size_t>(i)]] = i;

    const std::vector<std::string> required = {
        "order_exact_output_certifications",
        "order_certified_local_order_accepts",
        "order_certified_exact_cache_order_accepts",
        "order_certified_oracle_order_fallbacks",
        "etree_local_update_attempts",
        "etree_local_update_accepts",
        "etree_local_update_fallbacks"
    };
    for (const auto& name : required) {
        if (col.find(name) == col.end()) {
            throw std::runtime_error("CSV is missing symbolic-audit column: " + name);
        }
    }

    long long rows = 0;
    long long order_exact_output_certifications = 0;
    long long order_certified_local_order_accepts = 0;
    long long order_certified_exact_cache_order_accepts = 0;
    long long order_certified_oracle_order_fallbacks = 0;
    long long etree_local_update_attempts = 0;
    long long etree_local_update_accepts = 0;
    long long etree_local_update_fallbacks = 0;
    double max_row_order_oracle_fallback_rate = 0.0;
    double max_row_etree_fallback_rate = 0.0;

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        const auto row = split_csv_line(line);
        ++rows;
        order_exact_output_certifications += get_ll(col, row, "order_exact_output_certifications");
        order_certified_local_order_accepts += get_ll(col, row, "order_certified_local_order_accepts");
        order_certified_exact_cache_order_accepts += get_ll(col, row, "order_certified_exact_cache_order_accepts");
        order_certified_oracle_order_fallbacks += get_ll(col, row, "order_certified_oracle_order_fallbacks");
        etree_local_update_attempts += get_ll(col, row, "etree_local_update_attempts");
        etree_local_update_accepts += get_ll(col, row, "etree_local_update_accepts");
        etree_local_update_fallbacks += get_ll(col, row, "etree_local_update_fallbacks");
        max_row_order_oracle_fallback_rate = std::max(max_row_order_oracle_fallback_rate,
                                                       get_double(col, row, "order_oracle_fallback_rate"));
        const long long etree_attempts_row = get_ll(col, row, "etree_local_update_attempts");
        const long long etree_fallbacks_row = get_ll(col, row, "etree_local_update_fallbacks");
        if (etree_attempts_row > 0) {
            max_row_etree_fallback_rate = std::max(max_row_etree_fallback_rate,
                                                   static_cast<double>(etree_fallbacks_row) / static_cast<double>(etree_attempts_row));
        }
    }

    const double order_oracle_fallback_rate =
        (order_exact_output_certifications > 0)
            ? static_cast<double>(order_certified_oracle_order_fallbacks) / static_cast<double>(order_exact_output_certifications)
            : 0.0;
    const double order_certified_non_oracle_rate =
        (order_exact_output_certifications > 0)
            ? static_cast<double>(order_certified_local_order_accepts + order_certified_exact_cache_order_accepts) / static_cast<double>(order_exact_output_certifications)
            : 0.0;
    const double etree_fallback_rate =
        (etree_local_update_attempts > 0)
            ? static_cast<double>(etree_local_update_fallbacks) / static_cast<double>(etree_local_update_attempts)
            : 0.0;
    const double etree_accept_rate =
        (etree_local_update_attempts > 0)
            ? static_cast<double>(etree_local_update_accepts) / static_cast<double>(etree_local_update_attempts)
            : 0.0;

    const bool pass = order_oracle_fallback_rate <= max_order_oracle_fallback_rate
                   && etree_fallback_rate <= max_etree_fallback_rate;

    std::cout << "{\n"
              << "  \"csv\": \"" << csv_path.string() << "\",\n"
              << "  \"rows\": " << rows << ",\n"
              << "  \"order_exact_output_certifications\": " << order_exact_output_certifications << ",\n"
              << "  \"order_certified_local_order_accepts\": " << order_certified_local_order_accepts << ",\n"
              << "  \"order_certified_exact_cache_order_accepts\": " << order_certified_exact_cache_order_accepts << ",\n"
              << "  \"order_certified_oracle_order_fallbacks\": " << order_certified_oracle_order_fallbacks << ",\n"
              << "  \"order_certified_non_oracle_rate\": " << order_certified_non_oracle_rate << ",\n"
              << "  \"order_oracle_fallback_rate\": " << order_oracle_fallback_rate << ",\n"
              << "  \"max_row_order_oracle_fallback_rate\": " << max_row_order_oracle_fallback_rate << ",\n"
              << "  \"etree_local_update_attempts\": " << etree_local_update_attempts << ",\n"
              << "  \"etree_local_update_accepts\": " << etree_local_update_accepts << ",\n"
              << "  \"etree_local_update_fallbacks\": " << etree_local_update_fallbacks << ",\n"
              << "  \"etree_accept_rate\": " << etree_accept_rate << ",\n"
              << "  \"etree_fallback_rate\": " << etree_fallback_rate << ",\n"
              << "  \"max_row_etree_fallback_rate\": " << max_row_etree_fallback_rate << ",\n"
              << "  \"max_allowed_order_oracle_fallback_rate\": " << max_order_oracle_fallback_rate << ",\n"
              << "  \"max_allowed_etree_fallback_rate\": " << max_etree_fallback_rate << ",\n"
              << "  \"pass\": " << (pass ? "true" : "false") << "\n"
              << "}\n";
    return pass ? 0 : 2;
}

int run_symbolic_replay_audit(const fs::path& csv_path) {
    std::ifstream ifs(csv_path);
    if (!ifs) throw std::runtime_error("Failed to open benchmark CSV for replay audit: " + csv_path.string());

    std::string header_line;
    if (!std::getline(ifs, header_line)) throw std::runtime_error("Empty benchmark CSV: " + csv_path.string());
    const auto header = split_csv_line(header_line);
    std::map<std::string, int> col;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[static_cast<size_t>(i)]] = i;

    const std::vector<std::string> required = {
        "order_certified_oracle_order_fallbacks", "order_no_oracle_cache_hits",
        "order_no_oracle_cache_misses", "etree_exact_recomputes",
        "etree_exact_cache_hits", "etree_exact_cache_misses",
        "etree_no_exact_cache_hits", "etree_no_exact_cache_misses",
        "etree_local_update_fallbacks"
    };
    for (const auto& name : required) {
        if (col.find(name) == col.end()) throw std::runtime_error("CSV is missing replay-audit column: " + name);
    }

    long long rows = 0;
    long long order_oracle_fallbacks = 0, order_no_oracle_cache_hits = 0, order_no_oracle_cache_misses = 0;
    long long etree_exact_recomputes = 0, etree_exact_cache_hits = 0, etree_exact_cache_misses = 0;
    long long etree_no_exact_cache_hits = 0, etree_no_exact_cache_misses = 0, etree_local_update_fallbacks = 0;

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        const auto row = split_csv_line(line);
        ++rows;
        order_oracle_fallbacks += get_ll(col, row, "order_certified_oracle_order_fallbacks");
        order_no_oracle_cache_hits += get_ll(col, row, "order_no_oracle_cache_hits");
        order_no_oracle_cache_misses += get_ll(col, row, "order_no_oracle_cache_misses");
        etree_exact_recomputes += get_ll(col, row, "etree_exact_recomputes");
        etree_exact_cache_hits += get_ll(col, row, "etree_exact_cache_hits");
        etree_exact_cache_misses += get_ll(col, row, "etree_exact_cache_misses");
        etree_no_exact_cache_hits += get_ll(col, row, "etree_no_exact_cache_hits");
        etree_no_exact_cache_misses += get_ll(col, row, "etree_no_exact_cache_misses");
        etree_local_update_fallbacks += get_ll(col, row, "etree_local_update_fallbacks");
    }

    const bool pass = order_oracle_fallbacks == 0 && order_no_oracle_cache_misses == 0 &&
                      etree_no_exact_cache_misses == 0 && etree_local_update_fallbacks == 0;

    std::cout << "{\n"
              << "  \"csv\": \"" << csv_path.string() << "\",\n"
              << "  \"rows\": " << rows << ",\n"
              << "  \"order_certified_oracle_order_fallbacks\": " << order_oracle_fallbacks << ",\n"
              << "  \"order_no_oracle_cache_hits\": " << order_no_oracle_cache_hits << ",\n"
              << "  \"order_no_oracle_cache_misses\": " << order_no_oracle_cache_misses << ",\n"
              << "  \"etree_exact_recomputes\": " << etree_exact_recomputes << ",\n"
              << "  \"etree_exact_cache_hits\": " << etree_exact_cache_hits << ",\n"
              << "  \"etree_exact_cache_misses\": " << etree_exact_cache_misses << ",\n"
              << "  \"etree_no_exact_cache_hits\": " << etree_no_exact_cache_hits << ",\n"
              << "  \"etree_no_exact_cache_misses\": " << etree_no_exact_cache_misses << ",\n"
              << "  \"etree_local_update_fallbacks\": " << etree_local_update_fallbacks << ",\n"
              << "  \"symbolic_replay_audit_passed\": " << (pass ? "true" : "false") << "\n"
              << "}\n";
    return pass ? 0 : 2;
}


} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cout
                << "Usage:\n"
                << "  islam_run_benchmark synthetic <num_poses> <csv_out> [loop_period] [max_batches] [max_gn_iter] [selective_alpha] [lc_gap] [eta_threshold] [gating] [use_spo] [backend] [dx_threshold]\n"
                << "  islam_run_benchmark <graph.{g2o|graph}> <csv_out> [max_batches] [max_gn_iter] [selective_alpha] [lc_gap] [eta_threshold] [gating] [use_spo] [backend] [dx_threshold]\n"
                << "  islam_run_benchmark audit-symbolic <benchmark.csv> [max_order_oracle_fallback_rate] [max_etree_fallback_rate]\n"
                << "  islam_run_benchmark audit-symbolic-replay <benchmark.csv>\n";
            return 0;
        }

        const std::string mode_or_path = argv[1];
        if (mode_or_path == "audit-symbolic") {
            if (argc < 3) throw std::runtime_error("audit-symbolic mode requires <benchmark.csv>");
            const double max_order_oracle_fallback_rate = (argc >= 4) ? std::stod(argv[3]) : 0.0;
            const double max_etree_fallback_rate = (argc >= 5) ? std::stod(argv[4]) : 0.0;
            return run_symbolic_audit(fs::path(argv[2]), max_order_oracle_fallback_rate, max_etree_fallback_rate);
        }
        if (mode_or_path == "audit-symbolic-replay") {
            if (argc < 3) throw std::runtime_error("audit-symbolic-replay mode requires <benchmark.csv>");
            return run_symbolic_replay_audit(fs::path(argv[2]));
        }

        islam::Graph full;
        fs::path csv_out;
        int max_batches = 200;
        int max_gn_iter = 10;
        double selective_alpha = 0.3;
        int lc_gap = 5;
        double eta_threshold = 1.0;
        double dx_threshold = 1e-3;
        islam::GatingMode gating_mode = islam::GatingMode::IGG;
        bool use_spo = true;
        islam::LinearSolverBackend backend = islam::LinearSolverBackend::Auto;

        if (mode_or_path == "synthetic") {
            if (argc < 4) throw std::runtime_error("synthetic mode requires <num_poses> <csv_out>");
            const int num_poses = std::max(2, std::stoi(argv[2]));
            csv_out = fs::path(argv[3]);
            const int loop_period = (argc >= 5) ? std::max(2, std::stoi(argv[4])) : 5;
            max_batches = (argc >= 6) ? std::max(1, std::stoi(argv[5])) : 200;
            max_gn_iter = (argc >= 7) ? std::max(1, std::stoi(argv[6])) : 10;
            selective_alpha = (argc >= 8) ? std::stod(argv[7]) : 0.3;
            lc_gap = (argc >= 9) ? std::max(0, std::stoi(argv[8])) : 5;
            eta_threshold = (argc >= 10) ? std::stod(argv[9]) : 1.0;
            gating_mode = (argc >= 11) ? islam::gating_mode_from_string(argv[10]) : islam::GatingMode::IGG;
            use_spo = (argc >= 12) ? (std::stoi(argv[11]) != 0) : true;
            backend = (argc >= 13) ? islam::linear_solver_backend_from_string(argv[12]) : islam::LinearSolverBackend::Auto;
            dx_threshold = (argc >= 14) ? std::stod(argv[13]) : 1e-3;
            full = islam::make_pose_chain_with_periodic_loops(num_poses, loop_period);
        } else {
            const fs::path path(mode_or_path);
            csv_out = (argc >= 3) ? fs::path(argv[2]) : fs::path("benchmark.csv");
            max_batches = (argc >= 4) ? std::max(1, std::stoi(argv[3])) : 200;
            max_gn_iter = (argc >= 5) ? std::max(1, std::stoi(argv[4])) : 10;
            selective_alpha = (argc >= 6) ? std::stod(argv[5]) : 0.3;
            lc_gap = (argc >= 7) ? std::max(0, std::stoi(argv[6])) : 5;
            eta_threshold = (argc >= 8) ? std::stod(argv[7]) : 1.0;
            gating_mode = (argc >= 9) ? islam::gating_mode_from_string(argv[8]) : islam::GatingMode::IGG;
            use_spo = (argc >= 10) ? (std::stoi(argv[9]) != 0) : true;
            backend = (argc >= 11) ? islam::linear_solver_backend_from_string(argv[10]) : islam::LinearSolverBackend::Auto;
            dx_threshold = (argc >= 12) ? std::stod(argv[11]) : 1e-3;
            const auto ext = path.extension().string();
            if (ext == ".g2o") full = islam::read_graph_g2o(path);
            else if (ext == ".graph") full = islam::read_graph_toro(path);
            else throw std::runtime_error("Unsupported file extension: " + ext);
            full = islam::reorder_edges(full);
        }

        const auto stats = islam::run_selective_benchmark(
            full, max_batches, max_gn_iter, dx_threshold, 1.0, selective_alpha, lc_gap, gating_mode, use_spo, eta_threshold, backend);
        islam::write_benchmark_csv(csv_out, stats);
        std::cout << "Wrote benchmark CSV: " << csv_out << " (" << stats.size() << " rows, requested backend=" << islam::to_string(backend) << ")\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
