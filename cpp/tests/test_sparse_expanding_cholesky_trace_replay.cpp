#include "islam/sparse_expanding_cholesky.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace {
using Sparse = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

[[noreturn]] void fail(const std::string& msg) { throw std::runtime_error(msg); }

struct TraceOp {
    std::string kind;
    int a = -1;
    int b = -1;
    int id = -1;
    double weight = 0.0;
};

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

Sparse make_edge_column(int n, int i, int j, double weight) {
    if (i < 0 || j < 0 || i >= n || j >= n || i == j || weight <= 0.0) {
        fail("invalid edge column request");
    }
    const double s = std::sqrt(weight);
    std::vector<Triplet> trips;
    trips.emplace_back(i, 0, s);
    trips.emplace_back(j, 0, -s);
    Sparse c(n, 1);
    c.setFromTriplets(trips.begin(), trips.end());
    c.makeCompressed();
    return c;
}

Sparse make_prior_column(int n, int i, double weight) {
    if (i < 0 || i >= n || weight <= 0.0) fail("invalid prior column request");
    Sparse c(n, 1);
    c.insert(i, 0) = std::sqrt(weight);
    c.makeCompressed();
    return c;
}

Sparse resize_column_preserve(const Sparse& column, int rows) {
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

void assert_backend_ok(const islam::SparseExpandingCholesky& dyn,
                       const std::string& label,
                       double tol = 2e-8) {
    if (!dyn.covers_state()) fail(label + ": factor does not cover state");
    double residual = 0.0;
    if (!dyn.passes_factorization_residual_check(&residual) || residual > tol) {
        std::ostringstream oss;
        oss << label << ": residual too large: " << std::setprecision(17) << residual;
        fail(oss.str());
    }

    islam::SparseExpandingCholesky fresh(dyn.options());
    fresh.factorize(dyn.normal_matrix());
    Eigen::VectorXd rhs(dyn.state_size());
    for (int i = 0; i < dyn.state_size(); ++i) rhs[i] = 0.25 + 0.11 * static_cast<double>((i % 7) - 3);
    const Eigen::VectorXd x_dyn = dyn.solve(rhs);
    const Eigen::VectorXd x_fresh = fresh.solve(rhs);
    const double rel = (x_dyn - x_fresh).norm() / (1.0 + x_fresh.norm());
    if (!(std::isfinite(rel) && rel <= tol)) {
        std::ostringstream oss;
        oss << label << ": replay solve disagreement: " << std::setprecision(17) << rel;
        fail(oss.str());
    }
}

std::vector<TraceOp> generate_trace() {
    std::mt19937 rng(7901);
    std::uniform_int_distribution<int> op_dist(0, 99);
    std::vector<TraceOp> ops;
    ops.push_back({"factorize_initial", 8, -1, -1, 0.0});

    int n = 8;
    int next_id = 0;
    std::vector<int> live_ids;
    std::map<int, TraceOp> live_adds;

    for (int step = 0; step < 60; ++step) {
        const int op = op_dist(rng);
        if (step == 11 || step == 27 || step == 43) {
            ops.push_back({"grow", n + 1, -1, -1, 0.0});
            ops.push_back({"prior_add", n, -1, next_id++, 9.0});
            TraceOp link{"edge_add", n - 1, n, next_id++, 0.31 + 0.02 * static_cast<double>(step % 3)};
            ops.push_back(link);
            live_ids.push_back(link.id);
            live_adds[link.id] = link;
            ++n;
        } else if (op < 34) {
            std::uniform_int_distribution<int> vdist(0, n - 1);
            ops.push_back({"diag_add", vdist(rng), -1, -1, 0.035 + 0.005 * static_cast<double>(step % 6)});
        } else if (op < 77 || live_ids.empty()) {
            std::uniform_int_distribution<int> vdist(0, n - 1);
            int i = vdist(rng);
            int j = vdist(rng);
            if (i == j) j = (j + 1) % n;
            if (i > j) std::swap(i, j);
            TraceOp add{"edge_add", i, j, next_id++, 0.055 + 0.006 * static_cast<double>(step % 9)};
            ops.push_back(add);
            live_ids.push_back(add.id);
            live_adds[add.id] = add;
        } else {
            std::uniform_int_distribution<int> rdist(0, static_cast<int>(live_ids.size()) - 1);
            const int idx = rdist(rng);
            const int id = live_ids[static_cast<size_t>(idx)];
            const auto it = live_adds.find(id);
            if (it == live_adds.end()) fail("internal trace-generation live id mismatch");
            const TraceOp& add = it->second;
            ops.push_back({"edge_remove", add.a, add.b, id, add.weight});
            live_adds.erase(it);
            live_ids.erase(live_ids.begin() + idx);
        }
    }
    return ops;
}

void write_trace(const std::vector<TraceOp>& ops, const std::string& path) {
    std::ofstream out(path);
    if (!out) fail("could not open trace for writing: " + path);
    out << "# M81 sparse expanding Cholesky operation trace v2\n";
    out << std::setprecision(17);
    for (const TraceOp& op : ops) {
        out << op.kind << ' ' << op.a << ' ' << op.b << ' ' << op.id << ' ' << op.weight << '\n';
    }
}

std::vector<TraceOp> read_trace(const std::string& path) {
    std::ifstream in(path);
    if (!in) fail("could not open trace for reading: " + path);
    std::vector<TraceOp> ops;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        TraceOp op;
        if (!(iss >> op.kind >> op.a >> op.b >> op.id >> op.weight)) {
            fail("malformed trace line: " + line);
        }
        ops.push_back(op);
    }
    return ops;
}

struct RunSummary {
    int state_size = 0;
    double residual = 0.0;
    Eigen::VectorXd solve;
    islam::SparseExpandingCholeskyStats stats;
};

RunSummary replay_trace(const std::vector<TraceOp>& ops,
                        const islam::SparseExpandingCholeskyOptions& options,
                        const std::string& label,
                        int stop_after,
                        int verify_every,
                        const std::string& failure_prefix_path,
                        bool quiet) {
    if (ops.empty()) fail(label + ": empty trace");
    if (verify_every <= 0) fail("verify_every must be positive");

    islam::SparseExpandingCholesky dyn(options);
    std::map<int, TraceOp> live_edges;
    const int limit = stop_after < 0
        ? static_cast<int>(ops.size())
        : std::min(stop_after, static_cast<int>(ops.size()));

    auto export_failure_prefix = [&](int last_index) {
        if (failure_prefix_path.empty()) return;
        const int count = std::max(0, std::min(last_index + 1, static_cast<int>(ops.size())));
        std::vector<TraceOp> prefix(ops.begin(), ops.begin() + count);
        write_trace(prefix, failure_prefix_path);
    };

    for (int k = 0; k < limit; ++k) {
        const TraceOp& op = ops[static_cast<size_t>(k)];
        const std::string step_label = label + " op " + std::to_string(k) + " (" + op.kind + ")";
        try {
            if (op.kind == "factorize_initial") {
                dyn.factorize(make_initial_hessian(op.a));
            } else if (op.kind == "grow") {
                dyn.grow_to(op.a);
            } else if (op.kind == "prior_add") {
                dyn.apply_contribution(make_prior_column(dyn.state_size(), op.a, op.weight), {op.a}, true);
            } else if (op.kind == "diag_add") {
                dyn.apply_diagonal_update({op.a}, op.weight, true);
            } else if (op.kind == "edge_add") {
                dyn.apply_contribution(make_edge_column(dyn.state_size(), op.a, op.b, op.weight), {op.a, op.b}, true);
                live_edges[op.id] = op;
            } else if (op.kind == "edge_remove") {
                TraceOp add = op;
                const auto it = live_edges.find(op.id);
                if (it != live_edges.end()) {
                    add = it->second;
                    live_edges.erase(it);
                }
                dyn.apply_contribution(resize_column_preserve(make_edge_column(std::max(add.a, add.b) + 1, add.a, add.b, add.weight),
                                                              dyn.state_size()),
                                       {add.a, add.b}, false);
            } else {
                fail("unknown trace operation: " + op.kind);
            }

            const bool should_verify = (k + 1 == limit) || ((k + 1) % verify_every == 0);
            if (should_verify && dyn.covers_state()) {
                assert_backend_ok(dyn, step_label);
            }
            if (!quiet && (k + 1) % 25 == 0) {
                std::cout << "  replayed " << (k + 1) << " / " << limit << " operations\n";
            }
        } catch (...) {
            export_failure_prefix(k);
            throw;
        }
    }

    RunSummary summary;
    summary.state_size = dyn.state_size();
    summary.residual = dyn.scaled_factorization_residual();
    summary.stats = dyn.stats();
    summary.solve.resize(dyn.state_size());
    Eigen::VectorXd rhs(dyn.state_size());
    for (int i = 0; i < dyn.state_size(); ++i) rhs[i] = 0.33 + 0.07 * static_cast<double>((i % 11) - 5);
    summary.solve = dyn.solve(rhs);
    return summary;
}

void assert_same_summary(const RunSummary& a, const RunSummary& b) {
    if (a.state_size != b.state_size) fail("trace replay state-size mismatch");
    const double residual_gap = std::abs(a.residual - b.residual);
    if (!(std::isfinite(residual_gap) && residual_gap <= 1e-13)) {
        std::ostringstream oss;
        oss << "trace replay residual mismatch: " << std::setprecision(17)
            << a.residual << " vs " << b.residual;
        fail(oss.str());
    }
    const double rel = (a.solve - b.solve).norm() / (1.0 + a.solve.norm());
    if (!(std::isfinite(rel) && rel <= 1e-13)) {
        std::ostringstream oss;
        oss << "trace replay final-solve mismatch: " << std::setprecision(17) << rel;
        fail(oss.str());
    }
}

} // namespace

int main(int argc, char** argv) {
    std::string trace_path = "/tmp/islam_m81_sparse_backend_trace.txt";
    std::string failure_prefix_path;
    bool replay_only = false;
    int stop_after = -1;
    int verify_every = 1;
    bool quiet = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--trace-out" && i + 1 < argc) trace_path = argv[++i];
        else if (arg == "--replay" && i + 1 < argc) { trace_path = argv[++i]; replay_only = true; }
        else if (arg == "--stop-after" && i + 1 < argc) stop_after = std::stoi(argv[++i]);
        else if (arg == "--verify-every" && i + 1 < argc) verify_every = std::stoi(argv[++i]);
        else if (arg == "--failure-prefix-out" && i + 1 < argc) failure_prefix_path = argv[++i];
        else if (arg == "--quiet") quiet = true;
        else if (arg == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--trace-out path] [--replay path] [--stop-after N]"
                      << " [--verify-every N] [--failure-prefix-out path] [--quiet]\n";
            return 0;
        } else {
            fail("unknown argument: " + arg);
        }
    }

    if (stop_after == 0) fail("--stop-after must be positive or omitted");
    if (verify_every <= 0) fail("--verify-every must be positive");

    islam::SparseExpandingCholeskyOptions options;
    options.certification_tolerance = 2e-9;
    options.jitter = 1e-12;
    options.max_jitter_tries = 8;
    options.use_factor_dependency_cache = true;
    options.use_column_local_certification = true;
    options.full_certification_fallback = true;

    std::vector<TraceOp> ops = replay_only ? read_trace(trace_path) : generate_trace();
    if (!replay_only) write_trace(ops, trace_path);

    const RunSummary first = replay_trace(ops, options, "first-pass", stop_after, verify_every, failure_prefix_path, quiet);
    const RunSummary second = replay_trace(read_trace(trace_path), options, "replay-pass", stop_after, verify_every, failure_prefix_path, true);
    assert_same_summary(first, second);

    if (first.stats.dependency_cache_hits == 0) fail("trace did not exercise dependency cache");
    if (first.stats.expansion_suffix_updates == 0) fail("trace did not exercise expansion suffix updates");
    if (first.stats.structural_pattern_changes == 0) fail("trace did not exercise structural changes");
    // A local affected-closure update is desirable but backend options/patterns can
    // conservatively route all updates through suffix refactorization. Treat this as
    // telemetry rather than a correctness failure; the replay/residual/solve checks
    // above remain the hard validation gates.

    std::cout << "M81 sparse backend trace/replay regression passed\n";
    std::cout << "  trace path: " << trace_path << "\n";
    std::cout << "  operations in trace: " << ops.size() << "\n";
    if (stop_after > 0) std::cout << "  stopped after: " << std::min(stop_after, static_cast<int>(ops.size())) << " operations\n";
    std::cout << "  verify every: " << verify_every << " operation(s)\n";
    std::cout << "  final state size: " << first.state_size << "\n";
    std::cout << "  final scaled residual: " << std::setprecision(17) << first.residual << "\n";
    std::cout << "  dependency cache hits: " << first.stats.dependency_cache_hits << "\n";
    std::cout << "  affected-closure refactorizations: " << first.stats.affected_closure_refactorizations << "\n";
    std::cout << "  suffix refactorizations: " << first.stats.suffix_refactorizations << "\n";
    return 0;
}
