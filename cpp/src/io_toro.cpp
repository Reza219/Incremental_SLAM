#include "islam/io_toro.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace islam {
namespace {

std::string strip_comment(const std::string& line) {
    const auto pos = line.find('#');
    return pos == std::string::npos ? line : line.substr(0, pos);
}

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    for (std::string token; iss >> token;) {
        tokens.push_back(token);
    }
    return tokens;
}

Eigen::MatrixXd symmetrize(const Eigen::MatrixXd& M) {
    return 0.5 * (M + M.transpose());
}

} // namespace

Graph read_graph_toro(const std::filesystem::path& path) {
    std::ifstream fin(path);
    if (!fin) {
        throw std::runtime_error("Failed to open TORO file: " + path.string());
    }

    struct RawNode {
        int id = -1;
        Eigen::Vector3d value = Eigen::Vector3d::Zero();
    };

    struct RawEdge {
        int from_id = -1;
        int to_id = -1;
        Eigen::Vector3d z = Eigen::Vector3d::Zero();
        Eigen::Matrix3d Omega = Eigen::Matrix3d::Zero();
    };

    std::vector<RawNode> raw_nodes;
    std::vector<RawEdge> raw_edges;
    std::unordered_set<int> seen_ids;

    std::string line;
    while (std::getline(fin, line)) {
        line = strip_comment(line);
        const auto tokens = split_ws(line);
        if (tokens.empty()) {
            continue;
        }

        const std::string& tag = tokens.front();

        if (tag == "VERTEX2") {
            if (tokens.size() < 5) {
                throw std::runtime_error("Malformed VERTEX2 line");
            }

            const int id = static_cast<int>(std::stod(tokens[1]));
            if (!seen_ids.insert(id).second) {
                throw std::runtime_error("Duplicate vertex id in TORO");
            }

            Eigen::Vector3d v;
            v << std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4]);
            raw_nodes.push_back(RawNode{id, v});
        } else if (tag == "EDGE2") {
            if (tokens.size() < 12) {
                throw std::runtime_error("Malformed EDGE2 line");
            }

            RawEdge e;
            e.from_id = static_cast<int>(std::stod(tokens[1]));
            e.to_id = static_cast<int>(std::stod(tokens[2]));
            e.z << std::stod(tokens[3]), std::stod(tokens[4]), std::stod(tokens[5]);

            e.Omega.setZero();
            e.Omega(0, 0) = std::stod(tokens[6]);
            e.Omega(0, 1) = std::stod(tokens[7]);
            e.Omega(1, 0) = e.Omega(0, 1);
            e.Omega(1, 1) = std::stod(tokens[8]);
            e.Omega(2, 2) = std::stod(tokens[9]);
            e.Omega(0, 2) = std::stod(tokens[10]);
            e.Omega(2, 0) = e.Omega(0, 2);
            e.Omega(1, 2) = std::stod(tokens[11]);
            e.Omega(2, 1) = e.Omega(1, 2);
            e.Omega = symmetrize(e.Omega);

            raw_edges.push_back(e);
        }
    }

    Graph g;
    g.x.resize(static_cast<int>(raw_nodes.size()) * 3);

    int offset = 0;
    for (const auto& node : raw_nodes) {
        g.id_lookup[node.id] = NodeRef{offset, 3};
        g.x.segment(offset, 3) = node.value;
        g.var2node.push_back(node.id);
        g.var2node.push_back(node.id);
        g.var2node.push_back(node.id);
        offset += 3;
    }

    g.edges.reserve(raw_edges.size());
    for (const auto& e : raw_edges) {
        const auto it_from = g.id_lookup.find(e.from_id);
        const auto it_to = g.id_lookup.find(e.to_id);
        if (it_from == g.id_lookup.end() || it_to == g.id_lookup.end()) {
            throw std::runtime_error("Edge references unknown node id in TORO");
        }

        Edge edge;
        edge.type = EdgeType::PosePose;
        edge.from_id = e.from_id;
        edge.to_id = e.to_id;
        edge.from_idx = it_from->second.offset;
        edge.to_idx = it_to->second.offset;
        edge.measurement = e.z;
        edge.information = e.Omega;
        g.edges.push_back(std::move(edge));
    }

    return g;
}

} // namespace islam
