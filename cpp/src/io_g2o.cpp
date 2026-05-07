#include "islam/io_g2o.hpp"
#include "islam/se3.hpp"

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

Eigen::MatrixXd parse_upper_information(const std::vector<std::string>& tokens, size_t start, int dim) {
    const size_t need = start + static_cast<size_t>(dim * (dim + 1) / 2);
    if (tokens.size() < need) {
        throw std::runtime_error("Malformed upper-triangular information matrix in g2o");
    }
    Eigen::MatrixXd Omega = Eigen::MatrixXd::Zero(dim, dim);
    size_t k = start;
    for (int r = 0; r < dim; ++r) {
        for (int c = r; c < dim; ++c) {
            const double v = std::stod(tokens[k++]);
            Omega(r, c) = v;
            Omega(c, r) = v;
        }
    }
    return symmetrize(Omega);
}

} // namespace

Graph read_graph_g2o(const std::filesystem::path& path) {
    std::ifstream fin(path);
    if (!fin) {
        throw std::runtime_error("Failed to open g2o file: " + path.string());
    }

    struct RawNode {
        int id = -1;
        int dim = 0;
        Eigen::VectorXd value;
    };

    struct RawEdge {
        EdgeType type = EdgeType::PosePose;
        int from_id = -1;
        int to_id = -1;
        Eigen::VectorXd z;
        Eigen::MatrixXd Omega;
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

        if (tag == "VERTEX_SE2") {
            if (tokens.size() < 5) {
                throw std::runtime_error("Malformed VERTEX_SE2 line");
            }
            const int id = std::stoi(tokens[1]);
            if (!seen_ids.insert(id).second) {
                throw std::runtime_error("Duplicate vertex id in g2o");
            }

            Eigen::Vector3d v;
            v << std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4]);
            raw_nodes.push_back(RawNode{id, 3, v});
        } else if (tag == "VERTEX_XY") {
            if (tokens.size() < 4) {
                throw std::runtime_error("Malformed VERTEX_XY line");
            }
            const int id = std::stoi(tokens[1]);
            if (!seen_ids.insert(id).second) {
                throw std::runtime_error("Duplicate vertex id in g2o");
            }

            Eigen::Vector2d v;
            v << std::stod(tokens[2]), std::stod(tokens[3]);
            raw_nodes.push_back(RawNode{id, 2, v});
        } else if (tag == "VERTEX_SE3:QUAT") {
            if (tokens.size() < 9) {
                throw std::runtime_error("Malformed VERTEX_SE3:QUAT line");
            }
            const int id = std::stoi(tokens[1]);
            if (!seen_ids.insert(id).second) {
                throw std::runtime_error("Duplicate vertex id in g2o");
            }
            const auto v = se3_xyz_quat_to_vec(std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4]),
                                               std::stod(tokens[5]), std::stod(tokens[6]), std::stod(tokens[7]),
                                               std::stod(tokens[8]));
            raw_nodes.push_back(RawNode{id, 6, v});
        } else if (tag == "EDGE_SE2") {
            if (tokens.size() < 12) {
                throw std::runtime_error("Malformed EDGE_SE2 line");
            }

            Eigen::Vector3d z;
            z << std::stod(tokens[3]), std::stod(tokens[4]), std::stod(tokens[5]);

            Eigen::Matrix3d Omega = Eigen::Matrix3d::Zero();
            const double Ixx = std::stod(tokens[6]);
            const double Ixy = std::stod(tokens[7]);
            const double Ixth = std::stod(tokens[8]);
            const double Iyy = std::stod(tokens[9]);
            const double Iyth = std::stod(tokens[10]);
            const double Ithth = std::stod(tokens[11]);

            Omega(0, 0) = Ixx;
            Omega(0, 1) = Ixy;
            Omega(1, 0) = Ixy;
            Omega(0, 2) = Ixth;
            Omega(2, 0) = Ixth;
            Omega(1, 1) = Iyy;
            Omega(1, 2) = Iyth;
            Omega(2, 1) = Iyth;
            Omega(2, 2) = Ithth;

            raw_edges.push_back(RawEdge{
                EdgeType::PosePose,
                std::stoi(tokens[1]),
                std::stoi(tokens[2]),
                z,
                symmetrize(Omega)
            });
        } else if (tag == "EDGE_SE2_XY") {
            if (tokens.size() < 8) {
                throw std::runtime_error("Malformed EDGE_SE2_XY line");
            }

            Eigen::Vector2d z;
            z << std::stod(tokens[3]), std::stod(tokens[4]);

            Eigen::Matrix2d Omega = Eigen::Matrix2d::Zero();
            const double Ixx = std::stod(tokens[5]);
            const double Ixy = std::stod(tokens[6]);
            const double Iyy = std::stod(tokens[7]);
            Omega(0, 0) = Ixx;
            Omega(0, 1) = Ixy;
            Omega(1, 0) = Ixy;
            Omega(1, 1) = Iyy;

            raw_edges.push_back(RawEdge{
                EdgeType::PoseLandmark,
                std::stoi(tokens[1]),
                std::stoi(tokens[2]),
                z,
                symmetrize(Omega)
            });
        } else if (tag == "EDGE_SE3:QUAT") {
            if (tokens.size() < 30) {
                throw std::runtime_error("Malformed EDGE_SE3:QUAT line");
            }
            const auto z = se3_xyz_quat_to_vec(std::stod(tokens[3]), std::stod(tokens[4]), std::stod(tokens[5]),
                                               std::stod(tokens[6]), std::stod(tokens[7]), std::stod(tokens[8]),
                                               std::stod(tokens[9]));
            raw_edges.push_back(RawEdge{
                EdgeType::PosePose3D,
                std::stoi(tokens[1]),
                std::stoi(tokens[2]),
                z,
                parse_upper_information(tokens, 10, 6)
            });
        } else if (tag == "EDGE_GPS3" || tag == "EDGE_POSITION3") {
            if (tokens.size() < 11) {
                throw std::runtime_error("Malformed EDGE_GPS3/EDGE_POSITION3 line");
            }
            Eigen::Vector3d z;
            z << std::stod(tokens[2]), std::stod(tokens[3]), std::stod(tokens[4]);
            raw_edges.push_back(RawEdge{
                EdgeType::GpsPrior3D,
                std::stoi(tokens[1]),
                -1,
                z,
                parse_upper_information(tokens, 5, 3)
            });
        }
    }

    Graph g;
    int offset = 0;
    int total_dim = 0;
    for (const auto& node : raw_nodes) {
        total_dim += node.dim;
    }
    g.x.resize(total_dim);

    for (const auto& node : raw_nodes) {
        g.id_lookup[node.id] = NodeRef{offset, node.dim};
        g.x.segment(offset, node.dim) = node.value;
        for (int k = 0; k < node.dim; ++k) {
            g.var2node.push_back(node.id);
        }
        offset += node.dim;
    }

    g.edges.reserve(raw_edges.size());
    for (const auto& e : raw_edges) {
        const auto it_from = g.id_lookup.find(e.from_id);
        if (it_from == g.id_lookup.end()) {
            throw std::runtime_error("Edge references unknown from-node id in g2o");
        }
        auto it_to = g.id_lookup.end();
        if (e.to_id >= 0) {
            it_to = g.id_lookup.find(e.to_id);
            if (it_to == g.id_lookup.end()) {
                throw std::runtime_error("Edge references unknown to-node id in g2o");
            }
        }

        Edge edge;
        edge.type = e.type;
        edge.from_id = e.from_id;
        edge.to_id = e.to_id;
        edge.from_idx = it_from->second.offset;
        edge.to_idx = (e.to_id >= 0) ? it_to->second.offset : -1;
        edge.measurement = e.z;
        edge.information = e.Omega;
        g.edges.push_back(std::move(edge));
    }

    return g;
}

} // namespace islam
