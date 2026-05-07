#pragma once

#include <Eigen/Dense>

#include <string>

namespace islam {

enum class EdgeType {
    PosePose,
    PoseLandmark,
    GpsPrior,
    PosePose3D,
    GpsPrior3D
};

inline std::string to_string(EdgeType type) {
    switch (type) {
    case EdgeType::PosePose: return "P";
    case EdgeType::PoseLandmark: return "L";
    case EdgeType::GpsPrior: return "G";
    case EdgeType::PosePose3D: return "P3";
    case EdgeType::GpsPrior3D: return "G3";
    default: return "?";
    }
}

struct NodeRef {
    int offset = 0;
    int dimension = 0;
};

struct Edge {
    EdgeType type = EdgeType::PosePose;

    int from_id = -1;
    int to_id = -1;

    int from_idx = -1;
    int to_idx = -1;

    Eigen::VectorXd measurement;
    Eigen::MatrixXd information;

    [[nodiscard]] bool is_unary() const noexcept { return to_id < 0 || to_idx < 0; }
};

} // namespace islam
