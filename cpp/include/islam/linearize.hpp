#pragma once

#include "islam/graph.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace islam {

struct EdgeLinearization {
    Eigen::SparseMatrix<double> J;
    Eigen::VectorXd r;
};

void linearize_pose_pose(
    const Eigen::Vector3d& x1,
    const Eigen::Vector3d& x2,
    const Eigen::Vector3d& z,
    Eigen::Vector3d& e,
    Eigen::Matrix<double, 3, 3>& A,
    Eigen::Matrix<double, 3, 3>& B);

void linearize_pose_landmark(
    const Eigen::Vector3d& x,
    const Eigen::Vector2d& l,
    const Eigen::Vector2d& z,
    Eigen::Vector2d& e,
    Eigen::Matrix<double, 2, 3>& A,
    Eigen::Matrix<double, 2, 2>& B);

void linearize_pose_gnss(
    const Eigen::Vector3d& x,
    const Eigen::Vector2d& z,
    Eigen::Vector2d& e,
    Eigen::Matrix<double, 2, 3>& A);

void linearize_pose_pose_3d(
    const Eigen::Matrix<double, 6, 1>& x1,
    const Eigen::Matrix<double, 6, 1>& x2,
    const Eigen::Matrix<double, 6, 1>& z,
    Eigen::Matrix<double, 6, 1>& e,
    Eigen::Matrix<double, 6, 6>& A,
    Eigen::Matrix<double, 6, 6>& B);

void linearize_pose_gnss_3d(
    const Eigen::Matrix<double, 6, 1>& x,
    const Eigen::Vector3d& z,
    Eigen::Vector3d& e,
    Eigen::Matrix<double, 3, 6>& A);

EdgeLinearization jacobian_edge_jr(const Edge& edge, const Graph& graph);

} // namespace islam
