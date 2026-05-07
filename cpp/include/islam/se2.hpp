#pragma once

#include <Eigen/Dense>

namespace islam {

double normalize_angle(double angle);

Eigen::Matrix3d v2t(const Eigen::Vector3d& v);
Eigen::Vector3d t2v(const Eigen::Matrix3d& T);
Eigen::Matrix3d invt(const Eigen::Matrix3d& T);

} // namespace islam
