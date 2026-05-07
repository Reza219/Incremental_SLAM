#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace islam {

Eigen::Matrix3d skew3(const Eigen::Vector3d& v);
Eigen::Matrix3d so3_exp(const Eigen::Vector3d& w);
Eigen::Vector3d so3_log(const Eigen::Matrix3d& R);
Eigen::Matrix3d so3_left_jacobian(const Eigen::Vector3d& w);
Eigen::Matrix3d so3_right_jacobian(const Eigen::Vector3d& w);
Eigen::Matrix3d so3_left_jacobian_inverse(const Eigen::Vector3d& w);
Eigen::Matrix3d so3_right_jacobian_inverse(const Eigen::Vector3d& w);
Eigen::Isometry3d se3_vec_to_iso(const Eigen::Matrix<double, 6, 1>& v);
Eigen::Matrix<double, 6, 1> se3_iso_to_vec(const Eigen::Isometry3d& T);
Eigen::Matrix<double, 6, 1> se3_xyz_quat_to_vec(double x,
                                                 double y,
                                                 double z,
                                                 double qx,
                                                 double qy,
                                                 double qz,
                                                 double qw);
Eigen::Matrix<double, 6, 1> se3_relative_error(const Eigen::Matrix<double, 6, 1>& x1,
                                                const Eigen::Matrix<double, 6, 1>& x2,
                                                const Eigen::Matrix<double, 6, 1>& z);

void se3_relative_error_jacobians(const Eigen::Matrix<double, 6, 1>& x1,
                                  const Eigen::Matrix<double, 6, 1>& x2,
                                  const Eigen::Matrix<double, 6, 1>& z,
                                  Eigen::Matrix<double, 6, 1>& e,
                                  Eigen::Matrix<double, 6, 6>& A,
                                  Eigen::Matrix<double, 6, 6>& B);

} // namespace islam
