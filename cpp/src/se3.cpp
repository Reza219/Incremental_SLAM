#include "islam/se3.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace islam {
namespace {

constexpr double kSmallAngle = 1e-8;
constexpr double kTinyAngle = 1e-12;

} // namespace

Eigen::Matrix3d skew3(const Eigen::Vector3d& v) {
    Eigen::Matrix3d S;
    S << 0.0, -v.z(), v.y(),
         v.z(), 0.0, -v.x(),
        -v.y(), v.x(), 0.0;
    return S;
}

Eigen::Matrix3d so3_exp(const Eigen::Vector3d& w) {
    const double theta = w.norm();
    const Eigen::Matrix3d W = skew3(w);
    if (theta < kTinyAngle) {
        return Eigen::Matrix3d::Identity() + W + 0.5 * W * W;
    }
    const double a = std::sin(theta) / theta;
    const double b = (1.0 - std::cos(theta)) / (theta * theta);
    return Eigen::Matrix3d::Identity() + a * W + b * W * W;
}

Eigen::Vector3d so3_log(const Eigen::Matrix3d& R_in) {
    Eigen::Matrix3d R = R_in;
    // Project gently onto SO(3) to reduce roundoff drift from numerical probes.
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    R = svd.matrixU() * svd.matrixV().transpose();
    if (R.determinant() < 0.0) {
        Eigen::Matrix3d U = svd.matrixU();
        U.col(2) *= -1.0;
        R = U * svd.matrixV().transpose();
    }

    const double cos_theta = std::clamp((R.trace() - 1.0) * 0.5, -1.0, 1.0);
    const double theta = std::acos(cos_theta);
    Eigen::Vector3d vee;
    vee << R(2,1) - R(1,2), R(0,2) - R(2,0), R(1,0) - R(0,1);
    if (theta < kTinyAngle) {
        return 0.5 * vee;
    }
    if (3.141592653589793238462643383279502884 - theta < 1e-7) {
        Eigen::AngleAxisd aa(R);
        return aa.angle() * aa.axis();
    }
    return (theta / (2.0 * std::sin(theta))) * vee;
}

Eigen::Matrix3d so3_left_jacobian(const Eigen::Vector3d& w) {
    const double theta = w.norm();
    const Eigen::Matrix3d W = skew3(w);
    const Eigen::Matrix3d W2 = W * W;
    if (theta < kSmallAngle) {
        return Eigen::Matrix3d::Identity() + 0.5 * W + (1.0 / 6.0) * W2;
    }
    const double theta2 = theta * theta;
    const double theta3 = theta2 * theta;
    return Eigen::Matrix3d::Identity()
         + ((1.0 - std::cos(theta)) / theta2) * W
         + ((theta - std::sin(theta)) / theta3) * W2;
}

Eigen::Matrix3d so3_right_jacobian(const Eigen::Vector3d& w) {
    return so3_left_jacobian(-w);
}

Eigen::Matrix3d so3_left_jacobian_inverse(const Eigen::Vector3d& w) {
    const double theta = w.norm();
    const Eigen::Matrix3d W = skew3(w);
    const Eigen::Matrix3d W2 = W * W;
    if (theta < kSmallAngle) {
        return Eigen::Matrix3d::Identity() - 0.5 * W + (1.0 / 12.0) * W2;
    }
    const double sin_theta = std::sin(theta);
    if (std::abs(sin_theta) < kTinyAngle) {
        // The closed form is singular exactly at integer multiples of pi. Fall back to
        // solving the well-conditioned Jacobian expression used in the forward map.
        return so3_left_jacobian(w).fullPivLu().inverse();
    }
    const double theta2 = theta * theta;
    const double coeff = (1.0 / theta2) - ((1.0 + std::cos(theta)) / (2.0 * theta * sin_theta));
    return Eigen::Matrix3d::Identity() - 0.5 * W + coeff * W2;
}

Eigen::Matrix3d so3_right_jacobian_inverse(const Eigen::Vector3d& w) {
    return so3_left_jacobian_inverse(-w);
}

Eigen::Isometry3d se3_vec_to_iso(const Eigen::Matrix<double, 6, 1>& v) {
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.translation() = v.head<3>();
    T.linear() = so3_exp(v.tail<3>());
    return T;
}

Eigen::Matrix<double, 6, 1> se3_iso_to_vec(const Eigen::Isometry3d& T) {
    Eigen::Matrix<double, 6, 1> v;
    v.head<3>() = T.translation();
    v.tail<3>() = so3_log(T.linear());
    return v;
}

Eigen::Matrix<double, 6, 1> se3_xyz_quat_to_vec(double x,
                                                 double y,
                                                 double z,
                                                 double qx,
                                                 double qy,
                                                 double qz,
                                                 double qw) {
    Eigen::Quaterniond q(qw, qx, qy, qz);
    const double n = q.norm();
    if (!(n > 0.0)) {
        throw std::runtime_error("Invalid zero-norm quaternion in SE3 input");
    }
    q.normalize();
    Eigen::Matrix<double, 6, 1> v;
    v.head<3>() << x, y, z;
    v.tail<3>() = so3_log(q.toRotationMatrix());
    return v;
}

Eigen::Matrix<double, 6, 1> se3_relative_error(const Eigen::Matrix<double, 6, 1>& x1,
                                                const Eigen::Matrix<double, 6, 1>& x2,
                                                const Eigen::Matrix<double, 6, 1>& z) {
    const Eigen::Isometry3d X1 = se3_vec_to_iso(x1);
    const Eigen::Isometry3d X2 = se3_vec_to_iso(x2);
    const Eigen::Isometry3d Z = se3_vec_to_iso(z);
    const Eigen::Isometry3d E = Z.inverse() * X1.inverse() * X2;
    Eigen::Matrix<double, 6, 1> out;
    out.head<3>() = E.translation();
    out.tail<3>() = so3_log(E.linear());
    return out;
}

void se3_relative_error_jacobians(const Eigen::Matrix<double, 6, 1>& x1,
                                  const Eigen::Matrix<double, 6, 1>& x2,
                                  const Eigen::Matrix<double, 6, 1>& z,
                                  Eigen::Matrix<double, 6, 1>& e,
                                  Eigen::Matrix<double, 6, 6>& A,
                                  Eigen::Matrix<double, 6, 6>& B) {
    const Eigen::Vector3d t1 = x1.head<3>();
    const Eigen::Vector3d w1 = x1.tail<3>();
    const Eigen::Vector3d t2 = x2.head<3>();
    const Eigen::Vector3d w2 = x2.tail<3>();
    const Eigen::Vector3d tz = z.head<3>();
    const Eigen::Vector3d wz = z.tail<3>();

    const Eigen::Matrix3d R1 = so3_exp(w1);
    const Eigen::Matrix3d R2 = so3_exp(w2);
    const Eigen::Matrix3d Rz = so3_exp(wz);
    const Eigen::Matrix3d R1T = R1.transpose();
    const Eigen::Matrix3d RzT = Rz.transpose();

    const Eigen::Vector3d dt = t2 - t1;
    const Eigen::Vector3d p = R1T * dt;
    const Eigen::Matrix3d Re = RzT * R1T * R2;
    const Eigen::Vector3d phi = so3_log(Re);

    e.head<3>() = RzT * (p - tz);
    e.tail<3>() = phi;

    A.setZero();
    B.setZero();

    // The state stores translation plus rotation-vector coordinates [t, w].
    // The expressions below are analytic derivatives with respect to additive
    // perturbations in those six coordinates. For a rotation-vector perturbation
    // dw, Exp(w + dw) = Exp(w) Exp(J_r(w) dw) + O(||dw||^2).
    const Eigen::Matrix3d Jr1 = so3_right_jacobian(w1);
    const Eigen::Matrix3d Jr2 = so3_right_jacobian(w2);
    const Eigen::Matrix3d Jlinv_phi = so3_left_jacobian_inverse(phi);
    const Eigen::Matrix3d Jrinv_phi = so3_right_jacobian_inverse(phi);

    A.block<3,3>(0,0) = -RzT * R1T;
    B.block<3,3>(0,0) =  RzT * R1T;
    A.block<3,3>(0,3) =  RzT * skew3(p) * Jr1;

    // Re = Rz^T R1^T R2. Perturbing R1 on the right creates a left perturbation
    // Exp(-Rz^T J_r(w1) dw1) Re; perturbing R2 creates the right perturbation
    // Re Exp(J_r(w2) dw2).
    A.block<3,3>(3,3) = -Jlinv_phi * RzT * Jr1;
    B.block<3,3>(3,3) =  Jrinv_phi * Jr2;
}

} // namespace islam
