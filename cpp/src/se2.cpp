#include "islam/se2.hpp"

#include <cmath>

namespace islam {

double normalize_angle(double angle) {
    constexpr double kPi = 3.14159265358979323846;
    const double wrapped = std::fmod(angle + kPi, 2.0 * kPi);
    const double positive = wrapped < 0.0 ? wrapped + 2.0 * kPi : wrapped;
    return positive - kPi;
}

Eigen::Matrix3d v2t(const Eigen::Vector3d& v) {
    const double c = std::cos(v.z());
    const double s = std::sin(v.z());

    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    T(0, 0) = c;  T(0, 1) = -s; T(0, 2) = v.x();
    T(1, 0) = s;  T(1, 1) =  c; T(1, 2) = v.y();
    return T;
}

Eigen::Vector3d t2v(const Eigen::Matrix3d& T) {
    Eigen::Vector3d v;
    v.x() = T(0, 2);
    v.y() = T(1, 2);
    v.z() = std::atan2(T(1, 0), T(0, 0));
    return v;
}

Eigen::Matrix3d invt(const Eigen::Matrix3d& T) {
    Eigen::Matrix3d out = Eigen::Matrix3d::Identity();
    const Eigen::Matrix2d R = T.topLeftCorner<2, 2>();
    const Eigen::Vector2d t = T.topRightCorner<2, 1>();
    out.topLeftCorner<2, 2>() = R.transpose();
    out.topRightCorner<2, 1>() = -R.transpose() * t;
    return out;
}

} // namespace islam
