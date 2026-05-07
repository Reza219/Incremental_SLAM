#include "islam/linearize.hpp"
#include "islam/se2.hpp"
#include "islam/se3.hpp"

#include <Eigen/Eigenvalues>
#include <stdexcept>
#include <vector>

namespace islam {
namespace {

Eigen::MatrixXd whitening_factor(const Eigen::MatrixXd& information) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(0.5 * (information + information.transpose()));
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Failed eigen decomposition for whitening factor");
    }
    Eigen::VectorXd evals = es.eigenvalues().cwiseMax(0.0).cwiseSqrt();
    return es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();
}

int from_dim(const Edge& edge) {
    switch (edge.type) {
    case EdgeType::PosePose:
    case EdgeType::PoseLandmark:
    case EdgeType::GpsPrior:
        return 3;
    case EdgeType::PosePose3D:
    case EdgeType::GpsPrior3D:
        return 6;
    default:
        throw std::runtime_error("Unsupported edge type");
    }
}

int to_dim(const Edge& edge) {
    switch (edge.type) {
    case EdgeType::PosePose: return 3;
    case EdgeType::PoseLandmark: return 2;
    case EdgeType::GpsPrior: return 0;
    case EdgeType::PosePose3D: return 6;
    case EdgeType::GpsPrior3D: return 0;
    default:
        throw std::runtime_error("Unsupported edge type");
    }
}

} // namespace

void linearize_pose_pose(const Eigen::Vector3d& x1,
                         const Eigen::Vector3d& x2,
                         const Eigen::Vector3d& z,
                         Eigen::Vector3d& e,
                         Eigen::Matrix<double, 3, 3>& A,
                         Eigen::Matrix<double, 3, 3>& B) {
    const double x_1 = x1.x();
    const double y_1 = x1.y();
    const double th1 = x1.z();
    const double x_2 = x2.x();
    const double y_2 = x2.y();
    const double th2 = x2.z();

    const double cz = std::cos(z.z());
    const double sz = std::sin(z.z());
    const double c1 = std::cos(th1);
    const double s1 = std::sin(th1);

    Eigen::Matrix2d R1;
    R1 << c1, -s1,
          s1,  c1;

    Eigen::Matrix2d Rz;
    Rz << cz, -sz,
          sz,  cz;

    const Eigen::Vector2d delta_t(x_2 - x_1, y_2 - y_1);
    const Eigen::Vector2d pred_xy = R1.transpose() * delta_t;
    const Eigen::Vector2d e_xy = Rz.transpose() * (pred_xy - z.head<2>());
    const double e_th = normalize_angle(th2 - th1 - z.z());

    e << e_xy.x(), e_xy.y(), e_th;

    Eigen::Matrix<double, 2, 3> Q1;
    Q1 << -c1, -s1, (x_1 - x_2) * s1 + (y_2 - y_1) * c1,
           s1, -c1, (x_1 - x_2) * c1 + (y_1 - y_2) * s1;

    Eigen::Matrix<double, 2, 3> Q2;
    Q2 <<  c1,  s1, 0.0,
          -s1,  c1, 0.0;

    A.setZero();
    B.setZero();
    A.topRows<2>() = Rz.transpose() * Q1;
    B.topRows<2>() = Rz.transpose() * Q2;
    A(2, 2) = -1.0;
    B(2, 2) = 1.0;
}

void linearize_pose_landmark(const Eigen::Vector3d& x,
                             const Eigen::Vector2d& l,
                             const Eigen::Vector2d& z,
                             Eigen::Vector2d& e,
                             Eigen::Matrix<double, 2, 3>& A,
                             Eigen::Matrix<double, 2, 2>& B) {
    const double xi = x.x();
    const double yi = x.y();
    const double th = x.z();
    const double lx = l.x();
    const double ly = l.y();

    const double c = std::cos(th);
    const double s = std::sin(th);

    Eigen::Matrix2d Rt;
    Rt <<  c, s,
          -s, c;

    e = Rt * (Eigen::Vector2d(lx, ly) - Eigen::Vector2d(xi, yi)) - z;

    A << -c, -s, (xi - lx) * s + (ly - yi) * c,
          s, -c, (xi - lx) * c + (yi - ly) * s;

    B <<  c, s,
         -s, c;
}

void linearize_pose_gnss(const Eigen::Vector3d& x,
                         const Eigen::Vector2d& z,
                         Eigen::Vector2d& e,
                         Eigen::Matrix<double, 2, 3>& A) {
    e << x.x() - z.x(), x.y() - z.y();
    A << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0;
}

void linearize_pose_pose_3d(const Eigen::Matrix<double, 6, 1>& x1,
                            const Eigen::Matrix<double, 6, 1>& x2,
                            const Eigen::Matrix<double, 6, 1>& z,
                            Eigen::Matrix<double, 6, 1>& e,
                            Eigen::Matrix<double, 6, 6>& A,
                            Eigen::Matrix<double, 6, 6>& B) {
    se3_relative_error_jacobians(x1, x2, z, e, A, B);
}

void linearize_pose_gnss_3d(const Eigen::Matrix<double, 6, 1>& x,
                            const Eigen::Vector3d& z,
                            Eigen::Vector3d& e,
                            Eigen::Matrix<double, 3, 6>& A) {
    e = x.head<3>() - z;
    A.setZero();
    A.block<3,3>(0,0).setIdentity();
}

EdgeLinearization jacobian_edge_jr(const Edge& edge, const Graph& graph) {
    const int n = graph.state_size();
    const int fd = from_dim(edge);
    const int td = to_dim(edge);

    Eigen::VectorXd residual;
    Eigen::MatrixXd Jlocal;

    if (edge.type == EdgeType::PosePose) {
        if (edge.from_idx < 0 || edge.to_idx < 0) {
            throw std::runtime_error("Pose-pose edge missing indices");
        }
        const Eigen::Vector3d xf = graph.block(edge.from_idx, 3);
        const Eigen::Vector3d xt = graph.block(edge.to_idx, 3);
        const Eigen::Vector3d z = edge.measurement;
        Eigen::Vector3d e;
        Eigen::Matrix3d A, B;
        linearize_pose_pose(xf, xt, z, e, A, B);
        residual = e;
        Jlocal.resize(3, fd + td);
        Jlocal.leftCols(fd) = A;
        Jlocal.rightCols(td) = B;
    } else if (edge.type == EdgeType::PosePose3D) {
        if (edge.from_idx < 0 || edge.to_idx < 0) {
            throw std::runtime_error("3D pose-pose edge missing indices");
        }
        const Eigen::Matrix<double, 6, 1> xf = graph.block(edge.from_idx, 6);
        const Eigen::Matrix<double, 6, 1> xt = graph.block(edge.to_idx, 6);
        const Eigen::Matrix<double, 6, 1> z = edge.measurement;
        Eigen::Matrix<double, 6, 1> e;
        Eigen::Matrix<double, 6, 6> A, B;
        linearize_pose_pose_3d(xf, xt, z, e, A, B);
        residual = e;
        Jlocal.resize(6, fd + td);
        Jlocal.leftCols(fd) = A;
        Jlocal.rightCols(td) = B;
    } else if (edge.type == EdgeType::PoseLandmark) {
        if (edge.from_idx < 0 || edge.to_idx < 0) {
            throw std::runtime_error("Pose-landmark edge missing indices");
        }
        const Eigen::Vector3d xf = graph.block(edge.from_idx, 3);
        const Eigen::Vector2d xl = graph.block(edge.to_idx, 2);
        const Eigen::Vector2d z = edge.measurement;
        Eigen::Vector2d e;
        Eigen::Matrix<double, 2, 3> A;
        Eigen::Matrix<double, 2, 2> B;
        linearize_pose_landmark(xf, xl, z, e, A, B);
        residual = e;
        Jlocal.resize(2, fd + td);
        Jlocal.leftCols(fd) = A;
        Jlocal.rightCols(td) = B;
    } else if (edge.type == EdgeType::GpsPrior) {
        if (edge.from_idx < 0) {
            throw std::runtime_error("GPS edge missing from index");
        }
        const Eigen::Vector3d xf = graph.block(edge.from_idx, 3);
        const Eigen::Vector2d z = edge.measurement;
        Eigen::Vector2d e;
        Eigen::Matrix<double, 2, 3> A;
        linearize_pose_gnss(xf, z, e, A);
        residual = e;
        Jlocal = A;
    } else if (edge.type == EdgeType::GpsPrior3D) {
        if (edge.from_idx < 0) {
            throw std::runtime_error("3D GPS edge missing from index");
        }
        const Eigen::Matrix<double, 6, 1> xf = graph.block(edge.from_idx, 6);
        const Eigen::Vector3d z = edge.measurement;
        Eigen::Vector3d e;
        Eigen::Matrix<double, 3, 6> A;
        linearize_pose_gnss_3d(xf, z, e, A);
        residual = e;
        Jlocal = A;
    } else {
        throw std::runtime_error("Unsupported edge type in jacobian_edge_jr");
    }

    const Eigen::MatrixXd L = whitening_factor(edge.information);
    const Eigen::MatrixXd Jw = L * Jlocal;
    const Eigen::VectorXd rw = L * residual;

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(Jw.rows() * Jw.cols()));

    for (int r = 0; r < Jw.rows(); ++r) {
        for (int c = 0; c < fd; ++c) {
            const double val = Jw(r, c);
            if (val != 0.0) triplets.emplace_back(r, edge.from_idx + c, val);
        }
        for (int c = 0; c < td; ++c) {
            const double val = Jw(r, fd + c);
            if (val != 0.0) triplets.emplace_back(r, edge.to_idx + c, val);
        }
    }

    EdgeLinearization out;
    out.J.resize(Jw.rows(), n);
    out.J.setFromTriplets(triplets.begin(), triplets.end());
    out.r = rw;
    return out;
}

} // namespace islam
