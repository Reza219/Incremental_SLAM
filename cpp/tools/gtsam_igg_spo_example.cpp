#include "islam/gtsam_igg_spo.hpp"

#ifdef ISLAM_HAS_GTSAM

#include <gtsam/geometry/Pose2.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/inference/Symbol.h>

#include <iostream>

int main() {
    using gtsam::symbol_shorthand::X;

    auto prior_noise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(3) << 1e-3, 1e-3, 1e-3).finished());
    auto odo_noise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(3) << 0.1, 0.1, 0.05).finished());

    gtsam::NonlinearFactorGraph initial_graph;
    initial_graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose2>>(X(0), gtsam::Pose2(0.0, 0.0, 0.0), prior_noise);

    gtsam::Values initial_values;
    initial_values.insert(X(0), gtsam::Pose2(0.0, 0.0, 0.0));

    islam::GtsamIggSpoParams params;
    params.eta_threshold = 0.25;
    params.dx_threshold = 1e-6;
    params.max_gn_iterations = 5;
    params.selective_alpha = 0.8;

    islam::GtsamIggSpoSolver solver(params);
    solver.reset(initial_graph, initial_values);

    for (int i = 1; i <= 4; ++i) {
        gtsam::NonlinearFactorGraph new_factors;
        new_factors.emplace_shared<gtsam::BetweenFactor<gtsam::Pose2>>(
            X(i - 1), X(i), gtsam::Pose2(1.0, 0.0, 0.0), odo_noise);

        gtsam::Values new_values;
        new_values.insert(X(i), gtsam::Pose2(static_cast<double>(i), 0.0, 0.0));

        const auto stats = solver.update(new_factors, new_values);
        std::cout << "increment=" << i
                  << " global_gate=" << stats.global_gate
                  << " delta_eta=" << stats.delta_eta
                  << " iterations=" << stats.gn_iterations
                  << " error=" << stats.final_error << '\n';
    }

    return 0;
}

#else
int main() { return 0; }
#endif
