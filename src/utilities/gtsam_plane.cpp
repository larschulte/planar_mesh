#include <iostream>
#include <Eigen/Dense>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/sam/BearingRangeFactor.h>

#include <matplot/matplot.h>
namespace plt = matplot;

class RangeFactor : public gtsam::NoiseModelFactor3<double, gtsam::Unit3, gtsam::Unit3>
{
private:

    // plane 
    Eigen::Vector3d po_plane_; // plane origin
    Eigen::Vector3d v_plane_; // plane direction

    // point
    Eigen::Vector3d po_point_; // point origin
    double t_point_; // point direction

public:
    RangeFactor(gtsam::Key KEY_pp_distance, gtsam::Key KEY_pp_normal, gtsam::Key KEY_v_point, Eigen::Vector3d po_plane, Eigen::Vector3d v_plane, Eigen::Vector3d po_point, double t_point, const gtsam::SharedNoiseModel &model)
        : 
        gtsam::NoiseModelFactor3<double, gtsam::Unit3, gtsam::Unit3>(model, KEY_pp_distance, KEY_pp_normal, KEY_v_point),
        po_plane_(po_plane),
        v_plane_(v_plane),
        po_point_(po_point),
        t_point_(t_point)
    {}

    gtsam::Vector evaluateError(const double& _t_plane, const gtsam::Unit3& _n_plane, const gtsam::Unit3& _v_point, boost::optional<gtsam::Matrix &> H_distance = boost::none, boost::optional<gtsam::Matrix &> H_normal = boost::none, boost::optional<gtsam::Matrix &> H_bearing = boost::none) const override
    {
        // Convert to Eigen
        Eigen::Vector3d n_plane = _n_plane.unitVector();
        Eigen::Vector3d v_point = _v_point.unitVector();

        // h(q)
        Eigen::Vector3d pp_plane = po_plane_ + v_plane_ * _t_plane;
        double h = n_plane.dot(pp_plane - po_point_) / n_plane.dot(v_point);

        // e
        double e = h - t_point_;

        // error
        gtsam::Vector error(1);
        error << e;

        // Compute Jacobians
        if (H_distance) 
        {
            // "e" = n_plane.dot(po_plane_ + v_plane_ * "_t_plane" - po_point_) / n_plane.dot(v_point) - t_point_
            double d_e__d_t_plane = n_plane.dot(v_plane_) / n_plane.dot(v_point);

            // fill in the Jacobian
            *H_distance = gtsam::Matrix(1, 1);
            (*H_distance)(0, 0) = d_e__d_t_plane;
        }

        if (H_normal) 
        {
            // "e" = "n_plane".dot(pp_plane - po_point_) / "n_plane".dot(v_point) - t_point_
            double v = n_plane.dot(v_point);
            double u = n_plane.dot(pp_plane - po_point_);
            Eigen::Vector3d d_v__d_n_plane = v_point;
            Eigen::Vector3d d_u__d_n_plane = pp_plane - po_point_;
            Eigen::Vector3d d_e__d_n_plane = (v * d_u__d_n_plane - u * d_v__d_n_plane) / (v * v);

            // fill in the Jacobian
            *H_normal = gtsam::Matrix(1, 2);
            (*H_normal).block<1, 2>(0, 0) = d_e__d_n_plane.transpose() * _n_plane.basis();
        }

        if (H_bearing)
        {
            // "e" = n_plane.dot(pp_plane - po_point_) / n_plane.dot("v_point") - t_point_
            Eigen::Vector3d d_e_range__d_v_point = -1 * n_plane.dot(pp_plane - po_point_) * n_plane / (n_plane.dot(v_point) * n_plane.dot(v_point));

            // fill in the Jacobian
            *H_bearing = gtsam::Matrix::Zero(1, 2);
            (*H_bearing).block<1, 2>(0, 0) = d_e_range__d_v_point.transpose() * _v_point.basis();
        }

        return error;
    }
};

int main()
{
    // DATASET
    // origin + position
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> dataset
    {
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1)},
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 1)},
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 1)},
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1)}
    };

    // MEASUREMENT
    std::vector<double> MEASUREMENT_RANGES;
    std::vector<gtsam::Unit3> MEASUREMENT_BEARINGS;
    for (auto &data : dataset)
    {
        MEASUREMENT_RANGES.push_back((data.second - data.first).norm());
        MEASUREMENT_BEARINGS.push_back(gtsam::Unit3((data.second - data.first).normalized()));
    }

    // VARIABLES
    gtsam::Symbol VARIABLE_t_plane('t', 1);
    gtsam::Symbol VARIABLE_n_plane('n', 1);
    std::vector<gtsam::Symbol> VARIABLES_BEARINGS;
    for (std::size_t i = 0; i < dataset.size(); i++) VARIABLES_BEARINGS.push_back(gtsam::Symbol('v', i + 1));

    // CONSTRAINT
    Eigen::Vector3d CONSTRAINT_po_plane = dataset[0].first;
    Eigen::Vector3d CONSTRAINT_v_plane = (dataset[0].second - dataset[0].first).normalized();
        
    // FACTORS
    gtsam::NonlinearFactorGraph graph;
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_bearing = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(2) << 0.01, 0.01).finished());
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_range = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << 0.01).finished());
    for (std::size_t i = 0; i < dataset.size(); i++)
    {
        auto FACTOR_range = boost::make_shared<RangeFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLES_BEARINGS[i], CONSTRAINT_po_plane, CONSTRAINT_v_plane, dataset[i].first, MEASUREMENT_RANGES[i], NOISE_range);
        auto FACTOR_bearing = gtsam::PriorFactor<gtsam::Unit3>(VARIABLES_BEARINGS[i], MEASUREMENT_BEARINGS[i], NOISE_bearing);
        graph.add(FACTOR_range);
        graph.add(FACTOR_bearing);
    }

    // INITIAL
    gtsam::Values initial;
    initial.insert(VARIABLE_t_plane, (dataset[0].second - dataset[0].first).norm());
    initial.insert(VARIABLE_n_plane, gtsam::Unit3((dataset[0].first - dataset[0].second).normalized()));
    for (std::size_t i = 0; i < dataset.size(); i++) initial.insert(VARIABLES_BEARINGS[i], MEASUREMENT_BEARINGS[i]);
    
    // OPTIMIZATION
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimize();

    // OUTPUT
    graph.print("Factor graph:\n");
    result.print("Final result:\n");
    gtsam::Marginals marginals(graph, result);
    std::cout.precision(2);
    std::cout << "distance covariance:\n" << marginals.marginalCovariance(VARIABLE_t_plane) << std::endl;
    std::cout << "normal covariance:\n" << marginals.marginalCovariance(VARIABLE_n_plane) << std::endl;

    // PLOT
    Eigen::Vector3d plane_position = CONSTRAINT_po_plane + CONSTRAINT_v_plane * result.at<double>(VARIABLE_t_plane);
    Eigen::Vector3d plane_normal = result.at<gtsam::Unit3>(VARIABLE_n_plane).unitVector();

    // compute intersection of point and plane
    std::vector<Eigen::Vector3d> point_positions;
    for (std::size_t i = 0; i < dataset.size(); i++)
    {
        Eigen::Vector3d v_point = result.at<gtsam::Unit3>(VARIABLES_BEARINGS[i]).unitVector();
        double t_point = plane_normal.dot(plane_position - dataset[i].first) / plane_normal.dot(v_point);
        Eigen::Vector3d point_position = dataset[i].first + v_point * t_point;
        point_positions.push_back(point_position);
    }

    // plot 
    plt::figure();
    plt::hold(true);
    for (std::size_t i = 0; i < dataset.size(); i++)
    {
        // original line
        std::vector<double> x = {dataset[i].first(0), dataset[i].second(0)};
        std::vector<double> y = {dataset[i].first(1), dataset[i].second(1)};
        std::vector<double> z = {dataset[i].first(2), dataset[i].second(2)};
        plt::plot3(x, y, z, "r-");

        // new line
        std::vector<double> x_point = {dataset[i].first(0), point_positions[i](0)};
        std::vector<double> y_point = {dataset[i].first(1), point_positions[i](1)};
        std::vector<double> z_point = {dataset[i].first(2), point_positions[i](2)};
        plt::plot3(x_point, y_point, z_point, "g-");
    }
    // plot fitted plane
    { 
        // plot surface specified by plane position, normal
        auto [X, Y] = plt::meshgrid(plt::iota(-1, 0.1, 1), plt::iota(-1, 0.1, 1));
        auto Z = plt::transform(X, Y, [&](double x, double y) { return plane_position(2) - (plane_normal(0) * (x - plane_position(0)) + plane_normal(1) * (y - plane_position(1))) / plane_normal(2); });
        plt::surf(X, Y, Z)->face_alpha(0.1);
    }
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::zlabel("Z");
    plt::grid(true);
    plt::axis("equal");
    plt::show();

    // END
    return 0;
}