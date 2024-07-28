#include "utilities/gtsam_plane_2d.hpp"
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>

#include <matplot/matplot.h>

RangeFactor::RangeFactor(gtsam::Key KEY_pp_distance, gtsam::Key KEY_pp_normal, gtsam::Key KEY_v_point, Eigen::Vector2d po_plane, Eigen::Vector2d v_plane, Eigen::Vector2d po_point, double t_point, const gtsam::SharedNoiseModel &model)
    : 
    gtsam::NoiseModelFactor3<double, gtsam::Rot2, gtsam::Rot2>(model, KEY_pp_distance, KEY_pp_normal, KEY_v_point),
    po_plane_(po_plane),
    v_plane_(v_plane),
    po_point_(po_point),
    t_point_(t_point)
{}

gtsam::Vector RangeFactor::evaluateError(const double& _t_plane, const gtsam::Rot2& _n_plane, const gtsam::Rot2& _v_point, boost::optional<gtsam::Matrix &> H_distance, boost::optional<gtsam::Matrix &> H_normal, boost::optional<gtsam::Matrix &> H_bearing) const
{
    // Convert to Eigen
    Eigen::Vector2d n_plane(std::cos(_n_plane.theta()), std::sin(_n_plane.theta()));
    Eigen::Vector2d v_point(std::cos(_v_point.theta()), std::sin(_v_point.theta()));

    // h(q)
    Eigen::Vector2d pp_plane = po_plane_ + v_plane_ * _t_plane;
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
        Eigen::Vector2d d_v__d_n_plane = v_point;
        Eigen::Vector2d d_u__d_n_plane = pp_plane - po_point_;
        Eigen::Vector2d d_e__d_n_plane = (v * d_u__d_n_plane - u * d_v__d_n_plane) / (v * v);

        // fill in the Jacobian
        // Convert d_e__d_n_plane to derivative with respect to the angle of Rot2
        double angle_n = _n_plane.theta();
        Eigen::Matrix2d R;
        R << -std::sin(angle_n), -std::cos(angle_n), std::cos(angle_n), -std::sin(angle_n);
        Eigen::Vector2d d_e__d_theta_n = R * d_e__d_n_plane;

        // fill in the Jacobian
        *H_normal = gtsam::Matrix(1, 1);
        (*H_normal)(0, 0) = d_e__d_theta_n.sum();
    }

    if (H_bearing)
    {
        // "e" = n_plane.dot(pp_plane - po_point_) / n_plane.dot("v_point") - t_point_
        Eigen::Vector2d d_e_range__d_v_point = -1 * n_plane.dot(pp_plane - po_point_) * n_plane / (n_plane.dot(v_point) * n_plane.dot(v_point));

        // fill in the Jacobian
        // Convert d_e__d_v_point to derivative with respect to the angle of Rot2
        double angle_v = _v_point.theta();
        Eigen::Matrix2d R;
        R << -std::sin(angle_v), -std::cos(angle_v), std::cos(angle_v), -std::sin(angle_v);
        Eigen::Vector2d d_e__d_theta_v = R * d_e_range__d_v_point;

        // fill in the Jacobian
        *H_bearing = gtsam::Matrix(1, 1);
        (*H_bearing)(0, 0) = d_e__d_theta_v.sum();
    }

    return error;
}

void fit_plane_to_points(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& dataset, Eigen::Vector2d &plane_position, Eigen::Vector2d &plane_normal, double bearing_noise, double range_noise, bool plot_graph)
{
    // CHECK NUMBER OF POINTS
    if (dataset.size() < 2) throw std::runtime_error("At least 3 points are required to fit a plane.");

    // MEASUREMENT
    std::vector<double> MEASUREMENT_RANGES;
    std::vector<gtsam::Rot2> MEASUREMENT_BEARINGS;
    for (auto &data : dataset)
    {
        MEASUREMENT_RANGES.push_back((data.second - data.first).norm());
        MEASUREMENT_BEARINGS.push_back(gtsam::Rot2::fromCosSin((data.second - data.first).normalized().x(), (data.second - data.first).normalized().y()));
    }

    // VARIABLES
    gtsam::Symbol VARIABLE_t_plane('t', 1);
    gtsam::Symbol VARIABLE_n_plane('n', 1);
    std::vector<gtsam::Symbol> VARIABLES_BEARINGS;
    for (std::size_t i = 0; i < dataset.size(); i++) VARIABLES_BEARINGS.push_back(gtsam::Symbol('v', i + 1));

    // CONSTRAINT
    Eigen::Vector2d CONSTRAINT_po_plane = dataset[0].first;
    Eigen::Vector2d CONSTRAINT_v_plane = (plane_position != Eigen::Vector2d::Zero()) ? (plane_position - dataset[0].first).normalized() : (dataset[0].second - dataset[0].first).normalized();
        
    // FACTORS
    gtsam::NonlinearFactorGraph graph;
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_bearing = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << bearing_noise).finished());
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_range = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << range_noise).finished());
    for (std::size_t i = 0; i < dataset.size(); i++)
    {
        auto FACTOR_range = boost::make_shared<RangeFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLES_BEARINGS[i], CONSTRAINT_po_plane, CONSTRAINT_v_plane, dataset[i].first, MEASUREMENT_RANGES[i], NOISE_range);
        auto FACTOR_bearing = gtsam::PriorFactor<gtsam::Rot2>(VARIABLES_BEARINGS[i], MEASUREMENT_BEARINGS[i], NOISE_bearing);
        graph.add(FACTOR_range);
        graph.add(FACTOR_bearing);
    }

    // INITIAL
    gtsam::Values initial;
    initial.insert(VARIABLE_n_plane, (plane_normal != Eigen::Vector2d::Zero()) ? gtsam::Rot2::fromCosSin(plane_normal.x(), plane_normal.y()) : gtsam::Rot2::fromCosSin((dataset[0].first - dataset[0].second).normalized().x(), (dataset[0].first - dataset[0].second).normalized().y()));
    initial.insert(VARIABLE_t_plane, (plane_position != Eigen::Vector2d::Zero()) ? (plane_position - dataset[0].first).norm() : (dataset[0].second - dataset[0].first).norm());
    for (std::size_t i = 0; i < dataset.size(); i++) initial.insert(VARIABLES_BEARINGS[i], MEASUREMENT_BEARINGS[i]);
    
    // OPTIMIZATION
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimize();
    plane_position = CONSTRAINT_po_plane + CONSTRAINT_v_plane * result.at<double>(VARIABLE_t_plane);
    plane_normal = Eigen::Vector2d(std::cos(result.at<gtsam::Rot2>(VARIABLE_n_plane).theta()), std::sin(result.at<gtsam::Rot2>(VARIABLE_n_plane).theta()));

    // PLOT GRAPH
    if (plot_graph)
    {
        // compute locations to plot
        std::vector<Eigen::Vector2d> points_origin;
        std::vector<Eigen::Vector2d> points_position;
        std::vector<Eigen::Vector2d> optimized_points_position;
        for (std::size_t i = 0; i < dataset.size(); i++)
        {
            // origin
            Eigen::Vector2d point_origin = dataset[i].first; 
            points_origin.push_back(point_origin);

            // position
            Eigen::Vector2d point_position = dataset[i].second;
            points_position.push_back(point_position);

            // optimized position
            Eigen::Vector2d optimized_bearing = Eigen::Vector2d(std::cos(result.at<gtsam::Rot2>(VARIABLES_BEARINGS[i]).theta()), std::sin(result.at<gtsam::Rot2>(VARIABLES_BEARINGS[i]).theta()));
            double optimized_range = (plane_position - point_origin).dot(plane_normal) / (plane_normal.dot(optimized_bearing)) ;
            Eigen::Vector2d optimized_position = point_origin + optimized_range * optimized_bearing;
            optimized_points_position.push_back(optimized_position);
        }

        // add figure
        matplot::figure();
        matplot::hold(matplot::on);
        matplot::axis("equal");

        // plot ray and position
        for (std::size_t i = 0; i < dataset.size(); i++)
        {
            // ray
            matplot::plot({points_origin[i].x(), points_position[i].x()}, {points_origin[i].y(), points_position[i].y()}, "r-");

            // position
            matplot::plot({points_position[i].x()}, {points_position[i].y()}, "ro");

            // optimized ray
            matplot::plot({points_origin[i].x(), optimized_points_position[i].x()}, {points_origin[i].y(), optimized_points_position[i].y()}, "g-");

            // optimized position
            matplot::plot({optimized_points_position[i].x()}, {optimized_points_position[i].y()}, "go");
        }
        
        // plot optimized plane
        Eigen::Vector2d plane_normal_rotated = Eigen::Vector2d(-plane_normal.y(), plane_normal.x());
        Eigen::Vector2d line_start = plane_position - 3 * plane_normal_rotated;
        Eigen::Vector2d line_end = plane_position + 3 * plane_normal_rotated;
        matplot::plot({line_start.x(), line_end.x()}, {line_start.y(), line_end.y()}, "g--");
        
        // show
        matplot::show();
    }

    // // LOG
    // graph.print("Factor graph:\n");
    // result.print("Final result:\n");
    // gtsam::Marginals marginals(graph, result);
    // std::cout.precision(2);
    // std::cout << "distance covariance:\n" << marginals.marginalCovariance(VARIABLE_t_plane) << std::endl;
    // std::cout << "normal covariance:\n" << marginals.marginalCovariance(VARIABLE_n_plane) << std::endl;
}