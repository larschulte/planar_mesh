#include "utilities/gtsam_plane_2d.hpp"
#include <Eigen/Dense>
#include <iostream>

#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/DoglegOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>

#include <matplot/matplot.h>

std::vector<Eigen::Vector2d> generate_ray_directions(Eigen::Vector2d center_ray_direction, int num_points, double min_angle_degree, double max_angle_degree)
{
    // convert to radian
    double min_angle = min_angle_degree / 180.0 * M_PI;
    double max_angle = max_angle_degree / 180.0 * M_PI;

    // normalized vectors
    center_ray_direction.normalize();

    std::vector<Eigen::Vector2d> ray_directions;
    double step_angle = (max_angle - min_angle) / (double)(num_points - 1.0);
    for (int i = 0; i < num_points; i++)
    {
        double ray_angle = min_angle + (double)i * step_angle;
        ray_directions.push_back(Eigen::Rotation2Dd(ray_angle) * center_ray_direction);
    }
    return ray_directions;
}

void add_dataset(
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& dataset_gt, 
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& dataset_noisy,
    Eigen::Vector2d gt_plane_position,
    Eigen::Vector2d gt_plane_normal,
    Eigen::Vector2d sensor_origin,
    std::vector<Eigen::Vector2d> ray_directions,
    double angle_std,
    double range_std
    )
{
    // noise generators
    std::mt19937 random_number_generator(std::random_device{}());
    std::normal_distribution<> nd_angle(0, angle_std);
    std::normal_distribution<> nd_range(0, range_std);

    // generate dataset
    for (std::size_t i = 0; i < ray_directions.size(); i++)
    {
        // ray direction
        double angle_noise = nd_angle(random_number_generator);
        Eigen::Vector2d ray_direction_gt = ray_directions[i];
        Eigen::Vector2d ray_direction_noisy = Eigen::Rotation2Dd(angle_noise) * ray_direction_gt;

        // range
        double range_noise = nd_range(random_number_generator);
        double range_gt = (gt_plane_normal.dot(gt_plane_position - sensor_origin)) / (gt_plane_normal.dot(ray_direction_noisy));
        double range_noisy = range_gt + range_noise;

        // intersection point
        Eigen::Vector2d point_gt = sensor_origin + ray_direction_noisy * range_gt;
        Eigen::Vector2d point_noisy = sensor_origin + ray_direction_gt * range_noisy;

        // add to dataset gt
        dataset_gt.push_back(std::make_pair(sensor_origin, point_gt));
        dataset_noisy.push_back(std::make_pair(sensor_origin, point_noisy));
    }
}

void fit_plane(
    const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& dataset_noisy, 
    const double& angle_std, 
    const double& range_std, 
    const Eigen::Vector2d& initial_plane_position,
    const Eigen::Vector2d& initial_plane_normal,
    Eigen::Vector2d& optimized_plane_position,
    Eigen::Vector2d& optimized_plane_normal,
    std::vector<gtsam::Rot2>& optimized_bearings
    )
{
    // CHECK NUMBER OF POINTS
    if (dataset_noisy.size() < 2) throw std::runtime_error("At least 2 points are required to fit a plane.");

    // MEASUREMENT
    std::vector<double> MEASUREMENT_RANGES;
    std::vector<gtsam::Rot2> MEASUREMENT_BEARINGS;
    for (auto &data_noisy : dataset_noisy)
    {
        MEASUREMENT_RANGES.push_back((data_noisy.second - data_noisy.first).norm());
        MEASUREMENT_BEARINGS.push_back(gtsam::Rot2::fromCosSin((data_noisy.second - data_noisy.first).normalized().x(), (data_noisy.second - data_noisy.first).normalized().y()));
    }

    // VARIABLES
    gtsam::Symbol VARIABLE_t_plane('t', 1);
    gtsam::Symbol VARIABLE_n_plane('n', 1);
    std::vector<gtsam::Symbol> VARIABLES_BEARINGS;
    for (std::size_t i = 0; i < dataset_noisy.size(); i++) VARIABLES_BEARINGS.push_back(gtsam::Symbol('v', i + 1));

    // CONSTRAINT
    Eigen::Vector2d CONSTRAINT_po_plane = dataset_noisy[0].first;
    Eigen::Vector2d CONSTRAINT_v_plane = (initial_plane_position != Eigen::Vector2d::Zero()) ? (initial_plane_position - dataset_noisy[0].first).normalized() : (dataset_noisy[0].second - dataset_noisy[0].first).normalized();
        
    // FACTORS
    gtsam::NonlinearFactorGraph graph;
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_bearing = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << angle_std).finished());
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_range = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << range_std).finished());
    for (std::size_t i = 0; i < dataset_noisy.size(); i++)
    {
        auto FACTOR_range = boost::make_shared<RangeFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLES_BEARINGS[i], CONSTRAINT_po_plane, CONSTRAINT_v_plane, dataset_noisy[i].first, MEASUREMENT_RANGES[i], NOISE_range);
        auto FACTOR_bearing = gtsam::PriorFactor<gtsam::Rot2>(VARIABLES_BEARINGS[i], MEASUREMENT_BEARINGS[i], NOISE_bearing);
        graph.add(FACTOR_range);
        graph.add(FACTOR_bearing);
    }

    // INITIAL
    gtsam::Values initial;
    initial.insert(VARIABLE_n_plane, (initial_plane_normal != Eigen::Vector2d::Zero()) ? gtsam::Rot2::fromCosSin(initial_plane_normal.x(), initial_plane_normal.y()) : gtsam::Rot2::fromCosSin((dataset_noisy[0].first - dataset_noisy[0].second).normalized().x(), (dataset_noisy[0].first - dataset_noisy[0].second).normalized().y()));
    initial.insert(VARIABLE_t_plane, (initial_plane_position != Eigen::Vector2d::Zero()) ? (initial_plane_position - dataset_noisy[0].first).norm() : (dataset_noisy[0].second - dataset_noisy[0].first).norm());
    for (std::size_t i = 0; i < dataset_noisy.size(); i++) initial.insert(VARIABLES_BEARINGS[i], MEASUREMENT_BEARINGS[i]);
    
    // OPTIMIZATION
    // gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimizeSafely();
    // gtsam::Values result = gtsam::GaussNewtonOptimizer(graph, initial).optimize();
    gtsam::Values result = gtsam::DoglegOptimizer(graph, initial).optimize();
    result.print();
    std::cout << "error: " << graph.error(result) << std::endl;

    // RESULTS
    double OPTIMIZED_t_plane = result.at<double>(VARIABLE_t_plane);
    gtsam::Rot2 OPTIMIZED_n_plane = result.at<gtsam::Rot2>(VARIABLE_n_plane);
    optimized_plane_position = CONSTRAINT_po_plane + CONSTRAINT_v_plane * OPTIMIZED_t_plane;
    optimized_plane_normal = Eigen::Vector2d(std::cos(OPTIMIZED_n_plane.theta()), std::sin(OPTIMIZED_n_plane.theta()));
    optimized_bearings.clear();
    for (std::size_t i = 0; i < dataset_noisy.size(); i++) optimized_bearings.push_back(result.at<gtsam::Rot2>(VARIABLES_BEARINGS[i]));

    // print
    std::cout << "initial plane_position: " << initial_plane_position.transpose() << std::endl;
    std::cout << "initial plane_normal: " << initial_plane_normal.transpose() << std::endl;
    std::cout << "optimized plane_position: " << optimized_plane_position.transpose() << std::endl;
    std::cout << "optimized plane_normal: " << optimized_plane_normal.transpose() << std::endl;
}

void plot_graph(
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> dataset_gt, 
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> dataset_noisy, 
    Eigen::Vector2d gt_plane_position, 
    Eigen::Vector2d gt_plane_normal, 
    Eigen::Vector2d optimized_plane_position, 
    Eigen::Vector2d optimized_plane_normal, 
    std::vector<gtsam::Rot2> optimized_bearings
    )
{
    // compute locations to plot
    std::vector<Eigen::Vector2d> points_origin;
    std::vector<Eigen::Vector2d> points_position;
    std::vector<Eigen::Vector2d> optimized_points_position;
    for (std::size_t i = 0; i < dataset_noisy.size(); i++)
    {
        // origin
        Eigen::Vector2d point_origin = dataset_noisy[i].first; 
        points_origin.push_back(point_origin);

        // position
        Eigen::Vector2d point_position = dataset_noisy[i].second;
        points_position.push_back(point_position);

        // optimized position
        Eigen::Vector2d optimized_bearing = Eigen::Vector2d(std::cos(optimized_bearings[i].theta()), std::sin(optimized_bearings[i].theta()));
        double optimized_range = (optimized_plane_position - point_origin).dot(optimized_plane_normal) / (optimized_plane_normal.dot(optimized_bearing)) ;
        Eigen::Vector2d optimized_position = point_origin + optimized_range * optimized_bearing;
        optimized_points_position.push_back(optimized_position);
    }

    // add figure
    matplot::figure();
    matplot::hold(matplot::on);
    matplot::axis("equal");
    matplot::axis("square");

    // plot ray and position
    for (std::size_t i = 0; i < dataset_noisy.size(); i++)
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
    Eigen::Vector2d optimized_plane_normal_rotated = Eigen::Vector2d(-optimized_plane_normal.y(), optimized_plane_normal.x());
    Eigen::Vector2d optimized_line_start = optimized_plane_position - 3 * optimized_plane_normal_rotated;
    Eigen::Vector2d optimized_line_end = optimized_plane_position + 3 * optimized_plane_normal_rotated;
    matplot::plot({optimized_line_start.x(), optimized_line_end.x()}, {optimized_line_start.y(), optimized_line_end.y()}, "g--");
    
    // plot gt ray and position
    for (std::size_t i = 0; i < dataset_gt.size(); i++)
    {
        // ray
        matplot::plot({dataset_gt[i].first.x(), dataset_gt[i].second.x()}, {dataset_gt[i].first.y(), dataset_gt[i].second.y()}, "k-");

        // position
        matplot::plot({dataset_gt[i].second.x()}, {dataset_gt[i].second.y()}, "ko");
    }
    
    // plot gt plane
    Eigen::Vector2d gt_plane_normal_rotated = Eigen::Vector2d(-gt_plane_normal.y(), gt_plane_normal.x());
    Eigen::Vector2d gt_line_start = gt_plane_position - 3 * gt_plane_normal_rotated;
    Eigen::Vector2d gt_line_end = gt_plane_position + 3 * gt_plane_normal_rotated;
    matplot::plot({gt_line_start.x(), gt_line_end.x()}, {gt_line_start.y(), gt_line_end.y()}, "k--");

    // show
    matplot::show();
}

int main()
{
    // settings
    Eigen::Vector2d gt_plane_position(0, 1);
    Eigen::Vector2d gt_plane_normal = Eigen::Vector2d(-1, -1).normalized();
    double angle_std = 1.0 / 180.0 * M_PI;
    double range_std = 0.01;

    // dataset
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> dataset_gt;
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> dataset_noisy;
    add_dataset(dataset_gt, dataset_noisy, gt_plane_position, gt_plane_normal, Eigen::Vector2d(0, 0), generate_ray_directions(Eigen::Vector2d(0, 1), 31, -30.0, 30.0), angle_std, range_std);
    add_dataset(dataset_gt, dataset_noisy, gt_plane_position, gt_plane_normal, Eigen::Vector2d(-1, 0), generate_ray_directions(Eigen::Vector2d(1, 1), 31, -30.0, 30.0), angle_std, range_std);

    // fit plane
    Eigen::Vector2d initial_plane_position(0, 0);
    Eigen::Vector2d initial_plane_normal(0, 1);
    Eigen::Vector2d optimized_plane_position;
    Eigen::Vector2d optimized_plane_normal;
    std::vector<gtsam::Rot2> optimized_bearings;
    fit_plane(dataset_noisy, angle_std, range_std, initial_plane_position, initial_plane_normal, optimized_plane_position, optimized_plane_normal, optimized_bearings);

    // plot graph
    plot_graph(dataset_gt, dataset_noisy, gt_plane_position, gt_plane_normal, optimized_plane_position, optimized_plane_normal, optimized_bearings);

    // END
    return 0;
}