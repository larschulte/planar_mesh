#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <gtsam/geometry/Rot2.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>

class RangeFactor : public gtsam::NoiseModelFactor3<double, gtsam::Rot2, gtsam::Rot2>
{
public:
    RangeFactor(gtsam::Key KEY_pp_distance, gtsam::Key KEY_pp_normal, gtsam::Key KEY_v_point, Eigen::Vector2d po_plane, Eigen::Vector2d v_plane, Eigen::Vector2d po_point, double t_point, const gtsam::SharedNoiseModel &model);
    gtsam::Vector evaluateError(const double& _t_plane, const gtsam::Rot2& _n_plane, const gtsam::Rot2& _v_point, boost::optional<gtsam::Matrix &> H_distance = boost::none, boost::optional<gtsam::Matrix &> H_normal = boost::none, boost::optional<gtsam::Matrix &> H_bearing = boost::none) const override;

private:
    // plane 
    Eigen::Vector2d po_plane_; // plane origin
    Eigen::Vector2d v_plane_; // plane direction

    // point
    Eigen::Vector2d po_point_; // point origin
    double t_point_; // point direction
};

void fit_plane_to_points(const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>& dataset, Eigen::Vector2d &plane_position, Eigen::Vector2d &plane_normal, double bearing_noise, double range_noise);