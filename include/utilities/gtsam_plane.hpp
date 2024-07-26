#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <gtsam/geometry/Unit3.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>

class RangeFactor : public gtsam::NoiseModelFactor3<double, gtsam::Unit3, gtsam::Unit3>
{
public:
    RangeFactor(gtsam::Key KEY_pp_distance, gtsam::Key KEY_pp_normal, gtsam::Key KEY_v_point, Eigen::Vector3d po_plane, Eigen::Vector3d v_plane, Eigen::Vector3d po_point, double t_point, const gtsam::SharedNoiseModel &model);
    gtsam::Vector evaluateError(const double& _t_plane, const gtsam::Unit3& _n_plane, const gtsam::Unit3& _v_point, boost::optional<gtsam::Matrix &> H_distance = boost::none, boost::optional<gtsam::Matrix &> H_normal = boost::none, boost::optional<gtsam::Matrix &> H_bearing = boost::none) const override;

private:
    // plane 
    Eigen::Vector3d po_plane_; // plane origin
    Eigen::Vector3d v_plane_; // plane direction

    // point
    Eigen::Vector3d po_point_; // point origin
    double t_point_; // point direction
};

void fit_plane_to_points(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> dataset, Eigen::Vector3d &plane_position, Eigen::Vector3d &plane_normal, double bearing_noise, double range_noise);