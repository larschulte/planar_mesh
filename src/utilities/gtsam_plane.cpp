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

class ProjectionFactor : public gtsam::NoiseModelFactor2<double, gtsam::Unit3>
{
private:

    // plane 
    Eigen::Vector3d o_plane_; // plane origin
    Eigen::Vector3d v_plane_; // plane direction

    // point
    Eigen::Vector3d o_point_; // point origin
    Eigen::Vector3d p_point_; // point position
    Eigen::Vector3d v_point_; // point direction
    double t_point_; // point distance

public:
    ProjectionFactor(gtsam::Key KEY_p_distance, gtsam::Key KEY_p_normal, Eigen::Vector3d o_plane, Eigen::Vector3d v_plane, Eigen::Vector3d o_point, Eigen::Vector3d p_point, const gtsam::SharedNoiseModel &model)
        : 
        gtsam::NoiseModelFactor2<double, gtsam::Unit3>(model, KEY_p_distance, KEY_p_normal),
        o_plane_(o_plane),
        v_plane_(v_plane),
        o_point_(o_point),
        p_point_(p_point)
    {
        v_point_ = (p_point - o_point).normalized();
        t_point_ = (p_point - o_point).norm();
    }

    gtsam::Vector evaluateError(const double& _t_plane, const gtsam::Unit3& _n_plane , boost::optional<gtsam::Matrix &> H_distance = boost::none, boost::optional<gtsam::Matrix &> H_normal = boost::none) const override
    {
        // h(q)
        Eigen::Vector3d n_plane = _n_plane.unitVector();
        Eigen::Vector3d p_position = o_plane_ + v_plane_ * _t_plane;
        double h = n_plane.dot(p_position - o_point_) / n_plane.dot(v_point_);

        // e
        double e = h - t_point_;

        // error
        gtsam::Vector error(1);
        error << e;

        // Compute Jacobians
        if (H_distance) 
        {
            // "e" = n_plane.dot(o_plane_ + v_plane_ * "_t_plane" - o_point_) / n_plane.dot(v_point_) - t_point_
            double d_e__d_t_plane = n_plane.dot(v_plane_) / n_plane.dot(v_point_);

            // fill in the Jacobian
            *H_distance = gtsam::Matrix(1, 1);
            (*H_distance)(0, 0) = d_e__d_t_plane;
        }

        if (H_normal) 
        {
            // "e" = "n_plane".dot(p_position - o_point_) / "n_plane".dot(v_point_) - t_point_
            double v = n_plane.dot(v_point_);
            double u = n_plane.dot(p_position - o_point_);
            Eigen::Vector3d d_v__d_n_plane = v_point_;
            Eigen::Vector3d d_u__d_n_plane = p_position - o_point_;
            Eigen::Vector3d d_e__d_n_plane = (v * d_u__d_n_plane - u * d_v__d_n_plane) / (v * v);

            // fill in the Jacobian
            *H_normal = gtsam::Matrix(1, 2);
            (*H_normal).block<1, 2>(0, 0) = d_e__d_n_plane.transpose() * _n_plane.basis();
        }

        return error;
    }
};

int main()
{
    // Create an empty nonlinear factor graph
    gtsam::NonlinearFactorGraph graph;
    
    // MEASUREMENT
    Eigen::Vector3d x1_origin(0, 0, 0); Eigen::Vector3d x1_position(0, 0, 1);
    Eigen::Vector3d x2_origin(0, 0, 0); Eigen::Vector3d x2_position(1, 0, 1);
    Eigen::Vector3d x3_origin(0, 0, 0); Eigen::Vector3d x3_position(0, 1, 1);

    // PLANE
    // key
    gtsam::Symbol KEY_p_distance('p', 1);
    gtsam::Symbol KEY_p_normal('p', 2);
    // sample point to use
    Eigen::Vector3d x0_origin = x1_origin;
    Eigen::Vector3d x0_position = x1_position;
    // contraint
    Eigen::Vector3d p_origin = x0_origin;
    Eigen::Vector3d p_direction = (x0_position - x0_origin).normalized();
    // initial guess
    gtsam::Values initial;
    gtsam::Unit3 initial_normal(1, 1, 1);
    double initial_distance = (x0_position - x0_origin).norm();
    initial.insert(KEY_p_distance, initial_distance);
    initial.insert(KEY_p_normal, initial_normal);

    // FACTORS
    gtsam::noiseModel::Diagonal::shared_ptr range_noise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << 0.01).finished());
    graph.add(boost::make_shared<ProjectionFactor>(KEY_p_distance, KEY_p_normal, p_origin, p_direction, x1_origin, x1_position, range_noise));
    graph.add(boost::make_shared<ProjectionFactor>(KEY_p_distance, KEY_p_normal, p_origin, p_direction, x2_origin, x2_position, range_noise));
    graph.add(boost::make_shared<ProjectionFactor>(KEY_p_distance, KEY_p_normal, p_origin, p_direction, x3_origin, x3_position, range_noise));

    // print graph
    graph.print("Factor graph:\n");

    // Optimize using Levenberg-Marquardt optimization
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimize();

    // print result
    result.print("Final result:\n");

    // Query the marginals
    std::cout.precision(2);
    gtsam::Marginals marginals(graph, result);
    std::cout << "distance covariance:\n" << marginals.marginalCovariance(KEY_p_distance) << std::endl;
    std::cout << "normal covariance:\n" << marginals.marginalCovariance(KEY_p_normal) << std::endl;

    return 0;
}