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
    Eigen::Vector3d po_plane_; // plane origin
    Eigen::Vector3d v_plane_; // plane direction

    // point
    Eigen::Vector3d po_point_; // point origin
    Eigen::Vector3d pp_point_; // point position
    Eigen::Vector3d v_point_; // point direction
    double t_point_; // point distance

public:
    ProjectionFactor(gtsam::Key KEY_pp_distance, gtsam::Key KEY_pp_normal, Eigen::Vector3d po_plane, Eigen::Vector3d v_plane, Eigen::Vector3d po_point, Eigen::Vector3d pp_point, const gtsam::SharedNoiseModel &model)
        : 
        gtsam::NoiseModelFactor2<double, gtsam::Unit3>(model, KEY_pp_distance, KEY_pp_normal),
        po_plane_(po_plane),
        v_plane_(v_plane),
        po_point_(po_point),
        pp_point_(pp_point)
    {
        v_point_ = (pp_point - po_point).normalized();
        t_point_ = (pp_point - po_point).norm();
    }

    gtsam::Vector evaluateError(const double& _t_plane, const gtsam::Unit3& _n_plane , boost::optional<gtsam::Matrix &> H_distance = boost::none, boost::optional<gtsam::Matrix &> H_normal = boost::none) const override
    {
        // h(q)
        Eigen::Vector3d n_plane = _n_plane.unitVector();
        Eigen::Vector3d pp_plane = po_plane_ + v_plane_ * _t_plane;
        double h = n_plane.dot(pp_plane - po_point_) / n_plane.dot(v_point_);

        // e
        double e = h - t_point_;

        // error
        gtsam::Vector error(1);
        error << e;

        // Compute Jacobians
        if (H_distance) 
        {
            // "e" = n_plane.dot(po_plane_ + v_plane_ * "_t_plane" - po_point_) / n_plane.dot(v_point_) - t_point_
            double d_e__d_t_plane = n_plane.dot(v_plane_) / n_plane.dot(v_point_);

            // fill in the Jacobian
            *H_distance = gtsam::Matrix(1, 1);
            (*H_distance)(0, 0) = d_e__d_t_plane;
        }

        if (H_normal) 
        {
            // "e" = "n_plane".dot(pp_plane - po_point_) / "n_plane".dot(v_point_) - t_point_
            double v = n_plane.dot(v_point_);
            double u = n_plane.dot(pp_plane - po_point_);
            Eigen::Vector3d d_v__d_n_plane = v_point_;
            Eigen::Vector3d d_u__d_n_plane = pp_plane - po_point_;
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
    Eigen::Vector3d po_point1(0, 0, 0); Eigen::Vector3d pp_point1(0, 0, 1);
    Eigen::Vector3d po_point2(0, 0, 0); Eigen::Vector3d pp_point2(1, 0, 1);
    Eigen::Vector3d po_point3(0, 0, 0); Eigen::Vector3d pp_point3(0, 1, 1);

    // PLANE
    // sample point to define plane location
    Eigen::Vector3d po_sample = po_point1;
    Eigen::Vector3d pp_sample = pp_point1;
    // key
    gtsam::Symbol KEY_t_plane('p', 1);
    gtsam::Symbol KEY_n_plane('p', 2);
    // contraint
    Eigen::Vector3d CONSTRAINT_po_plane = po_sample;
    Eigen::Vector3d CONSTRAINT_v_plane = (pp_sample - po_sample).normalized();
    // initial guess
    gtsam::Values initial;
    gtsam::Unit3 INITIAL_n_plane(1, 1, 1);
    double INITIAL_t_plane = (pp_sample - po_sample).norm();
    initial.insert(KEY_t_plane, INITIAL_t_plane);
    initial.insert(KEY_n_plane, INITIAL_n_plane);

    // FACTORS
    gtsam::noiseModel::Diagonal::shared_ptr range_noise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << 0.01).finished());
    graph.add(boost::make_shared<ProjectionFactor>(KEY_t_plane, KEY_n_plane, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point1, pp_point1, range_noise));
    graph.add(boost::make_shared<ProjectionFactor>(KEY_t_plane, KEY_n_plane, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point2, pp_point2, range_noise));
    graph.add(boost::make_shared<ProjectionFactor>(KEY_t_plane, KEY_n_plane, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point3, pp_point3, range_noise));

    // Optimize using Levenberg-Marquardt optimization
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimize();

    // PRINT
    // print graph
    graph.print("Factor graph:\n");
    // print result
    result.print("Final result:\n");
    // print covariance
    gtsam::Marginals marginals(graph, result);
    std::cout.precision(2);
    std::cout << "distance covariance:\n" << marginals.marginalCovariance(KEY_t_plane) << std::endl;
    std::cout << "normal covariance:\n" << marginals.marginalCovariance(KEY_n_plane) << std::endl;

    return 0;
}