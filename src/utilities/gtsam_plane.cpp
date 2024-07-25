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

class ProjectionFactor : public gtsam::NoiseModelFactor3<double, gtsam::Unit3, gtsam::Unit3>
{
private:

    // plane 
    Eigen::Vector3d po_plane_; // plane origin
    Eigen::Vector3d v_plane_; // plane direction

    // point
    Eigen::Vector3d po_point_; // point origin
    double t_point_; // point direction

public:
    ProjectionFactor(gtsam::Key KEY_pp_distance, gtsam::Key KEY_pp_normal, gtsam::Key KEY_v_point, Eigen::Vector3d po_plane, Eigen::Vector3d v_plane, Eigen::Vector3d po_point, double t_point, const gtsam::SharedNoiseModel &model)
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
    // MEASUREMENT
    Eigen::Vector3d po_point1(0, 0, 0); Eigen::Vector3d pp_point1(0, 0, 1); 
    Eigen::Vector3d po_point2(0, 0, 0); Eigen::Vector3d pp_point2(1, 0, 1); 
    Eigen::Vector3d po_point3(0, 0, 0); Eigen::Vector3d pp_point3(0, 1, 1); 
    Eigen::Vector3d po_point4(0, 0, 0); Eigen::Vector3d pp_point4(1, 1, 1); 
    double MEASUREMENT_range1 = (pp_point1 - po_point1).norm();
    double MEASUREMENT_range2 = (pp_point2 - po_point2).norm();
    double MEASUREMENT_range3 = (pp_point3 - po_point3).norm();
    double MEASUREMENT_range4 = (pp_point4 - po_point4).norm();
    gtsam::Unit3 MEASUREMENT_bearing1((pp_point1 - po_point1).normalized()); 
    gtsam::Unit3 MEASUREMENT_bearing2((pp_point2 - po_point2).normalized());
    gtsam::Unit3 MEASUREMENT_bearing3((pp_point3 - po_point3).normalized());
    gtsam::Unit3 MEASUREMENT_bearing4((pp_point4 - po_point4).normalized());

    // CONSTRAINT
    // sample point to use
    Eigen::Vector3d po_sample = po_point1;
    Eigen::Vector3d pp_sample = pp_point1;
    Eigen::Vector3d CONSTRAINT_po_plane = po_sample;
    Eigen::Vector3d CONSTRAINT_v_plane = (pp_sample - po_sample).normalized();

    // VARIABLES
    gtsam::Symbol VARIABLE_t_plane('t', 1);
    gtsam::Symbol VARIABLE_n_plane('n', 1);
    gtsam::Symbol VARIABLE_bearing1('v', 1);
    gtsam::Symbol VARIABLE_bearing2('v', 2);
    gtsam::Symbol VARIABLE_bearing3('v', 3);
    gtsam::Symbol VARIABLE_bearing4('v', 4);

    // INITIAL GUESS
    gtsam::Values initial;
    double INITIAL_t_plane = (pp_sample - po_sample).norm();
    gtsam::Unit3 INITIAL_n_plane((po_sample - pp_sample).normalized());
    initial.insert(VARIABLE_t_plane, INITIAL_t_plane);
    initial.insert(VARIABLE_n_plane, INITIAL_n_plane);
    initial.insert(VARIABLE_bearing1, MEASUREMENT_bearing1);
    initial.insert(VARIABLE_bearing2, MEASUREMENT_bearing2);
    initial.insert(VARIABLE_bearing3, MEASUREMENT_bearing3);
    initial.insert(VARIABLE_bearing4, MEASUREMENT_bearing4);

    // FACTORS
    gtsam::NonlinearFactorGraph graph;
    // bearing factor
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_bearing = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(2) << 0.01, 0.01).finished());
    graph.add(gtsam::PriorFactor<gtsam::Unit3>(VARIABLE_bearing1, MEASUREMENT_bearing1, NOISE_bearing));
    graph.add(gtsam::PriorFactor<gtsam::Unit3>(VARIABLE_bearing2, MEASUREMENT_bearing2, NOISE_bearing));
    graph.add(gtsam::PriorFactor<gtsam::Unit3>(VARIABLE_bearing3, MEASUREMENT_bearing3, NOISE_bearing));
    graph.add(gtsam::PriorFactor<gtsam::Unit3>(VARIABLE_bearing4, MEASUREMENT_bearing4, NOISE_bearing));
    // range factor
    gtsam::noiseModel::Diagonal::shared_ptr NOISE_range = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << 0.01).finished());
    graph.add(boost::make_shared<ProjectionFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLE_bearing1, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point1, MEASUREMENT_range1, NOISE_range));
    graph.add(boost::make_shared<ProjectionFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLE_bearing2, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point2, MEASUREMENT_range2, NOISE_range));
    graph.add(boost::make_shared<ProjectionFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLE_bearing3, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point3, MEASUREMENT_range3, NOISE_range));
    graph.add(boost::make_shared<ProjectionFactor>(VARIABLE_t_plane, VARIABLE_n_plane, VARIABLE_bearing4, CONSTRAINT_po_plane, CONSTRAINT_v_plane, po_point4, MEASUREMENT_range4, NOISE_range));
    
    // OPTIMIZATION
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimize();

    // OUTPUT
    // print graph
    graph.print("Factor graph:\n");
    // print result
    result.print("Final result:\n");
    // print covariance
    gtsam::Marginals marginals(graph, result);
    std::cout.precision(2);
    std::cout << "distance covariance:\n" << marginals.marginalCovariance(VARIABLE_t_plane) << std::endl;
    std::cout << "normal covariance:\n" << marginals.marginalCovariance(VARIABLE_n_plane) << std::endl;

    // RETURN
    return 0;
}