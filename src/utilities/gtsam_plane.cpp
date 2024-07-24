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

class ProjectionFactor : public gtsam::NoiseModelFactor1<gtsam::Pose3>
{
private:
    Eigen::Vector3d o_; // origin
    Eigen::Vector3d p_; // point
    Eigen::Vector3d n_; // ray direction
    double m_; // distance

public:
    ProjectionFactor(gtsam::Key j, Eigen::Vector3d o, Eigen::Vector3d p, const gtsam::SharedNoiseModel &model)
        : 
        gtsam::NoiseModelFactor1<gtsam::Pose3>(model, j),
        o_(o),
        p_(p)
    {
        n_ = (p - o).normalized();
        m_ = (p - o).norm();
    }

    gtsam::Vector evaluateError(const gtsam::Pose3 &q, boost::optional<gtsam::Matrix &> H = boost::none) const override
    {
        // position and normal of the unkown pose
        Eigen::Vector3d p0 = q.translation();
        Eigen::Vector3d n0 = q.rotation().matrix().col(2);

        // h(q)
        double h = n0.dot(p0 - o_) / n0.dot(n_);

        // Optional: compute the Jacobian
        if (H)
        {
            Eigen::MatrixXd H_p(1, 3); // Partial derivative w.r.t. translation
            Eigen::MatrixXd H_R(1, 3); // Partial derivative w.r.t. rotation

            // Compute partial derivative w.r.t. translation
            H_p = n0.transpose() / n0.dot(n_);

            // Compute partial derivative w.r.t. rotation
            Eigen::Matrix3d skew_n0 = gtsam::skewSymmetric(n0);
            H_R = ((p0 - o_).transpose() * skew_n0) / n0.dot(n_);

            // Combine into a single Jacobian matrix
            H->resize(1, 6);
            H->block<1, 3>(0, 0) = H_p;
            H->block<1, 3>(0, 3) = H_R;
        }

        // e
        double e = h - m_;

        // error
        gtsam::Vector error(1);
        error << e;
        
        return error;
    }
};

int main()
{
    // Create an empty nonlinear factor graph
    gtsam::NonlinearFactorGraph graph;

    // // Prior factor
    // gtsam::Pose3 priorMean(gtsam::Rot3::RzRyRx(0.0, 0.0, 0.0), gtsam::Point3(0.0, 0.0, 0.0));
    // gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1).finished());
    // graph.add(boost::make_shared<gtsam::PriorFactor<gtsam::Pose3>>(1, priorMean, priorNoise));

    // // Odometry factors
    // gtsam::noiseModel::Diagonal::shared_ptr odometryNoise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.2, 0.2, 0.2, 0.1, 0.1, 0.1).finished());
    // graph.add(boost::make_shared<gtsam::BetweenFactor<gtsam::Pose3>>(1, 2, gtsam::Pose3(gtsam::Rot3(), gtsam::Point3(2.0, 0.0, 0.0)), odometryNoise));
    // graph.add(boost::make_shared<gtsam::BetweenFactor<gtsam::Pose3>>(2, 3, gtsam::Pose3(gtsam::Rot3(), gtsam::Point3(2.0, 0.0, 0.0)), odometryNoise));

    Eigen::Vector3d ray_origin(0, 0, 0);
    Eigen::Vector3d ray_direction(0, 0, 1);
    gtsam::Pose3 priorMean(gtsam::Rot3::RzRyRx(0.0, 0.0, 0.0), gtsam::Point3(ray_origin + ray_direction));
    
    // Custom noise model: High confidence along the ray, less confidence perpendicular to the ray
    gtsam::Vector6 priorSigmas;
    priorSigmas << 0.1, 0.1, 0.1, 0.01, 0.01, 0.01; // Small sigma for position along the ray direction
    gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::Sigmas(priorSigmas);
    
    graph.add(boost::make_shared<gtsam::PriorFactor<gtsam::Pose3>>(1, priorMean, priorNoise));



    // Measurement factors
    gtsam::noiseModel::Diagonal::shared_ptr projection_noise = gtsam::noiseModel::Diagonal::Sigmas((gtsam::Vector(1) << 0.01).finished());
    graph.add(boost::make_shared<ProjectionFactor>(1, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1), projection_noise));
    graph.add(boost::make_shared<ProjectionFactor>(1, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 1), projection_noise));
    graph.add(boost::make_shared<ProjectionFactor>(1, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 1), projection_noise));

    // Create (deliberately inaccurate) initial estimate
    gtsam::Values initial;
    initial.insert(1, gtsam::Pose3(gtsam::Rot3::RzRyRx(0.1, 0.1, 0.1), gtsam::Point3(0.5, 0.0, 0.2)));

    // print graph
    graph.print("Factor graph:\n");

    // Optimize using Levenberg-Marquardt optimization
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initial).optimize();

    // print result
    result.print("Final result:\n");

    // Query the marginals
    std::cout.precision(2);
    gtsam::Marginals marginals(graph, result);
    std::cout << "plane covariance:\n" << marginals.marginalCovariance(1) << std::endl;

    return 0;
}
