#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

// cost functor -> a function that computes residuals
// cost function -> differentiable cost functor
// loss function -> weight scaling of residual, like huber, soft_l1, cauchy, etc.
// residual -> sum of squared residuals is minimized

struct PointData 
{
    Eigen::Vector3d point;
    Eigen::Vector3d direction;
    double sigma;
};

struct PlaneFittingCostFunctor 
{
    PlaneFittingCostFunctor(PointData point_data) : point_(point_data.point), direction_(point_data.direction.normalized()), sigma_(point_data.sigma) {}

    template <typename T>
    bool operator()(const T* const plane_normal, const T* const plane_position, T* residual) const 
    {
        // plane
        Eigen::Matrix<T, 3, 1> p0(plane_position[0], plane_position[1], plane_position[2]);
        Eigen::Matrix<T, 3, 1> n0(plane_normal[0], plane_normal[1], plane_normal[2]);
        n0.normalize();

        // point
        Eigen::Matrix<T, 3, 1> p(point_.cast<T>());
        Eigen::Matrix<T, 3, 1> u(direction_.cast<T>());
        
        // intersection distance 
        T d = (n0.dot(p0 - p)) / (n0.dot(u));

        // scaled intersection distance
        T d_effective = d / T(sigma_);

        // return residual
        residual[0] = d_effective;
        return true;
    }

private:
    const Eigen::Vector3d point_;
    const Eigen::Vector3d direction_;
    const double sigma_;
};


int main() 
{
    // data
    Eigen::Vector3d plane_normal(1, 1, 1);
    Eigen::Vector3d plane_position(0, 0, 0);
    std::vector<PointData> point_data_list;
    point_data_list.push_back({{0, 0, 0}, {1, 1, 1}, 1.0});
    point_data_list.push_back({{1, 0, 0}, {1, 1, 1}, 1.0});
    point_data_list.push_back({{0, 1, 0}, {1, 1, 1}, 1.0});
    

    // option
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

    // problem
    ceres::Problem problem;
    for (const PointData& point_data : point_data_list) 
    {
        PlaneFittingCostFunctor* cost_functor = new PlaneFittingCostFunctor(point_data);
        ceres::CostFunction* differentiable_cost_functor = new ceres::AutoDiffCostFunction<PlaneFittingCostFunctor, 1, 3, 3>(cost_functor); // problem non analytical, thus autodiff
        problem.AddResidualBlock(differentiable_cost_functor, nullptr, plane_normal.data(), plane_position.data());
    }

    // summary
    ceres::Solver::Summary summary;

    // solve
    ceres::Solve(options, &problem, &summary);

    // output
    std::cout << summary.FullReport() << std::endl;
    std::cout << "Optimal normal vector: " << plane_normal.normalized().transpose() << std::endl;
    std::cout << "Optimal point on the plane: " << plane_position.transpose() << std::endl;


    return 0;
}
