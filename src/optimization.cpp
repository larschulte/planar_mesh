#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <ceres/rotation.h>

struct PointData 
{
    Eigen::Vector3d point;
    Eigen::Vector3d direction;
    double sigma;
};

struct PlaneFittingCostFunctor 
{
    PlaneFittingCostFunctor(const Eigen::Vector3d& point, const Eigen::Vector3d& direction, double sigma)
        : point_(point), direction_(direction), sigma_(sigma) {}

    template <typename T>
    bool operator()(const T* const plane_normal, const T* const plane_point, T* residual) const 
    {
        Eigen::Matrix<T, 3, 1> n(plane_normal[0], plane_normal[1], plane_normal[2]);
        Eigen::Matrix<T, 3, 1> p0(plane_point[0], plane_point[1], plane_point[2]);

        n.normalize();

        Eigen::Matrix<T, 3, 1> p(point_.cast<T>());
        Eigen::Matrix<T, 3, 1> u(direction_.cast<T>());

        T t = (n.dot(p0 - p)) / (n.dot(u));
        Eigen::Matrix<T, 3, 1> intersection_point = p + t * u;

        T d_eff = (intersection_point - p).norm();

        residual[0] = d_eff / T(sigma_);

        return true;
    }

private:
    const Eigen::Vector3d point_;
    const Eigen::Vector3d direction_;
    const double sigma_;
};


int main() 
{
    std::vector<PointData> points = 
    {
        {{1, 2, 3}, {1, 0, 0}, 0.1},
        {{4, 5, 6}, {0, 1, 0}, 0.2},
        {{7, 8, 9}, {0, 0, 1}, 0.3}
    };

    Eigen::Vector3d initial_normal(1, 1, 1);
    Eigen::Vector3d initial_point_on_plane(0, 0, 0);

    initial_normal.normalize();

    ceres::Problem problem;
    for (const auto& point_data : points) 
    {
        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<PlaneFittingCostFunctor, 1, 3, 3>(
                new PlaneFittingCostFunctor(point_data.point, point_data.direction, point_data.sigma)),
            nullptr,
            initial_normal.data(),
            initial_point_on_plane.data()
        );
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;

    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << std::endl;
    std::cout << "Optimal normal vector: " << initial_normal.transpose() << std::endl;
    std::cout << "Optimal point on the plane: " << initial_point_on_plane.transpose() << std::endl;

    return 0;
}
