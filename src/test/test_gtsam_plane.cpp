#include "utilities/gtsam_plane.hpp"
#include <Eigen/Dense>
#include <iostream>

int main()
{
    // DATASET
    // origin + position
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> dataset
    {
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1)},
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 1)},
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 1)},
        {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 2)}
    };
    
    // PLANE
    Eigen::Vector3d plane_position;
    Eigen::Vector3d plane_normal;

    // FIT PLANE
    double bearing_noise = 0.01;
    double range_noise = 0.01;
    fit_plane_to_points(dataset, plane_position, plane_normal, bearing_noise, range_noise);

    // OUTPUT
    std::cout << "plane_position: " << plane_position.transpose() << std::endl;
    std::cout << "plane_normal: " << plane_normal.transpose() << std::endl;

    // END
    return 0;
}