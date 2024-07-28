#include "utilities/gtsam_plane_2d.hpp"
#include <Eigen/Dense>
#include <iostream>

int main()
{
    // DATASET
    // origin + position
    std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> dataset
    {
        {Eigen::Vector2d(0, 0), Eigen::Vector2d(0, 1)},
        {Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1)},
        {Eigen::Vector2d(0, 0), Eigen::Vector2d(2, 1.5)},
    };
    
    // PLANE
    Eigen::Vector2d plane_position;
    Eigen::Vector2d plane_normal;

    // FIT PLANE
    double bearing_noise = 0.01;
    double range_noise = 0.01;
    bool plot_graph = true;
    fit_plane_to_points(dataset, plane_position, plane_normal, bearing_noise, range_noise, plot_graph);

    // OUTPUT
    std::cout << "plane_position: " << plane_position.transpose() << std::endl;
    std::cout << "plane_normal: " << plane_normal.transpose() << std::endl;

    // END
    return 0;
}