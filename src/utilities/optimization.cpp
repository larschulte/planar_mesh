#include "utilities/optimization.hpp"

int main() 
{
    // data
    std::set<int> pointIDs = {0, 1, 2};
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    point_to_vector3d_map[0] = Eigen::Vector3d(0, 0, 1);
    point_to_vector3d_map[1] = Eigen::Vector3d(1, 0, 1);
    point_to_vector3d_map[2] = Eigen::Vector3d(0, 1, 1);
    std::map<int, Eigen::Vector3d> point_to_origin_vector3d_map;
    point_to_origin_vector3d_map[0] = Eigen::Vector3d(0, 0, 0);
    point_to_origin_vector3d_map[1] = Eigen::Vector3d(0, 0, 0);
    point_to_origin_vector3d_map[2] = Eigen::Vector3d(0, 0, 0);

    // plane
    Eigen::Vector3d plane_normal(0, 1, 1);
    Eigen::Vector3d plane_position(0, 0, 0);

    // fit
    fit_plane_to_lidar_points(pointIDs, point_to_vector3d_map, point_to_origin_vector3d_map, plane_normal, plane_position);

    return 0;
}
