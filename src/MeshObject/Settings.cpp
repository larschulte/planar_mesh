#include "MeshObject/Settings.hpp"
#include <map>
#include <math.h>

Settings::Settings() 
{
    // application settings
    std::map<std::string, std::pair<std::string, std::string>> dataset_map;
    dataset_map["room"] = std::make_pair(
        "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/",
        "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam"
    );
    dataset_map["osney"] = std::make_pair(
        "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/",
        "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam"
    );
    dataset_map["blenheim"] = std::make_pair(
        "/home/jiahao/datasets/2024-03-14-09-09-02-lenord-walk-for-lintong/individual_clouds/",
        "/home/jiahao/datasets/2024-03-14-09-09-02-lenord-walk-for-lintong/slam_pose_graph.g2o"
    );
    dataset_map["christchurch"] = std::make_pair(
        "/home/jiahao/datasets/christ church spires/2024-03-18-09-43-37/individual_clouds/",
        "/home/jiahao/datasets/christ church spires/2024-03-18-09-43-37/slam_pose_graph.g2o"
    );
    
    std::string dataset = "room";
    cloud_path = dataset_map[dataset].first;
    pose_path = dataset_map[dataset].second;

    use_sim_data = true;
    sim_object = 2;
    
    start_cloud = 50;
    start_point = 0;
    distance_threshold = 0.05;
    fit_plane_threshold = 3;
    remove_low_confidence_threshold = 10;
    projective_std_threshold = 0.1;
    merged_eigenvalue_threshold = 15e-5;
    shuffle_pointcloud = false;
    use_radius_value = true;
    pointcloud_fraction = 1;
    radius_value = 2;
    radius_ratio = tan(4 * M_PI / 180);
    range_noise_std = 0.001;
    
    // interactive viewer settings
    show_generic_points = true;
    show_interior_points = true;
    show_pointcloud = true;
    show_triangle = true;
    show_edge = true;
    show_projected_point = false;
    show_confirmed_only = false;
    show_keycode = false;
    show_singular_edge = false;
    show_singular_vertex = false;
    color_mode = 0;
    surface_denominator = 10.0;
    siblings_denominator = 5.0;
    radius_denominator = 2.0;
    show_wireframe = true;
    show_sphere = false;
    number_of_spheres_to_display = 60;
}
