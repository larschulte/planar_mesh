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
    dataset_map["math"] = std::make_pair(
        "/home/jiahao/datasets/math/individual_clouds/",
        "/home/jiahao/datasets/math/slam_pose_graph.g2o"
    );
    
    std::string dataset = "math";
    cloud_path = dataset_map[dataset].first;
    pose_path = dataset_map[dataset].second;

    use_sim_data = false;
    sim_object = 0;
    std::string lidar = "64";
    if (lidar == "32")
    {
        range_precision = 0.005;
        range_accuracy = 0.01; 
        radius_ratio = tan(4 * M_PI / 180);
    }
    else if (lidar == "64")
    {
        range_precision = 0.02;
        range_accuracy = 0.03; 
        radius_ratio = tan(4 * M_PI / 180);
    }
    process_every_n_points = 10;
    
    start_cloud = 50;
    start_point = 0;
    fit_plane_threshold = 4;
    shuffle_pointcloud = false;
    use_radius_value = false;
    pointcloud_fraction = 1;
    radius_value = 2;

    abnormal_size = 1.5;
    envelope_size = 3.5;

    // min_face_angle = 20;

    // log settings
    log.add_point_by_radius_search = false;
    log.load_point_cloud = true;
    log.step = true;
    log.refine_surfaces = false;
    log.process_point = false;
    log.initialize = false;
    log.deletion = false;
    log.review_surfaces = false;
    log.connect_by_edges_and_faces = false;
    log.can_merge = false;
    log.merge_surface = false;

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
    siblings_denominator = 3.0;
    radius_denominator = 0.5;
    positional_uncertainty_denominator = 0.002;
    show_wireframe = true;
    show_sphere = false;
    number_of_spheres_to_display = 60;
}
