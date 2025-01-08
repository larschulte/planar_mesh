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
    dataset_map["nottinghill"] = std::make_pair(
        "/home/jiahao/datasets/2024-10-23_14-40-02_rec004_nottinghill-rad101/slam_clouds/",
        "/home/jiahao/datasets/2024-10-23_14-40-02_rec004_nottinghill-rad101/slam_pose_graph.slam"
    );
    dataset_map["office"] = std::make_pair(
        "/home/jiahao/datasets/2024-11-19_14-05-08_rec013/slam_clouds/",
        "/home/jiahao/datasets/2024-11-19_14-05-08_rec013/slam_pose_graph.slam"
    );
    dataset_map["abingdon"] = std::make_pair(
        "/home/jiahao/datasets/abingdon logs/2024-12-03_10-25-53_rec001_rad301_run2/slam_clouds/",
        "/home/jiahao/datasets/abingdon logs/2024-12-03_10-25-53_rec001_rad301_run2/slam_pose_graph.slam"
    );

    std::string dataset = "nottinghill";
    data_loader_settings.pcd_file_folder = dataset_map[dataset].first;
    data_loader_settings.pose_file_path = dataset_map[dataset].second;

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
        radius_ratio = tan(6 * M_PI / 180);
        data_loader_settings.remove_double_return_flag = true;
        data_loader_settings.filter_low_intensity_flag = true;
        data_loader_settings.azimuth_resolution = 0.6;
        data_loader_settings.altitude_resolution = 1.5;
    }
    data_loader_settings.start_cloud = 0;
    data_loader_settings.start_point = 0;

    process_every_n_points = 1;
    
    fit_plane_threshold = 3;
    shuffle_pointcloud = true;
    use_radius_value = true;
    pointcloud_fraction = 1;
    radius_value = 1;
    duplicated_point_distance_threshold = 0.0; // if two points are closer than this distance, they are considered the same point

    abnormal_size = 1.5;
    envelope_size = 3.5;

    num_of_delete_before_put_to_repeated_queue = 2;
    
    num_threads = 4;
    record_countent_surface_count = false;

    use_queue = true;
    retry_threshold = 10;
    num_iterations = 10;

    // min_face_angle = 20;

    // log settings
    bool show_all = false;
    if (show_all)
    {
        log.add_point_by_radius_search = true;
        log.load_point_cloud = true;
        log.step = true;
        log.refine_surfaces = true;
        log.process_point = true;
        log.initialize = true;
        log.deletion = true;
        log.review_surfaces = true;
        log.connect_by_edges_and_faces = true;
        log.can_merge = true;
        log.merge_surface = true;    
        log.duplicated_point = true;
        log.num_of_concurrent_processes = true;
        log.total_processed_points = true;
        log.show_contented_surface = true;
    }
    else
    {
        log.add_point_by_radius_search = false;
        log.load_point_cloud = true;
        log.step = false;
        log.refine_surfaces = false;
        log.process_point = false;
        log.initialize = false;
        log.deletion = false;
        log.review_surfaces = false;
        log.connect_by_edges_and_faces = false;
        log.can_merge = false;
        log.merge_surface = false;
        log.duplicated_point = false;
        log.num_of_concurrent_processes = false;
        log.total_processed_points = false;
        log.show_contented_surface = false;
    }

    // output time
    output_time = false;
    output_file_name = "with_sah_threads_more_iteration_" + std::to_string(num_threads) + ".txt";
    turn_off_sah = false;

    // interactive viewer settings
    update_display = true;
    flip_color = false;
    show_internal_vertices = true;
    show_generic_points = false;
    show_interior_points = true;
    show_pointcloud = true;
    show_triangle = true;
    show_edge = true;
    show_confirmed_only = false;
    show_keycode = false;
    show_singular_edge = false;
    show_singular_vertex = false;
    color_mode = ColorMode::ID;
    point_mode = PointMode::USED;
    surface_denominator = 10.0;
    siblings_denominator = 3.0;
    radius_denominator = 0.3;
    positional_uncertainty_denominator = 0.002;
    contention_denominator = 500;
    show_wireframe = true;
    show_sphere = false;
    number_of_spheres_to_display = 60;
}

bool Settings::edge_is_short_enough(const double& edge_length, const double& radius0, const double& radius1) const
{
    const double minimum_edge_length = 0.0f;
    return edge_length < radius0 + minimum_edge_length && edge_length < radius1 + minimum_edge_length;
}

double Settings::compute_rrs_half_size(const double& radius) const
{
    const double minimum_rrs_half_size = 0.0f;
    return std::max(minimum_rrs_half_size, radius);
}