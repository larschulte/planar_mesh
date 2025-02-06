#pragma once

#include <string>
#include <ostream>

struct LOG
{
    bool add_point_by_radius_search;
    bool load_point_cloud;
    bool step;
    bool refine_surfaces;
    bool process_point;
    bool initialize;
    bool deletion;
    bool review_surfaces;
    bool connect_by_edges_and_faces;
    bool can_merge;
    bool merge_surface;
    bool duplicated_point;
    bool num_of_concurrent_processes;
    bool total_processed_points;
    bool show_contented_surface;
};

enum class ColorMode
{
    ID,
    POSITIONAL_UNCERTAINTY,
    POSITIONAL_UNCERTAINTY_NORMALIZED,
    PROJECTED_UNCERTAINTY,
    RADIUS,
    SIBLINGS,
    SURFACE_UNCERTAINTY,
    MAX_DISTANCE_TRAVELLED,
    DISTANCE_TRAVELLED,
    WEIGHT
};

enum class PointMode
{
    USED,
    ORIGINAL,
    PROJECTED
};

inline std::ostream& operator<<(std::ostream& os, const PointMode& mode)
{
    switch (mode)
    {
        case PointMode::USED:
            os << "   USED   |          |           |";
            break;
        case PointMode::ORIGINAL:
            os << "          | ORIGINAL |           |";
            break;
        case PointMode::PROJECTED:
            os << "          |          | PROJECTED |";
            break;
        default:
            os << "not set";
            break;
    }
    return os;
}

struct DataLoader_Settings
{
    std::string pcd_file_folder; 
    std::string pose_file_path; 
    double azimuth_resolution; 
    double altitude_resolution;
    bool remove_double_return_flag;
    bool filter_low_intensity_flag;
    int start_cloud;
    int start_point;
};

struct Settings 
{
    Settings();

    bool edge_is_short_enough(const double& edge_length, const double& radius0, const double& radius1) const;
    double compute_rrs_half_size(const double& radius) const;

    // data loader settings
    DataLoader_Settings data_loader_settings;

    // application settings
    bool headless_mode;
    int num_scans;
    std::string save_folder;
    bool use_sim_data;
    int sim_object;
    double range_precision;
    double range_accuracy;    
    std::size_t fit_plane_threshold;
    bool shuffle_pointcloud;
    bool use_radius_value;
    double pointcloud_fraction;
    double radius_value;
    double extra_radius;
    double radius_ratio; // distance to radius ratio
    std::size_t process_every_n_points;
    double duplicated_point_distance_threshold;
    
    double odometry_position_uncertainty_rate;
    double odometry_angular_uncertainty_rate;

    double confidence_interval_multiplier;

    double high_incident_angle_threshold_std;

    double abnormal_size; // number of std
    double envelope_size; // number of std

    double seed_surface_area_threshold;

    unsigned int num_of_delete_before_put_to_repeated_queue;

    unsigned int num_threads;
    bool use_queue;

    bool record_countent_surface_count;

    unsigned int retry_threshold;
    unsigned int num_iterations;

    // double min_face_angle;

    // log settings
    LOG log;

    // output time
    bool output_time;
    std::string output_file_name;
    bool turn_off_sah;

    // interactive viewer settings
    bool update_display;
    bool flip_color;
    bool show_internal_vertices;
    bool show_generic_points;
    bool show_interior_points;
    bool show_pointcloud;
    bool show_triangle;
    bool show_edge;
    bool show_keycode;
    bool show_seed_surface;
    ColorMode color_mode;
    PointMode point_mode;
    double surface_denominator;
    double siblings_denominator;
    double radius_denominator;
    double positional_uncertainty_denominator;
    double contention_denominator;
    bool show_wireframe;
    bool show_sphere;
    int number_of_spheres_to_display;  
};