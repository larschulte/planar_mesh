#pragma once

#include <string>

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
    RADIUS,
    SIBLINGS,
    SURFACE_UNCERTAINTY,
    CONTENTION
};

struct Settings 
{
    Settings();

    // application settings
    bool use_sim_data;
    int sim_object;
    double range_precision;
    double range_accuracy;
    std::string cloud_path;
    std::string pose_path;
    int start_cloud;
    std::size_t start_point;
    std::size_t fit_plane_threshold;
    bool shuffle_pointcloud;
    bool use_radius_value;
    double pointcloud_fraction;
    double radius_value;
    double radius_ratio; // distance to radius ratio
    std::size_t process_every_n_points;
    double duplicated_point_distance_threshold;

    double abnormal_size; // number of std
    double envelope_size; // number of std

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
    bool flip_color;
    bool show_generic_points;
    bool show_interior_points;
    bool show_pointcloud;
    bool show_triangle;
    bool show_edge;
    bool show_projected_point;
    bool show_confirmed_only;
    bool show_keycode;
    bool show_singular_edge;
    bool show_singular_vertex;
    ColorMode color_mode;
    double surface_denominator;
    double siblings_denominator;
    double radius_denominator;
    double positional_uncertainty_denominator;
    double contention_denominator;
    bool show_wireframe;
    bool show_sphere;
    int number_of_spheres_to_display;  
};