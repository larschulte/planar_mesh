#pragma once

#include <string>

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

    double abnormal_size; // number of std
    double envelope_size; // number of std

    // interactive viewer settings
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
    int color_mode;
    double surface_denominator;
    double siblings_denominator;
    double radius_denominator;
    double positional_uncertainty_denominator;
    bool show_wireframe;
    bool show_sphere;
    int number_of_spheres_to_display;  
};