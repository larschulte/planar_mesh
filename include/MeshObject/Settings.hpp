#pragma once

#include <string>

struct Settings 
{
    Settings();

    // application settings
    std::string cloud_path;
    std::string pose_path;
    int start_cloud;
    std::size_t start_point;
    double distance_threshold;
    std::size_t fit_plane_threshold;
    std::size_t remove_low_confidence_threshold;
    double average_projective_distance_threshold;
    double merged_eigenvalue_threshold;
    bool shuffle_pointcloud;
    bool use_radius_value;
    double pointcloud_fraction;
    double radius_value;
    double radius_ratio; // distance to radius ratio

    // interactive viewer settings
    bool show_generic_points;
    bool show_interior_points;
    bool show_pointcloud;
    bool show_triangle;
    bool show_edge;
    bool show_projected_point;
    int color_mode;
    double surface_denominator;
    bool show_wireframe;
    bool show_sphere;
    int number_of_spheres_to_display;  
};