#pragma once

#include "MeshObject/Application.hpp"
#include <pcl/visualization/pcl_visualizer.h>

template <typename PointT>
class InteractiveViewer 
{
public:
    InteractiveViewer(Application<PointT>& app);

private:
    // functions
    void update_display();
    void keyboard_callback(const pcl::visualization::KeyboardEvent &event, void*);

    // objects
    Application<PointT>& app_;
    pcl::visualization::PCLVisualizer::Ptr viewer_;

    // settings
    bool show_generic_points = true;
    bool show_interior_points = true;
    bool show_pointcloud = true;
    bool show_triangle = true;
    bool show_edge = true;
    bool show_projected_point = false;
    bool show_error_color = false;
    bool show_wireframe = true;
    bool show_sphere = false;
    int number_of_spheres_to_display;

    // temp storage
    std::vector<std::string> sphere_name_list;
};
