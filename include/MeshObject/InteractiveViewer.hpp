#pragma once

#include "MeshObject/Application.hpp"
#include <pcl/visualization/pcl_visualizer.h>
#include "MeshObject/Settings.hpp"

template <typename PointT>
class InteractiveViewer 
{
public:
    InteractiveViewer(Application<PointT>& app);

private:
    // functions
    void update_display(bool export_ply = false);
    void save_simplified_surfaces();
    void keyboard_callback(const pcl::visualization::KeyboardEvent &event, void*);

    // objects
    Application<PointT>& app_;
    std::shared_ptr<Storage> storage_;
    pcl::visualization::PCLVisualizer::Ptr viewer_;

    // settings
    Settings settings_;

    // temp storage
    std::vector<std::string> sphere_name_list;
};
