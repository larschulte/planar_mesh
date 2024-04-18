#pragma once
#define PCL_NO_PRECOMPILE

#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "header_library/matplotlibcpp.h"

namespace plt = matplotlibcpp;

// define color tuple
typedef std::tuple<int, int, int> color_tuple;

template <typename PointT> void 
add_to_viewer(pcl::visualization::PCLVisualizer::Ptr viewer, int viewport, typename pcl::PointCloud<PointT>::Ptr cloud, std::string cloud_name, color_tuple color, float point_size)
{
    // color
    pcl::visualization::PointCloudColorHandlerCustom<PointT> color_handler(cloud, std::get<0>(color), std::get<1>(color), std::get<2>(color));

    // add to viewer
    viewer->addPointCloud<PointT> (cloud, color_handler, cloud_name, viewport);

    // set point size
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, point_size, cloud_name);
}

// add triangles to viewer
template<typename PointT>
void viewer_add_triangles(pcl::visualization::PCLVisualizer::Ptr viewer, int viewport, typename pcl::PointCloud<PointT>::Ptr cloud, typename pcl::PointCloud<PointT>::Ptr transformed_cloud, delaunator::Delaunator d)
{
    // plot triangles
    for(std::size_t i = 0; i < d.triangles.size(); i+=30) {
        // indices
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        // get local points
        PointT p0, p1, p2;
        p0 = cloud->points[i1];
        p1 = cloud->points[i2];
        p2 = cloud->points[i3];

        // filter triangles using normal
        // edge
        Eigen::Vector3f e1(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
        Eigen::Vector3f e2(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
        Eigen::Vector3f e3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
        // view direction
        Eigen::Vector3f observation_direction = Eigen::Vector3f(p0.x, p0.y, p0.z);

        // comptue edge to view direction angle
        double angle1 = acos(std::abs(e1.dot(observation_direction)) / (e1.norm() * observation_direction.norm()));
        double angle2 = acos(std::abs(e2.dot(observation_direction)) / (e2.norm() * observation_direction.norm()));
        double angle3 = acos(std::abs(e3.dot(observation_direction)) / (e3.norm() * observation_direction.norm()));

        // convert angle to degree
        angle1 = angle1 * 180 / M_PI;
        angle2 = angle2 * 180 / M_PI;
        angle3 = angle3 * 180 / M_PI;

        // if any angle is less than 20 degree, skip
        double angle_threshold = 20;
        if (angle1 < angle_threshold || angle2 < angle_threshold || angle3 < angle_threshold)
        {
            continue;
        }

        // get global points
        PointT p0_global, p1_global, p2_global;
        p0_global = transformed_cloud->points[i1];
        p1_global = transformed_cloud->points[i2];
        p2_global = transformed_cloud->points[i3];

        // pcl - plot the triangle
        viewer->addLine(p0_global, p1_global, 1, 1, 1, std::to_string(i) + "_line1", viewport);
        viewer->addLine(p1_global, p2_global, 1, 1, 1, std::to_string(i) + "_line2", viewport);
        viewer->addLine(p2_global, p0_global, 1, 1, 1, std::to_string(i) + "_line3", viewport);

        // message
        std::cout << "adding triange " << i <<  " out of " << d.triangles.size() << std::endl;
    }
}

void plt_plot_black_background()
{
    std::vector<double> background_x_list = {-180, 180, 180, -180, -180};
    std::vector<double> background_y_list = {-90, -90, 90, 90, -90};
    std::map<std::string, std::string> keywords = {{"color", "black"}};
    plt::fill(background_x_list, background_y_list, keywords);

    plt::xlabel("azimuth");
    plt::ylabel("altitude");   
}

// plot triangles in plt
void plt_plot_triangles(delaunator::Delaunator d)
{
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        // indices
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        // plt - plot the triangle
        float tx0 = d.coords[2 * i1];
        float ty0 = d.coords[2 * i1 + 1];
        float tx1 = d.coords[2 * i2];
        float ty1 = d.coords[2 * i2 + 1];
        float tx2 = d.coords[2 * i3];
        float ty2 = d.coords[2 * i3 + 1];
        plt::plot({tx0, tx1, tx2, tx0}, {ty0, ty1, ty2, ty0}, std::map<std::string, std::string>{{"linewidth", "0.1"}, {"color", "white"}});
        
        // message
        std::cout << "adding triange " << i <<  " out of " << d.triangles.size() << std::endl;
    }   

    // reverse x axis order
    plt::xlim(180, -180);
}

// plot coordinates in plt
void plt_scatter_plot_coords(delaunator::Delaunator d, std::vector<double> colors)
{
    // extract x, y coordinates
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    for (std::size_t i = 0; i < d.coords.size(); i+=2)
    {
        x_coords.push_back(d.coords[i]);
        y_coords.push_back(d.coords[i + 1]);
    }

    // scatter plots
    plt::scatter_colored(x_coords, y_coords, colors, 3);

    // reverse x axis order
    plt::xlim(180, -180);
}