#pragma once

#include <tuple>
#include <Eigen/Dense>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>

std::tuple<int, int, int> valueToJet(float value) 
{
    // Ensure value is within [0, 1]
    if (value < 0.0f) value = 0.0f;
    if (value > 1.0f) value = 1.0f;

    float r = 0, g = 0, b = 0;

    if (value < 0.25f) 
    {
        r = 0;
        g = 4 * value;
        b = 1;
    } 
    else if (value < 0.5f) 
    {
        r = 0;
        g = 1;
        b = 1 - 4 * (value - 0.25f);
    } 
    else if (value < 0.75f) 
    {
        r = 4 * (value - 0.5f);
        g = 1;
        b = 0;
    } 
    else 
    {
        r = 1;
        g = 1 - 4 * (value - 0.75f);
        b = 0;
    }

    return std::make_tuple(static_cast<int>(r * 255), static_cast<int>(g * 255), static_cast<int>(b * 255));
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
{
    // transform point cloud
    typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
    pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
    return transformed_cloud;
}
