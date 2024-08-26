#pragma once

#include "utilities/utilities.hpp"

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
{
    // transform point cloud
    typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
    pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
    return transformed_cloud;
}

template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}