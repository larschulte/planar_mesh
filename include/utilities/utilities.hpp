#pragma once

#include <tuple>
#include <Eigen/Dense>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>

std::tuple<int, int, int> valueToJet(float value);

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen);

// is_point_in_triangle(surface_coordinate0, surface_coordinate1, surface_coordinate2, surface_coordinate)
bool is_point_in_triangle(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c, const Eigen::Vector2d& p);

#include "utilities/utilities.tpp"