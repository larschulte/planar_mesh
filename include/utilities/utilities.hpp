#pragma once

#include <tuple>
#include <Eigen/Dense>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>

std::tuple<int, int, int> valueToJet(double value);

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen);

// is_point_in_triangle(surface_coordinate0, surface_coordinate1, surface_coordinate2, surface_coordinate)
bool is_point_in_triangle(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c, const Eigen::Vector2d& p);

Eigen::Matrix3d generate_orthogonal_basis(const Eigen::Vector3d& unit_vector);
Eigen::Matrix3d generate_unit_vector_covariance(const Eigen::Vector3d& unit_vector, double angular_uncertainty_in_radian, double epsilon);

double shortest_distance_to_line_segment(const Eigen::Vector3d& rayOrigin, const Eigen::Vector3d& rayEnd, const Eigen::Vector3d& targetPoint);

#include "utilities/utilities.tpp"