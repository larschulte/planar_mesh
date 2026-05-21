#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <map>
#include "MeshObject/Settings.hpp"

std::vector<std::string> read_under_folder(std::string pcd_file_folder);
Eigen::Affine3d find_pose(const std::string& pcd_file, const std::string& pose_file);
std::map<std::string, Eigen::Affine3d> create_file_to_pose_map(std::vector<std::string> pcd_file_list, std::string pose_file_path);

template<typename PointT> 
typename pcl::PointCloud<PointT>::Ptr load_pointcloud(std::string pcd_file);

template <typename PointT>
class DataLoader
{
public:
    DataLoader();

    void load_dataset(DataLoader_Settings settings);

    typename pcl::PointCloud<PointT>::Ptr get_cloud(int i);
    Eigen::Affine3d get_pose(int i);
    int size();

private:
    std::vector<std::string> pcd_file_list_;
    std::map<std::string, Eigen::Affine3d> file_to_pose_map_;

    // settings
    double azimuth_resolution_;
    double altitude_resolution_;
    bool remove_double_return_flag_;
    bool filter_low_intensity_flag_;

    typename pcl::PointCloud<PointT>::Ptr remove_double_return(typename pcl::PointCloud<PointT>::Ptr input_pointcloud);
    typename pcl::PointCloud<PointT>::Ptr remove_double_return_2(typename pcl::PointCloud<PointT>::Ptr input_pointcloud);
    typename pcl::PointCloud<PointT>::Ptr filter_low_intensity(typename pcl::PointCloud<PointT>::Ptr input_pointcloud);
};