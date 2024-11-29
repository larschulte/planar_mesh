#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <map>

std::vector<std::string> read_under_folder(std::string pcd_file_folder);
Eigen::Affine3d find_pose(std::string pcd_file, std::string pose_file);
std::map<std::string, Eigen::Affine3d> create_file_to_pose_map(std::vector<std::string> pcd_file_list, std::string pose_file_path);

template<typename PointT> 
typename pcl::PointCloud<PointT>::Ptr load_pointcloud(std::string pcd_file);

template <typename PointT>
class DataLoader
{
public:
    DataLoader();

    void load_dataset(std::string pcd_file_folder, std::string pose_file_path, double azimuth_resolution, double altitude_resolution);
    DataLoader(std::string pcd_file_folder, std::string pose_file_path);

    typename pcl::PointCloud<PointT>::Ptr get_cloud(int i, bool remove_double_return_flag, bool filter_low_intensity_flag);
    Eigen::Affine3d get_pose(int i);
    int size();

private:
    std::vector<std::string> pcd_file_list_;
    std::map<std::string, Eigen::Affine3d> file_to_pose_map_;

    double azimuth_resolution_;
    double altitude_resolution_;

    typename pcl::PointCloud<PointT>::Ptr remove_double_return(typename pcl::PointCloud<PointT>::Ptr input_pointcloud);
    typename pcl::PointCloud<PointT>::Ptr filter_low_intensity(typename pcl::PointCloud<PointT>::Ptr input_pointcloud);
};