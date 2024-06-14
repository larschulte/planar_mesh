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

    void load_dataset(std::string pcd_file_folder, std::string pose_file_path);
    DataLoader(std::string pcd_file_folder, std::string pose_file_path);

    typename pcl::PointCloud<PointT>::Ptr get_cloud(int i);
    Eigen::Affine3d get_pose(int i);
    int size();

private:
    std::vector<std::string> pcd_file_list_;
    std::map<std::string, Eigen::Affine3d> file_to_pose_map_;
};