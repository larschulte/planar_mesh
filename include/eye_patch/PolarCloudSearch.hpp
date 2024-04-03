#pragma once

#include <Eigen/Geometry>
#include <pcl_conversions/pcl_conversions.h>
#include "eye_patch/EyePointT.hpp"
#include <pcl/common/common.h>

struct GridIndexBounds 
{
    int start_azimuth;
    int end_azimuth;
    int start_altitude;
    int end_altitude;
};


// template class
template <typename PointT>
class PolarCloudSearch
{
public:
    PolarCloudSearch();
    void setCloudPointer(pcl::PointCloud<PointT>::Ptr old_pointcloud);
    void setResolution(float res_azimuth, float res_altitude);
    void prepareSearch(const Eigen::Isometry3f &T_world2lidar);
    void addPoint(const PointT& lidar_point, const Eigen::Isometry3f &T_world2lidar);

    pcl::Indices extractNearbyIndices(const PointT& lidar_point,
                                      const Eigen::Isometry3f &T_world2lidar);

private:
    pcl::PointCloud<PointT>::Ptr old_pointcloud_;
    std::vector<std::vector<pcl::Indices>> grid_;

    void clearGrid();
    Eigen::Vector3f computeSphericalCoordinate(const Eigen::Vector3f &XL);
    GridIndexBounds calculateGridIndexBounds(const PointT& surfel, const Eigen::Isometry3f& T_world2lidar);
    void addIndexToBounds(const GridIndexBounds& bounds, int point_index);
    pcl::Indices extractIndicesFromBounds(const GridIndexBounds& bounds);

    float res_azimuth_;
    float res_altitude_;
    int number_of_azimuth_; // first dimension of polar_image_
    int number_of_altitude_; // second dimension of polar_image_
    int index_of_azimuth(float azimuth_angle);
    int index_of_altitude(float altitude_angle);
};
