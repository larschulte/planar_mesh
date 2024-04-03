#pragma once
#define PCL_NO_PRECOMPILE // must define PCL_NO_PRECOMPILE before including any PCL templates when using custom point type

#include "eye_patch/VilensPointT.hpp"
#include "eye_patch/EyePointT.hpp"

#include <pcl/filters/passthrough.h>
#include <pcl/common/transforms.h>
#include <geometry_msgs/PointStamped.h>

namespace eye_patch
{

    class EyePatchAlgo
    {
    public:
        EyePatchAlgo();
        ~EyePatchAlgo();

        void add_pointcloud(const pcl::PointCloud<VilensPointT> &callback_lidar_pointcloud, const Eigen::Isometry3f &tf_lidar2world_eigen_f);
        void getResults(pcl::PointCloud<EyePointT>::Ptr output_surfel_pointcloud_ptr_);

    };
}