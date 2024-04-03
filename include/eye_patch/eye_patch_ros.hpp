#pragma once
#define PCL_NO_PRECOMPILE // must define PCL_NO_PRECOMPILE before including any PCL templates when using custom point type

// #include "surfel_map/LidarPointT.hpp"
// #include "surfel_map/EyePointT.hpp"
// #include "surfel_map/ros_param_tools.hpp"

#include <ros/ros.h>
#include <tf2_ros/transform_listener.h>
// #include "eye_patch/eye_patch_algo.hpp"

#include "eye_patch/VilensPointT.hpp"
#include "eye_patch/EyePointT.hpp"
#include "eye_patch/BagPointT.hpp"

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/common/transforms.h>
#include <eigen_conversions/eigen_msg.h>
#include <Eigen/Geometry>
#include <iostream>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

namespace eye_patch
{
    class EyePatchRos
    {
    public:
        EyePatchRos(ros::NodeHandle nh);
        ~EyePatchRos();

    private:
        void callback(const sensor_msgs::PointCloud2::ConstPtr &callback_pointcloud_msg_constptr);
        void publish(const pcl::PointCloud<EyePointT>::Ptr input_pointcloud_ptr, ros::Time stamp, std::string frame_id);

        // node handle for input and output
        ros::NodeHandle nh_;

        // input
        tf2_ros::Buffer tfBuffer_;
        tf2_ros::TransformListener tfListener_;
        ros::Subscriber sub_vilens_pointcloud_;
        

        // // algorithm
        // std::unique_ptr<EyePatchAlgo> algo_;

        // output
        ros::Publisher pub_pointcloud_;
    };
}


int main(int argc, char **argv)
{
    // node
    ros::init(argc, argv, "eye_patch_ros_node");
    ros::NodeHandle nh("~");

    // start object
    eye_patch::EyePatchRos eye_patch_ros_object(nh);

    // keep node running
    ros::spin();
    return 0;
}