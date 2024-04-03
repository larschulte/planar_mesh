#include "eye_patch/eye_patch_ros.hpp"

using namespace eye_patch;

EyePatchRos::EyePatchRos(ros::NodeHandle nh) : nh_(nh), tfListener_(tfBuffer_)
{
    // subcriber and publisher
    sub_vilens_pointcloud_ = nh_.subscribe("/vilens_pointcloud", 1, &EyePatchRos::EyePatchRos::callback, this);
    pub_pointcloud_ = nh_.advertise<sensor_msgs::PointCloud2>("/eye_patch_pointcloud", 0, true);

    // // algorithm 
    // algo_ = std::make_unique<EyePatchAlgo>();
};

EyePatchRos::~EyePatchRos() {};

void EyePatchRos::callback(const sensor_msgs::PointCloud2::ConstPtr &callback_pointcloud_msg_constptr)
{
    
    // DATA
    
    // header
    std_msgs::Header header = callback_pointcloud_msg_constptr->header; // header
    
    // pointcloud
    pcl::PointCloud<BagPointT> callback_pointcloud;
    pcl::fromROSMsg(*callback_pointcloud_msg_constptr, callback_pointcloud);
    
    // transform
    std::string source = header.frame_id;
    std::string target = "map";
    ros::Time time_stamp = header.stamp;
    if (!tfBuffer_.canTransform(target, source, time_stamp))
    {
        return;
    }
    geometry_msgs::TransformStamped tf_lidar2world_msg = tfBuffer_.lookupTransform(target, source, time_stamp);
    
    // convert to eigen transformation matrix
    Eigen::Isometry3d tf_lidar2map_eigen_d;
    tf::transformMsgToEigen(tf_lidar2world_msg.transform, tf_lidar2map_eigen_d);
    
    // // transform pointcloud to global
    // pcl::PointCloud<VilensPointT> global_pointcloud;
    // pcl::transformPointCloudWithNormals(callback_pointcloud, global_pointcloud, tf_lidar2map_eigen_d.cast<float>());


    // convert to eyepoint pointcloud with id
    pcl::PointCloud<EyePointT> eyepoint_pointcloud;
    eyepoint_pointcloud.resize(callback_pointcloud.size());

    // // grid to store points 
    // std::vector<std::vector<pcl::Indices>> grid_;
    // res_azimuth_ = 0.6;
    // res_altitude_ = 1;
    // number_of_azimuth_ = std::ceil(360.f/res_azimuth_);
    // number_of_altitude_ = std::ceil(180.f/res_altitude_);
    // grid_.resize(number_of_azimuth_, std::vector<pcl::Indices>(number_of_altitude_));

    std::vector<double> x_list;
    std::vector<double> y_list;
    std::vector<double> z_list;
    x_list.resize(callback_pointcloud.size());
    y_list.resize(callback_pointcloud.size());
    z_list.resize(callback_pointcloud.size());

    for (int i = 0; i < callback_pointcloud.size(); i++)
    {
        EyePointT eye_point;
        eye_point.x = callback_pointcloud.points[i].x;
        eye_point.y = callback_pointcloud.points[i].y;
        eye_point.z = callback_pointcloud.points[i].z;
        eye_point.id = i;
        eyepoint_pointcloud.points[i] = eye_point;

        // compute the azimuth and altitude
        float x = callback_pointcloud.points[i].x;
        float y = callback_pointcloud.points[i].y;
        float z = callback_pointcloud.points[i].z;
        float r = sqrt(x*x + y*y + z*z);
        double azimuth = atan2(y, x);
        double altitude = asin(z/r);

        // convert to degree
        azimuth = azimuth * 180 / M_PI;
        altitude = altitude * 180 / M_PI;

        // scatter plot azimuth and altitude
        // matplotlib
        // Sample data
        x_list.at(i) = azimuth;
        y_list.at(i) = altitude;
        z_list.at(i) = callback_pointcloud.points[i].intensity;

        // output to cout
        std::cout << "point " << i << " has azimuth " << azimuth << " and altitude " << altitude << std::endl;
    }

    // Create a scatter plot
    // todo - plot the bag pcd as well, and also to read pcd from file directly
    std::map<std::string, std::string> keywords;
    plt::scatter_colored(x_list, y_list, z_list, 1, keywords);
    plt::title("Sample Scatter Plot");
    plt::xlabel("x axis");
    plt::ylabel("y axis");
    plt::show();
    
    // publish world frame
    sensor_msgs::PointCloud2 msg;
    pcl::toROSMsg(eyepoint_pointcloud, msg);
    msg.header.stamp = header.stamp;
    msg.header.frame_id = "lidar";
    pub_pointcloud_.publish(msg);

    // ALGORITHM
    // pointcloud: callback_pointcloud, global_pointcloud
    // transform: tf_lidar2map_eigen_d



    // // run algorithm
    // algo_->add_pointcloud(callback_pointcloud, tf_lidar2map_eigen_d.cast<float>());
    // ROS_INFO("run algorithm complete");

    // // publish
    // callback_counts_ += 1;
    // if (callback_counts_ % 1 == 0)
    // {
    //     pcl::PointCloud<EyePointT>::Ptr publish_surfel_pointcloud_ptr(new pcl::PointCloud<EyePointT>);
    //     algo_->getResults(publish_surfel_pointcloud_ptr);

    //     publish(publish_surfel_pointcloud_ptr, header.stamp, p_ros_.frame_world);
    //     ROS_INFO("publish surfel cloud with size [%li]", publish_surfel_pointcloud_ptr->size());
    // }


    

    // // test publish
    // sensor_msgs::PointCloud2 msg = *callback_pointcloud_msg_constptr;
    // pub_pointcloud_.publish(msg);
};


// void EyePatchRos::publish(const pcl::PointCloud<EyePointT>::Ptr input_pointcloud_ptr, ros::Time stamp, std::string frame_id)
// {
//     sensor_msgs::PointCloud2 msg;
//     pcl::toROSMsg(*input_pointcloud_ptr, msg);
//     msg.header.stamp = stamp;
//     msg.header.frame_id = frame_id;
//     pub_pointcloud_.publish(msg);
// }