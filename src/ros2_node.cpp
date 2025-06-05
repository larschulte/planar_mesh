// cpp
#include <iostream>

// third party
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/string.hpp"
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <tf2_ros/transform_listener.h>
#include <tf2_ros/buffer.h>

// project
#include "point_type/VilensPointT.hpp"
#include "MeshObject/Application.hpp"

class PlanarMeshNode : public rclcpp::Node
{
public:
    PlanarMeshNode(): 
        Node("planarmesh_node")
    {
        RCLCPP_INFO(this->get_logger(), "PlanarMeshNode has been started");

        // settings
        const std::string subscription_topic_name = "/vilens/point_cloud_transformed_processed";
        const std::string publisher_topic_name = "planarmesh_output";
        const std::string base_frame_name = "odom_vilens";
        
        // create subscription and publisher
        subscription_ = this->create_subscription<sensor_msgs::msg::PointCloud2>(subscription_topic_name, 1, std::bind(&PlanarMeshNode::pointcloud_callback, this, std::placeholders::_1));
        publisher_ = this->create_publisher<sensor_msgs::msg::PointCloud2>(publisher_topic_name, 10);

        // create tf buffer and listener
        target_frame_ = base_frame_name;
        tf_buffer_ = std::make_unique<tf2_ros::Buffer>(this->get_clock());
        tf_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf_buffer_);        
    }

private:
    void pointcloud_callback(const sensor_msgs::msg::PointCloud2::SharedPtr msg)
    {
        // Process the incoming point cloud message
        unsigned int num_points = msg->width * msg->height;
        RCLCPP_INFO(this->get_logger(), "Received point cloud with %u points", num_points);

        // get frame id of the msg
        std::string frame_id = msg->header.frame_id;
        RCLCPP_INFO(this->get_logger(), "Point cloud frame id: %s", frame_id.c_str());
        
        // get pose of the msg by querying the tf buffer
        geometry_msgs::msg::TransformStamped transform;
        try {
            transform = tf_buffer_->lookupTransform(target_frame_, frame_id, msg->header.stamp);
            RCLCPP_INFO(this->get_logger(), "Transform from %s to %s found", frame_id.c_str(), target_frame_.c_str());
        } catch (const tf2::TransformException &ex) {
            RCLCPP_ERROR(this->get_logger(), "Could not transform %s to %s: %s", frame_id.c_str(), target_frame_.c_str(), ex.what());
            return;
        }

        // convert to pcl point cloud and eigen affine3d pose
        pcl::PointCloud<VilensPointT>::Ptr cloud(new pcl::PointCloud<VilensPointT>());
        pcl::fromROSMsg(*msg, *cloud);

        // convert geometry_msgs::msg::TransformStamped to Eigen::Affine3d
        Eigen::Affine3d eigen_transform;
        eigen_transform.translation() = Eigen::Vector3d(
            transform.transform.translation.x,
            transform.transform.translation.y,
            transform.transform.translation.z);
        eigen_transform.linear() = Eigen::Quaterniond(
            transform.transform.rotation.w,
            transform.transform.rotation.x,
            transform.transform.rotation.y,
            transform.transform.rotation.z).toRotationMatrix();
        
        // feed application with the point cloud and transform
        app_.load_pointcloud(cloud, eigen_transform);
        app_.process_pointcloud();

        // get the processed point cloud
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr output_pointcloud;
        app_.get_output_pointcloud(output_pointcloud);

        // publish the processed point cloud
        sensor_msgs::msg::PointCloud2 output_msg;
        pcl::toROSMsg(*output_pointcloud, output_msg);
        output_msg.header.frame_id = target_frame_;
        output_msg.header.stamp = msg->header.stamp;
        publisher_->publish(output_msg);

        // info
        RCLCPP_INFO(this->get_logger(), "Published processed point cloud with %zu points", output_pointcloud->size());
    }

    // ros related
    rclcpp::Subscription<sensor_msgs::msg::PointCloud2>::SharedPtr subscription_;
    rclcpp::Publisher<sensor_msgs::msg::PointCloud2>::SharedPtr publisher_;
    std::unique_ptr<tf2_ros::Buffer> tf_buffer_;
    std::shared_ptr<tf2_ros::TransformListener> tf_listener_{nullptr};
    std::string target_frame_;

    // algortihm related
    Application<VilensPointT> app_;
};

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<PlanarMeshNode>());
    rclcpp::shutdown();
    return 0;
}