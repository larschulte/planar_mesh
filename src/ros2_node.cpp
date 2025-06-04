// cpp
#include <iostream>

// third party
#include <pcl/point_types.h>
#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/string.hpp"

// project
#include "point_type/VilensPointT.hpp"
#include "MeshObject/Application.hpp"

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    
    // print hello world
    std::cout << "Hello, ROS2 World!" << std::endl;
    
    return 0;
}