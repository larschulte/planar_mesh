#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"


using InputPointT = VilensPointT;
int main()
{
    // parameters
    std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
    std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";
    std::size_t number_of_files = 150;
    double sensor_range_std = 0.01;

    // data loader
    DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);
    
    // algorithm
    Algorithm<InputPointT> algorithm(sensor_range_std);

    // control cloud (for comparing with updated old cloud)
    pcl::PointCloud<InputPointT>::Ptr control_cloud (new pcl::PointCloud<InputPointT>);

    for (std::size_t i = 0; i < number_of_files; i++)
    {
        // load cloud and pose
        pcl::PointCloud<InputPointT>::Ptr new_cloud = data_loader.get_cloud(i);
        Eigen::Affine3d new_pose = data_loader.get_pose(i);

        // algorithms
        *control_cloud += *transform_to_global<InputPointT>(new_cloud, new_pose);
        algorithm.add_pointcloud_and_pose(new_cloud, new_pose);

        // message
        std::cout << i << " old cloud number of points: " << algorithm.get_old_cloud()->size() << std::endl;
    }

    // the current update assume planar surface within each triangle, and does not filter the planar surface even if the triangle is very large
    // this will be solved when introducing eye patch

    // ------------------------------------------------------ pclvisuliazer    
    // set up viewer
    pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    
    // set up viewports
    int port1(0);
    viewer->createViewPort (0.0, 0.0, 0.5, 1.0, port1);
    viewer->setBackgroundColor (0, 0, 0, port1);
    int port2(0);
    viewer->createViewPort (0.5, 0.0, 1.0, 1.0, port2);
    viewer->setBackgroundColor (0, 0, 0, port2);

    // // set up viewports
    // int port1(0);
    // viewer->createViewPort (0.0, 0.0, 1, 1.0, port1);
    // viewer->setBackgroundColor (0, 0, 0, port1);
    
    // set up coordinate system
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1);

    // add_to_viewer<InputPointT>(viewer, port1, old_cloud_near, "oldcloud near updated", color_tuple(0, 255, 0), 2); // y
    // add_to_viewer<InputPointT>(viewer, port1, old_cloud_far, "oldcloud far updated", color_tuple(255, 0, 0), 2); // g

    // // add to viewer
    // add_to_viewer<InputPointT>(viewer, port2, control_cloud, "control cloud", color_tuple(255, 0, 0), 3); // b
    // add_to_viewer<InputPointT>(viewer, port2, old_cloud_mean, "oldcloud mean", color_tuple(0, 255, 0), 3); // b

    // add to viewer
    add_to_viewer<InputPointT>(viewer, port1, control_cloud, "control cloud", color_tuple(0, 255, 0), 3); // b
    add_to_viewer<InputPointT>(viewer, port2, algorithm.get_old_cloud(), "oldcloud mean", color_tuple(0, 255, 0), 3); // b
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "control cloud");
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "oldcloud mean");
    
    
    
    // display
    viewer->spin();

    return (0);
}