#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"


template <typename PointT>
class InteractiveViewer 
{
public:
    InteractiveViewer(DataLoader<PointT>& data_loader, Algorithm<PointT>& algorithm) : 
        viewer(new pcl::visualization::PCLVisualizer ("3D Viewer")), 
        control_cloud(new typename pcl::PointCloud<PointT>()),
        data_loader(data_loader),
        algorithm(algorithm)
    {   
        // turn off warning
        viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
        
        // set up viewports
        viewer->createViewPort (0.0, 0.0, 0.5, 1.0, port1);
        viewer->setBackgroundColor (0, 0, 0, port1);
        viewer->createViewPort (0.5, 0.0, 1.0, 1.0, port2);
        viewer->setBackgroundColor (0, 0, 0, port2);
   
        // set up coordinate system
        viewer->initCameraParameters();
        viewer->addCoordinateSystem(1);
        
        // register keyboard callback
        viewer->registerKeyboardCallback(&InteractiveViewer::callback, *this, nullptr);

        // spin
        viewer->spin();
    }

private:
    pcl::visualization::PCLVisualizer::Ptr viewer;
    int port1;
    int port2;
    int i = 0;
    typename pcl::PointCloud<PointT>::Ptr control_cloud;
    DataLoader<PointT>& data_loader;
    Algorithm<PointT>& algorithm;  


    void callback(const pcl::visualization::KeyboardEvent &event, void*) 
    {
        if (event.getKeySym() == "space" && event.keyDown()) {
            // load cloud and pose
            typename pcl::PointCloud<PointT>::Ptr new_cloud = data_loader.get_cloud(i);
            Eigen::Affine3d new_pose = data_loader.get_pose(i);

            // algorithms
            *control_cloud += *transform_to_global<PointT>(new_cloud, new_pose);
            algorithm.add_pointcloud_and_pose(new_cloud, new_pose);

            // update viewer
            if (viewer->contains("control cloud"))
            {
                viewer->removePointCloud("control cloud");
            }
            if (viewer->contains("oldcloud mean"))
            {
                viewer->removePointCloud("oldcloud mean");
            }
            add_to_viewer<PointT>(viewer, port1, control_cloud, "control cloud", color_tuple(0, 255, 0), 3); // b
            add_to_viewer<PointT>(viewer, port2, algorithm.get_old_cloud(), "oldcloud mean", color_tuple(0, 255, 0), 3); // b

            // update i
            i++;
        }
    }  
};



using InputPointT = VilensPointT;
int main()
{
    // parameters
    std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
    std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";
    double sensor_range_std = 0.01;

    // data loader
    DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);
    
    // algorithm
    Algorithm<InputPointT> algorithm(sensor_range_std);

    // viewer
    InteractiveViewer<InputPointT> viewer(data_loader, algorithm);

    return (0);
}

// int main()
// {
//     // add new points to old cloud if, the new point has lower uncertainty than the old counter part
    
//     // use the maximum allowed curvature on each new triangle, before sphere couldn't form (this surface is connected, but we have no information about how the surface is connected)
//     // if the point inside have a lower curvature, then use this lower curvature sphere to update the old point

//     // or, start each scan with a small curvature, assuming the surface is flat
//     // when the old point

//     // lower curvature surface needs fewer points
//     // higher curvature surface needs more points

//     // curvature update
//     // 

//     // we want large curvature to capture details
//     // we want low curvature to achieve lower number of points on flat surface

//     // need a way to determine what curvature to use during update

//     // before that, need to first have to way to add more points into the scan
    
//     // the per scan update philosophy allows simple creation of delaunay triangulation, since points within each scan can be assumed to originate from the same origin, making 
//     // delaunaytion easier
    

//     return 0;
// }