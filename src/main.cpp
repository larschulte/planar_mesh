#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"


// application class
template <typename PointT>
class Application
{
public:
    Application(std::string pcd_file_folder, std::string pose_file_path, double sensor_range_std) 
        : 
        data_loader_(pcd_file_folder, pose_file_path), 
        algorithm_(sensor_range_std)
        {}

    void step()
    {
        // load cloud and pose
        typename pcl::PointCloud<PointT>::Ptr new_cloud = data_loader_.get_cloud(i_);
        Eigen::Affine3d new_pose = data_loader_.get_pose(i_);

        // algorithms
        algorithm_.input(new_cloud, new_pose);

        // increment
        i_++;
    }

    typename pcl::PointCloud<PointT>::Ptr get_control_cloud() { return algorithm_.get_control_cloud(); }
    typename pcl::PointCloud<PointT>::Ptr get_old_cloud() { return algorithm_.get_old_cloud(); }
    typename pcl::PointCloud<PointT>::Ptr get_near_cloud() { return algorithm_.get_near_cloud(); }
    typename pcl::PointCloud<PointT>::Ptr get_far_cloud() { return algorithm_.get_far_cloud(); }
    
private:
    DataLoader<PointT> data_loader_;
    Algorithm<PointT> algorithm_;

    int i_ = 0;
};


// interactive viewer class
template <typename PointT>
class InteractiveViewer 
{
public:
    InteractiveViewer(Application<PointT>& app) 
        : 
        app_(app),
        viewer_(new pcl::visualization::PCLVisualizer ("3D Viewer"))
    {   
        // turn off warning
        viewer_->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
        
        // set up viewports
        viewer_->createViewPort (0.0, 0.0, 0.5, 1.0, port1);
        viewer_->setBackgroundColor (0, 0, 0, port1);
        viewer_->createViewPort (0.5, 0.0, 1.0, 1.0, port2);
        viewer_->setBackgroundColor (0, 0, 0, port2);
   
        // set up coordinate system
        viewer_->initCameraParameters();
        viewer_->addCoordinateSystem(1);

        // set up initial pointcloud
        // color
        pcl::visualization::PointCloudColorHandlerCustom<PointT> color_control_cloud(app_.get_control_cloud(), 0, 255, 0);
        pcl::visualization::PointCloudColorHandlerCustom<PointT> color_old_cloud(app_.get_old_cloud(), 0, 255, 0);
        // add to viewer
        viewer_->addPointCloud<PointT> (app_.get_control_cloud(), color_control_cloud, "control cloud", port1);
        viewer_->addPointCloud<PointT> (app_.get_old_cloud(), color_old_cloud, "old cloud", port2);
        // set point size
        viewer_->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "control cloud");
        viewer_->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "old cloud");

        // register keyboard callback
        viewer_->registerKeyboardCallback(&InteractiveViewer::keyboard_callback, *this, nullptr);

        // spin
        viewer_->spin();
    }

private:
    Application<PointT>& app_;

    pcl::visualization::PCLVisualizer::Ptr viewer_;
    int port1;
    int port2;
    
    void keyboard_callback(const pcl::visualization::KeyboardEvent &event, void*) 
    {
        bool space_down = event.getKeySym() == "space" && event.keyDown();
        if (space_down)
        {
            // step application
            app_.step();

            // update pointcloud
            pcl::visualization::PointCloudColorHandlerCustom<PointT> color_control_cloud(app_.get_control_cloud(), 0, 255, 0);
            pcl::visualization::PointCloudColorHandlerCustom<PointT> color_old_cloud(app_.get_old_cloud(), 0, 255, 0);
            viewer_->updatePointCloud<PointT>(app_.get_control_cloud(), color_control_cloud, "control cloud");
            viewer_->updatePointCloud<PointT>(app_.get_old_cloud(), color_old_cloud, "old cloud");
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

    // application
    Application<InputPointT> app(pcd_file_folder, pose_file_path, sensor_range_std);

    // viewer
    InteractiveViewer<InputPointT> viewer(app);

    return (0);
}


// add new points to old cloud if, the new point has lower uncertainty than the old counter part

// use the maximum allowed curvature on each new triangle, before sphere couldn't form (this surface is connected, but we have no information about how the surface is connected)
// if the point inside have a lower curvature, then use this lower curvature sphere to update the old point

// or, start each scan with a small curvature, assuming the surface is flat
// when the old point

// lower curvature surface needs fewer points
// higher curvature surface needs more points

// curvature update

// we want large curvature to capture details
// we want low curvature to achieve lower number of points on flat surface

// need a way to determine what curvature to use during update

// before that, need to first have to way to add more points into the scan

// the per scan update philosophy allows simple creation of delaunay triangulation, since points within each scan can be assumed to originate from the same origin, making 
// delaunaytion easier