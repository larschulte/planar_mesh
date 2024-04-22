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

        // store output
        control_cloud_ = algorithm_.get_control_cloud();
        old_cloud_ = algorithm_.get_old_cloud();
        near_cloud_ = algorithm_.get_near_cloud();
        far_cloud_ = algorithm_.get_far_cloud();

        // output 
        i_++;
    }
    
    typename pcl::PointCloud<PointT>::Ptr control_cloud_;
    typename pcl::PointCloud<PointT>::Ptr old_cloud_;
    typename pcl::PointCloud<PointT>::Ptr near_cloud_;
    typename pcl::PointCloud<PointT>::Ptr far_cloud_;

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

            // remove old visual
            if (viewer_->contains("control cloud"))
            {
                viewer_->removePointCloud("control cloud");
            }
            if (viewer_->contains("oldcloud mean"))
            {
                viewer_->removePointCloud("oldcloud mean");
            }

            // add new visual
            add_to_viewer<PointT>(viewer_, port1, app_.control_cloud_, "control cloud", color_tuple(0, 255, 0), 3); // b
            add_to_viewer<PointT>(viewer_, port2, app_.old_cloud_, "oldcloud mean", color_tuple(0, 255, 0), 3); // b
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