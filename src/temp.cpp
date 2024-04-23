#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"

template <typename PointT>
pcl::PolygonMesh obtain_mesh(const typename pcl::PointCloud<PointT>::Ptr cloud, const delaunator::Delaunator& d)
{
    // initialize
    pcl::PolygonMesh mesh;

    // process
    // store triangles
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        pcl::Vertices v;
        v.vertices.push_back(d.triangles[i]);
        v.vertices.push_back(d.triangles[i + 1]);
        v.vertices.push_back(d.triangles[i + 2]);
        mesh.polygons.push_back(v);
    }
    // store pointcloud
    // pcl::PCLPointCloud2 pcl_pc2;
    // pcl::toPCLPointCloud2(*cloud, pcl_pc2);
    // mesh.cloud = pcl_pc2;
    pcl::toPCLPointCloud2(*cloud, mesh.cloud);

    // return
    return mesh;
}


// using InputPointT = VilensPointT;
// int main()
// {
//     std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
//     std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";

//     DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);

//     // // ------------------------------------------------------- PLT
//     // // viewer
//     // pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     // viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
//     // viewer->initCameraParameters();
//     // viewer->addCoordinateSystem(1);

//     // // load cloud and pose
//     // int i1 = 0;
//     // typename pcl::PointCloud<InputPointT>::Ptr cloud1 = data_loader.get_cloud(i1);
//     // Eigen::Affine3d pose1 = data_loader.get_pose(i1);
//     // typename pcl::PointCloud<InputPointT>::Ptr cloud1_transformed = transform_to_global<InputPointT> (cloud1, pose1);
//     // pcl::visualization::PointCloudColorHandlerCustom<InputPointT> color1(cloud1_transformed, 0, 255, 0);
//     // viewer->addPointCloud<InputPointT> (cloud1_transformed, color1, "cloud1");
//     // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud1");

//     // int i2 = 60;
//     // typename pcl::PointCloud<InputPointT>::Ptr cloud2 = data_loader.get_cloud(i2);
//     // Eigen::Affine3d pose2 = data_loader.get_pose(i2);
//     // typename pcl::PointCloud<InputPointT>::Ptr cloud2_transformed = transform_to_global<InputPointT> (cloud2, pose2);
//     // pcl::visualization::PointCloudColorHandlerCustom<InputPointT> color2(cloud2_transformed, 255, 0, 0);
//     // viewer->addPointCloud<InputPointT> (cloud2_transformed, color2, "cloud2");
//     // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud2");

//     // // spin
//     // viewer->spin();


//     // // // ------------------------------------------------------- PLT
//     // // plt triangulation
//     // plt_plot_black_background();

//     // // plt plot triangles
//     // for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
//     //     // indices
//     //     int i1 = d.triangles[i];
//     //     int i2 = d.triangles[i + 1];
//     //     int i3 = d.triangles[i + 2];

//     //     // plt - plot the triangle
//     //     float tx0 = d.coords[2 * i1];
//     //     float ty0 = d.coords[2 * i1 + 1];
//     //     float tx1 = d.coords[2 * i2];
//     //     float ty1 = d.coords[2 * i2 + 1];
//     //     float tx2 = d.coords[2 * i3];
//     //     float ty2 = d.coords[2 * i3 + 1];
//     //     plt::plot({tx0, tx1, tx2, tx0}, {ty0, ty1, ty2, ty0}, std::map<std::string, std::string>{{"linewidth", "0.1"}, {"color", "white"}});
        
//     //     // message
//     //     std::cout << "adding triange " << i <<  " out of " << d.triangles.size() << std::endl;
//     // }

//     // // reverse x axis order
//     // plt::xlim(180, -180);

//     // // show
//     // plt::show();

    
//     int i1 = 0;
//     typename pcl::PointCloud<InputPointT>::Ptr old_cloud = data_loader.get_cloud(i1);
//     Eigen::Affine3d old_pose = data_loader.get_pose(i1);

//     int i2 = 50;
//     typename pcl::PointCloud<InputPointT>::Ptr new_cloud = data_loader.get_cloud(i2);
//     Eigen::Affine3d new_pose = data_loader.get_pose(i2);

//     // obtain_triangulation
//     delaunator::Delaunator old_d = obtain_triangulation<InputPointT> (old_cloud);
//     pcl::PolygonMesh old_mesh = obtain_mesh<InputPointT>(old_cloud, old_d);
    
    
//     float sensor_range_std = 0.1;
//     std::vector<Eigen::Vector3f> old_cloud_direction = compute_point_directions<InputPointT>(old_cloud);
//     std::vector<float> old_cloud_variance(old_cloud->size(), std::pow(sensor_range_std, 2));

//     // add cloud 2 to cloud 1 mesh
//     typename pcl::PointCloud<InputPointT>::Ptr old_cloud_local = transform_to_frame<PointT>(old_cloud, Eigen::Isometry3d::Identity(), new_pose);
//     update_pointcloud<InputPointT>(old_cloud_local, old_cloud_direction, old_cloud_variance, new_cloud, new_cloud_direction, new_cloud_variance, sensor_range_std);
//     typename pcl::PointCloud<InputPointT>::Ptr old_cloud_global = transform_to_frame<PointT>(old_cloud_local, new_pose, Eigen::Isometry3d::Identity());
//     *old_cloud = *old_cloud_global;







//     // end add cloud 2 to cloud 1 mesh




//     // viewer
//     pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
//     viewer->initCameraParameters();
//     viewer->addCoordinateSystem(1);
//     // // add mesh
//     // viewer->addPolygonMesh(old_mesh, "mesh");
//     // add polyline
//     viewer->addPolylineFromPolygonMesh(old_mesh, "polyline");

//     // spin
//     viewer->spin();

//     return 0;
// }



using InputPointT = VilensPointT;
int main()
{
    std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
    std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";

    DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);

    // obtain pointcloud and direction
    int i1 = 100;
    typename pcl::PointCloud<InputPointT>::Ptr old_cloud = data_loader.get_cloud(i1);
    Eigen::Affine3d old_pose = data_loader.get_pose(i1);
    std::vector<Eigen::Vector3f> old_cloud_direction = compute_point_directions<InputPointT>(old_cloud);

    // transform
    typename pcl::PointCloud<InputPointT>::Ptr old_cloud_transformed = transform_cloud_to_global<InputPointT> (old_cloud, old_pose);
    std::vector<Eigen::Vector3f> old_cloud_direction_transformed = transform_direction_to_global(old_cloud_direction, old_pose);
    

    // transform origin
    Eigen::Vector3f origin = Eigen::Vector3f::Zero();
    Eigen::Vector3f origin_transformed = old_pose.cast<float>() * origin;

    // output old_pose
    std::cout << "old_pose: " << old_pose.matrix() << std::endl;



    // visualize point and direction

    // viewer
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    viewer->initCameraParameters();

    // // add point with color and size
    // pcl::visualization::PointCloudColorHandlerCustom<InputPointT> old_color(old_cloud_transformed, 0, 255, 0);
    // viewer->addPointCloud<InputPointT> (old_cloud_transformed, old_color, "cloud1");
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud1");


    // convert to viewer cloud
    typename pcl::PointCloud<InputPointT>::Ptr cloud_to_use = old_cloud_transformed;
    std::vector<Eigen::Vector3f> direction_to_use = old_cloud_direction_transformed;
    Eigen::Vector3f origin_to_use = origin_transformed;
    // typename pcl::PointCloud<InputPointT>::Ptr cloud_to_use = old_cloud;
    // std::vector<Eigen::Vector3f> direction_to_use = old_cloud_direction;
    // Eigen::Vector3f origin_to_use = origin;

    pcl::PointCloud<pcl::PointXYZINormal>::Ptr viewer_cloud(new pcl::PointCloud<pcl::PointXYZINormal>);
    viewer_cloud->resize(cloud_to_use->size());
    for (std::size_t i = 0; i < cloud_to_use->size(); i++)
    {
        pcl::PointXYZINormal viewer_point;
        viewer_point.x = cloud_to_use->points[i].x;
        viewer_point.y = cloud_to_use->points[i].y;
        viewer_point.z = cloud_to_use->points[i].z;
        viewer_point.intensity = 0;
        viewer_point.normal_x = -direction_to_use[i][0];
        viewer_point.normal_y = -direction_to_use[i][1];
        viewer_point.normal_z = -direction_to_use[i][2];
        viewer_cloud->points[i] = viewer_point;
    }
    // add to viewer
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZINormal> viewer_cloud_color(viewer_cloud, 255, 0, 0);
    viewer->addPointCloud<pcl::PointXYZINormal> (viewer_cloud, viewer_cloud_color, "pointcloud");
    viewer->addPointCloudNormals<pcl::PointXYZINormal> (viewer_cloud, 1, 0.05, "normals");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "pointcloud");

    viewer->addCoordinateSystem(0.5);
    viewer->addCoordinateSystem(1, origin_to_use[0], origin_to_use[1], origin_to_use[2], "pose");

    // spin
    viewer->spin();

    return 0;
}