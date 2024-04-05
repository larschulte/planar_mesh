#include <iostream>

#include "eye_patch/BagPointT.hpp"
#include "eye_patch/VilensPointT.hpp"


#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

#include <delaunator.hpp>
#include <pcl/visualization/pcl_visualizer.h>
#include <string>

#include "eye_patch/matplotlibcpp.h"
namespace plt = matplotlibcpp;

#include <pcl/common/transforms.h>
#include <pcl/kdtree/kdtree_flann.h>


// define color tuple
typedef std::tuple<int, int, int> color_tuple;


template <typename PointT> void 
add_to_viewer(pcl::visualization::PCLVisualizer::Ptr viewer, int viewport, typename pcl::PointCloud<PointT>::Ptr cloud, std::string cloud_name, color_tuple color, float point_size)
{
    // convert to xyz cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr xyz_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    xyz_cloud->resize(cloud->size());
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        pcl::PointXYZ xyz_point;
        xyz_point.x = cloud->points[i].x;
        xyz_point.y = cloud->points[i].y;
        xyz_point.z = cloud->points[i].z;
        xyz_cloud->points[i] = xyz_point;
    }

    // color
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> color_handler(xyz_cloud, std::get<0>(color), std::get<1>(color), std::get<2>(color));

    // add to viewer
    viewer->addPointCloud<pcl::PointXYZ> (xyz_cloud, color_handler, cloud_name, viewport);

    // set point size
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, point_size, cloud_name);
}


template<typename PointT> typename pcl::PointCloud<PointT>::Ptr 
load_pointcloud(std::string pcd_file)
{
    // load the pcd file
    typename pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
    if (pcl::io::loadPCDFile<PointT> (pcd_file, *cloud) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
    }
    std::cout << "Loaded " << cloud->size() << " points" << std::endl;
    return cloud;
}

Eigen::Affine3d 
find_pose(std::string pcd_file, std::string pose_file)
{
    // get sec and nsec from pcd_file, using split
    std::string pcd_file_name = pcd_file.substr(pcd_file.find_last_of("/\\") + 1);
    std::string sec_str = pcd_file_name.substr(6, 10);
    std::string nsec_str = pcd_file_name.substr(17, 9);
    std::cout << "sec: " << sec_str << std::endl;
    std::cout << "nsec: " << nsec_str << std::endl;

    // find pose 
    std::ifstream pose_stream(pose_file);
    std::string line;
    double pose_x, pose_y, pose_z, pose_qx, pose_qy, pose_qz, pose_qw;
    bool pose_found = false;
    while (std::getline(pose_stream, line))
    {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "VERTEX_SE3:QUAT_TIME")
        {
            int vertex_id;
            double x, y, z, qx, qy, qz, qw;
            int timestamp_sec, timestamp_nsec;
            iss >> vertex_id >> x >> y >> z >> qx >> qy >> qz >> qw >> timestamp_sec >> timestamp_nsec;

            if (timestamp_sec == std::stoi(sec_str) && timestamp_nsec == std::stoi(nsec_str))
            {
                std::cout << "found pose" << std::endl;
                pose_x = x;
                pose_y = y;
                pose_z = z;
                pose_qx = qx;
                pose_qy = qy;
                pose_qz = qz;
                pose_qw = qw;
                pose_found = true;
                break;
            }
        }
    }
    if (!pose_found)
    {
        std::cout << "pose not found" << std::endl;
        return Eigen::Isometry3d::Identity();
    }

    // convert pose to eigen::isometry3d
    Eigen::Affine3d pose_eigen = Eigen::Isometry3d::Identity();
    pose_eigen.translation() << pose_x, pose_y, pose_z;
    Eigen::Quaterniond q(pose_qw, pose_qx, pose_qy, pose_qz);
    pose_eigen.rotate(q);

    return pose_eigen;
}

template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
{
    // transform point cloud
    typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
    pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
    return transformed_cloud;
}

// triangulation
// triangles list stores [trig1-vertex1, trig1-vertex2, trig1-vertex3, trig2-vertex1, trig2-vertex2, trig2-vertex3, ...]
template<typename PointT>
delaunator::Delaunator obtain_triangulation(typename pcl::PointCloud<PointT>::Ptr cloud)
{

    // compute azimuth and altitude coordinates
    std::vector<double>* coords = new std::vector<double>();
    coords->reserve(cloud->size() * 2);
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        float x = cloud->points[i].x;
        float y = cloud->points[i].y;
        float z = cloud->points[i].z;
        
        double r = sqrt(x * x + y * y + z * z);
        double azimuth = atan2(y, x) * 180 / M_PI;
        double altitude = asin(z / r) * 180 / M_PI;

        coords->push_back(azimuth);
        coords->push_back(altitude);
    }

    // create delaunay triangles given x and y
    delaunator::Delaunator d(*coords);
    std::cout << "Number of triangles: " << d.triangles.size() << std::endl;

    return d;
}


// add triangles to viewer
template<typename PointT>
void add_triangles_to_viewer(pcl::visualization::PCLVisualizer::Ptr viewer, int viewport, typename pcl::PointCloud<PointT>::Ptr cloud, typename pcl::PointCloud<PointT>::Ptr transformed_cloud, delaunator::Delaunator d)
{
    // plot triangles
    for(std::size_t i = 0; i < d.triangles.size(); i+=30) {
        // indices
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        // get local points
        PointT p0, p1, p2;
        p0 = cloud->points[i1];
        p1 = cloud->points[i2];
        p2 = cloud->points[i3];

        // filter triangles using normal
        // edge
        Eigen::Vector3f e1(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
        Eigen::Vector3f e2(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
        Eigen::Vector3f e3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
        // view direction
        Eigen::Vector3f observation_direction = Eigen::Vector3f(p0.x, p0.y, p0.z);

        // comptue edge to view direction angle
        double angle1 = acos(std::abs(e1.dot(observation_direction)) / (e1.norm() * observation_direction.norm()));
        double angle2 = acos(std::abs(e2.dot(observation_direction)) / (e2.norm() * observation_direction.norm()));
        double angle3 = acos(std::abs(e3.dot(observation_direction)) / (e3.norm() * observation_direction.norm()));

        // convert angle to degree
        angle1 = angle1 * 180 / M_PI;
        angle2 = angle2 * 180 / M_PI;
        angle3 = angle3 * 180 / M_PI;

        // if any angle is less than 20 degree, skip
        double angle_threshold = 20;
        if (angle1 < angle_threshold || angle2 < angle_threshold || angle3 < angle_threshold)
        {
            continue;
        }

        // get global points
        PointT p0_global, p1_global, p2_global;
        p0_global = transformed_cloud->points[i1];
        p1_global = transformed_cloud->points[i2];
        p2_global = transformed_cloud->points[i3];

        // pcl - plot the triangle
        viewer->addLine(p0_global, p1_global, 1, 1, 1, std::to_string(i) + "_line1", viewport);
        viewer->addLine(p1_global, p2_global, 1, 1, 1, std::to_string(i) + "_line2", viewport);
        viewer->addLine(p2_global, p0_global, 1, 1, 1, std::to_string(i) + "_line3", viewport);

        // message
        std::cout << "adding triange " << i <<  " out of " << d.triangles.size() << std::endl;
    }
}

void plt_plot_black_background()
{
    std::vector<double> background_x_list = {-180, 180, 180, -180, -180};
    std::vector<double> background_y_list = {-90, -90, 90, 90, -90};
    std::map<std::string, std::string> keywords = {{"color", "black"}};
    plt::fill(background_x_list, background_y_list, keywords);

    plt::xlabel("azimuth");
    plt::ylabel("altitude");   
}

// plot triangles in plt
void plt_plot_triangles(delaunator::Delaunator d)
{
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        // indices
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        // plt - plot the triangle
        float tx0 = d.coords[2 * i1];
        float ty0 = d.coords[2 * i1 + 1];
        float tx1 = d.coords[2 * i2];
        float ty1 = d.coords[2 * i2 + 1];
        float tx2 = d.coords[2 * i3];
        float ty2 = d.coords[2 * i3 + 1];
        plt::plot({tx0, tx1, tx2, tx0}, {ty0, ty1, ty2, ty0}, std::map<std::string, std::string>{{"linewidth", "0.1"}, {"color", "white"}});
        
        // message
        std::cout << "adding triange " << i <<  " out of " << d.triangles.size() << std::endl;
    }   
}



using InputPointT = VilensPointT;

// int main()
// {
//     // given index number, add pointcloud to display
//     std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";

//     // std::string pcd_file1 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460854_847538000.pcd";
//     // pcl::PointCloud<InputPointT>::Ptr cloud1 = load_pointcloud<InputPointT>(pcd_file1);
//     // Eigen::Affine3d pose_eigen1 = find_pose(pcd_file1, pose_file);
//     // pcl::PointCloud<InputPointT>::Ptr transformed_cloud1 = transform_to_global<InputPointT>(cloud1, pose_eigen1);

//     std::string pcd_file2 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd";
//     pcl::PointCloud<InputPointT>::Ptr cloud2 = load_pointcloud<InputPointT>(pcd_file2);
//     Eigen::Affine3d pose_eigen2 = find_pose(pcd_file2, pose_file);
//     pcl::PointCloud<InputPointT>::Ptr transformed_cloud2 = transform_to_global<InputPointT>(cloud2, pose_eigen2);


//     // update old points using triangles

//     // // transform cloud1 to cloud2 frame 
//     // Eigen::Affine3d pose1_to_pose2 = pose_eigen2.inverse() * pose_eigen1;
//     // pcl::PointCloud<InputPointT>::Ptr cloud1_in_cloud2_frame (new pcl::PointCloud<InputPointT> ());
//     // pcl::transformPointCloud (*cloud1, *cloud1_in_cloud2_frame, pose1_to_pose2);

//     // triangulate cloud2
//     delaunator::Delaunator d = obtain_triangulation<InputPointT>(cloud2);


//     // ------------------------------------------------------ plt
//     // plot black background
//     plt_plot_black_background();
//     // plot triangles
//     plt_plot_triangles(d);
//     // display
//     plt::show();



//     // // ------------------------------------------------------ pclvisuliazer    
//     // // set up viewer
//     // pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     // viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    
//     // // set up viewports
//     // int port1(0);
//     // viewer->createViewPort (0.0, 0.0, 0.5, 1.0, port1);
//     // viewer->setBackgroundColor (0, 0, 0, port1);
//     // int port2(0);
//     // viewer->createViewPort (0.5, 0.0, 1.0, 1.0, port2);
//     // viewer->setBackgroundColor (0, 0, 0, port2);

//     // // set up coordinate system
//     // viewer->initCameraParameters();
//     // viewer->addCoordinateSystem(1);

//     // // // add triangles to viewer
//     // // add_triangles_to_viewer<InputPointT>(viewer, port1, cloud2, transformed_cloud2, d2);
//     // // display pointclouds
//     // add_to_viewer<InputPointT>(viewer, port1, transformed_cloud1, "transformed cloud", color_tuple(0, 255, 0), 1);
//     // add_to_viewer<InputPointT>(viewer, port1, transformed_cloud2, "transformed cloud2", color_tuple(255, 0, 0), 1);
//     // add_to_viewer<InputPointT>(viewer, port2, cloud1_in_cloud2_frame, "cloud1 in cloud2 frame", color_tuple(0, 255, 0), 1);
//     // add_to_viewer<InputPointT>(viewer, port2, cloud2, "cloud2", color_tuple(255, 0, 0), 1);
    
//     // // display
//     // viewer->spin();
    
//     return (0);
// }


int main()
{
    // given index number, add pointcloud to display
    std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";
    std::string pcd_file2 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd";
    pcl::PointCloud<InputPointT>::Ptr cloud2 = load_pointcloud<InputPointT>(pcd_file2);

    // triangulate cloud2
    delaunator::Delaunator d = obtain_triangulation<InputPointT>(cloud2);

    // update old points using triangles


    // 2d cloud from coords
    pcl::PointCloud<pcl::PointXY>::Ptr cloud2d (new pcl::PointCloud<pcl::PointXY>);
    cloud2d->resize(d.coords.size() / 2);
    for (std::size_t i = 0; i < d.coords.size(); i+=2)
    {
        pcl::PointXY point;
        point.x = d.coords[i];
        point.y = d.coords[i + 1];
        cloud2d->points[i / 2] = point;
    }

    // kd tree search
    pcl::KdTreeFLANN<pcl::PointXY> kdtree;
    kdtree.setInputCloud (cloud2d);
    pcl::PointXY searchPoint;
    searchPoint.x = 0;
    searchPoint.y = 0;

    // K nearest neighbor search
    int K = 1;
    std::vector<int> pointIdxKNNSearch(K);
    std::vector<float> pointKNNSquaredDistance(K);
    if ( kdtree.nearestKSearch (searchPoint, K, pointIdxKNNSearch, pointKNNSquaredDistance) > 0 )
    {
        for (std::size_t i = 0; i < pointIdxKNNSearch.size (); ++i)
        std::cout << "    "  <<   cloud2d->points[pointIdxKNNSearch[i]].x 
                    << " " << cloud2d->points[pointIdxKNNSearch[i]].y 
                    << " (squared distance: " << pointKNNSquaredDistance[i] << ")" << std::endl;
    }

    // list triangles that contains the search point
    std::vector<int> triangles_containing_search_point;
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        if (i1 == pointIdxKNNSearch[0] || i2 == pointIdxKNNSearch[0] || i3 == pointIdxKNNSearch[0])
        {
            triangles_containing_search_point.push_back(i);
        }
    }

    // plot those triangles
    for (std::size_t i = 0; i < triangles_containing_search_point.size(); i++)
    {
        int i1 = d.triangles[triangles_containing_search_point[i]];
        int i2 = d.triangles[triangles_containing_search_point[i] + 1];
        int i3 = d.triangles[triangles_containing_search_point[i] + 2];

        // plt - plot the triangle
        float tx0 = d.coords[2 * i1];
        float ty0 = d.coords[2 * i1 + 1];
        float tx1 = d.coords[2 * i2];
        float ty1 = d.coords[2 * i2 + 1];
        float tx2 = d.coords[2 * i3];
        float ty2 = d.coords[2 * i3 + 1];
        plt::plot({tx0, tx1, tx2, tx0}, {ty0, ty1, ty2, ty0}, std::map<std::string, std::string>{{"linewidth", "1"}, {"color", "black"}});
    }

    // plt - plot points
    std::vector<double> x_list, y_list;
    for (std::size_t i = 0; i < cloud2d->size(); i++)
    {
        x_list.push_back((*cloud2d)[i].x);
        y_list.push_back((*cloud2d)[i].y);
    }
    plt::scatter(x_list, y_list, 10);

    // plt - plot search point
    std::vector<double> search_x_list = {searchPoint.x};
    std::vector<double> search_y_list = {searchPoint.y};
    plt::scatter(search_x_list, search_y_list, 10, std::map<std::string, std::string>{{"color", "red"}});

    // plt - plot K nearest neighbors
    std::vector<double> knn_x_list, knn_y_list;
    for (std::size_t i = 0; i < pointIdxKNNSearch.size(); i++)
    {
        knn_x_list.push_back((*cloud2d)[pointIdxKNNSearch[i]].x);
        knn_y_list.push_back((*cloud2d)[pointIdxKNNSearch[i]].y);
    }
    plt::scatter(knn_x_list, knn_y_list, 10, std::map<std::string, std::string>{{"color", "green"}});

    plt::show();

    return (0);
}

// int main()
// {
//     pcl::PointCloud<pcl::PointXY>::Ptr cloud (new pcl::PointCloud<pcl::PointXY>);

//     // Generate pointcloud data
//     cloud->width = 1000;
//     cloud->height = 1;
//     cloud->points.resize (cloud->width * cloud->height);

//     for (std::size_t i = 0; i < cloud->size (); ++i)
//     {
//         (*cloud)[i].x = 1024.0f * rand () / (RAND_MAX + 1.0f);
//         (*cloud)[i].y = 1024.0f * rand () / (RAND_MAX + 1.0f);
//     }

    


//     pcl::KdTreeFLANN<pcl::PointXY> kdtree;
//     kdtree.setInputCloud (cloud);
//     pcl::PointXY searchPoint;
//     searchPoint.x = 1024.0f * rand () / (RAND_MAX + 1.0f);
//     searchPoint.y = 1024.0f * rand () / (RAND_MAX + 1.0f);

//     // K nearest neighbor search

//     int K = 10;

//     std::vector<int> pointIdxKNNSearch(K);
//     std::vector<float> pointKNNSquaredDistance(K);

//     std::cout << "K nearest neighbor search at (" << searchPoint.x 
//                 << " " << searchPoint.y 
//                 << ") with K=" << K << std::endl;

//     if ( kdtree.nearestKSearch (searchPoint, K, pointIdxKNNSearch, pointKNNSquaredDistance) > 0 )
//     {
//         for (std::size_t i = 0; i < pointIdxKNNSearch.size (); ++i)
//         std::cout << "    "  <<   (*cloud)[ pointIdxKNNSearch[i] ].x 
//                     << " " << (*cloud)[ pointIdxKNNSearch[i] ].y 
//                     << " (squared distance: " << pointKNNSquaredDistance[i] << ")" << std::endl;
//     }

//     // plt - plot points
//     std::vector<double> x_list, y_list;
//     for (std::size_t i = 0; i < cloud->size(); i++)
//     {
//         x_list.push_back((*cloud)[i].x);
//         y_list.push_back((*cloud)[i].y);
//     }
//     plt::scatter(x_list, y_list, 10);

//     // plt - plot search point
//     std::vector<double> search_x_list = {searchPoint.x};
//     std::vector<double> search_y_list = {searchPoint.y};
//     plt::scatter(search_x_list, search_y_list, 10, std::map<std::string, std::string>{{"color", "red"}});

//     // plt - plot K nearest neighbors
//     std::vector<double> knn_x_list, knn_y_list;
//     for (std::size_t i = 0; i < pointIdxKNNSearch.size(); i++)
//     {
//         knn_x_list.push_back((*cloud)[pointIdxKNNSearch[i]].x);
//         knn_y_list.push_back((*cloud)[pointIdxKNNSearch[i]].y);
//     }
//     plt::scatter(knn_x_list, knn_y_list, 10, std::map<std::string, std::string>{{"color", "green"}});

//     plt::show();

// }
