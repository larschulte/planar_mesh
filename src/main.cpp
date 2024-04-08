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
void viewer_add_triangles(pcl::visualization::PCLVisualizer::Ptr viewer, int viewport, typename pcl::PointCloud<PointT>::Ptr cloud, typename pcl::PointCloud<PointT>::Ptr transformed_cloud, delaunator::Delaunator d)
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

// point to triangle map for d2
std::map<int, std::vector<int>> point_to_triangle_map(delaunator::Delaunator d)
{
    std::map<int, std::vector<int>> pt_map;
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        // indices
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        // store
        pt_map[i1].push_back(i);
        pt_map[i2].push_back(i);
        pt_map[i3].push_back(i);
    }

    return pt_map;
}

// convert to 2d polar cloud
template<typename PointT>
pcl::PointCloud<pcl::PointXY>::Ptr obtain_2d_polar_cloud(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    pcl::PointCloud<pcl::PointXY>::Ptr cloud_polar (new pcl::PointCloud<pcl::PointXY>);
    cloud_polar->resize(cloud->size());
    for (std::size_t i = 0; i < cloud->size(); i+=1)
    {
        pcl::PointXY point;

        // compute azimuth and altitude
        float x = cloud->points[i].x;
        float y = cloud->points[i].y;
        float z = cloud->points[i].z;
        double r = sqrt(x * x + y * y + z * z);
        double azimuth = atan2(y, x) * 180 / M_PI;
        double altitude = asin(z / r) * 180 / M_PI;

        point.x = azimuth;
        point.y = altitude;
        cloud_polar->points[i] = point;
    }

    return cloud_polar;
}


// convert to target frame using pose
template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr transform_to_frame(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d current_pose, Eigen::Affine3d target_pose)
{
    Eigen::Affine3d pose_current2target = target_pose.inverse() * current_pose; // the cloud will first be transformed to current pose, then to target pose
    return transform_to_global<PointT>(cloud, pose_current2target);
}

using InputPointT = VilensPointT;

int main()
{
    // given index number, add pointcloud to display
    std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";
    std::string pcd_file1 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd";
    std::string pcd_file2 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460870_532271000.pcd";
    pcl::PointCloud<InputPointT>::Ptr cloud1 = load_pointcloud<InputPointT>(pcd_file1);
    pcl::PointCloud<InputPointT>::Ptr cloud2 = load_pointcloud<InputPointT>(pcd_file2);

    Eigen::Affine3d pose1 = find_pose(pcd_file1, pose_file);
    Eigen::Affine3d pose2 = find_pose(pcd_file2, pose_file);

    // 3std 
    double range_std = 0.01;
    double range_std_x3 = range_std * 3;


    // compute near cloud and far for cloud1
    pcl::PointCloud<InputPointT>::Ptr cloud1_near (new pcl::PointCloud<InputPointT>);
    pcl::PointCloud<InputPointT>::Ptr cloud1_far (new pcl::PointCloud<InputPointT>);
    // resize
    cloud1_near->resize(cloud1->size());
    cloud1_far->resize(cloud1->size());
    // vector to store eigen direction
    std::vector<Eigen::Vector3f> directions1;
    // compute azimuth and altitude coordinates
    for (std::size_t i = 0; i < cloud1->size(); i++)
    {
        float x = cloud1->points[i].x;
        float y = cloud1->points[i].y;
        float z = cloud1->points[i].z;

        Eigen::Vector3f point(x, y, z);
        Eigen::Vector3f direction = point.normalized();
        Eigen::Vector3f point_near = point - direction * range_std_x3;
        Eigen::Vector3f point_far = point + direction * range_std_x3;

        // store
        cloud1_near->points[i].x = point_near(0);
        cloud1_near->points[i].y = point_near(1);
        cloud1_near->points[i].z = point_near(2);

        cloud1_far->points[i].x = point_far(0);
        cloud1_far->points[i].y = point_far(1);
        cloud1_far->points[i].z = point_far(2);

        directions1.push_back(direction);
    }


    // compute near cloud and far cloud for cloud 2
    pcl::PointCloud<InputPointT>::Ptr cloud2_near (new pcl::PointCloud<InputPointT>);
    pcl::PointCloud<InputPointT>::Ptr cloud2_far (new pcl::PointCloud<InputPointT>);
    // resize
    cloud2_near->resize(cloud2->size());
    cloud2_far->resize(cloud2->size());
    // vector to store eigen direction
    std::vector<Eigen::Vector3f> directions2;
    // compute azimuth and altitude coordinates
    for (std::size_t i = 0; i < cloud2->size(); i++)
    {
        float x = cloud2->points[i].x;
        float y = cloud2->points[i].y;
        float z = cloud2->points[i].z;

        Eigen::Vector3f point(x, y, z);
        Eigen::Vector3f direction = point.normalized();
        Eigen::Vector3f point_near = point - direction * range_std_x3;
        Eigen::Vector3f point_far = point + direction * range_std_x3;

        // store
        cloud2_near->points[i].x = point_near(0);
        cloud2_near->points[i].y = point_near(1);
        cloud2_near->points[i].z = point_near(2);

        cloud2_far->points[i].x = point_far(0);
        cloud2_far->points[i].y = point_far(1);
        cloud2_far->points[i].z = point_far(2);

        directions2.push_back(direction);
    }


    // update cloud1 using cloud2

    // triangulate cloud 2
    delaunator::Delaunator d2 = obtain_triangulation<InputPointT>(cloud2);

    // point to triangle map for d2
    std::map<int, std::vector<int>> pt_map2 = point_to_triangle_map(d2);

    
    // prepare kd tree search for cloud 2
    pcl::PointCloud<pcl::PointXY>::Ptr cloud2_polar = obtain_2d_polar_cloud<InputPointT>(cloud2);
    pcl::KdTreeFLANN<pcl::PointXY> kdtree;
    kdtree.setInputCloud(cloud2_polar);

    // convert cloud1_near to cloud2 coordinate
    pcl::PointCloud<InputPointT>::Ptr cloud1_near_in_cloud2_frame = transform_to_frame<InputPointT>(cloud1_near, pose1, pose2);
    pcl::PointCloud<InputPointT>::Ptr cloud1_far_cloud2_frame = transform_to_frame<InputPointT>(cloud1_far, pose1, pose2);


    // cloud to store updated point
    pcl::PointCloud<InputPointT>::Ptr cloud1_near_in_cloud2_frame_updated (new pcl::PointCloud<InputPointT>);
    
    // iterate through polar of cloud1_near_in_cloud2_frame
    pcl::PointCloud<pcl::PointXY>::Ptr cloud1_in_cloud2_frame_polar = obtain_2d_polar_cloud<InputPointT>(cloud1_near_in_cloud2_frame);
    for (std::size_t i = 0; i < cloud1_in_cloud2_frame_polar->size(); i++)
    {
        InputPointT current_point = cloud1_near_in_cloud2_frame->points[i];
        pcl::PointXY current_point_polar = cloud1_in_cloud2_frame_polar->points[i];
        Eigen::Vector3f current_point_direction = directions1[i];

        // search the closest point in cloud2 to current point in cloud 1 in polar coordinate
        int K = 1;
        pcl::Indices pointIdxKNNSearch(K);
        std::vector<float> pointKNNSquaredDistance(K);
        kdtree.nearestKSearch(current_point_polar, K, pointIdxKNNSearch, pointKNNSquaredDistance);

        // from the closest point, find all possible triangles that contains the current point
        std::vector<int> triangles_containing_search_point = pt_map2[pointIdxKNNSearch[0]];

        // find triangle that intersects with the search point in given direction
        for (int triangle_index : triangles_containing_search_point)
        {
            // indices
            int i1 = d2.triangles[triangle_index];
            int i2 = d2.triangles[triangle_index + 1];
            int i3 = d2.triangles[triangle_index + 2];

            // get local points
            InputPointT p0, p1, p2;
            p0 = cloud2_near->points[i1];
            p1 = cloud2_near->points[i2];
            p2 = cloud2_near->points[i3];

            // compute point to plane intersection
            Eigen::Vector3f normal = (p1.getVector3fMap() - p0.getVector3fMap()).cross(p2.getVector3fMap() - p0.getVector3fMap());
            normal.normalize();
            // make normal point towards 0, 0, 0 of cloud2
            if (normal.dot(p0.getVector3fMap()) < 0)
            {
                normal = -normal;
            }
            double distance = (p0.getVector3fMap() - current_point).dot(normal) / current_point_direction.dot(normal);
            Eigen::Vector3f intersection = current_point + current_point_direction * distance;

            // compute barycentric coordinates of intersection
            Eigen::Vector3f v0 = p1.getVector3fMap() - p0.getVector3fMap();
            Eigen::Vector3f v1 = p2.getVector3fMap() - p0.getVector3fMap();
            Eigen::Vector3f v2 = intersection - p0.getVector3fMap();
            float dot00 = v0.dot(v0);
            float dot01 = v0.dot(v1);
            float dot02 = v0.dot(v2);
            float dot11 = v1.dot(v1);
            float dot12 = v1.dot(v2);
            float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
            float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
            // check if point is inside the triangle
            bool inside = (u >= 0) && (v >= 0) && (u + v < 1);

            if (inside)
            {
                // check if intersection is farther than the point
                // update cloud1_near_in_cloud2_frame point if new triangle is farther than the point
                if (distance > 0 && distance < 2 * range_std_x3)
                {
                    InputPointT updated_point;
                    updated_point.x = intersection(0);
                    updated_point.y = intersection(1);
                    updated_point.z = intersection(2);
                    cloud1_near_in_cloud2_frame_updated->push_back(updated_point);
                }

                break;
            }
        }
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

    // display pointclouds
    // add_to_viewer<InputPointT>(viewer, port1, cloud1_near_in_cloud2_frame, "cloud1 near", color_tuple(255, 0, 0), 2); // r
    add_to_viewer<InputPointT>(viewer, port1, cloud2_near, "cloud2 near", color_tuple(0, 255, 0), 2); // g

    // add_to_viewer<InputPointT>(viewer, port2, cloud1_near_in_cloud2_frame, "cloud1 near2", color_tuple(255, 0, 0), 2); // r
    add_to_viewer<InputPointT>(viewer, port2, cloud1_near_in_cloud2_frame_updated, "cloud1 near updated", color_tuple(255, 255, 0), 2); // y
    
    // display
    viewer->spin();

    return (0);
}