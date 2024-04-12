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
#include <pcl/filters/extract_indices.h>

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

// check if within triangle
bool is_inside_triangle(Eigen::Vector3f p0, Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f point)
{
    Eigen::Vector3f v0 = p1 - p0;
    Eigen::Vector3f v1 = p2 - p0;
    Eigen::Vector3f v2 = point - p0;
    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);
    float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return (u >= 0) && (v >= 0) && (u + v < 1);
}

// compute ray to triangle intersection
Eigen::Vector3f ray_triangle_intersection(Eigen::Vector3f ray_origin, Eigen::Vector3f ray_direction, Eigen::Vector3f p0, Eigen::Vector3f p1, Eigen::Vector3f p2)
{
    // compute triangle normal
    Eigen::Vector3f normal = (p1 - p0).cross(p2 - p0);
    
    // if parallel, return NaN
    if (normal.dot(ray_direction) == 0)
    {
        return Eigen::Vector3f(NAN, NAN, NAN);
    }

    // compute intersection
    double distance = (p0 - ray_origin).dot(normal) / ray_direction.dot(normal);
    Eigen::Vector3f intersection = ray_origin + ray_direction * distance;

    // return
    return intersection;
}

// generate near and far pointcloud and direction vector
template<typename PointT>
void compute_near_far_cloud(typename pcl::PointCloud<PointT>::Ptr cloud, double range_std, typename pcl::PointCloud<PointT>::Ptr cloud_near, typename pcl::PointCloud<PointT>::Ptr cloud_far, std::vector<Eigen::Vector3f>& directions)
{
    // 3std
    double range_std_x3 = range_std * 3;
    // resize
    cloud_near->resize(cloud->size());
    cloud_far->resize(cloud->size());
    // vector to store eigen direction
    directions.clear();
    // compute azimuth and altitude coordinates
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        float x = cloud->points[i].x;
        float y = cloud->points[i].y;
        float z = cloud->points[i].z;

        Eigen::Vector3f point(x, y, z);
        Eigen::Vector3f direction = point.normalized();
        Eigen::Vector3f point_near = point - direction * range_std_x3;
        Eigen::Vector3f point_far = point + direction * range_std_x3;

        // store
        cloud_near->points[i].x = point_near(0);
        cloud_near->points[i].y = point_near(1);
        cloud_near->points[i].z = point_near(2);

        cloud_far->points[i].x = point_far(0);
        cloud_far->points[i].y = point_far(1);
        cloud_far->points[i].z = point_far(2);

        directions.push_back(direction);
    }
}

// compute direction vector
template<typename PointT>
std::vector<Eigen::Vector3f> compute_point_directions(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    // initialize
    std::vector<Eigen::Vector3f> directions;
    directions.resize(cloud->size());

    // compute
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        Eigen::Vector3f direction = cloud->points[i].getVector3fMap().normalized();
        directions[i] = direction;
    }

    // return
    return directions;
}


// get k nearest neighbor in tree
template<typename PointT>
pcl::Indices get_k_nearest_neighbor(pcl::KdTreeFLANN<PointT> kdtree, PointT point, int k)
{
    pcl::Indices nearest_neighbor(k);
    std::vector<float> nearest_neighbor_distance(k);
    kdtree.nearestKSearch(point, k, nearest_neighbor, nearest_neighbor_distance);
    return nearest_neighbor;
}

// compute triangle center
template<typename PointT>
typename pcl::PointCloud<pcl::PointXYZ>::Ptr computer_triangle_center(typename pcl::PointCloud<PointT>::Ptr vertex_cloud, delaunator::Delaunator d)
{
    // initialize
    typename pcl::PointCloud<pcl::PointXYZ>::Ptr center_cloud (new typename pcl::PointCloud<pcl::PointXYZ>);
    center_cloud->resize(d.triangles.size() / 3);

    // compute triangle centers
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        // vertcies index
        int v1_index = d.triangles[i];
        int v2_index = d.triangles[i + 1];
        int v3_index = d.triangles[i + 2];

        // vertices
        Eigen::Vector3f v1 = vertex_cloud->points[v1_index].getVector3fMap();
        Eigen::Vector3f v2 = vertex_cloud->points[v2_index].getVector3fMap();
        Eigen::Vector3f v3 = vertex_cloud->points[v3_index].getVector3fMap();

        // center
        Eigen::Vector3f center = (v1 + v2 + v3) / 3;

        // store center
        center_cloud->points[i / 3].x = center(0);
        center_cloud->points[i / 3].y = center(1);
        center_cloud->points[i / 3].z = center(2);
    }

    // return
    return center_cloud;
}

// compute triangle center to vertices index map
// map[triangle_index in d] = [v1_index, v2_index, v3_index in vertex_cloud]
std::map<int, std::vector<int>> compute_center_to_vertices_index_map(delaunator::Delaunator d)
{
    // initialize
    std::map<int, std::vector<int>> map;
    
    // compute
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        // vertcies index
        int v1_index = d.triangles[i];
        int v2_index = d.triangles[i + 1];
        int v3_index = d.triangles[i + 2];

        // store center to vertex indices map
        map[i / 3].push_back(v1_index);
        map[i / 3].push_back(v2_index);
        map[i / 3].push_back(v3_index);
    }

    // return
    return map;
}

// compute updated pointcloud
// define macros type NEAR and FAR
#define NEAR 0
#define FAR 1
// function
template <typename PointT>
void update_pointcloud(
    typename pcl::PointCloud<PointT>::Ptr old_cloud, 
    const std::vector<Eigen::Vector3f> old_cloud_direction, 
    std::vector<float>& old_cloud_variance, 
    const typename pcl::PointCloud<PointT>::Ptr new_cloud, 
    std::vector<int>& used_points,
    double range_std
    )
{
    // triangulate new cloud
    delaunator::Delaunator d = obtain_triangulation<PointT>(new_cloud);
    
    // compute triangle centers 
    pcl::PointCloud<pcl::PointXYZ>::Ptr center_cloud = computer_triangle_center<PointT>(new_cloud, d);

    // compute triangle center to vertices index map
    std::map<int, std::vector<int>> center_to_vertices_index_map = compute_center_to_vertices_index_map(d);
    
    // prepare kd tree search
    pcl::KdTreeFLANN<pcl::PointXY> new_cloud_triangle_center_polar_kdtree;
    new_cloud_triangle_center_polar_kdtree.setInputCloud(obtain_2d_polar_cloud<pcl::PointXYZ>(center_cloud));
    
    // iterate through polar of old cloud
    std::vector<int> used_triangles;
    pcl::PointCloud<pcl::PointXY>::Ptr old_cloud_polar = obtain_2d_polar_cloud<PointT>(old_cloud);
    for (std::size_t i = 0; i < old_cloud_polar->size(); i++)
    {
        // current point
        Eigen::Vector3f current_point = old_cloud->points[i].getVector3fMap();
        pcl::PointXY current_point_polar = old_cloud_polar->points[i];
        Eigen::Vector3f current_point_direction = old_cloud_direction[i];
        
        // 1. search closest triangle center
        int K = 4;
        pcl::Indices center_indices_searched = get_k_nearest_neighbor(new_cloud_triangle_center_polar_kdtree, current_point_polar, K);
        
        // 2. find intersections to those triangles
        std::vector<int> all_intersected_triangle_indices;
        std::vector<Eigen::Vector3f> all_intersections;
        for (int center_index : center_indices_searched)
        {
            // get triangle vertices xyz
            Eigen::Vector3f v1 = new_cloud->points[center_to_vertices_index_map[center_index][0]].getVector3fMap();
            Eigen::Vector3f v2 = new_cloud->points[center_to_vertices_index_map[center_index][1]].getVector3fMap();
            Eigen::Vector3f v3 = new_cloud->points[center_to_vertices_index_map[center_index][2]].getVector3fMap();

            // compute intersection
            Eigen::Vector3f intersection = ray_triangle_intersection(current_point, current_point_direction, v1, v2, v3);

            // check if intersection is inside the triangle
            bool inside = is_inside_triangle(v1, v2, v3, intersection);
            if (inside)
            {
                all_intersected_triangle_indices.push_back(center_index);
                all_intersections.push_back(intersection);
            }
        }
        bool no_intersections = all_intersections.size() == 0;
        if (no_intersections) continue;

        // 3. find cloest intersection
        auto iterator = std::min_element(all_intersections.begin(), all_intersections.end(),
            [current_point](const Eigen::Vector3f& a, const Eigen::Vector3f& b) {return (a - current_point).norm() < (b - current_point).norm();});
        int closest_intersected_triangle_index = all_intersected_triangle_indices[std::distance(all_intersections.begin(), iterator)];
        Eigen::Vector3f closest_intersection = *iterator;
        
        // 4. compute posterior
        Eigen::Vector3f prior_point = current_point;
        float prior_variance = old_cloud_variance[i];

        Eigen::Vector3f likelihood_point = closest_intersection;
        float likelihood_variance = std::pow(range_std, 2);

        Eigen::Vector3f posterior_point = (prior_point / prior_variance + likelihood_point / likelihood_variance) / (1 / prior_variance + 1 / likelihood_variance);
        float posterior_variance = 1 / (1 / prior_variance + 1 / likelihood_variance);

        float update_threshold_distance = 3 * std::pow(prior_variance, 0.5) + 3 * std::pow(likelihood_variance, 0.5);
        bool too_far = (prior_point - likelihood_point).norm() > update_threshold_distance;
        if (too_far) continue;


        // 5. update old cloud
        old_cloud->points[i].x = posterior_point(0);
        old_cloud->points[i].y = posterior_point(1);
        old_cloud->points[i].z = posterior_point(2);
        old_cloud_variance[i] = posterior_variance;
        used_triangles.push_back(closest_intersected_triangle_index); // record used triangles
    }

    // record used_points
    for (int used_triangle : used_triangles)
    {
        for (int vertex_index : center_to_vertices_index_map[used_triangle])
        {
            used_points.push_back(vertex_index);
        }
    }
}

// unique a vector
template<typename T>
void unique_vector(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


// remove used points
template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr remove_used_points(typename pcl::PointCloud<PointT>::Ptr cloud, std::vector<int> used_points)
{
    // convert into indices
    pcl::PointIndices::Ptr used_points_indices (new pcl::PointIndices);
    used_points_indices->indices = used_points;

    // remove used points
    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(used_points_indices);
    extract.setNegative(true);
    typename pcl::PointCloud<PointT>::Ptr cloud_copied (new typename pcl::PointCloud<PointT>);
    extract.filter(*cloud_copied);
    
    return cloud_copied;
}

// remove used direction
std::vector<Eigen::Vector3f> remove_used_directions(std::vector<Eigen::Vector3f> directions, std::vector<int> used_points)
{
    std::vector<Eigen::Vector3f> directions_copied;
    for (std::size_t i = 0; i < directions.size(); i++)
    {
        // if i not in used points
        if (std::find(used_points.begin(), used_points.end(), i) == used_points.end())
        {
            directions_copied.push_back(directions[i]);
        }
    }
    return directions_copied;
}

// remove used variance
std::vector<float> remove_used_variance(std::vector<float> variances, std::vector<int> used_points)
{
    std::vector<float> variances_copied;
    for (std::size_t i = 0; i < variances.size(); i++)
    {
        // if i not in used points
        if (std::find(used_points.begin(), used_points.end(), i) == used_points.end())
        {
            variances_copied.push_back(variances[i]);
        }
    }
    return variances_copied;
}

using InputPointT = VilensPointT;


int main()
{
    // given index number, add pointcloud to display
    // std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";
    std::string pose_file = "/home/jiahao/datasets/bag2pcd_output/slam_poses/slam_poss_graph.slam";

    // pcd files
    // need to test on broadstreet dataset
    std::vector<std::string> pcd_file_list;
    // pcd_file_list.push_back("/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd");
    // pcd_file_list.push_back("/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460870_532271000.pcd");
    // pcd_file_list.push_back("/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460871_630561000.pcd");
    // pcd_file_list.push_back("/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460873_030014000.pcd");
    // pcd_file_list.push_back("/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460876_026726000.pcd");
    pcd_file_list.push_back("/home/jiahao/datasets/bag2pcd_output/slam_clouds/cloud_1712763500_887906000.pcd");
    pcd_file_list.push_back("/home/jiahao/datasets/bag2pcd_output/slam_clouds/cloud_1712763501_388037000.pcd");
    pcd_file_list.push_back("/home/jiahao/datasets/bag2pcd_output/slam_clouds/cloud_1712763501_888041000.pcd");

    // algorithm parameters
    double range_std = 0.01;

    // control cloud (for comparing with updated old cloud)
    pcl::PointCloud<InputPointT>::Ptr control_cloud (new pcl::PointCloud<InputPointT>);

    // old cloud
    pcl::PointCloud<InputPointT>::Ptr old_cloud (new pcl::PointCloud<InputPointT>);
    std::vector<Eigen::Vector3f> old_cloud_direction;
    std::vector<float> old_cloud_variance;

    for (std::string pcd_file : pcd_file_list)
    {
        // load cloud and pose
        pcl::PointCloud<InputPointT>::Ptr new_cloud = load_pointcloud<InputPointT>(pcd_file);
        Eigen::Affine3d new_pose = find_pose(pcd_file, pose_file);

        // add its global to control cloud
        *control_cloud += *transform_to_global<InputPointT>(new_cloud, new_pose);

        // compute near cloud and far for new cloud
        std::vector<Eigen::Vector3f> new_cloud_direction = compute_point_directions<InputPointT>(new_cloud);
        std::vector<float> new_cloud_variance(new_cloud->size(), std::pow(range_std, 2));

        // transform old cloud to local frame
        pcl::PointCloud<InputPointT>::Ptr old_cloud_local = transform_to_frame<InputPointT>(old_cloud, Eigen::Isometry3d::Identity(), new_pose);

        // collect used point, update old cloud
        std::vector<int> used_points;
        update_pointcloud<InputPointT>(old_cloud_local, old_cloud_direction, old_cloud_variance, new_cloud, used_points,range_std);
        unique_vector<int>(used_points);
        
        // add ununsed points / directions to old cloud
        pcl::PointCloud<InputPointT>::Ptr filtered_new_cloud = remove_used_points<InputPointT>(new_cloud, used_points);
        *old_cloud_local += *filtered_new_cloud;

        std::vector<Eigen::Vector3f> filtered_directions = remove_used_directions(new_cloud_direction, used_points);
        old_cloud_direction.insert(old_cloud_direction.end(), filtered_directions.begin(), filtered_directions.end());

        std::vector<float> filtered_varaince = remove_used_variance(new_cloud_variance, used_points);
        old_cloud_variance.insert(old_cloud_variance.end(), filtered_varaince.begin(), filtered_varaince.end());

        // transform old cloud to global frame
        pcl::PointCloud<InputPointT>::Ptr old_cloud_global = transform_to_frame<InputPointT>(old_cloud_local, new_pose, Eigen::Isometry3d::Identity());

        // update old cloud
        *old_cloud = *old_cloud_global;
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
    add_to_viewer<InputPointT>(viewer, port2, old_cloud, "oldcloud mean", color_tuple(0, 255, 0), 3); // b
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "control cloud");
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "oldcloud mean");
    
    
    
    // display
    viewer->spin();

    return (0);
}









// int main()
// {
//     // given index number, add pointcloud to display
//     std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";


//     std::string pcd_file1 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd";
//     std::string pcd_file2 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460870_532271000.pcd";
//     std::string pcd_file3 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460871_630561000.pcd";

//     pcl::PointCloud<InputPointT>::Ptr cloud1 = load_pointcloud<InputPointT>(pcd_file1);
//     pcl::PointCloud<InputPointT>::Ptr cloud2 = load_pointcloud<InputPointT>(pcd_file2);
//     pcl::PointCloud<InputPointT>::Ptr cloud3 = load_pointcloud<InputPointT>(pcd_file3);

//     Eigen::Affine3d pose1 = find_pose(pcd_file1, pose_file);
//     Eigen::Affine3d pose2 = find_pose(pcd_file2, pose_file);
//     Eigen::Affine3d pose3 = find_pose(pcd_file3, pose_file);

//     // 3std 
//     double range_std = 0.01;

//     // compute near cloud and far for cloud1
//     pcl::PointCloud<InputPointT>::Ptr cloud1_near (new pcl::PointCloud<InputPointT>);
//     pcl::PointCloud<InputPointT>::Ptr cloud1_far (new pcl::PointCloud<InputPointT>);
//     std::vector<Eigen::Vector3f> cloud1_direction;
//     generate_near_far_cloud<InputPointT>(cloud1, range_std, cloud1_near, cloud1_far, cloud1_direction);
    
//     // compute near cloud and far cloud for cloud 2
//     pcl::PointCloud<InputPointT>::Ptr cloud2_near (new pcl::PointCloud<InputPointT>);
//     pcl::PointCloud<InputPointT>::Ptr cloud2_far (new pcl::PointCloud<InputPointT>);
//     std::vector<Eigen::Vector3f> cloud2_direction;
//     generate_near_far_cloud<InputPointT>(cloud2, range_std, cloud2_near, cloud2_far, cloud2_direction);

//     // compute near cloud and far cloud for cloud 3
//     pcl::PointCloud<InputPointT>::Ptr cloud3_near (new pcl::PointCloud<InputPointT>);
//     pcl::PointCloud<InputPointT>::Ptr cloud3_far (new pcl::PointCloud<InputPointT>);
//     std::vector<Eigen::Vector3f> cloud3_direction;
//     generate_near_far_cloud<InputPointT>(cloud3, range_std, cloud3_near, cloud3_far, cloud3_direction);
    
//     // update pointcloud
//     pcl::PointCloud<InputPointT>::Ptr cloud1_near_in_cloud2_frame = transform_to_frame<InputPointT>(cloud1_near, pose1, pose2);
//     pcl::PointCloud<InputPointT>::Ptr cloud1_far_in_cloud2_frame = transform_to_frame<InputPointT>(cloud1_far, pose1, pose2);
//     update_pointcloud<InputPointT>(cloud1_near_in_cloud2_frame, cloud1_direction, cloud2_near, cloud2_direction, range_std, NEAR);
//     update_pointcloud<InputPointT>(cloud1_far_in_cloud2_frame, cloud1_direction, cloud2_far, cloud2_direction, range_std, FAR);

//     // update pointcloud
//     pcl::PointCloud<InputPointT>::Ptr old_cloud_near = transform_to_frame<InputPointT>(cloud1_near_in_cloud2_frame, pose2, pose3);
//     pcl::PointCloud<InputPointT>::Ptr old_cloud_far = transform_to_frame<InputPointT>(cloud1_far_in_cloud2_frame, pose2, pose3);
//     update_pointcloud<InputPointT>(old_cloud_near, cloud1_direction, cloud3_near, cloud3_direction, range_std, NEAR);
//     update_pointcloud<InputPointT>(old_cloud_far, cloud1_direction, cloud3_far, cloud3_direction, range_std, FAR);








//     // the current update assume planar surface within each triangle, and does not filter the planar surface even if the triangle is very large
//     // this will be solved when introducing eye patch

//     // ------------------------------------------------------ pclvisuliazer    
//     // set up viewer
//     pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    
//     // // set up viewports
//     // int port1(0);
//     // viewer->createViewPort (0.0, 0.0, 0.5, 1.0, port1);
//     // viewer->setBackgroundColor (0, 0, 0, port1);
//     // int port2(0);
//     // viewer->createViewPort (0.5, 0.0, 1.0, 1.0, port2);
//     // viewer->setBackgroundColor (0, 0, 0, port2);

//     // set up viewports
//     int port1(0);
//     viewer->createViewPort (0.0, 0.0, 1, 1.0, port1);
//     viewer->setBackgroundColor (0, 0, 0, port1);
    
//     // set up coordinate system
//     viewer->initCameraParameters();
//     viewer->addCoordinateSystem(1);

//     // // display pointclouds
//     // // add_to_viewer<InputPointT>(viewer, port1, cloud1_near_in_cloud2_frame, "cloud1 near", color_tuple(255, 0, 0), 2); // r
//     // add_to_viewer<InputPointT>(viewer, port1, cloud1_near_in_cloud2_frame, "cloud1 near original", color_tuple(0, 255, 0), 2); // y
//     // add_to_viewer<InputPointT>(viewer, port1, cloud1_far_in_cloud2_frame, "cloud1 far original", color_tuple(255, 0, 0), 2); // g
//     // add_to_viewer<InputPointT>(viewer, port1, cloud2_near, "cloud2 near original", color_tuple(0, 255, 0), 2); // y
//     // add_to_viewer<InputPointT>(viewer, port1, cloud2_far, "cloud2 far original", color_tuple(255, 0, 0), 2); // g
//     // add_to_viewer<InputPointT>(viewer, port1, cloud3_near, "cloud3 near original", color_tuple(0, 255, 0), 2); // y
//     // add_to_viewer<InputPointT>(viewer, port1, cloud3_far, "cloud3 far original", color_tuple(255, 0, 0), 2); // g


//     add_to_viewer<InputPointT>(viewer, port1, old_cloud_near, "oldcloud near updated", color_tuple(0, 255, 0), 2); // y
//     add_to_viewer<InputPointT>(viewer, port1, old_cloud_far, "oldcloud far updated", color_tuple(255, 0, 0), 2); // g
    
//     // display
//     viewer->spin();

//     return (0);
// }