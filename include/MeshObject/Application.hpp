#pragma once

#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "point_type/VilensPointT.hpp"


#include "utilities/DataLoader.hpp"
#include <unordered_set>
#include "MeshObject/MeshObject.hpp"
#include "MeshObject/Settings.hpp"

#include <mutex>

class Vertex;
class Edge;
class Face;
class Surface;
class GenericPoint;
class InteriorPoint;
class Storage;

template <typename PointT>
class Application
{
public:
    Application();

    // helper function
    Eigen::Matrix3d merge_covariances_of_surfaces(std::shared_ptr<Surface> surface1, std::shared_ptr<Surface> surface2);
    double compute_eigenvalue_of_merged_surfaces(std::shared_ptr<Surface> surface1, std::shared_ptr<Surface> surface2);
    // void try_merge_surfaces(std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& surfaces_to_merge);
    
    // algorithm
    void process_point(const std::shared_ptr<GenericPoint>& generic_point);
    bool add_point_by_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, double& radius, std::vector<std::shared_ptr<Face>>& searched_faces);
    bool add_point_by_radius_search(const std::shared_ptr<GenericPoint>& generic_point, double& radius, std::vector<std::shared_ptr<Vertex>>& neighboring_vertices_vector);
    void add_point_by_new_surface(const std::shared_ptr<GenericPoint>& generic_point, double& radius);
    
    // interaction
    void refine_surfaces();
    void change_color();
    void get_lidar_data(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_sim_data(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void step();
    void loop();
    void restart();
    void rebuild_tree();

    // getter
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr compute_vertex_point_pointcloud(const Settings& settings);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr compute_interior_point_pointcloud(const Settings& settings);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr compute_generic_point_pointcloud();
    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map();
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_faces();
    const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& get_edges();
    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices();
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> get_boundary_edges();
    
    // data
    void load_point_cloud();
    int ith_cloud;
    std::size_t ith_point = 0;
    std::size_t ith_size = 0;

private:
    // objects
    DataLoader<VilensPointT> data_loader;
    std::shared_ptr<Storage> storage_;
    
    // dataset
    std::string dataset;
    typename pcl::PointCloud<VilensPointT>::Ptr pointcloud;
    Eigen::Vector3d origin;

    // settings
    Settings settings_;

    // viewer related
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map;

    std::atomic<unsigned int> num_of_concurrent_processes = 0;
    std::atomic<unsigned int> accumulated_points = 0;
    std::mutex process_point_mutex;

    std::chrono::time_point<std::chrono::high_resolution_clock> t_init;
    std::chrono::time_point<std::chrono::high_resolution_clock> t_last;
    std::mutex t_last_mutex;
};