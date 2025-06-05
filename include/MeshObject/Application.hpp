#pragma once

#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>


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
    void add_point_to_map(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>> bvh_results, std::vector<std::shared_ptr<Vertex>> rrs_results);
    
    // interaction
    void refine_surfaces();
    void change_color();
    void get_lidar_data(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_sim_data(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void step();
    void process_pointcloud();
    void process_the_rest();
    void restart();
    std::shared_ptr<Storage> get_storage();

    // getter
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr compute_vertex_point_pointcloud(const Settings& settings);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr compute_interior_point_pointcloud(const Settings& settings);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr compute_generic_point_pointcloud();
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_faces();
    const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& get_edges();
    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices();
    
    // data
    void load_pointcloud_from_dataloader();
    void load_pointcloud(typename pcl::PointCloud<PointT>::Ptr pointcloud_local, Eigen::Affine3d& pose, bool already_in_global_frame = false);
    void get_output_pointcloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& pointcloud_out);

    void write_mesh();

private:
    // objects
    DataLoader<PointT> data_loader;
    std::shared_ptr<Storage> storage_;
        
    // state
    bool first_cloud_;
    int ith_cloud;
    std::size_t ith_point;
    std::size_t ith_size;
    Eigen::Vector3d origin;
    double distance_travelled_;
    typename pcl::PointCloud<PointT>::Ptr pointcloud;

    // settings
    Settings settings_;

    // viewer related
    std::unordered_map<std::shared_ptr<Vertex>, int, MeshObjectHash> vertex_to_cloud_indices_map;

    std::atomic<unsigned int> num_of_concurrent_processes = 0;
    std::atomic<unsigned int> accumulated_points = 0;
    std::mutex process_point_mutex;

    std::chrono::time_point<std::chrono::high_resolution_clock> t_init;
    std::chrono::time_point<std::chrono::high_resolution_clock> t_last;
    std::mutex t_last_mutex;

    double total_duration;
    std::vector<double> duration_list;
    std::vector<double> rrs_search_duration_list;
    std::vector<double> rrs_update_duration_list;
    std::vector<double> bvh_search_duration_list;
    std::vector<double> bvh_update_duration_list;
    std::vector<double> add_to_map_duration_list;
    std::vector<double> delete_from_map_duration_list;
    std::vector<double> relative_position_duration_list;

    std::vector<double> rrs_search_duration_per_thread;
    std::vector<double> rrs_update_duration_per_thread;
    std::vector<double> bvh_search_duration_per_thread;
    std::vector<double> bvh_update_duration_per_thread;
    std::vector<double> add_to_map_duration_per_thread;
    std::vector<double> delete_from_map_duration_per_thread;
    std::vector<double> relative_position_duration_per_thread;

    unsigned int total_loops;

    std::unordered_map<std::shared_ptr<Surface>, unsigned int, MeshObjectHash> surface_to_contention_count;

};