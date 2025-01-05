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
    bool add_point_by_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& searched_faces, std::vector<std::shared_ptr<Vertex>>& rrs_results, std::shared_ptr<Surface>& added_surface, std::shared_ptr<Face>& face_added_by_intersection_search);
    bool add_point_by_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& neighboring_vertices_vector, std::vector<std::shared_ptr<Face>>& bvh_results, std::shared_ptr<Surface>& added_surface, std::shared_ptr<Vertex>& vertex_added_by_radius_search);
    
    // interaction
    void refine_surfaces();
    void change_color();
    void get_lidar_data(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_sim_data(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void step();
    void loop();
    void process_the_rest();
    void restart();
    void rebuild_tree();
    void cleanup_surfaces();
    void remove_non_manifold_edges();

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

private:
    // objects
    DataLoader<PointT> data_loader;
    std::shared_ptr<Storage> storage_;
        
    // state
    int ith_cloud;
    std::size_t ith_point;
    std::size_t ith_size;
    Eigen::Vector3d origin;
    double distance_travelled_;
    typename pcl::PointCloud<PointT>::Ptr pointcloud;

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
    double total_duration;
    unsigned int total_loops;

    std::unordered_map<std::shared_ptr<Surface>, unsigned int, MeshObjectHash> surface_to_contention_count;
};