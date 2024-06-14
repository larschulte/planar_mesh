#pragma once

#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "point_type/VilensPointT.hpp"


#include "utilities/DataLoader.hpp"

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
    void try_merge_surfaces(std::set<std::shared_ptr<Surface>>& surfaces_to_merge);
    
    // algorithm
    void process_point(Eigen::Vector3d thisPointOriginVEC, Eigen::Vector3d thisPointVEC);
    void add_point_by_radius_search(const Eigen::Vector3d& thisPointVEC, const Eigen::Vector3d& thisPointOriginVEC);
    
    // interaction
    void refine_surfaces();
    void change_color();
    void add_back_generic_points();
    void step();
    void loop();

    // getter
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_colored_cloud();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_distance_cloud();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_colored_cloud();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_distance_cloud();
    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map();
    const std::set<std::shared_ptr<Face>>& get_faces();
    const std::set<std::shared_ptr<Edge>>& get_edges();
    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices();
    std::set<std::shared_ptr<Edge>> get_boundary_edges();
    
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
    double distance_threshold;
    std::size_t fit_plane_threshold;
    double merged_eigenvalue_threshold;
    bool shuffle_pointcloud;
    double pointcloud_fraction;
    double distance_to_radius_ratio;

    // viewer related
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map;
};