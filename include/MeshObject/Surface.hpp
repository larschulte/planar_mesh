#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include "MeshObject/EdgeBVH.hpp"

// forward declarations
class Vertex;
class Edge;
class Face;
class InteriorPoint;
class Storage;

class Surface : public std::enable_shared_from_this<Surface> 
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage);
    void delete_();

public:
    double compute_point_to_surface_distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const;
    double compute_point_to_surface_distance_with_improved_covariance(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const;
    Eigen::Vector3d compute_point_to_surface_position(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const;
    
    void merge_surface(const std::shared_ptr<Surface>& surface);

    const int& get_id() const;
    const Eigen::Vector3d& get_mean() const;
    const Eigen::Matrix3d& get_covariance() const;
    const Eigen::Matrix3d& get_eigenvectors() const;
    const Eigen::Vector3d& get_eigenvalues() const;
    const Eigen::Vector3d& get_normal() const;
    std::size_t get_total_point_size() const;
    const std::tuple<int, int, int>& get_color() const;
    bool is_expired() const;

    void connect(const std::shared_ptr<Vertex>& vertex);
    void connect(const std::shared_ptr<Vertex>& vertex, const std::set<std::shared_ptr<Vertex>>& nearby_vertices);
    void connect(const std::shared_ptr<Edge>& edge);
    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<InteriorPoint>& interior_point);
    void disconnect(const std::shared_ptr<Vertex>& vertex);
    void disconnect(const std::shared_ptr<Edge>& edge);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<InteriorPoint>& interior_point);
    void set_random_color();

    void refine_surface();

    void add_searchable_edge(const std::shared_ptr<Edge>& edge);
    void remove_searchable_edge(const std::shared_ptr<Edge>& edge);
    
private:
    bool deleting_ = false;
    bool is_expired_ = true;

    EdgeBVH edge_bvh_;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::set<std::shared_ptr<Vertex>> vertices_;
    std::set<std::shared_ptr<Edge>> edges_;
    std::set<std::shared_ptr<Face>> faces_;
    std::set<std::shared_ptr<InteriorPoint>> interior_points_;

    void add_point_to_surface_fitting(const Eigen::Vector3d& point, const Eigen::Vector3d& origin);
    void remove_point_from_surface_fitting(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    Eigen::Vector3d mean_;
    Eigen::Matrix3d covariance_;
    Eigen::Matrix3d eigenvectors_;
    Eigen::Vector3d eigenvalues_;
    Eigen::Vector3d normal_;

    std::tuple<int, int, int> color_;
};

bool operator<(const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs);
bool operator==(const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs);
bool operator>=(const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs);
bool operator!=(const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs);