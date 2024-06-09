#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

// forward declarations
class Vertex;
class Edge;
class Face;
class Storage;

class Surface : public std::enable_shared_from_this<Surface> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage);
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Surface> surface1, std::weak_ptr<Surface> surface2);
    void delete_();

public:
    int get_id() const;
    double compute_point_to_surface_distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const;
    Eigen::Vector3d compute_point_to_surface_position(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const;
    
    Eigen::Vector3d get_mean() const;
    Eigen::Matrix3d get_covariance() const;
    Eigen::Matrix3d get_eigenvectors() const;
    Eigen::Vector3d get_eigenvalues() const;
    Eigen::Vector3d get_normal() const;
    int get_total_point_size() const;

    void connect(std::weak_ptr<Vertex> vertex);
    void connect(std::weak_ptr<Edge> edge);
    void connect(std::weak_ptr<Face> face);
    void connect(std::weak_ptr<InteriorPoint> interior_point);
    void disconnect(std::weak_ptr<Vertex> vertex);
    void disconnect(std::weak_ptr<Edge> edge);
    void disconnect(std::weak_ptr<Face> face);
    void disconnect(std::weak_ptr<InteriorPoint> interior_point);
    
private:
    bool deleting_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Vertex>> vertices_;
    std::set<std::weak_ptr<Edge>> edges_;
    std::set<std::weak_ptr<Face>> faces_;
    std::set<std::weak_ptr<InteriorPoint>> interior_points_;

    void add_point_to_surface_fitting(Eigen::Vector3d point, Eigen::Vector3d origin);
    Eigen::Vector3d mean_;
    Eigen::Matrix3d covariance_;
    Eigen::Matrix3d eigenvectors_;
    Eigen::Vector3d eigenvalues_;
    Eigen::Vector3d normal_;

    void set_random_color_();
    std::tuple<int, int, int> color_;
};

bool operator<(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs);
bool operator==(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs);
bool operator>=(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs);