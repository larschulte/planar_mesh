#pragma once

#include <memory>
#include <Eigen/Dense>
#include <set>

// Forward declarations
class Edge;
class Face;
class Storage;
class Surface;

class Vertex : public std::enable_shared_from_this<Vertex> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin);
    void initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin, double radius);
    void delete_();

public:
    int get_id() const;
    Eigen::Vector3d get_position() const;
    Eigen::Vector3d get_projected_position() const;
    Eigen::Vector3d get_origin() const;
    std::weak_ptr<Surface> get_surface() const;

    void connect(std::weak_ptr<Edge> edge);
    void connect(std::weak_ptr<Face> face);
    void connect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<Edge> edge);
    void disconnect(std::weak_ptr<Face> face);
    void disconnect(std::weak_ptr<Surface> surface);

    void update_boundary_state();

public: // for reverse radius search
    void set_reverse_radius_search_radius(double radius);
    Eigen::Vector3d get_min() const;
    Eigen::Vector3d get_max() const;
    double get_radius() const;
    bool contains(const Eigen::Vector3d& point) const;
    bool approx_contains(const Eigen::Vector3d& point) const;

private: // for reverse radius search
    double reverse_search_radius_;
    Eigen::Vector3d min_;
    Eigen::Vector3d max_;

private:
    bool deleting_ = false;
    bool is_boundary_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Edge>> edges_;
    std::set<std::weak_ptr<Face>> faces_;
    std::set<std::weak_ptr<Surface>> surfaces_;

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
};

bool operator<(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs);
bool operator==(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs);