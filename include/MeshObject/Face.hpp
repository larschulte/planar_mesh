#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

// Forward declarations
class Vertex;
class Edge;
class Storage;
class Surface;
class InteriorPoint;

class Face : public std::enable_shared_from_this<Face> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vertex> vertex0, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2);
    void delete_();

public:
    int get_id() const;
    Eigen::Vector3d get_center() const;
    std::set<std::weak_ptr<Vertex>> get_vertices() const;
    std::weak_ptr<Vertex> get_vertex(int index) const;
    std::weak_ptr<Surface> get_surface() const;

    bool intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);
    Eigen::Vector3d compute_intersection_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);

    void connect(std::weak_ptr<Vertex> vertex);
    void connect(std::weak_ptr<Edge> edge);
    void connect(std::weak_ptr<Surface> surface);
    void connect(std::weak_ptr<InteriorPoint> interior_point);
    void disconnect(std::weak_ptr<Vertex> vertex);
    void disconnect(std::weak_ptr<Edge> edge);
    void disconnect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<InteriorPoint> interior_point);

private:
    Eigen::Vector3d center_;

    bool deleting_ = false;
    bool is_searchable_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Vertex>> vertices_;
    std::set<std::weak_ptr<Edge>> edges_;
    std::set<std::weak_ptr<Surface>> surfaces_;
    std::set<std::weak_ptr<InteriorPoint>> interior_points_;
};

bool operator<(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs);
bool operator==(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs);