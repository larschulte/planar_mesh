#pragma once

#include <memory>
#include <set>

// Forward declarations
class Vertex;
class Face;
class Storage;
class Surface;

class Edge : public std::enable_shared_from_this<Edge> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2);
    void delete_();

public:
    int get_id() const;
    std::weak_ptr<Vertex> get_vertex(int index) const;

    void connect(std::weak_ptr<Vertex> vertex);
    void connect(std::weak_ptr<Face> face);
    void connect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<Vertex> vertex);
    void disconnect(std::weak_ptr<Face> face);
    void disconnect(std::weak_ptr<Surface> surface);

    bool has_vertex(std::weak_ptr<Vertex> vertex) const;
    bool is_boundary() const;
    void update_boundary_state();

    Eigen::Vector3d get_center() const;
    Eigen::Vector3d get_max() const;
    Eigen::Vector3d get_min() const;

    bool intersects_edge(std::weak_ptr<Vertex> vertex0, std::weak_ptr<Vertex> vertex1);

private:
    bool deleting_ = false;
    bool is_boundary_ = false;

    Eigen::Vector3d center_;
    Eigen::Vector3d max_;
    Eigen::Vector3d min_;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Vertex>> vertices_;
    std::set<std::weak_ptr<Face>> faces_;
    std::set<std::weak_ptr<Surface>> surfaces_;
};

bool operator<(const std::weak_ptr<Edge>& lhs, const std::weak_ptr<Edge>& rhs);
bool operator==(const std::weak_ptr<Edge>& lhs, const std::weak_ptr<Edge>& rhs);

bool segments_intersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);