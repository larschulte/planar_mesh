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
    void initialize_(std::shared_ptr<Storage> storage, std::shared_ptr<Vertex> vertex0, std::shared_ptr<Vertex> vertex1, std::shared_ptr<Vertex> vertex2);
    void delete_();

public:
    const int& get_id() const;
    const Eigen::Vector3d& get_center() const;
    const std::set<std::shared_ptr<Vertex>>& get_vertices() const;
    const std::shared_ptr<Vertex>& get_vertex(int index) const;
    const std::shared_ptr<Surface>& get_surface() const;
    bool is_expired() const;

    bool intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);
    Eigen::Vector3d compute_intersection_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);

    void connect(const std::shared_ptr<Vertex>& vertex);
    void connect(const std::shared_ptr<Edge>& edge);
    void connect(const std::shared_ptr<Surface>& surface);
    void connect(const std::shared_ptr<InteriorPoint>& interior_point);
    void disconnect(const std::shared_ptr<Vertex>& vertex);
    void disconnect(const std::shared_ptr<Edge>& edge);
    void disconnect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<InteriorPoint>& interior_point);

private:
    Eigen::Vector3d center_;

    bool deleting_ = false;
    bool is_searchable_ = false;
    bool is_expired_ = true;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::set<std::shared_ptr<Vertex>> vertices_;
    std::set<std::shared_ptr<Edge>> edges_;
    std::set<std::shared_ptr<Surface>> surfaces_;
    std::set<std::shared_ptr<InteriorPoint>> interior_points_;
};

bool operator<(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs);
bool operator==(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs);