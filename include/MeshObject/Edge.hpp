#pragma once

#include <memory>
#include <unordered_set>
#include <map>

#include "MeshObject/MeshObject.hpp"
#include "MeshObject/Settings.hpp"
// #include "MeshObject/EdgeBVH.hpp"

// Forward declarations
class Vertex;
class Face;
class Storage;
class Surface;

class Edge : public std::enable_shared_from_this<Edge>, public MeshObject
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    void delete_();

public:

    ~Edge();

    // std::shared_ptr<EdgeBVH::Node> node;

    // read and write lock
    mutable std::shared_mutex rwlock_faces_;

    mutable std::shared_mutex rwlock_lifecycle_;

    const int& get_id() const;
    std::shared_ptr<Vertex> get_vertex(int index) const;
    std::weak_ptr<Vertex> get_vertex_weak_ptr(int index) const;
    std::shared_ptr<Surface> get_surface() const;
    std::vector<std::weak_ptr<Face>> get_faces() const;
    const Eigen::Vector3d& get_center() const;
    const Eigen::Vector3d& get_max() const;
    const Eigen::Vector3d& get_min() const;
    const double& get_length() const;
    bool is_expired() const;

    void connect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Face>& face);

    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);

    void set_can_self_destruct(bool can_self_destruct);

    bool is_connected_to_boundary_edges(std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& all_connected_faces, std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& all_connected_edges) const;

    bool has_vertex(const std::shared_ptr<Vertex>& vertex) const;
    bool is_boundary() const;
    bool is_singular() const;
    bool is_deleting() const;
    void update_boundary_state();

    bool is_non_manifold() const;

    bool intersects_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1);

    void set_under_surface_deletion(bool flag);

private:
    static Settings settings_;

    bool deleting_ = false;
    bool is_expired_ = true;
    bool is_boundary_;
    bool is_singular_;
    bool is_searchable_;
    bool under_surface_deletion_ = false;

    bool can_self_destruct_ = true;

    Eigen::Vector3d center_;
    Eigen::Vector3d max_;
    Eigen::Vector3d min_;
    double length_;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::vector<std::weak_ptr<Vertex>> vertices_;
    std::vector<std::weak_ptr<Face>> faces_;
    std::weak_ptr<Surface> surface_;
};

bool operator<(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs);
bool operator==(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs);

bool segments_intersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);