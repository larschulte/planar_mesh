#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include "MeshObject/MeshObject.hpp"
#include "MeshObject/Settings.hpp"
#include "MeshObject/TriangleBVH.hpp"

// Forward declarations
class Vertex;
class Edge;
class Storage;
class Surface;
class InteriorPoint;
class GenericPoint;

class Face : public std::enable_shared_from_this<Face>, public MeshObject
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface> surface, const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    void initialize_(
        const std::shared_ptr<Storage>& storage, 
        const std::shared_ptr<Surface> surface, 
        const std::shared_ptr<Vertex>& vertex0, 
        const std::shared_ptr<Vertex>& vertex1, 
        const std::shared_ptr<Vertex>& vertex2,
        const std::shared_ptr<Edge>& edge0,
        const std::shared_ptr<Edge>& edge1,
        const std::shared_ptr<Edge>& edge2);
    void delete_();

public:
    // read write lock
    mutable std::shared_mutex rwlock_interior_points_;

    mutable std::shared_mutex rwlock_lifecycle_;

    void temp_initialize(const Eigen::Vector3d& end_point);

    void un_add_face();

    std::shared_ptr<Node> node;

    const int& get_id() const;
    const Eigen::Vector3d& get_center() const;
    const std::vector<std::shared_ptr<Vertex>>& get_vertices() const;
    std::vector<std::shared_ptr<InteriorPoint>> get_interior_points() const;
    const std::vector<std::shared_ptr<Edge>>& get_edges() const;
    const std::shared_ptr<Vertex>& get_vertex(int index) const;
    const std::shared_ptr<Surface>& get_surface() const;
    const std::shared_ptr<Vertex>& get_first_vertex() const;
    const Eigen::Vector3d& get_min() const;
    const Eigen::Vector3d& get_max() const;
    bool is_expired() const;
    bool is_deleting() const;
    bool is_searchable() const;
    bool has_vertex(const std::shared_ptr<Vertex>& vertex) const;

    bool intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction) const;
    bool intersects_point(const std::shared_ptr<GenericPoint>& generic_point);
    Eigen::Vector3d compute_intersection_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);

    void connect(const std::shared_ptr<InteriorPoint>& interior_point);
    void disconnect(const std::shared_ptr<InteriorPoint>& interior_point);

    bool is_connected_to_boundary_edges(std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& all_connected_faces, std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& all_connected_edges) const;
    bool is_non_manifold() const;

    void update_radius(const std::shared_ptr<GenericPoint>& generic_point);

private:
    static Settings settings_;

    Eigen::Vector3d center_;

    bool deleting_ = false;
    bool is_expired_ = true;

    bool can_self_destruct_ = true;

    Eigen::Vector3d min_;
    Eigen::Vector3d max_;
    Eigen::Vector3d v0_;
    Eigen::Vector3d v1_;
    Eigen::Vector3d v2_;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::vector<std::shared_ptr<Vertex>> vertices_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<InteriorPoint>> interior_points_;
    std::shared_ptr<Surface> surface_;

    std::shared_ptr<Vertex> first_vertex_;
};

bool operator<(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs);
bool operator==(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs);