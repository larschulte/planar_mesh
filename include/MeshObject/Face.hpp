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
    void delete_();

public:
    void temp_initialize(const Eigen::Vector3d& end_point);

    std::shared_ptr<Node> node;

    const int& get_id() const;
    const Eigen::Vector3d& get_center() const;
    const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& get_vertices() const;
    const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& get_interior_points() const;
    const std::shared_ptr<Vertex>& get_vertex(int index) const;
    const std::shared_ptr<Vertex>& get_first_vertex() const;
    const std::shared_ptr<Surface>& get_surface() const;
    const Eigen::Vector3d& get_min() const;
    const Eigen::Vector3d& get_max() const;
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_sibling_faces() const;
    bool is_expired() const;
    bool has_vertex(const std::shared_ptr<Vertex>& vertex) const;

    bool intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);
    bool intersects_point(const std::shared_ptr<GenericPoint>& generic_point);
    Eigen::Vector3d compute_intersection_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction);

    void connect(const std::shared_ptr<Vertex>& vertex);
    void connect(const std::shared_ptr<Edge>& edge);
    void connect(const std::shared_ptr<Surface>& surface);
    void connect(const std::shared_ptr<InteriorPoint>& interior_point);
    void connect(const std::shared_ptr<Face>& sibling_face);
    void disconnect(const std::shared_ptr<Vertex>& vertex);
    void disconnect(const std::shared_ptr<Edge>& edge);
    void disconnect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<InteriorPoint>& interior_point);
    void disconnect(const std::shared_ptr<Face>& sibling_face);

    void swap(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);

    void update_confirmed_status();
    bool is_confirmed() const;

    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);

    void update_radius(const std::shared_ptr<GenericPoint>& generic_point);

private:
    static Settings settings_;

    Eigen::Vector3d center_;

    bool deleting_ = false;
    bool is_searchable_ = false;
    bool is_expired_ = true;
    bool is_confirmed_ = false;

    bool can_self_destruct_ = true;

    Eigen::Vector3d min_;
    Eigen::Vector3d max_;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points_;
    std::shared_ptr<Surface> surface_;

    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> sibling_faces_;

    std::shared_ptr<Vertex> first_vertex_;
};

bool operator<(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs);
bool operator==(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs);