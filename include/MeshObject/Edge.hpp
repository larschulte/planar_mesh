#pragma once

#include <memory>
#include <unordered_set>
#include <map>

#include "MeshObject/MeshObject.hpp"

// Forward declarations
class Vertex;
class Face;
class Storage;
class Surface;

class Edge : public std::enable_shared_from_this<Edge>, public MeshObject
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    void delete_();

public:
    const int& get_id() const;
    const std::shared_ptr<Vertex>& get_vertex(int index) const;
    const std::shared_ptr<Surface>& get_surface() const;
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_faces() const;
    const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& get_sibling_edges() const;
    const Eigen::Vector3d& get_center() const;
    const Eigen::Vector3d& get_max() const;
    const Eigen::Vector3d& get_min() const;
    const double& get_length() const;
    bool is_expired() const;

    void connect(const std::shared_ptr<Vertex>& vertex);
    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void connect(const std::shared_ptr<Edge>& sibling_edge);
    void disconnect(const std::shared_ptr<Vertex>& vertex);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<Edge>& sibling_edge);

    void swap(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);

    void update_confirmed_status();
    bool is_confirmed() const;

    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);

    bool has_vertex(const std::shared_ptr<Vertex>& vertex) const;
    bool is_boundary() const;
    bool is_singular() const;
    void update_boundary_state();
    void update_singular_state();
    void update_searchable_state();
    void remove_searchable_state();

    bool intersects_edge(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1);

private:
    bool deleting_ = false;
    bool is_expired_ = true;
    bool is_boundary_;
    bool is_singular_;
    bool is_searchable_;

    bool can_self_destruct_ = true;

    std::size_t num_confirmed_faces = 0;
    bool is_confirmed_ = false;

    Eigen::Vector3d center_;
    Eigen::Vector3d max_;
    Eigen::Vector3d min_;
    double length_;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::shared_ptr<Surface> surface_;
    
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> sibling_edges_;
};

bool operator<(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs);
bool operator==(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs);

bool segments_intersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);