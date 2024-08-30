#pragma once

#include <memory>
#include <unordered_set>
#include <Eigen/Dense> 

#include "MeshObject/RRSTree.hpp"
#include "MeshObject/TriangleBVH.hpp"

// forward declarations
class Vertex;
class Edge;
class Face;
class Surface;
class GenericPoint;
class InteriorPoint;

enum class DeletedPointStorage
{
    NONE,
    GENERIC,
    PENETRATED,
    RADIUS_CHANGE
};

class Storage : public std::enable_shared_from_this<Storage> 
{
public: // to user
    Storage();
    ~Storage();

    const std::shared_ptr<Vertex>& add_vertex(const Eigen::Vector3d& origin, const Eigen::Vector3d& position);
    const std::shared_ptr<Vertex>& add_vertex(const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius);
    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<GenericPoint>& generic_point);
    const std::shared_ptr<Edge>& add_edge(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    const std::shared_ptr<Face>& add_face(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3);
    const std::shared_ptr<Surface>& add_surface();
    const std::shared_ptr<GenericPoint>& add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    const std::shared_ptr<InteriorPoint>& add_interior_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    const std::shared_ptr<InteriorPoint>& add_interior_point(const std::shared_ptr<GenericPoint>& generic_point);
    
    void delete_vertex(const std::shared_ptr<Vertex>& vertex);
    void delete_edge(const std::shared_ptr<Edge>& edge);
    void delete_face(const std::shared_ptr<Face>& face);
    void delete_surface(const std::shared_ptr<Surface>& surface);
    void delete_interior_point(const std::shared_ptr<InteriorPoint>& interior_point);

    void set_deleted_points_storage_name(const DeletedPointStorage& storage_name);
    DeletedPointStorage get_deleted_points_storage_name() const;
    void add_deleted_point(const std::shared_ptr<Vertex>& vertex);
    void add_deleted_point(const std::shared_ptr<InteriorPoint>& interior_point);
    const std::shared_ptr<GenericPoint>& add_generic_point(const std::shared_ptr<Vertex>& vertex);
    const std::shared_ptr<GenericPoint>& add_generic_point(const std::shared_ptr<InteriorPoint>& interior_point);
    const std::shared_ptr<GenericPoint>& add_penetrated_point(const std::shared_ptr<Vertex>& vertex);
    const std::shared_ptr<GenericPoint>& add_penetrated_point(const std::shared_ptr<InteriorPoint>& interior_point);
    const std::shared_ptr<GenericPoint>& add_radius_point(const std::shared_ptr<Vertex>& vertex);
    const std::shared_ptr<GenericPoint>& add_radius_point(const std::shared_ptr<InteriorPoint>& interior_point);
    void delete_generic_point(const std::shared_ptr<GenericPoint>& genertic_point);
    void delete_penetrated_point(const std::shared_ptr<GenericPoint>& penetrated_point);
    void delete_radius_point(const std::shared_ptr<GenericPoint>& radius_point);
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> pop_generic_points();
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> pop_penetrated_points();
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> pop_radius_points();

    bool can_reverse_radius_search();
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> reverse_radius_search(const Eigen::Vector3d& point);
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point);
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> reverse_radius_search(const std::shared_ptr<Vertex>& vertex);
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> face_intersection_search(const Eigen::Vector3d& origin, const Eigen::Vector3d& point);
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> face_intersection_search(const std::shared_ptr<GenericPoint>& generic_point);

    const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& get_vertices() const;
    const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& get_edges() const;
    const std::shared_ptr<Edge>& get_edge(std::shared_ptr<Vertex> vertex1, std::shared_ptr<Vertex> vertex2) const;
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_faces() const;
    const std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& get_surfaces() const;
    const std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash>& get_generic_points() const;
    void clear_generic_points();
    const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& get_interior_points() const;
    const std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash>& get_penetrated_points() const;
    void clear_penetrated_points();

    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices();
    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map() const;
    bool is_expired() const;

    void set_penetrating_point(const std::shared_ptr<GenericPoint>& generic_point);
    const Eigen::Vector3d& get_penetrating_point();
    void clear_penetrating_point();
    bool has_penetrating_point() const;

    void print_rrs() const;
    void print_bvh() const;
    void rebuild_tree();

private: // to Vertex and Face class
    friend class Vertex;
    friend class Face;
    void add_searchable_vertex(const std::shared_ptr<Vertex>& vertex);
    void remove_searchable_vertex(const std::shared_ptr<Vertex>& vertex);

    void add_searchable_face(const std::shared_ptr<Face>& face);
    void remove_searchable_face(const std::shared_ptr<Face>& face);

public: // to MeshObject class
    int get_next_vertex_id();
    int get_next_edge_id();
    int get_next_face_id();
    int get_next_surface_id();
    int get_next_generic_point_id();
    int get_next_interior_point_id();

private:
    bool is_expired_ = true;

    Eigen::Vector3d penetrating_point_;
    bool has_penetrating_point_ = false;

    RRSTree rrs_tree_;
    TriangleBVH triangle_bvh_;

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points_;

    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> genertic_points_;
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> penetrated_points_;
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> radius_points_;

    DeletedPointStorage deleted_points_storage_name_ = DeletedPointStorage::GENERIC;

    int next_vertex_id_ = 0;
    int next_edge_id_ = 0;
    int next_face_id_ = 0;
    int next_surface_id_ = 0;
    int next_genertic_point_id_ = 0;
    int next_interior_point_id_ = 0;
};