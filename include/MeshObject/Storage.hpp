#pragma once

#include <memory>
#include <set>
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

class Storage : public std::enable_shared_from_this<Storage> 
{
public: // to user
    Storage();
    ~Storage();

    const std::shared_ptr<Vertex>& add_vertex(const Eigen::Vector3d& origin, const Eigen::Vector3d& position);
    const std::shared_ptr<Vertex>& add_vertex(const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius);
    const std::shared_ptr<Edge>& add_edge(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    const std::shared_ptr<Face>& add_face(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3);
    const std::shared_ptr<Surface>& add_surface();
    const std::shared_ptr<GenericPoint>& add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    const std::shared_ptr<InteriorPoint>& add_interior_point(const std::shared_ptr<Face>& face, const Eigen::Vector3d& position, const Eigen::Vector3d& origin);

    void delete_vertex(const std::shared_ptr<Vertex>& vertex);
    void delete_edge(const std::shared_ptr<Edge>& edge);
    void delete_face(const std::shared_ptr<Face>& face);
    void delete_surface(const std::shared_ptr<Surface>& surface);
    void delete_genertic_point(const std::shared_ptr<GenericPoint>& genertic_point);
    void delete_interior_point(const std::shared_ptr<InteriorPoint>& interior_point);

    bool can_reverse_radius_search();
    std::set<std::shared_ptr<Vertex>> reverse_radius_search(const Eigen::Vector3d& point);
    std::set<std::shared_ptr<Face>> face_intersection_search(const Eigen::Vector3d& origin, const Eigen::Vector3d& point);

    const std::set<std::shared_ptr<Vertex>>& get_vertices() const;
    const std::set<std::shared_ptr<Edge>>& get_edges() const;
    const std::shared_ptr<Edge>& get_edge(std::shared_ptr<Vertex> vertex1, std::shared_ptr<Vertex> vertex2) const;
    const std::set<std::shared_ptr<Face>>& get_faces() const;
    const std::set<std::shared_ptr<Surface>>& get_surfaces() const;

    const std::vector<std::shared_ptr<Vertex>>& get_rrs_vertices();
    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map() const;
    bool is_expired() const;

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

    RRSTree rrs_tree_;
    TriangleBVH triangle_bvh_;

    std::set<std::shared_ptr<Vertex>> vertices_;
    std::set<std::shared_ptr<Edge>> edges_;
    std::set<std::shared_ptr<Face>> faces_;
    std::set<std::shared_ptr<Surface>> surfaces_;
    std::set<std::shared_ptr<GenericPoint>> genertic_points_;
    std::set<std::shared_ptr<InteriorPoint>> interior_points_;

    int next_vertex_id_ = 0;
    int next_edge_id_ = 0;
    int next_face_id_ = 0;
    int next_surface_id_ = 0;
    int next_genertic_point_id_ = 0;
    int next_interior_point_id_ = 0;
};