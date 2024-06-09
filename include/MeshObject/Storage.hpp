#pragma once

#include <memory>
#include <set>
#include <Eigen/Dense> 

#include "eye_patch/RRSTree.hpp"
#include "eye_patch/TriangleBVH.hpp"

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
    std::weak_ptr<Vertex> add_vertex(Eigen::Vector3d pos, Eigen::Vector3d origin);
    std::weak_ptr<Edge> add_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2);
    std::weak_ptr<Face> add_face(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2, std::weak_ptr<Vertex> vertex3);
    std::weak_ptr<Surface> add_surface();
    std::weak_ptr<Surface> add_surface(std::weak_ptr<Surface> surface1, std::weak_ptr<Surface> surface2);
    std::weak_ptr<GenericPoint> add_generic_point(Eigen::Vector3d pos, Eigen::Vector3d origin);
    std::weak_ptr<InteriorPoint> add_interior_point(std::weak_ptr<Face> face, Eigen::Vector3d pos, Eigen::Vector3d origin);

    void delete_vertex(std::weak_ptr<Vertex> vertex);
    void delete_edge(std::weak_ptr<Edge> edge);
    void delete_face(std::weak_ptr<Face> face);
    void delete_surface(std::weak_ptr<Surface> surface);
    void delete_genertic_point(std::weak_ptr<GenericPoint> genertic_point);
    void delete_interior_point(std::weak_ptr<InteriorPoint> interior_point);

    bool can_reverse_radius_search();
    std::set<std::weak_ptr<Vertex>> reverse_radius_search(Eigen::Vector3d point);
    std::set<std::weak_ptr<Face>> face_intersection_search(Eigen::Vector3d origin, Eigen::Vector3d point);

    std::set<std::weak_ptr<Vertex>> get_vertices() const;
    std::set<std::weak_ptr<Edge>> get_edges() const;
    std::set<std::weak_ptr<Face>> get_faces() const;
    std::map<std::weak_ptr<Vertex>, int> get_vertex_to_cloud_indices_map() const;

    std::weak_ptr<Edge> get_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2) const;

public: // to Vertex and Face class
    void add_searchable_vertex(std::weak_ptr<Vertex> vertex);
    void remove_searchable_vertex(std::weak_ptr<Vertex> vertex);

    void add_searchable_face(std::weak_ptr<Face> face);
    void remove_searchable_face(std::weak_ptr<Face> face);

public: // to MeshObject class
    int get_next_vertex_id();
    int get_next_edge_id();
    int get_next_face_id();
    int get_next_surface_id();
    int get_next_generic_point_id();
    int get_next_interior_point_id();

private:
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