#pragma once

#include <memory>
#include <set>
#include <Eigen/Dense> 

// forward declarations
class Vertex;
class Edge;
class Face;
class Surface;

class Storage : public std::enable_shared_from_this<Storage> 
{
public:
    std::weak_ptr<Vertex> add_vertex(Eigen::Vector3d pos);
    std::weak_ptr<Edge> add_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2);
    std::weak_ptr<Face> add_face(std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3);
    std::weak_ptr<Surface> add_surface();

    void delete_vertex(std::weak_ptr<Vertex> vertex);
    void delete_edge(std::weak_ptr<Edge> edge);
    void delete_face(std::weak_ptr<Face> face);
    void delete_surface(std::weak_ptr<Surface> surface);

    int get_next_vertex_id();
    int get_next_edge_id();
    int get_next_face_id();
    int get_next_surface_id();

    std::set<std::shared_ptr<Vertex>> vertices_;
    std::set<std::shared_ptr<Edge>> edges_;
    std::set<std::shared_ptr<Face>> faces_;
    std::set<std::shared_ptr<Surface>> surfaces_;

private:
    int next_vertex_id_ = 0;
    int next_edge_id_ = 0;
    int next_face_id_ = 0;
    int next_surface_id_ = 0;
};