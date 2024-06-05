#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"

std::weak_ptr<Vertex> Storage::add_vertex(Eigen::Vector3d pos) 
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), pos);

    // store
    vertices_.insert(vertex);

    // return
    return vertex;
}

std::weak_ptr<Edge> Storage::add_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2) 
{    
    // create
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), vertex1, vertex2);

    // store
    edges_.insert(edge);

    // return
    return edge;

}

std::weak_ptr<Face> Storage::add_face(std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3) 
{
    // create
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), edge1, edge2, edge3);

    // store
    faces_.insert(face);

    // return
    return face;
}

std::weak_ptr<Surface> Storage::add_surface() 
{
    // create
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this());

    // store
    surfaces_.insert(surface);

    // return
    return surface;
}

// need to ensure the vertex/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vertex(std::weak_ptr<Vertex> vertex) 
{
    // check input
    if (vertex.expired()) throw std::runtime_error("Attempts to delete expired vertex.");

    // member delete
    vertex.lock()->delete_();

    // storage delete
    vertices_.erase(vertex.lock());
}

void Storage::delete_edge(std::weak_ptr<Edge> edge) 
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to delete expired edge.");

    // member delete
    edge.lock()->delete_();

    // storage delete
    edges_.erase(edge.lock());
}

void Storage::delete_face(std::weak_ptr<Face> face) 
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to delete expired face.");

    // member delete
    face.lock()->delete_();

    // storage delete
    faces_.erase(face.lock());
}

void Storage::delete_surface(std::weak_ptr<Surface> surface) 
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to delete expired surface.");

    // member delete
    surface.lock()->delete_();

    // storage delete
    surfaces_.erase(surface.lock());
}

// get id
int Storage::get_next_vertex_id() { return next_vertex_id_++; }
int Storage::get_next_edge_id() { return next_edge_id_++; }
int Storage::get_next_face_id() { return next_face_id_++; }
int Storage::get_next_surface_id() { return next_surface_id_++; }