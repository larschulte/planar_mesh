#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"

void Storage::add_vertex(Eigen::Vector3d pos) 
{
    // initialize
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), pos);

    // put into list
    vertices_.push_back(vertex);
}

void Storage::add_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2) 
{    
    // initialize
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), vertex1, vertex2);

    // put into list
    edges_.push_back(edge);

}

void Storage::add_face(std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3) 
{
    // initialize
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), edge1, edge2, edge3);

    // put into list
    faces_.push_back(face);
}

void Storage::add_surface() 
{
    // initialize
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this());

    // put into list
    surfaces_.push_back(surface);
}

// need to ensure the vertex/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vertex(std::weak_ptr<Vertex> vertex) 
{
    // check if valid
    if (vertex.expired()) throw std::runtime_error("Attempts to delete expired vertex.");
    auto vertex_valid = vertex.lock();

    // member delete
    vertex_valid->delete_();

    // storage delete
    vertices_.erase(std::remove(vertices_.begin(), vertices_.end(), vertex_valid), vertices_.end());
}

void Storage::delete_edge(std::weak_ptr<Edge> edge) 
{
    // check if valid
    if (edge.expired()) throw std::runtime_error("Attempts to delete expired edge.");
    auto edge_valid = edge.lock();

    // member delete
    edge_valid->delete_();

    // storage delete
    edges_.erase(std::remove(edges_.begin(), edges_.end(), edge_valid), edges_.end());
}

void Storage::delete_face(std::weak_ptr<Face> face) 
{
    // check if valid
    if (face.expired()) throw std::runtime_error("Attempts to delete expired face.");
    auto face_valid = face.lock();

    // member delete
    face_valid->delete_();

    // storage delete
    faces_.erase(std::remove(faces_.begin(), faces_.end(), face_valid), faces_.end());
}

void Storage::delete_surface(std::weak_ptr<Surface> surface) 
{
    // check if valid
    if (surface.expired()) throw std::runtime_error("Attempts to delete expired surface.");
    auto surface_valid = surface.lock();

    // member delete
    surface_valid->delete_();

    // storage delete
    surfaces_.erase(std::remove(surfaces_.begin(), surfaces_.end(), surface_valid), surfaces_.end());
}

// get id
int Storage::get_next_vertex_id() { return next_vertex_id_++; }
int Storage::get_next_edge_id() { return next_edge_id_++; }
int Storage::get_next_face_id() { return next_face_id_++; }
int Storage::get_next_surface_id() { return next_surface_id_++; }