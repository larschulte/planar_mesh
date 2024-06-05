#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Storage.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

void Surface::initialize_(std::weak_ptr<Storage> storage)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create surface with invalid storage.");
    auto storage_valid = storage.lock();

    // get id
    id_ = storage_valid->get_next_surface_id();

    // store
    storage_ = storage_valid;

    // log
    std::cout << "Surface " << id_ << " created.\n";
}

void Surface::delete_()
{
    // log
    std::cout << "Destroying surface " << id_ << std::endl;

    // disconnect all vertices
    for (auto vertex : vertices_)
    {
        if (vertex.expired()) throw std::runtime_error("Surface holds pointer to deleted vertex");
        vertex.lock()->disconnect_surface(shared_from_this());
    }

    // disconnect all edges
    for (auto edge : edges_)
    {
        if (edge.expired()) throw std::runtime_error("Surface holds pointer to deleted edge");
        edge.lock()->disconnect_surface(shared_from_this());
    }

    // disconnect all faces
    for (auto face : faces_)
    {
        if (face.expired()) throw std::runtime_error("Surface holds pointer to deleted face");
        face.lock()->disconnect_surface(shared_from_this());
    }

    // log
    std::cout << "Surface " << id_ << " destroyed" << std::endl;
}

int Surface::get_id() const
{
    return id_;
}

void Surface::connect_vertex(std::weak_ptr<Vertex> vertex)
{
    // check pointer validity
    if (vertex.expired()) throw std::runtime_error("Attempts to connect surface with invalid vertex.");
    auto vertex_valid = vertex.lock();

    // store
    vertices_.push_back(vertex_valid);
}

void Surface::disconnect_vertex(std::weak_ptr<Vertex> vertex)
{
    // check pointer validity
    if (vertex.expired()) throw std::runtime_error("Attempts to disconnect surface from invalid vertex.");
    auto vertex_valid = vertex.lock();

    // delete
    vertices_.erase(std::remove_if(vertices_.begin(), vertices_.end(), [&](const std::weak_ptr<Vertex> &v){return v.lock() == vertex_valid;}), vertices_.end());
}

void Surface::connect_edge(std::weak_ptr<Edge> edge)
{
    // check pointer validity
    if (edge.expired()) throw std::runtime_error("Attempts to connect surface with invalid edge.");
    auto edge_valid = edge.lock();

    // store
    edges_.push_back(edge_valid);
}

void Surface::disconnect_edge(std::weak_ptr<Edge> edge)
{
    // check pointer validity
    if (edge.expired()) throw std::runtime_error("Attempts to disconnect surface from invalid edge.");
    auto edge_valid = edge.lock();

    // delete
    edges_.erase(std::remove_if(edges_.begin(), edges_.end(), [&](const std::weak_ptr<Edge> &e){return e.lock() == edge_valid;}), edges_.end());
}

void Surface::connect_face(std::weak_ptr<Face> face)
{
    // check pointer validity
    if (face.expired()) throw std::runtime_error("Attempts to connect surface with invalid face.");
    auto face_valid = face.lock();

    // store
    faces_.push_back(face_valid);
}

void Surface::disconnect_face(std::weak_ptr<Face> face)
{
    // check pointer validity
    if (face.expired()) throw std::runtime_error("Attempts to disconnect surface from invalid face.");
    auto face_valid = face.lock();

    // delete
    faces_.erase(std::remove_if(faces_.begin(), faces_.end(), [&](const std::weak_ptr<Face> &f){return f.lock() == face_valid;}), faces_.end());
}