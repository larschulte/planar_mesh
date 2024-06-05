#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

void Vertex::initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d pos)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create vertex with invalid storage.");
    auto storage_valid = storage.lock();

    // get id
    id_ = storage_valid->get_next_vertex_id();

    // store
    pos_ = pos;

    // log
    std::cout << "Vertex " << id_ << " created.\n";
}

void Vertex::delete_()
{
    // log
    std::cout << "Destroying vertex " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    while (!edges_.empty())
    {
        disconnect(*edges_.begin());
    }
    while (!faces_.empty())
    {
        disconnect(*faces_.begin());
    }
    while (!surfaces_.empty())
    {
        disconnect(*surfaces_.begin());
    }

    // log
    std::cout << "---------- vertex " << id_ << " destroyed" << std::endl;
}

int Vertex::get_id() const 
{ 
    return id_; 
}

Eigen::Vector3d Vertex::get_pos() const 
{ 
    return pos_; 
}

void Vertex::connect(std::weak_ptr<Edge> edge) 
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to connect vertex with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge.lock()->connect(shared_from_this());
}

void Vertex::connect(std::weak_ptr<Face> face) 
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to connect vertex with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face.lock()->connect(shared_from_this());
}

void Vertex::connect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to connect vertex with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface.lock()->connect(shared_from_this());
}

void Vertex::disconnect(std::weak_ptr<Edge> edge) 
{
    // check input
    if (edge.expired()) return;

    // disconnect
    bool erased = edges_.erase(edge);
    if (erased) edge.lock()->disconnect(shared_from_this());
}

void Vertex::disconnect(std::weak_ptr<Face> face)
{
    // check pointer validity
    if (face.expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face.lock()->disconnect(shared_from_this());
}

void Vertex::disconnect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface.lock()->disconnect(shared_from_this());
}

bool operator<(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}