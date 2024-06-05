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

    // set deletion flag
    deleting_ = true;

    // disconnect
    while (!vertices_.empty())
    {
        disconnect(*vertices_.begin());
    }
    while (!edges_.empty())
    {
        disconnect(*edges_.begin());
    }
    while (!faces_.empty())
    {
        disconnect(*faces_.begin());
    }

    // log
    std::cout << "---------- surface " << id_ << " destroyed" << std::endl;
}

int Surface::get_id() const
{
    return id_;
}

void Surface::connect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) throw std::runtime_error("Attempts to connect surface with invalid vertex.");

    // connect
    bool inserted = vertices_.insert(vertex.lock()).second;
    if (inserted) vertex.lock()->connect(shared_from_this());
}

void Surface::connect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to connect surface with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge.lock()->connect(shared_from_this());
}

void Surface::connect(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to connect surface with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face.lock()->connect(shared_from_this());
}

void Surface::disconnect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) return;

    // disconnect
    bool erased = vertices_.erase(vertex.lock());
    if (erased) vertex.lock()->disconnect(shared_from_this());
}

void Surface::disconnect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) return;

    // disconnect
    bool erased = edges_.erase(edge);
    if (erased) edge.lock()->disconnect(shared_from_this());
}

void Surface::disconnect(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face.lock()->disconnect(shared_from_this());
}

bool operator<(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs)
{
    // check pointer validity
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs)
{
    // check pointer validity
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}