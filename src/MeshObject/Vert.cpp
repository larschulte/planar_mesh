#include "MeshObject/Storage.hpp"
#include "MeshObject/Vert.hpp"
#include "MeshObject/Edge.hpp"

void Vert::initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d pos)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create vert with invalid storage.");
    auto storage_valid = storage.lock();

    // get id
    id_ = storage_valid->get_next_vert_id();

    // store
    pos_ = pos;

    // log
    std::cout << "Vert " << id_ << " created.\n";
}

void Vert::delete_()
{
    // log
    std::cout << "Destroying vert " << id_ << std::endl;

    // cascade delete edges
    for (std::weak_ptr<Edge>& edge : edges_)
    {
        if (!edge.expired()) edge.lock()->cascade_delete_from_vert(shared_from_this()); // cascade delete is depth first
    }

    // log
    std::cout << "Vert " << id_ << " destroyed" << std::endl;
}

int Vert::get_id() const 
{ 
    return id_; 
}

void Vert::connect_edge(std::weak_ptr<Edge> edge) 
{
    // check pointer validity
    if (edge.expired()) throw std::runtime_error("Attempts to connect vert with invalid edge.");
    auto edge_valid = edge.lock();

    // store
    edges_.push_back(edge_valid);
}

void Vert::disconnect_edge(std::weak_ptr<Edge> edge) 
{
    // check pointer validity
    if (edge.expired()) throw std::runtime_error("Attempts to disconnect vert from invalid edge.");
    auto edge_valid = edge.lock();

    // delete
    edges_.erase(std::remove_if(edges_.begin(), edges_.end(), [&](const std::weak_ptr<Edge> &e){return e.lock() == edge_valid;}), edges_.end());
}