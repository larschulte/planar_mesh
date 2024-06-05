#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"

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

    // cascade delete edges
    for (std::weak_ptr<Edge>& edge : edges_)
    {
        if (!edge.expired()) edge.lock()->cascade_delete_from_vertex(shared_from_this()); // cascade delete is depth first
    }

    // log
    std::cout << "Vertex " << id_ << " destroyed" << std::endl;
}

int Vertex::get_id() const 
{ 
    return id_; 
}

void Vertex::connect_edge(std::weak_ptr<Edge> edge) 
{
    // check pointer validity
    if (edge.expired()) throw std::runtime_error("Attempts to connect vertex with invalid edge.");
    auto edge_valid = edge.lock();

    // store
    edges_.push_back(edge_valid);
}

void Vertex::disconnect_edge(std::weak_ptr<Edge> edge) 
{
    // check pointer validity
    if (edge.expired()) throw std::runtime_error("Attempts to disconnect vertex from invalid edge.");
    auto edge_valid = edge.lock();

    // delete
    edges_.erase(std::remove_if(edges_.begin(), edges_.end(), [&](const std::weak_ptr<Edge> &e){return e.lock() == edge_valid;}), edges_.end());
}