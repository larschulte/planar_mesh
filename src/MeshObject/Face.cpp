#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

void Face::initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create face with invalid storage.");
    if (edge1.expired()) throw std::runtime_error("Attempts to create face with invalid edge1.");
    if (edge2.expired()) throw std::runtime_error("Attempts to create face with invalid edge2.");
    if (edge3.expired()) throw std::runtime_error("Attempts to create face with invalid edge3.");
    auto storage_valid = storage.lock();
    auto edge0_valid = edge1.lock();
    auto edge1_valid = edge2.lock();
    auto edge2_valid = edge3.lock();
    
    // get id
    id_ = storage_valid->get_next_face_id();

    // sort edges by edge id
    std::vector<std::shared_ptr<Edge>> edges = {edge0_valid, edge1_valid, edge2_valid};
    std::sort(edges.begin(), edges.end(), [](const std::shared_ptr<Edge> &a, const std::shared_ptr<Edge> &b){return a->get_id() < b->get_id();});

    // store
    storage_ = storage_valid;

    // connect
    connect(edges[0]);
    connect(edges[1]);
    connect(edges[2]);

    // log
    std::cout << "Face " << id_ << " created between edge " << edge0_valid->get_id() << ", edge " << edge1_valid->get_id() << " and edge " << edge2_valid->get_id() << std::endl;
}

void Face::delete_()
{
    // log
    std::cout << "Destroying face " << id_ << std::endl;

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
    while (!surfaces_.empty())
    {
        disconnect(*surfaces_.begin());
    }
    
    // log
    std::cout << "---------- face " << id_ << " destroyed" << std::endl;
}

int Face::get_id() const
{
    return id_;
}

void Face::connect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) throw std::runtime_error("Attempts to connect face with invalid vertex.");

    // connect
    bool inserted = vertices_.insert(vertex).second;
    if (inserted) vertex.lock()->connect(shared_from_this());

    // check size
    if (vertices_.size() > 3) throw std::runtime_error("Face connected to more than 3 vertices.");
}

void Face::connect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to connect face with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge.lock()->connect(shared_from_this());

    // check size
    if (edges_.size() > 3) throw std::runtime_error("Face connected to more than 3 edges.");
}

void Face::connect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to connect face with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface.lock()->connect(shared_from_this());
}

void Face::disconnect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) return;

    // delete
    bool erased = vertices_.erase(vertex);
    if (erased) vertex.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_.lock()->delete_face(shared_from_this());
}

void Face::disconnect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) return;

    // delete
    bool erased = edges_.erase(edge);
    if (erased) edge.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_.lock()->delete_face(shared_from_this());
}

void Face::disconnect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) return;

    // delete
    bool erased = surfaces_.erase(surface);
    if (erased) surface.lock()->disconnect(shared_from_this());
}

bool operator<(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}