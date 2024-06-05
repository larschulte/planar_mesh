#include "MeshObject/Storage.hpp"
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
    auto edge1_valid = edge1.lock();
    auto edge2_valid = edge2.lock();
    auto edge3_valid = edge3.lock();
    
    // get id
    id_ = storage_valid->get_next_face_id();

    // sort edges by edge id
    std::vector<std::shared_ptr<Edge>> edges = {edge1_valid, edge2_valid, edge3_valid};
    std::sort(edges.begin(), edges.end(), [](const std::shared_ptr<Edge> &a, const std::shared_ptr<Edge> &b){return a->get_id() < b->get_id();});
    edge1_valid = edges[0];
    edge2_valid = edges[1];
    edge3_valid = edges[2];

    // store
    edge1_ = edge1_valid;
    edge2_ = edge2_valid;
    edge3_ = edge3_valid;
    storage_ = storage_valid;

    // register face with edges
    edge1_valid->connect_face(shared_from_this());
    edge2_valid->connect_face(shared_from_this());
    edge3_valid->connect_face(shared_from_this());

    // log
    std::cout << "Face " << id_ << " created between edge " << edge1_valid->get_id() << ", edge " << edge2_valid->get_id() << " and edge " << edge3_valid->get_id() << std::endl;
}

void Face::delete_()
{
    // log
    std::cout << "Destroying face " << id_ << std::endl;

    // disconnect edges
    if (!edge1_.expired()) edge1_.lock()->disconnect_face(shared_from_this());
    if (!edge2_.expired()) edge2_.lock()->disconnect_face(shared_from_this());
    if (!edge3_.expired()) edge3_.lock()->disconnect_face(shared_from_this());

    // log
    std::cout << "Face " << id_ << " destroyed" << std::endl;
}

void Face::cascade_delete_from_edge(std::weak_ptr<Edge> edge)
{
    // check if edge is valid
    if (edge.expired()) throw std::runtime_error("Attempts to cascade delete face from invalid edge.");
    if (edge1_.expired()) throw std::runtime_error("Face holds pointer to deleted edge1");
    if (edge2_.expired()) throw std::runtime_error("Face holds pointer to deleted edge2");
    if (edge3_.expired()) throw std::runtime_error("Face holds pointer to deleted edge3");
    auto edge_valid = edge.lock();
    auto edge1_valid = edge1_.lock();
    auto edge2_valid = edge2_.lock();
    auto edge3_valid = edge3_.lock();

    // check if the edge is one of the edges of the face    
    if (edge_valid != edge1_valid && edge_valid != edge2_valid && edge_valid != edge3_valid) throw std::runtime_error("Edge is not part of the face.");

    // check if storage is valid
    if (storage_.expired()) throw std::runtime_error("Face holds pointer to deleted storage.");
    auto storage_valid = storage_.lock();
    
    // delete face
    storage_valid->delete_face(shared_from_this());
}

void Face::connect_surface(std::weak_ptr<Surface> surface)
{
    // check if surface is valid
    if (surface.expired()) throw std::runtime_error("Attempts to connect face with invalid surface.");
    auto surface_valid = surface.lock();

    // store
    surfaces_.push_back(surface_valid);
}

void Face::disconnect_surface(std::weak_ptr<Surface> surface)
{
    // check if surface is valid
    if (surface.expired()) throw std::runtime_error("Attempts to disconnect face from invalid surface.");
    auto surface_valid = surface.lock();

    // delete
    surfaces_.erase(std::remove_if(surfaces_.begin(), surfaces_.end(), [&](const std::weak_ptr<Surface> &s){return s.lock() == surface_valid;}), surfaces_.end());
}