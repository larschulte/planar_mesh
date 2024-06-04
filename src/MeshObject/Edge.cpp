#include "MeshObject/Storage.hpp"
#include "MeshObject/Vert.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

void Edge::initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vert> vert1, std::weak_ptr<Vert> vert2)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create edge with invalid storage.");
    if (vert1.expired()) throw std::runtime_error("Attempts to create edge with invalid vert1.");
    if (vert2.expired()) throw std::runtime_error("Attempts to create edge with invalid vert2.");
    auto storage_valid = storage.lock();
    auto vert1_valid = vert1.lock();
    auto vert2_valid = vert2.lock();

    // get id
    id_ = storage_valid->get_next_edge_id();

    // sort
    if (vert1_valid->get_id() > vert2_valid->get_id()) std::swap(vert1_valid, vert2_valid);

    // store
    vert1_ = vert1_valid;
    vert2_ = vert2_valid;
    storage_ = storage_valid;
    
    // register edge with verts
    vert1_valid->connect_edge(shared_from_this());
    vert2_valid->connect_edge(shared_from_this());

    // log
    std::cout << "Edge " << id_ << " created between vertex " << vert1_valid->get_id() << " and vertex " << vert2_valid->get_id() << std::endl;
}

void Edge::delete_()
{
    // log
    std::cout << "Destroying edge " << id_ << std::endl;

    // disconnect verts
    if (!vert1_.expired()) vert1_.lock()->disconnect_edge(shared_from_this());
    if (!vert2_.expired()) vert2_.lock()->disconnect_edge(shared_from_this());
    
    // cascade delete faces
    for (auto &face : faces_)
    {
        if (!face.expired()) face.lock()->cascade_delete_from_edge(shared_from_this());
    }

    // log
    std::cout << "Edge " << id_ << " destroyed" << std::endl;
}

void Edge::cascade_delete_from_vert(std::weak_ptr<Vert> vert)
{
    // get valid pointers
    if (vert.expired()) throw std::runtime_error("Try to cascade delete from expired verts");
    if (vert1_.expired()) throw std::runtime_error("Edge holds pointer to deleted vert1");
    if (vert2_.expired()) throw std::runtime_error("Edge holds pointer to deleted vert2");
    if (storage_.expired()) throw std::runtime_error("Edge holds pointer to deleted storage");
    auto vert_valid = vert.lock();
    auto vert1_valid = vert1_.lock();
    auto vert2_valid = vert2_.lock();
    auto storage_valid = storage_.lock();

    // check if vert is connected to edge
    if (vert_valid != vert1_valid && vert_valid != vert2_valid) throw std::runtime_error("Try to cascade delete from unconnected vert");

    // delete edge
    storage_valid->delete_edge(shared_from_this());
}

void Edge::connect_face(std::weak_ptr<Face> face) 
{
    // check if face is valid
    if (face.expired()) throw std::runtime_error("Attempts to connect edge with invalid face.");
    auto face_valid = face.lock();

    // store
    faces_.push_back(face_valid);
}

void Edge::disconnect_face(std::weak_ptr<Face> face)
{
    // check if face is valid
    if (face.expired()) throw std::runtime_error("Attempts to disconnect edge from invalid face.");
    auto face_valid = face.lock();

    // remove face
    faces_.erase(std::remove_if(faces_.begin(), faces_.end(), [&](const std::weak_ptr<Face> &f){ return f.lock() == face_valid; }), faces_.end());

    // check if edge should be deleted
    check_self_destruction();
}

void Edge::check_self_destruction()
{
    // get valid pointer
    if (storage_.expired()) throw std::runtime_error("Edge holds pointer to deleted storage.");
    auto storage_valid = storage_.lock();

    // check self destruction
    if (faces_.empty())
    {
        storage_valid->delete_edge(shared_from_this());
    }
}

std::weak_ptr<Vert> Edge::get_vert1() const 
{
    // get valid pointer
    if (vert1_.expired()) throw std::runtime_error("Edge holds pointer to deleted vert.");
    auto vert1_valid = vert1_.lock();

    // return
    return vert1_valid; 
}
std::weak_ptr<Vert> Edge::get_vert2() const 
{ 
    // get valid pointer
    if (vert2_.expired()) throw std::runtime_error("Edge holds pointer to deleted vert.");
    auto vert2_valid = vert2_.lock();

    // return
    return vert2_valid; 
}

int Edge::get_id() const 
{ 
    return id_; 
}