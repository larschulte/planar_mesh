#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

void Edge::initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create edge with invalid storage.");
    if (vertex1.expired()) throw std::runtime_error("Attempts to create edge with invalid vertex1.");
    if (vertex2.expired()) throw std::runtime_error("Attempts to create edge with invalid vertex2.");
    auto storage_valid = storage.lock();
    auto vertex1_valid = vertex1.lock();
    auto vertex2_valid = vertex2.lock();

    // get id
    id_ = storage_valid->get_next_edge_id();

    // sort
    if (vertex1_valid->get_id() > vertex2_valid->get_id()) std::swap(vertex1_valid, vertex2_valid);

    // store
    vertex1_ = vertex1_valid;
    vertex2_ = vertex2_valid;
    storage_ = storage_valid;
    
    // register edge with verts
    vertex1_valid->connect_edge(shared_from_this());
    vertex2_valid->connect_edge(shared_from_this());

    // log
    std::cout << "Edge " << id_ << " created between vertex " << vertex1_valid->get_id() << " and vertex " << vertex2_valid->get_id() << std::endl;
}

void Edge::delete_()
{
    // log
    std::cout << "Destroying edge " << id_ << std::endl;

    // disconnect verts
    if (!vertex1_.expired()) vertex1_.lock()->disconnect_edge(shared_from_this());
    if (!vertex2_.expired()) vertex2_.lock()->disconnect_edge(shared_from_this());
    
    // cascade delete faces
    for (auto &face : faces_)
    {
        if (!face.expired()) face.lock()->cascade_delete_from_edge(shared_from_this());
    }

    // log
    std::cout << "Edge " << id_ << " destroyed" << std::endl;
}

void Edge::cascade_delete_from_vertex(std::weak_ptr<Vertex> vertex)
{
    // get valid pointers
    if (vertex.expired()) throw std::runtime_error("Try to cascade delete from expired verts");
    if (vertex1_.expired()) throw std::runtime_error("Edge holds pointer to deleted vertex1");
    if (vertex2_.expired()) throw std::runtime_error("Edge holds pointer to deleted vertex2");
    if (storage_.expired()) throw std::runtime_error("Edge holds pointer to deleted storage");
    auto vertex_valid = vertex.lock();
    auto vertex1_valid = vertex1_.lock();
    auto vertex2_valid = vertex2_.lock();
    auto storage_valid = storage_.lock();

    // check if vertex is connected to edge
    if (vertex_valid != vertex1_valid && vertex_valid != vertex2_valid) throw std::runtime_error("Try to cascade delete from unconnected vertex");

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

std::weak_ptr<Vertex> Edge::get_vertex1() const 
{
    // get valid pointer
    if (vertex1_.expired()) throw std::runtime_error("Edge holds pointer to deleted vertex.");
    auto vertex1_valid = vertex1_.lock();

    // return
    return vertex1_valid; 
}
std::weak_ptr<Vertex> Edge::get_vertex2() const 
{ 
    // get valid pointer
    if (vertex2_.expired()) throw std::runtime_error("Edge holds pointer to deleted vertex.");
    auto vertex2_valid = vertex2_.lock();

    // return
    return vertex2_valid; 
}

int Edge::get_id() const 
{ 
    return id_; 
}