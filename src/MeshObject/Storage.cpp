#include "MeshObject/Storage.hpp"
#include "MeshObject/Vert.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

void Storage::add_vert(Eigen::Vector3d pos) 
{
    // initialize
    std::shared_ptr<Vert> vert = std::make_shared<Vert>();
    vert->initialize_(shared_from_this(), pos);

    // put into list
    verts_.push_back(vert);    
}

void Storage::add_edge(std::weak_ptr<Vert> vert1, std::weak_ptr<Vert> vert2) 
{    
    // initialize
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), vert1, vert2);

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


// need to ensure the vert/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vert(std::weak_ptr<Vert> vert) 
{
    // check if valid
    if (vert.expired()) throw std::runtime_error("Attempts to delete expired vert.");
    auto vert_valid = vert.lock();

    // member delete
    vert_valid->delete_();

    // storage delete
    verts_.erase(std::remove(verts_.begin(), verts_.end(), vert_valid), verts_.end());
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

// get id
int Storage::get_next_vert_id() { return next_vert_id_++; }
int Storage::get_next_edge_id() { return next_edge_id_++; }
int Storage::get_next_face_id() { return next_face_id_++; }