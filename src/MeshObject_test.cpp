#include "MeshObject/Storage.hpp"
#include "MeshObject/Vert.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

#include <iostream>
#include <memory>


int main()
{
    std::shared_ptr storage = std::make_shared<Storage>();

    // Creating verts
    storage->add_vert(Eigen::Vector3d(0, 0, 0));
    storage->add_vert(Eigen::Vector3d(1, 0, 0));
    storage->add_vert(Eigen::Vector3d(0, 1, 0));

    // Get weak_ptr to verts from storage
    std::weak_ptr<Vert> vert0 = storage->verts_[0];
    std::weak_ptr<Vert> vert1 = storage->verts_[1];
    std::weak_ptr<Vert> vert2 = storage->verts_[2];

    // Creating edges
    storage->add_edge(vert0, vert1);
    storage->add_edge(vert1, vert2);
    storage->add_edge(vert2, vert0);

    // Get weak_ptr to edges from storage
    std::weak_ptr<Edge> edge0 = storage->edges_[0];
    std::weak_ptr<Edge> edge1 = storage->edges_[1];
    std::weak_ptr<Edge> edge2 = storage->edges_[2];

    // Creating faces
    storage->add_face(edge0, edge1, edge2);

    // Get weak_ptr to faces from storage
    std::weak_ptr<Face> face0 = storage->faces_[0];

    // Deleting a vert to see the cascading deletions
    // storage->delete_vert(vert0.lock());
    // storage->delete_face(face0.lock());
    storage->delete_edge(edge0.lock());

    // Check if weak pointers are expired
    std::cout << "Vert 0 expired: " << vert0.expired() << std::endl;
    std::cout << "Edge 0 expired: " << edge0.expired() << std::endl;
    std::cout << "Edge 1 expired: " << edge1.expired() << std::endl;
    std::cout << "Edge 2 expired: " << edge2.expired() << std::endl;
    std::cout << "Face 0 expired: " << face0.expired() << std::endl;
    
    return 0;
}