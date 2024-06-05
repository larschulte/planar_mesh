#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

#include <iostream>
#include <memory>


int main()
{
    std::shared_ptr storage = std::make_shared<Storage>();

    // Creating verts
    storage->add_vertex(Eigen::Vector3d(0, 0, 0));
    storage->add_vertex(Eigen::Vector3d(1, 0, 0));
    storage->add_vertex(Eigen::Vector3d(0, 1, 0));

    // Get weak_ptr to verts from storage
    std::weak_ptr<Vertex> vertex0 = storage->vertices_[0];
    std::weak_ptr<Vertex> vertex1 = storage->vertices_[1];
    std::weak_ptr<Vertex> vertex2 = storage->vertices_[2];

    // Creating edges
    storage->add_edge(vertex0, vertex1);
    storage->add_edge(vertex1, vertex2);
    storage->add_edge(vertex2, vertex0);

    // Get weak_ptr to edges from storage
    std::weak_ptr<Edge> edge0 = storage->edges_[0];
    std::weak_ptr<Edge> edge1 = storage->edges_[1];
    std::weak_ptr<Edge> edge2 = storage->edges_[2];

    // Creating faces
    storage->add_face(edge0, edge1, edge2);

    // Get weak_ptr to faces from storage
    std::weak_ptr<Face> face0 = storage->faces_[0];

    // Deleting a vertex to see the cascading deletions
    storage->delete_vertex(vertex0.lock());
    // storage->delete_face(face0.lock());
    // storage->delete_edge(edge0.lock());

    // Check if weak pointers are expired
    std::cout << "Vertex 0 expired: " << vertex0.expired() << std::endl;
    std::cout << "Edge 0 expired: " << edge0.expired() << std::endl;
    std::cout << "Edge 1 expired: " << edge1.expired() << std::endl;
    std::cout << "Edge 2 expired: " << edge2.expired() << std::endl;
    std::cout << "Face 0 expired: " << face0.expired() << std::endl;
    
    return 0;
}