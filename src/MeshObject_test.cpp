#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

#include <iostream>
#include <memory>


int main()
{
    std::shared_ptr storage = std::make_shared<Storage>();

    // Create mesh
    std::weak_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, 0));
    std::weak_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(1, 0, 0));
    std::weak_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 1, 0));
    std::weak_ptr<Edge> edge0 = storage->add_edge(vertex0, vertex1);
    std::weak_ptr<Edge> edge1 = storage->add_edge(vertex1, vertex2);
    std::weak_ptr<Edge> edge2 = storage->add_edge(vertex2, vertex0);
    std::weak_ptr<Face> face0 = storage->add_face(edge0, edge1, edge2);

    // Delete to test cascade deletions
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