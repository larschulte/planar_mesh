#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

#include <iostream>

int main()
{
    std::shared_ptr storage = std::make_shared<Storage>();

    // Create mesh
    std::shared_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0.9));
    std::shared_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 0, 1));
    std::shared_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 1, 1));
    std::shared_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(-1, 0, 1));
    std::shared_ptr<Vertex> vertex4 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, -1, 1));
    std::shared_ptr<Edge> edge0 = storage->add_edge(vertex0, vertex1);
    std::shared_ptr<Edge> edge1 = storage->add_edge(vertex0, vertex2);
    std::shared_ptr<Edge> edge2 = storage->add_edge(vertex0, vertex3);
    std::shared_ptr<Edge> edge3 = storage->add_edge(vertex0, vertex4);
    std::shared_ptr<Edge> edge4 = storage->add_edge(vertex1, vertex2);
    std::shared_ptr<Edge> edge5 = storage->add_edge(vertex2, vertex3);
    std::shared_ptr<Edge> edge6 = storage->add_edge(vertex3, vertex4);
    std::shared_ptr<Edge> edge7 = storage->add_edge(vertex4, vertex1);
    std::shared_ptr<Face> face0 = storage->add_face(vertex0, vertex1, vertex2);
    std::shared_ptr<Face> face1 = storage->add_face(vertex0, vertex2, vertex3);
    std::shared_ptr<Face> face2 = storage->add_face(vertex0, vertex3, vertex4);
    std::shared_ptr<Face> face3 = storage->add_face(vertex0, vertex4, vertex1);

    std::shared_ptr<Edge> edge8 = storage->add_edge(vertex1, vertex3);
    std::shared_ptr<Face> face4 = storage->add_face(vertex1, vertex2, vertex3);

    std::cout << "\nPrint tree" << std::endl;
    storage->print_bvh();

    std::set<std::shared_ptr<Face>> searched_results = storage->face_intersection_search(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0.95));
    std::cout << "\nPrint searched results" << std::endl;
    for (auto face : searched_results)
    {
        std::cout << "Intersects face "<< face->get_id();
        std::cout << " (vertex " << face->get_vertex(0)->get_id();
        std::cout << ", vertex " << face->get_vertex(1)->get_id();
        std::cout << ", vertex " << face->get_vertex(2)->get_id() << ")" << std::endl;
    }
    std::cout << "\n" << std::endl;


    storage->delete_vertex(vertex1);
    std::cout << "\nPrint tree after delete" << std::endl;
    storage->print_bvh();

    searched_results = storage->face_intersection_search(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0.95));
    std::cout << "\nPrint searched results" << std::endl;
    for (auto face : searched_results)
    {
        std::cout << "Intersects face "<< face->get_id();
        std::cout << " (vertex " << face->get_vertex(0)->get_id();
        std::cout << ", vertex " << face->get_vertex(1)->get_id();
        std::cout << ", vertex " << face->get_vertex(2)->get_id() << ")" << std::endl;
    }
}