#include <iostream>

#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"

int main()
{
    std::shared_ptr storage = std::make_shared<Storage>();

    std::weak_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(1, 1, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(3, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex4 = storage->add_vertex(Eigen::Vector3d(4, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex5 = storage->add_vertex(Eigen::Vector3d(5, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex6 = storage->add_vertex(Eigen::Vector3d(6, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::weak_ptr<Vertex> vertex7 = storage->add_vertex(Eigen::Vector3d(7, 0, 0), Eigen::Vector3d(0, 0, 0), 2);

    storage->print();

    std::set<std::weak_ptr<Vertex>> search_results = storage->reverse_radius_search(Eigen::Vector3d(2.1, 0, 0));
    for (auto vertex : search_results)
    {
        std::cout << vertex.lock()->get_position().transpose() << std::endl;
    }
    

    vertex0.lock()->set_reverse_radius_search_radius(3);
    std::cout << "After adjusting radius" << std::endl;

    storage->print();

    search_results = storage->reverse_radius_search(Eigen::Vector3d(2.1, 0, 0));
    for (auto vertex : search_results)
    {
        std::cout << vertex.lock()->get_position().transpose() << std::endl;
    }


    storage->delete_vertex(vertex4);
    std::cout << "After deleting" << std::endl;
    storage->print();


    storage->rebuild_tree();
    std::cout << "After rebuilding" << std::endl;
    storage->print();
    
    return 0;
}