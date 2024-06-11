#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Surface.hpp"

// int main()
// {
//     std::shared_ptr<Storage> storage = std::make_shared<Storage>();

//     std::shared_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 0, 0));
//     std::shared_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 1, 0));
//     std::shared_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 0, 0));
//     std::shared_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 1, 0));

//     std::shared_ptr<Surface> surface = storage->add_surface();

//     std::cout << "at start, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
//     surface->connect(vertex0);
//     std::cout << "after vertex 0, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
//     surface->connect(vertex1);
//     std::cout << "after vertex 1, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
//     surface->connect(vertex2);
//     std::cout << "after vertex 2, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
//     // surface->connect(vertex3);
//     // std::cout << "after vertex 3, normal = [" << surface->get_normal().transpose() << "]" << std::endl;

//     double distance = surface->compute_point_to_surface_distance_with_improved_covariance(Eigen::Vector3d(0, 0, -10), Eigen::Vector3d(0, 0, -0.5));
//     std::cout << "distance = " << distance << std::endl;

//     return 0;
// }

int main()
{
    std::shared_ptr<Storage> storage = std::make_shared<Storage>();

    std::shared_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 0, 0));
    std::shared_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 1, 0));
    std::shared_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 0, 0));
    std::shared_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 1, 0.1));
    std::shared_ptr<Vertex> vertex4 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(2, 1, 0.1));
    std::shared_ptr<Vertex> vertex5 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 2, 0.1));
    std::shared_ptr<Vertex> vertex6 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(2, 2, 0.1));

    std::shared_ptr<Surface> surface = storage->add_surface();

    std::cout << "at start, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex0);
    std::cout << "after vertex 0, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex1);
    std::cout << "after vertex 1, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex2);
    std::cout << "after vertex 2, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex3);
    std::cout << "after vertex 3, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex4);
    std::cout << "after vertex 4, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex5);
    std::cout << "after vertex 5, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    surface->connect(vertex6);
    std::cout << "after vertex 6, normal = [" << surface->get_normal().transpose() << "]" << std::endl;
    
    double distance = surface->compute_point_to_surface_distance_with_improved_covariance(Eigen::Vector3d(0, 0, -10), Eigen::Vector3d(0, 0, 0));
    std::cout << "distance = " << distance << std::endl;

    return 0;
}