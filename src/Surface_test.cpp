#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Surface.hpp"

int main()
{
    std::shared_ptr<Storage> storage = std::make_shared<Storage>();

    std::weak_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 0, 0));
    std::weak_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 1, 0));
    std::weak_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 0, 0));
    std::weak_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 1, 0));

    std::weak_ptr<Surface> surface = storage->add_surface();

    std::cout << "at start, normal = [" << surface.lock()->get_normal().transpose() << "]" << std::endl;
    surface.lock()->connect(vertex0);
    std::cout << "after vertex 0, normal = [" << surface.lock()->get_normal().transpose() << "]" << std::endl;
    surface.lock()->connect(vertex1);
    std::cout << "after vertex 1, normal = [" << surface.lock()->get_normal().transpose() << "]" << std::endl;
    surface.lock()->connect(vertex2);
    std::cout << "after vertex 2, normal = [" << surface.lock()->get_normal().transpose() << "]" << std::endl;
    // surface.lock()->connect(vertex3);
    // std::cout << "after vertex 3, normal = [" << surface.lock()->get_normal().transpose() << "]" << std::endl;

    double distance = surface.lock()->compute_point_to_surface_distance_with_improved_covariance(Eigen::Vector3d(0, 0, -10), Eigen::Vector3d(0, 0, -0.5));
    std::cout << "distance = " << distance << std::endl;

    return 0;
}