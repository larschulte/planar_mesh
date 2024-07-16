#include "MeshObject/Surface.hpp"
#include "MeshObject/Storage.hpp"

void test1()
{
    std::shared_ptr<Storage> storage = std::make_shared<Storage>();
    std::shared_ptr<Surface> surface0 = storage->add_surface();
    std::shared_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 2, 3));
    std::shared_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2, 3, 4));

    surface0->connect(vertex1);

    surface0->print_info();
    surface0->connect(vertex0);
    surface0->print_info();
    surface0->disconnect(vertex0);
    surface0->print_info();
}

int main()
{
    test1();
    return 0;
}