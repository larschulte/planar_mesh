#include "MeshObject/UnionFind.hpp"
#include "MeshObject/Storage.hpp"

void test1()
{
    std::cout << "============================== Test 1 ==============================\n";

    std::shared_ptr<Storage> storage = std::make_shared<Storage>();

    std::shared_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(0, 0, 0));
    std::shared_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 0, 0));
    std::shared_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(1, 1, 0));
    std::shared_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(2, 1, 0));
    std::shared_ptr<Vertex> vertex4 = storage->add_vertex(Eigen::Vector3d(0, 0, -1), Eigen::Vector3d(2, 2, 0));

    std::shared_ptr<Edge> edge0 = storage->add_edge(vertex0, vertex1);
    std::shared_ptr<Edge> edge1 = storage->add_edge(vertex0, vertex2);
    std::shared_ptr<Edge> edge2 = storage->add_edge(vertex1, vertex2);
    
    std::shared_ptr<Edge> edge3 = storage->add_edge(vertex2, vertex3);
    std::shared_ptr<Edge> edge4 = storage->add_edge(vertex2, vertex4);
    std::shared_ptr<Edge> edge5 = storage->add_edge(vertex3, vertex4);

    UnionFind uf;
    uf.add_vertices(storage->get_vertices());
    uf.add_edges(storage->get_edges());
    uf.print_sorted_grouped_vertices();

    storage->delete_vertex(vertex2);
    UnionFind uf2;
    uf2.add_vertices(storage->get_vertices());
    uf2.add_edges(storage->get_edges());
    uf2.print_sorted_grouped_vertices();
}

int main()
{
    test1();
    return 0;
}