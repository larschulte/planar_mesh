#include "MeshObject/RRSTree.hpp"
#include "MeshObject/Storage.hpp"

void test1()
{
    RRSTree rrs_tree;

    std::shared_ptr<Storage> storage = std::make_shared<Storage>();
    std::shared_ptr<Vertex> vertex0 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 2);
    std::shared_ptr<Vertex> vertex1 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 0), 2);
    std::shared_ptr<Vertex> vertex2 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2, 0, 0), 2);
    std::shared_ptr<Vertex> vertex3 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(3, 0, 0), 2);
    std::shared_ptr<Vertex> vertex4 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(4, 0, 0), 2);
    std::shared_ptr<Vertex> vertex5 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(5, 0, 0), 2);
    std::shared_ptr<Vertex> vertex6 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(6, 0, 0), 2);
    std::shared_ptr<Vertex> vertex7 = storage->add_vertex(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(7, 0, 0), 2);


    rrs_tree.tree_add_vertex(vertex0);
    rrs_tree.tree_add_vertex(vertex1);
    rrs_tree.tree_add_vertex(vertex2);
    rrs_tree.tree_add_vertex(vertex3);
    rrs_tree.tree_add_vertex(vertex4);
    rrs_tree.tree_add_vertex(vertex5);
    rrs_tree.tree_add_vertex(vertex6);

    rrs_tree.tree_print();   
    rrs_tree.rebuild();
    rrs_tree.tree_print();   

    rrs_tree.tree_delete_vertex(vertex3);
    rrs_tree.tree_print(); 
}

int main()
{
    test1();
    return 0;
}