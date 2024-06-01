#include "eye_patch/RRSTree.hpp"
#include <iostream>

int main()
{
    RRSTree tree;

    // tree.addBoundaryPoint(0, Eigen::Vector3d(0, 0, 0), 2);
    // // tree.addBoundaryPoint(1, Eigen::Vector3d(1, 0, 0), 2);
    // tree.adjustRadius(0, Eigen::Vector3d(0, 0, 0), 3);

    tree.addBoundaryPoint(0, Eigen::Vector3d(0, 0, 0), 2);
    tree.addBoundaryPoint(1, Eigen::Vector3d(1, 1, 0), 2);
    tree.addBoundaryPoint(2, Eigen::Vector3d(2, 0, 0), 2);
    tree.addBoundaryPoint(3, Eigen::Vector3d(3, 1, 0), 2);
    tree.addBoundaryPoint(4, Eigen::Vector3d(4, 0, 0), 2);
    tree.addBoundaryPoint(5, Eigen::Vector3d(5, 0, 0), 2);
    tree.addBoundaryPoint(6, Eigen::Vector3d(6, 0, 0), 2);
    tree.addBoundaryPoint(7, Eigen::Vector3d(7, 0, 0), 2);

    tree.printTree();

    std::set<int> search_results = tree.reverseRadiusSearch(Eigen::Vector3d(2.1, 0, 0));
    for (auto it = search_results.begin(); it != search_results.end(); ++it)
    {
        std::cout << *it << std::endl;
    }

    tree.adjustRadius(0, Eigen::Vector3d(0, 0, 0), 3);

    std::cout << "After adjusting radius" << std::endl;

    tree.printTree();

    std::set<int> search_results2 = tree.reverseRadiusSearch(Eigen::Vector3d(2.1, 0, 0));
    for (auto it = search_results2.begin(); it != search_results2.end(); ++it)
    {
        std::cout << *it << std::endl;
    }


    tree.deleteBoundaryPoint(4, Eigen::Vector3d(4, 0, 0));

    tree.printTree();

    tree.rebuildTree();

    tree.printTree();
    
    return 0;
}