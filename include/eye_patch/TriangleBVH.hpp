#ifndef TRIANGLEBVH_H
#define TRIANGLEBVH_H

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <array>

class TriangleBVH
{
private:
    struct BoundingBox 
    {
        Eigen::Vector3d min;
        Eigen::Vector3d max;

        BoundingBox();

        void expand(const Eigen::Vector3d& point);

        bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const;

        bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) const;

        int get_longest_axis();
    };

    struct Node 
    {
        BoundingBox box;
        double split_value;
        int split_axis;
        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;
        bool isLeaf() const;
        std::set<int> triangleIDs;
    };

    Eigen::Vector3d compute_triangle_center(int triangleID);
    double sort_triangle_list_in_axis(std::vector<int>& triangle_list, int axis, int start, int mid, int end);
    void expand_node_box(std::shared_ptr<Node> node, int triangle_id);
    bool rayTriangleIntersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& outIntersection);
    std::set<int> intersectHierarchy(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir);
    void convert_leaf_to_branch(std::shared_ptr<Node> node);
    std::shared_ptr<Node> build_node(std::vector<int> triangle_list);
    void addTriangleToNode(std::shared_ptr<Node> node, int triangleID);
    void deleteTriangleFromNode(std::shared_ptr<Node> node, int triangleID);

    double rebuild_threshold;
    std::vector<int> triangle_list;
    std::map<int, std::array<int, 3>> triangle_to_indices_map;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, Eigen::Vector3d> triangle_to_center_vector3d_map;
    std::shared_ptr<Node> root;
    int size_at_last_rebuild;

public:
    TriangleBVH();
    void addData(std::vector<int> _triangle_list, std::map<int, std::array<int, 3>> _triangle_to_indices_map, std::map<int, Eigen::Vector3d> _point_to_vector3d_map);
    void rebuild();
    void addTriangle(int triangleID, std::array<int, 3> indices, Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2);
    void deleteTriangle(int triangleID);
    std::set<int> intersectionSearch(Eigen::Vector3d origin, Eigen::Vector3d endPoint);
};

#endif // TRIANGLEBVH_H
