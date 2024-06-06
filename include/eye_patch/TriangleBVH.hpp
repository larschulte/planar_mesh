#ifndef TRIANGLEBVH_H
#define TRIANGLEBVH_H

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <array>

#include "MeshObject/Face.hpp"
#include "MeshObject/Vertex.hpp"

bool ray_triangle_intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& outIntersection);

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
        std::set<std::weak_ptr<Face>> faces;
    };

    double sort_face_list_in_axis(std::vector<std::weak_ptr<Face>>& face_list, int axis, int start, int mid, int end);
    void expand_node_box(std::shared_ptr<Node> node, std::weak_ptr<Face> face);
    std::set<std::weak_ptr<Face>> intersectHierarchy(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir);
    void convert_leaf_to_branch(std::shared_ptr<Node> node);
    std::shared_ptr<Node> build_node(std::vector<std::weak_ptr<Face>> face_list);
    void add_face_to_node(std::shared_ptr<Node> node, std::weak_ptr<Face> face);
    void delete_face_from_node(std::shared_ptr<Node> node, std::weak_ptr<Face> face);

    double rebuild_threshold;
    std::vector<std::weak_ptr<Face>> face_list;
    std::map<int, std::array<int, 3>> face_to_indices_map;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, Eigen::Vector3d> face_to_center_vector3d_map;
    std::shared_ptr<Node> root;
    int size_at_last_rebuild;

    std::set<std::weak_ptr<Face>> face_set;

public:
    TriangleBVH();
    void rebuild();
    void add_face(std::weak_ptr<Face> face);
    void delete_face(std::weak_ptr<Face> face);
    std::set<std::weak_ptr<Face>> intersectionSearch(Eigen::Vector3d origin, Eigen::Vector3d endPoint);
};

#endif // TRIANGLEBVH_H
