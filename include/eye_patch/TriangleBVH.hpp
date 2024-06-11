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
        std::vector<std::weak_ptr<Face>> faces;
    };

    double sort_face_list_in_axis(std::vector<std::weak_ptr<Face>>& face_list, int axis, int start, int mid, int end);
    void expand_node_box(std::shared_ptr<Node> node, std::weak_ptr<Face> face);
    
    std::shared_ptr<Node> build_node(std::vector<std::weak_ptr<Face>> face_list);
    void convert_leaf_to_branch(std::shared_ptr<Node> node);

    void node_intersection_search(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, std::set<std::weak_ptr<Face>>& faces_intersected) const;
    void node_add_face(std::shared_ptr<Node> node, std::weak_ptr<Face> face);
    void node_delete_face(std::shared_ptr<Node> node, std::weak_ptr<Face> face);
    void node_print(const std::shared_ptr<Node>& node, int level) const;
    void node_flatten(std::shared_ptr<TriangleBVH::Node> node, std::vector<std::weak_ptr<Face>>& face_list) const;

    std::vector<std::weak_ptr<Face>> get_face_list() const;

    double rebuild_threshold;
    int size_at_last_rebuild;

    
    std::shared_ptr<Node> root;

    std::set<std::weak_ptr<Face>> face_set;

    int face_size;

public:
    TriangleBVH();
    void rebuild();

    void add_face(std::weak_ptr<Face> face);
    void delete_face(std::weak_ptr<Face> face);
    std::set<std::weak_ptr<Face>> intersection_search(Eigen::Vector3d origin, Eigen::Vector3d endPoint);
    void print() const;
};

#endif // TRIANGLEBVH_H
