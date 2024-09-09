#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <unordered_set>
#include "MeshObject/MeshObject.hpp"
#include <memory>
#include <array>

#include "utilities/omp_utilities.hpp"

class Face;
class Vertex;

bool ray_triangle_intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& outIntersection);

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
    std::vector<std::shared_ptr<Face>> faces;

    omp_lock_t lock;
};

};

class TriangleBVH
{

private:
    std::shared_ptr<Node> root;
    double rebuild_threshold;
    int size_at_last_rebuild;
    int face_size;
    unsigned int leaf_size;

    double sort_face_list_in_axis(std::vector<std::shared_ptr<Face>>& face_list, int axis, int start, int mid, int end);
    void expand_node_box(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face);
    
    std::shared_ptr<Node> build_node(const std::vector<std::shared_ptr<Face>>& face_list, const int& start, const int& end);
    void convert_leaf_to_branch(const std::shared_ptr<Node>& node);

    void node_intersection_search(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, std::vector<std::shared_ptr<Face>>& faces_intersected) const;
    void node_add_face(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face);
    bool node_delete_face(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face);
    void node_print(const std::shared_ptr<Node>& node, int level) const;
    void node_flatten(const std::shared_ptr<Node>& node, std::vector<std::shared_ptr<Face>>& face_list) const;

public:
    TriangleBVH();
    void check_rebuild();
    void rebuild();
    std::vector<std::shared_ptr<Face>> get_face_list() const;

    void tree_add_face(std::shared_ptr<Face> face);
    void tree_delete_face(std::shared_ptr<Face> face);
    void tree_intersection_search(Eigen::Vector3d origin, Eigen::Vector3d endPoint, std::vector<std::shared_ptr<Face>>& faces_intersected) const;
    void tree_print() const;
};