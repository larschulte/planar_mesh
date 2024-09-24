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

#include <set>

class Face;
class Vertex;
class Surface;
class GenericPoint;

bool ray_triangle_intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& outIntersection);

struct BoundingBox 
{
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    BoundingBox();
    BoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max);
    bool expand(const Eigen::Vector3d& point);
    void expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
    void expand_box_no_return(const BoundingBox& box);
    bool expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
    bool expand_box(const BoundingBox& box);
    bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const;
    bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) const;
    int get_longest_axis();
};

struct Node 
{
    BoundingBox box;
    double split_value;
    int split_axis;
    std::shared_ptr<Node> parent;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
    std::shared_ptr<Node> sibling;
    std::atomic<bool> isLeaf = true;
    std::vector<std::shared_ptr<Face>> faces;

    custom::custom_lock custom_lock;
    omp_nest_lock_t omp_lock;

    bool locked_children = false;

    void recursive_unlock();
    void recursive_expand_parent_box();
    void recursive_shrink_parent_box();
};

enum class BVHReturnType
{
    INTERSECTED,
    SKIP,
    ABORT
};

// overload <<
std::ostream& operator<<(std::ostream& os, const BVHReturnType& type);

class TriangleBVH
{

private:
    std::shared_ptr<Node> root;
    double rebuild_threshold;
    int size_at_last_rebuild;
    int face_size;
    unsigned int leaf_size;

    double sort_face_list_in_axis(std::vector<std::shared_ptr<Face>>& face_list, int axis, int start, int mid, int end);
    
    std::shared_ptr<Node> build_node(const std::vector<std::shared_ptr<Face>>& face_list, const int& start, const int& end);
    void convert_leaf_to_branch(const std::shared_ptr<Node>& node);

    BVHReturnType node_intersection_search(const std::shared_ptr<Node>& node, const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const;
    BVHReturnType node_find_leaf_node(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& endPoint, std::shared_ptr<Node>& return_node);
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
    BVHReturnType tree_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const;
    BVHReturnType tree_find_leaf_node(const Eigen::Vector3d& origin, const Eigen::Vector3d& endPoint, std::shared_ptr<Node>& return_node);
    void tree_print() const;

    unsigned int get_size() const;
};