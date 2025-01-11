#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <unordered_set>
#include "MeshObject/MeshObject.hpp"
#include <memory>
#include <array>

#include <shared_mutex>
#include "MeshObject/Settings.hpp"

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
    Eigen::Vector3d min_used_for_surface_area;
    Eigen::Vector3d max_used_for_surface_area;
    double surface_area;
    BoundingBox();
    BoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max);
    bool expand(const Eigen::Vector3d& point);
    void expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
    void expand_box_no_return(const BoundingBox& box);
    bool expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
    bool expand_box(const BoundingBox& box);
    bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& invdir, double& tMin, double& tMax) const;
    bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& invdir) const;
    int get_longest_axis();
    double compute_surface_area() const;
    const double& get_surface_area();
};

enum class BVHReturnType
{
    INTERSECTED,
    SKIP,
    ABORT
};
class Node : public std::enable_shared_from_this<Node>
{
    public:
    BoundingBox box;
    double split_value;
    int split_axis;
    std::shared_ptr<Node> parent;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
    std::shared_ptr<Node> sibling;
    std::atomic<bool> isLeaf = true;
    std::vector<std::shared_ptr<Face>> faces;

    bool locked_children = false;

    void recursive_expand_parent_box();
    void recursive_shrink_parent_box();

    BVHReturnType node_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const;
    void node_add_face(const std::shared_ptr<Face>& face);
    bool node_delete_face(const std::shared_ptr<Face>& face);
    void node_print(int level) const;
    void node_flatten(std::vector<std::shared_ptr<Face>>& face_list) const;
};

// overload <<
std::ostream& operator<<(std::ostream& os, const BVHReturnType& type);

class TriangleBVH
{

private:
    static Settings settings_;
    int face_size;
    unsigned int leaf_size;
    std::shared_ptr<Node> root;
    std::shared_ptr<Node> find_best_node(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face);

public:
    TriangleBVH();
    std::vector<std::shared_ptr<Face>> get_face_list() const;

    void tree_add_face(std::shared_ptr<Face> face);
    void tree_delete_face(std::shared_ptr<Face> face);
    BVHReturnType tree_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const;
    void tree_print() const;

    unsigned int get_size() const;
};