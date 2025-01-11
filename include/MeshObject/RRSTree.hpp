#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <unordered_set>
#include "MeshObject/MeshObject.hpp"
#include <memory>
#include <array>
#include <iostream> 

#include <shared_mutex>
#include <set>
#include "MeshObject/Settings.hpp"

class Vertex;
class Surface;
class GenericPoint;

struct RRSBoundingBox 
{
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    Eigen::Vector3d min_used_for_surface_area;
    Eigen::Vector3d max_used_for_surface_area;
    double surface_area;
    RRSBoundingBox();
    RRSBoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max);
    bool expand(const Eigen::Vector3d& point);
    void expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
    void expand_box_no_return(const RRSBoundingBox& box);
    bool expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
    bool expand_box(const RRSBoundingBox& box);
    bool contains(const Eigen::Vector3d& point);
    int get_longest_axis();
    double compute_surface_area() const;
    const double& get_surface_area();
};


enum class RRSReturnType
{
    INTERSECTED,
    SKIP,
    ABORT
};

class RRSNode : public std::enable_shared_from_this<RRSNode>
{
    public:
    RRSBoundingBox box;
    double split_value;
    int split_axis;
    std::shared_ptr<RRSNode> parent;
    std::shared_ptr<RRSNode> left;
    std::shared_ptr<RRSNode> right;
    std::shared_ptr<RRSNode> sibling;
    std::atomic<bool> isLeaf = true;
    std::vector<std::shared_ptr<Vertex>> boundary_vertices;

    bool locked_children = false;

    void recursive_expand_parent_box();
    void recursive_shrink_parent_box();

    void node_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    bool node_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    RRSReturnType node_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results);
    void node_print(int level) const;
    void node_flattern(std::vector<std::shared_ptr<Vertex>>& flatten_list);
};


// Overload the << operator for RRSReturnType
std::ostream& operator<<(std::ostream& os, const RRSReturnType& type);

class RRSTree
{    
private:
    static Settings settings_;

    std::shared_ptr<RRSNode> root;
    int tree_size;
    std::size_t leaf_size;
    std::shared_ptr<RRSNode> find_best_node(const std::shared_ptr<RRSNode>& root, const std::shared_ptr<Vertex>& boundary_vertex);

public:
    RRSTree();
    std::vector<std::shared_ptr<Vertex>> compute_vertices_list();
    void print_size();
    unsigned int get_size() const;
    
    void tree_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    void tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    RRSReturnType tree_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>> &search_results);
    void tree_print() const;
};
