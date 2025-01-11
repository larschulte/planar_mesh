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

struct RRSNode 
{
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
};

enum class RRSReturnType
{
    INTERSECTED,
    SKIP,
    ABORT
};

// Overload the << operator for RRSReturnType
std::ostream& operator<<(std::ostream& os, const RRSReturnType& type);

class RRSTree
{    
private:
    static Settings settings_;

    std::shared_ptr<RRSNode> root;
    double rebuild_threshold;
    int size_at_last_rebuild;
    int tree_size;
    std::size_t leaf_size;

    double sort_boundary_vertex_list_in_axis(std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, int axis, int start, int mid, int end);
    void sort_boundary_vertex_list_in_axis(std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, int axis, int start, int end);
    // void expand_node_box(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex);

    std::shared_ptr<RRSNode> build_node(const std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, const int& start, const int& end);
    void convert_leaf_to_branch(const std::shared_ptr<RRSNode>& node);

    std::shared_ptr<RRSNode> find_best_node(const std::shared_ptr<RRSNode>& root, const std::shared_ptr<Vertex>& boundary_vertex);
    void node_add_vertex(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex);
    double calculate_sah(RRSBoundingBox& parent_box, RRSBoundingBox& left_box, RRSBoundingBox& right_box, int left_count, int right_count);
    void node_increase_radius(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex);
    bool node_delete_vertex(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex);
    RRSReturnType node_reverse_radius_search(RRSNode* node, const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results);
    void node_print(const std::shared_ptr<RRSNode>& node, int level) const;
    void node_flattern(const std::shared_ptr<RRSNode>& node, std::vector<std::shared_ptr<Vertex>>& flatten_list);

public:
    RRSTree();
    void check_rebuild();
    void rebuild();
    bool can_reverse_radius_search();
    std::vector<std::shared_ptr<Vertex>> compute_vertices_list();
    void print_size();
    unsigned int get_size() const;
    
    void tree_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    void tree_increase_radius(std::shared_ptr<Vertex> boundary_vertex);
    void tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    RRSReturnType tree_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>> &search_results);
    void tree_print() const;
};
