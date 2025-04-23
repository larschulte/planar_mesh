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
#include <mutex>
#include <set>
#include "MeshObject/Settings.hpp"
#include <shared_mutex>

class Vertex;
class Surface;
class GenericPoint;

struct RRSBoundingBox 
{
    // // include a read write lock
    // mutable std::shared_mutex mutex_; // Read-write lock

    // copy assigmnet operator
    RRSBoundingBox& operator=(const RRSBoundingBox& other);
    RRSBoundingBox(const RRSBoundingBox& other);

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
    double get_surface_area();
};


enum class RRSReturnType
{
    INTERSECTED,
    SKIP
};

class RRSNode : public std::enable_shared_from_this<RRSNode>
{
    static std::atomic<unsigned int> counter_;

public:
    unsigned int id_;
    RRSBoundingBox box_;
    std::weak_ptr<RRSNode> parent_;
    std::shared_ptr<RRSNode> left_;
    std::shared_ptr<RRSNode> right_;
    std::atomic<bool> isLeaf_ = true;
    std::weak_ptr<Vertex> boundary_vertex_;

    bool locked_children = false;

    // read write lock
    mutable std::shared_mutex rwlock_node_;
    mutable std::shared_mutex rwlock_box_;

    RRSNode();
    void recursive_expand_parent_box();
    void recursive_shrink_parent_box();

    void node_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    void node_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    void node_delete_self();
    void node_delete_child(const std::shared_ptr<RRSNode> node);
    void node_update_vertex_box(const std::shared_ptr<Vertex>& boundary_vertex);
    RRSReturnType node_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results);
    void node_print(int level) const;
    void node_flattern(std::vector<std::shared_ptr<Vertex>>& flatten_list);
    void node_count(unsigned int& count) const;
};

struct RRSNodeHasher {
    std::size_t operator()(const std::shared_ptr<RRSNode>& node) const {
        return std::hash<unsigned int>()(node->id_);
    }
};

struct RRSNodeEqual {
    bool operator()(const std::shared_ptr<RRSNode>& a, const std::shared_ptr<RRSNode>& b) const {
        return a->id_ == b->id_;
    }
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
    unsigned int get_node_size() const;
    
    void tree_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    void tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex);
    void tree_update_vertex_box(const std::shared_ptr<Vertex>& boundary_vertex);
    RRSReturnType tree_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>> &search_results);
    void tree_print() const;
};
