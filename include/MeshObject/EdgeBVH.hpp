#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <unordered_set>
#include <memory>
#include <array>

#include <shared_mutex>
#include <mutex>

// forward declarations
class Vertex;
class Edge;
class Surface;

class EdgeBVH
{
public:
    struct BoundingBox 
    {
        Eigen::Vector3d min;
        Eigen::Vector3d max;
        Eigen::Vector3d min_used_for_surface_area;
        Eigen::Vector3d max_used_for_surface_area;
        double surface_area;
        BoundingBox();
        BoundingBox(const std::shared_ptr<Edge>& edge);
        void expand(const Eigen::Vector3d& point);
        void expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
        void expand_box_no_return(const std::shared_ptr<Edge>& edge);
        void expand_box_no_return(const BoundingBox& box);
        bool expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max);
        bool expand_box(const EdgeBVH::BoundingBox& box);
        bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir_inv, double& tMin, double& tMax) const;
        bool intersect(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const;
        int get_longest_axis();
        double compute_surface_area() const;
        const double& get_surface_area();
    };

    class Node : public std::enable_shared_from_this<EdgeBVH::Node>
    {
        public:
        BoundingBox box_;
        double split_value_;
        int split_axis_;
        std::shared_ptr<EdgeBVH::Node> parent_;
        std::shared_ptr<EdgeBVH::Node> left_;
        std::shared_ptr<EdgeBVH::Node> right_;
        bool isLeaf() const;
        std::vector<std::shared_ptr<Edge>> edges_;

        void recursive_expand_parent_box();
        void recursive_shrink_parent_box();

        // read write lock
        mutable std::shared_mutex rwlock_node_;

        void node_add_edge(const std::shared_ptr<Edge>& edge);
        void node_delete_edge(const std::shared_ptr<Edge>& edge);
        bool node_intersect_edge(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
        void node_print(int level) const;    
        void node_flatten(std::vector<std::shared_ptr<Edge>>& flat_vector) const;
    };

private:
    std::shared_ptr<EdgeBVH::Node> root_;
    int edge_size_;
    
    // use SAH for EdgeBVH
    std::shared_ptr<EdgeBVH::Node> find_best_node(const std::shared_ptr<EdgeBVH::Node>& root, const std::shared_ptr<Edge>& edge);

public:
    EdgeBVH();
    std::vector<std::shared_ptr<Edge>> get_edge_list() const;

    void tree_add_edge(const std::shared_ptr<Edge>& edge);
    void tree_delete_edge(const std::shared_ptr<Edge>& edge);
    bool tree_intersect_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1);
    void tree_print() const;
};