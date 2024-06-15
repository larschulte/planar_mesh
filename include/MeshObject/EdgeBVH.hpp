#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <unordered_set>
#include <memory>
#include <array>

// forward declarations
class Vertex;
class Edge;

class EdgeBVH
{
private:
    struct BoundingBox 
    {
        Eigen::Vector3d min;
        Eigen::Vector3d max;
        BoundingBox();
        void expand(const Eigen::Vector3d& point);
        bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const;
        bool intersect(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const;
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
        std::vector<std::shared_ptr<Edge>> edges;
    };

private:
    std::shared_ptr<Node> root;
    double rebuild_threshold;
    int size_at_last_rebuild;
    int edge_size;

    double sort_edge_list_in_axis(std::vector<std::shared_ptr<Edge>>& edge_list, int axis, int start, int mid, int end);
    void expand_node_box(const std::shared_ptr<Node>& node, const std::shared_ptr<Edge>& edge);
    
    std::shared_ptr<EdgeBVH::Node> build_node(const std::vector<std::shared_ptr<Edge>>& edge_list, const int& start, const int& end);
    void convert_leaf_to_branch(const std::shared_ptr<Node>& node);

    void node_add_edge(const std::shared_ptr<Node>& node, const std::shared_ptr<Edge>& edge);
    bool node_delete_edge(const std::shared_ptr<Node>& node, const std::shared_ptr<Edge>& edge);
    bool node_intersect_edge(const std::shared_ptr<Node>& node, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    void node_print(const std::shared_ptr<Node>& node, int level) const;    
    void node_flatten(const std::shared_ptr<Node>& node, std::vector<std::shared_ptr<Edge>>& flat_vector) const;

public:
    EdgeBVH();
    void rebuild();
    std::vector<std::shared_ptr<Edge>> get_edge_list() const;

    void tree_add_edge(const std::shared_ptr<Edge>& edge);
    void tree_delete_edge(const std::shared_ptr<Edge>& edge);
    bool tree_intersect_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1);
    void tree_print() const;
};