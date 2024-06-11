#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <set>
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
        std::set<std::weak_ptr<Edge>> edges;
    };

    double sort_edge_list_in_axis(std::vector<std::weak_ptr<Edge>>& edge_list, int axis, int start, int mid, int end);
    void expand_node_box(std::shared_ptr<Node> node, std::weak_ptr<Edge> edge);
    
    std::shared_ptr<EdgeBVH::Node> build_node(std::vector<std::weak_ptr<Edge>> edge_list);
    void convert_leaf_to_branch(std::shared_ptr<Node> node);

    void node_add_edge(std::shared_ptr<Node> node, std::weak_ptr<Edge> edge);
    void node_delete_edge(std::shared_ptr<Node> node, std::weak_ptr<Edge> edge);
    bool node_intersect_edge(const std::shared_ptr<Node>& node, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2);
    void node_print(const std::shared_ptr<Node>& node, int level) const;    
    void node_flatten(std::shared_ptr<Node> node, std::vector<std::weak_ptr<Edge>>& flat_vector) const;

    std::vector<std::weak_ptr<Edge>> get_edge_list() const;
    
    std::shared_ptr<Node> root;

    double rebuild_threshold;
    int size_at_last_rebuild;

    int edge_size;
    
    std::set<std::weak_ptr<Edge>> edge_set;

public:
    EdgeBVH();
    void rebuild();

    void add_edge(std::weak_ptr<Edge> edge);
    void delete_edge(std::weak_ptr<Edge> edge);
    bool intersect_edges(std::weak_ptr<Vertex> vertex0, std::weak_ptr<Vertex> vertex1);
    void print() const;
};