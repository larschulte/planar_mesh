#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <array>
#include <iostream> 

#include "MeshObject/Vertex.hpp"

class RRSTree
{
private:
    struct BoundingBox 
    {
        Eigen::Vector3d min;
        Eigen::Vector3d max;

        BoundingBox() : min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
                        max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

        void expand(const Eigen::Vector3d& point)
        {
            min = min.cwiseMin(point);
            max = max.cwiseMax(point);
        }

        bool contains(const Eigen::Vector3d& point)
        {
            return (point.array() >= min.array()).all() && (point.array() <= max.array()).all();
        }

        int get_longest_axis()
        {
            Eigen::Vector3d diagonal_line = max - min;
            int axis = 0;
            if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
            if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
            return axis;
        }
    };

    struct Node 
    {
        BoundingBox box;
        double split_value;
        int split_axis;
        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;
        std::vector<std::weak_ptr<Vertex>> boundary_vertices;

        bool isLeaf() 
        {
            return !left && !right;
        }
    };

    double sort_boundary_vertex_list_in_axis(std::vector<std::weak_ptr<Vertex>>& boundary_vertex_list, int axis, int start, int mid, int end)
    {
        std::nth_element(boundary_vertex_list.begin() + start, boundary_vertex_list.begin() + mid, boundary_vertex_list.begin() + end, 
            [&](const std::weak_ptr<Vertex>& boundary_vertex_a, const std::weak_ptr<Vertex>& boundary_vertex_b) 
            {
                return boundary_vertex_a.lock()->get_position()[axis] < boundary_vertex_b.lock()->get_position()[axis];
            });
        return boundary_vertex_list[mid].lock()->get_position()[axis];
    }

    void expand_node_box(std::shared_ptr<Node> node, std::weak_ptr<Vertex> boundary_vertex)
    {
        node->box.expand(boundary_vertex.lock()->get_min());
        node->box.expand(boundary_vertex.lock()->get_max());
    }

    void convert_leaf_to_branch(std::shared_ptr<Node> node)
    {
        std::vector<std::weak_ptr<Vertex>> boundary_vertex_list = node->boundary_vertices;
        int start = 0;
        int end = boundary_vertex_list.size();
        int mid = (start + end) / 2;
        int axis = node->box.get_longest_axis();
        double split_value = sort_boundary_vertex_list_in_axis(boundary_vertex_list, axis, start, mid, end);
        
        node->boundary_vertices.clear();
        node->split_axis = axis;
        node->split_value = split_value;
        node->left = build_node(std::vector<std::weak_ptr<Vertex>>(boundary_vertex_list.begin(), boundary_vertex_list.begin() + mid));
        node->right = build_node(std::vector<std::weak_ptr<Vertex>>(boundary_vertex_list.begin() + mid, boundary_vertex_list.end()));
    }

    std::shared_ptr<Node> build_node(std::vector<std::weak_ptr<Vertex>> boundary_vertex_list)
    {
        auto node = std::make_shared<Node>();
        for (std::weak_ptr<Vertex> boundary_vertex : boundary_vertex_list) 
        {
            expand_node_box(node, boundary_vertex);
            node->boundary_vertices.push_back(boundary_vertex);
        }

        if (node->boundary_vertices.size() > 4) convert_leaf_to_branch(node);

        return node;
    }

    void add_vertex_to_node(std::shared_ptr<Node> node, std::weak_ptr<Vertex> boundary_vertex)
    {
        if (node->isLeaf())
        {
            expand_node_box(node, boundary_vertex);
            node->boundary_vertices.push_back(boundary_vertex);
            if (node->boundary_vertices.size() > 4) convert_leaf_to_branch(node);
        }
        else
        {
            expand_node_box(node, boundary_vertex);

            if (boundary_vertex.lock()->get_position()[node->split_axis] < node->split_value)
            {
                add_vertex_to_node(node->left, boundary_vertex);
            }
            else 
            {
                add_vertex_to_node(node->right, boundary_vertex);
            }
        }
    }

    void increase_vertex_radius_to_node(std::shared_ptr<Node> node, std::weak_ptr<Vertex> boundary_vertex)
    {
        expand_node_box(node, boundary_vertex);

        if (!node->isLeaf())
        {
            if (boundary_vertex.lock()->get_position()[node->split_axis] < node->split_value)
            {
                increase_vertex_radius_to_node(node->left, boundary_vertex);
            }
            else 
            {
                increase_vertex_radius_to_node(node->right, boundary_vertex);
            }
        }
    }

    void delete_vertex_from_node(std::shared_ptr<Node> node, std::weak_ptr<Vertex> boundary_vertex)
    {
        if (node->isLeaf())
        {
            node->boundary_vertices.erase(std::remove(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex), node->boundary_vertices.end());
        }
        else
        {
            if (boundary_vertex.lock()->get_position()[node->split_axis] < node->split_value)
            {
                delete_vertex_from_node(node->left, boundary_vertex);
            }
            else
            {
                delete_vertex_from_node(node->right, boundary_vertex);
            }
        }
    }

    std::set<std::weak_ptr<Vertex>> reverse_radius_search_node(const std::shared_ptr<Node>& node, const Eigen::Vector3d& point) 
    {
        bool contained = node->box.contains(point);
        if (!contained) return std::set<std::weak_ptr<Vertex>>();
        
        if (node->isLeaf())
        {
            std::set<std::weak_ptr<Vertex>> contained_pointIDs;
            for (std::weak_ptr<Vertex> boundary_vertex : node->boundary_vertices)
            {
                if (boundary_vertex.lock()->approx_contains(point)) 
                {
                    contained_pointIDs.insert(boundary_vertex);
                }
            }
            return contained_pointIDs;
        }
        else
        {
            std::set<std::weak_ptr<Vertex>> contained_pointIDs;

            if (node->left->box.contains(point))
            {
                std::set<std::weak_ptr<Vertex>> contained_pointIDs_left = reverse_radius_search_node(node->left, point);
                contained_pointIDs.insert(contained_pointIDs_left.begin(), contained_pointIDs_left.end());
            }
            if (node->right->box.contains(point))
            {
                std::set<std::weak_ptr<Vertex>> contained_pointIDs_right = reverse_radius_search_node(node->right, point);
                contained_pointIDs.insert(contained_pointIDs_right.begin(), contained_pointIDs_right.end());
            }

            return contained_pointIDs;
        }
    }

    std::vector<std::weak_ptr<Vertex>> flatten_node(std::shared_ptr<Node> node)
    {
        std::vector<std::weak_ptr<Vertex>> boundary_vertex_list;
        if (node->isLeaf())
        {
            for (std::weak_ptr<Vertex> boundary_vertex : node->boundary_vertices)
            {
                boundary_vertex_list.push_back(boundary_vertex);
            }
        }
        else
        {
            std::vector<std::weak_ptr<Vertex>> boundary_vertex_list_left = flatten_node(node->left);
            std::vector<std::weak_ptr<Vertex>> boundary_vertex_list_right = flatten_node(node->right);
            boundary_vertex_list.insert(boundary_vertex_list.end(), boundary_vertex_list_left.begin(), boundary_vertex_list_left.end());
            boundary_vertex_list.insert(boundary_vertex_list.end(), boundary_vertex_list_right.begin(), boundary_vertex_list_right.end());
        }
        return boundary_vertex_list;
    }

    void print_node(std::shared_ptr<Node> node, int level) const
    {
        if (node->isLeaf())
        {
            for (std::weak_ptr<Vertex> boundary_vertex : node->boundary_vertices)
            {
                std::cout << "Level: " <<  level << " | ID: " << boundary_vertex.lock()->get_id() << " | Position: " << boundary_vertex.lock()->get_position().transpose() << " | Radius: " << boundary_vertex.lock()->get_radius() << std::endl;
            }
            std::cout << std::endl;
        }
        else
        {
            print_node(node->left, level+1);
            print_node(node->right, level+1);
        }
    }

    double rebuild_threshold;
    int size_at_last_rebuild;

    std::shared_ptr<Node> root;
    int tree_size;

    std::set<std::weak_ptr<Vertex>> vertex_set;

public:

    RRSTree() : rebuild_threshold(2), size_at_last_rebuild(0), tree_size(0)
    {
        rebuild();
        vertex_set.clear();
    }

    void rebuild()
    {
        if (tree_size == 0)
        {
            std::vector<std::weak_ptr<Vertex>> boundary_vertex_list = std::vector<std::weak_ptr<Vertex>>();
            root = build_node(boundary_vertex_list);
        }
        else
        {
            std::vector<std::weak_ptr<Vertex>> boundary_vertex_list = flatten_node(root);
            root = build_node(boundary_vertex_list);
        }
    }

    bool can_reverse_radius_search()
    {
        return (tree_size > 0);
    }

    void add_vertex(std::weak_ptr<Vertex> boundary_vertex)
    {   
        // check if vertex already exists
        if (vertex_set.find(boundary_vertex) != vertex_set.end()) return;

        // add to vertex set
        vertex_set.insert(boundary_vertex);

        // increase size
        tree_size++;

        // add to tree
        if (tree_size > size_at_last_rebuild * rebuild_threshold)
        {
            add_vertex_to_node(root, boundary_vertex);
            rebuild();
            size_at_last_rebuild = tree_size;
        }
        else
        {
            add_vertex_to_node(root, boundary_vertex);
        }
    }

    void delete_vertex(std::weak_ptr<Vertex> boundary_vertex)
    {
        // check if vertex exists
        if (vertex_set.find(boundary_vertex) == vertex_set.end()) return;

        // delete from vertex set
        vertex_set.erase(boundary_vertex);

        // decrease size
        tree_size--;

        // delete from BVH
        delete_vertex_from_node(root, boundary_vertex);
    }

    std::set<std::weak_ptr<Vertex>> reverse_radius_search(const Eigen::Vector3d& point)
    {
        return reverse_radius_search_node(root, point);
    }

    void increase_vertex_radius(std::weak_ptr<Vertex> boundary_vertex)
    {   
        increase_vertex_radius_to_node(root, boundary_vertex);
    }

    void print() const
    {
        print_node(root, 0);
    }

    void print_size()
    {
        std::cout << "Size: " << tree_size << std::endl;
    }

    std::vector<std::weak_ptr<Vertex>> get_vertices()
    {
        return flatten_node(root);
    }
};