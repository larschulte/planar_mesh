#pragma once

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <array>
#include <iostream> 

struct BoundaryPoint 
{
    int pointID;
    Eigen::Vector3d position;
    double radius;
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    
    BoundaryPoint(int pointID, Eigen::Vector3d position) : pointID(pointID), position(position), radius(0) {}

    BoundaryPoint(int pointID, Eigen::Vector3d position, double radius)
        : pointID(pointID), position(position), radius(radius)
    {
        min = position - Eigen::Vector3d::Constant(radius);
        max = position + Eigen::Vector3d::Constant(radius);
    }

    void adjustRadius(double radius)
    {
        this->radius = radius;
        min = position - Eigen::Vector3d::Constant(radius);
        max = position + Eigen::Vector3d::Constant(radius);
    }

    bool contains(const Eigen::Vector3d& point)
    {
        return (point.array() >= min.array()).all() && (point.array() <= max.array()).all();
    }

    // overload == operator for std::remove
    bool operator==(const BoundaryPoint& other) const
    {
        return pointID == other.pointID;
    }

    // overload < operator for std::set
    bool operator<(const BoundaryPoint& other) const
    {
        return pointID < other.pointID;
    }
};

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
        std::vector<BoundaryPoint> boundaryPoints;

        bool isLeaf() 
        {
            return !left && !right;
        }
    };

    double sort_boundary_point_list_in_axis(std::vector<BoundaryPoint>& boundary_point_list, int axis, int start, int mid, int end)
    {
        std::nth_element(boundary_point_list.begin() + start, boundary_point_list.begin() + mid, boundary_point_list.begin() + end, 
            [&](const BoundaryPoint& boundary_point_a, const BoundaryPoint& boundary_point_b) 
            {
                return boundary_point_a.position[axis] < boundary_point_b.position[axis];
            });
        return boundary_point_list[mid].position[axis];
    }

    void expand_node_box(std::shared_ptr<Node> node, BoundaryPoint boundary_point)
    {
        node->box.expand(boundary_point.min);
        node->box.expand(boundary_point.max);
    }

    void convert_leaf_to_branch(std::shared_ptr<Node> node)
    {
        std::vector<BoundaryPoint> boundary_point_list = node->boundaryPoints;
        int start = 0;
        int end = boundary_point_list.size();
        int mid = (start + end) / 2;
        int axis = node->box.get_longest_axis();
        double split_value = sort_boundary_point_list_in_axis(boundary_point_list, axis, start, mid, end);
        
        node->boundaryPoints.clear();
        node->split_axis = axis;
        node->split_value = split_value;
        node->left = build_node(std::vector<BoundaryPoint>(boundary_point_list.begin(), boundary_point_list.begin() + mid));
        node->right = build_node(std::vector<BoundaryPoint>(boundary_point_list.begin() + mid, boundary_point_list.end()));
    }

    std::shared_ptr<Node> build_node(std::vector<BoundaryPoint> boundary_point_list)
    {
        auto node = std::make_shared<Node>();
        for (BoundaryPoint boundary_point : boundary_point_list) 
        {
            expand_node_box(node, boundary_point);
            node->boundaryPoints.push_back(boundary_point);
        }

        if (node->boundaryPoints.size() > 4) convert_leaf_to_branch(node);

        return node;
    }

    void addBoundaryPointToNode(std::shared_ptr<Node> node, BoundaryPoint boundary_point)
    {
        if (node->isLeaf())
        {
            expand_node_box(node, boundary_point);
            node->boundaryPoints.push_back(boundary_point);
            if (node->boundaryPoints.size() > 4) convert_leaf_to_branch(node);
        }
        else
        {
            expand_node_box(node, boundary_point);

            if (boundary_point.position[node->split_axis] < node->split_value)
            {
                addBoundaryPointToNode(node->left, boundary_point);
            }
            else 
            {
                addBoundaryPointToNode(node->right, boundary_point);
            }
        }
    }

    void adjustBoundaryPointRadiusToNode(std::shared_ptr<Node> node, BoundaryPoint boundary_point)
    {
        if (node->isLeaf())
        {
            expand_node_box(node, boundary_point);

            for (BoundaryPoint& node_boundary_point : node->boundaryPoints)
            {
                if (node_boundary_point.pointID == boundary_point.pointID)
                {
                    node_boundary_point.adjustRadius(boundary_point.radius);
                }
            }
        }
        else
        {
            expand_node_box(node, boundary_point);

            if (boundary_point.position[node->split_axis] < node->split_value)
            {
                adjustBoundaryPointRadiusToNode(node->left, boundary_point);
            }
            else 
            {
                adjustBoundaryPointRadiusToNode(node->right, boundary_point);
            }
        }
    }


    void deleteBoundaryPointFromNode(std::shared_ptr<Node> node, BoundaryPoint boundary_point)
    {
        if (node->isLeaf())
        {
            node->boundaryPoints.erase(std::remove(node->boundaryPoints.begin(), node->boundaryPoints.end(), boundary_point), node->boundaryPoints.end());
        }
        else
        {
            if (boundary_point.position[node->split_axis] < node->split_value)
            {
                deleteBoundaryPointFromNode(node->left, boundary_point);
            }
            else
            {
                deleteBoundaryPointFromNode(node->right, boundary_point);
            }
        }
    }


    std::set<int> reverseRadiusSearchNode(const std::shared_ptr<Node>& node, const Eigen::Vector3d& point, std::map<int, double>& pointID_radius_map) 
    {
        bool contained = node->box.contains(point);
        if (!contained) return std::set<int>();
        
        if (node->isLeaf())
        {
            std::set<int> contained_pointIDs;
            for (BoundaryPoint boundary_point : node->boundaryPoints)
            {
                if (boundary_point.contains(point)) 
                {
                    contained_pointIDs.insert(boundary_point.pointID);
                    pointID_radius_map[boundary_point.pointID] = boundary_point.radius;
                }
            }
            return contained_pointIDs;
        }
        else
        {
            std::set<int> contained_pointIDs;

            if (node->left->box.contains(point))
            {
                std::set<int> contained_pointIDs_left = reverseRadiusSearchNode(node->left, point, pointID_radius_map);
                contained_pointIDs.insert(contained_pointIDs_left.begin(), contained_pointIDs_left.end());
            }
            if (node->right->box.contains(point))
            {
                std::set<int> contained_pointIDs_right = reverseRadiusSearchNode(node->right, point, pointID_radius_map);
                contained_pointIDs.insert(contained_pointIDs_right.begin(), contained_pointIDs_right.end());
            }

            return contained_pointIDs;
        }
    }

    std::vector<BoundaryPoint> flatten_node(std::shared_ptr<Node> node)
    {
        std::vector<BoundaryPoint> boundary_point_list;
        if (node->isLeaf())
        {
            for (BoundaryPoint boundary_point : node->boundaryPoints)
            {
                boundary_point_list.push_back(boundary_point);
            }
        }
        else
        {
            std::vector<BoundaryPoint> boundary_point_list_left = flatten_node(node->left);
            std::vector<BoundaryPoint> boundary_point_list_right = flatten_node(node->right);
            boundary_point_list.insert(boundary_point_list.end(), boundary_point_list_left.begin(), boundary_point_list_left.end());
            boundary_point_list.insert(boundary_point_list.end(), boundary_point_list_right.begin(), boundary_point_list_right.end());
        }
        return boundary_point_list;
    }

    void printNode(std::shared_ptr<Node> node, int level)
    {
        if (node->isLeaf())
        {
            for (BoundaryPoint boundary_point : node->boundaryPoints)
            {
                std::cout << "Level: " <<  level << " | ID: " << boundary_point.pointID << " | Position: " << boundary_point.position.transpose() << " | Radius: " << boundary_point.radius << std::endl;
            }
            std::cout << std::endl;
        }
        else
        {
            printNode(node->left, level+1);
            printNode(node->right, level+1);
        }
    }

    double rebuild_threshold;
    int size_at_last_rebuild;

    std::shared_ptr<Node> root;
    int tree_size;

public:

    RRSTree() : rebuild_threshold(2), size_at_last_rebuild(0), tree_size(0)
    {
        rebuildTree();
    }

    void rebuildTree()
    {
        if (tree_size == 0)
        {
            std::vector<BoundaryPoint> boundary_point_list = std::vector<BoundaryPoint>();
            root = build_node(boundary_point_list);
        }
        else
        {
            std::vector<BoundaryPoint> boundary_point_list = flatten_node(root);
            root = build_node(boundary_point_list);
        }
    }

    void addBoundaryPoint(int pointID, const Eigen::Vector3d& position, double radius)
    {
        BoundaryPoint boundary_point(pointID, position, radius);
        
        tree_size++;

        // add to tree
        if (tree_size > size_at_last_rebuild * rebuild_threshold)
        {
            addBoundaryPointToNode(root, boundary_point);
            rebuildTree();
            size_at_last_rebuild = tree_size;
        }
        else
        {
            addBoundaryPointToNode(root, boundary_point);
        }
    }

    void deleteBoundaryPoint(int pointID, const Eigen::Vector3d& position)
    {
        BoundaryPoint boundary_point(pointID, position);

        tree_size--;

        // delete from BVH
        deleteBoundaryPointFromNode(root, boundary_point);
    }

    std::set<int> reverseRadiusSearch(Eigen::Vector3d point, std::map<int, double>& pointID_radius_map)
    {
        return reverseRadiusSearchNode(root, point, pointID_radius_map);
    }

    void adjustRadius(int pointID, Eigen::Vector3d position, double radius)
    {   
        BoundaryPoint boundary_point(pointID, position, radius);
        adjustBoundaryPointRadiusToNode(root, boundary_point);
    }

    void printTree()
    {
        printNode(root, 0);
    }

    void printSize()
    {
        std::cout << "Size: " << tree_size << std::endl;
    }
};