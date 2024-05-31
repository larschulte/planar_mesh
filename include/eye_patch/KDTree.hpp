#pragma once

#include <set>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <memory>

struct Point
{
    int pointID;
    Eigen::Vector3d position;

    // overload == operator for std::remove
    bool operator==(const Point& other) const
    {
        return pointID == other.pointID;
    }

    // overload < operator for std::set
    bool operator<(const Point& other) const
    {
        return pointID < other.pointID;
    }
};

struct Node {
    Node(){}
    Point point;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
    bool deleted = false;
};


class KDTree {
public:
    KDTree()
        : size_at_last_rebuild(0), rebuild_threshold(2)
    {
        rebuild();
    }

    void rebuild()
    {
        root = build_node(point_list, 0);
    }

    void sort_points_in_axis(std::vector<Point>& point_list, int axis, int start, int mid, int end)
    {
        std::nth_element(point_list.begin() + start, point_list.begin() + mid, point_list.begin() + end, 
            [&](const Point& point_a, const Point& point_b) 
            {
                return point_a.position[axis] < point_b.position[axis];
            });
    }

    std::shared_ptr<Node> build_node(std::vector<Point> point_list, int depth)
    {
        auto node = std::make_shared<Node>();
        
        int size = point_list.size();
        if (size == 0)
        {
            return nullptr;
        }
        else if (size == 1)
        {
            node->point = point_list[0];
        }
        else
        {
            int axis = depth % 3;
            int start = 0;
            int end = point_list.size();
            int mid = end / 2;
            sort_points_in_axis(point_list, axis, start, mid, end);

            node->point = point_list[mid];
            node->left = build_node(std::vector<Point>(point_list.begin(), point_list.begin() + mid), depth + 1);
            node->right = build_node(std::vector<Point>(point_list.begin() + mid + 1, point_list.end()), depth + 1);
        }

        return node;
    }

    void addPointToNode(std::shared_ptr<Node>& node, Point point, int depth)
    {
        // null node
        if (node == nullptr)
        {
            node = std::make_shared<Node>();
            node->point = point;
            return;
        }
        // branch/leaf node
        else if (node->point == point)
        {
            node->deleted = false;
            return;
        }
        // children nodes
        else
        {
            int axis = depth % 3;
            if (point.position[axis] < node->point.position[axis])
            {
                addPointToNode(node->left, point, depth + 1);
            }
            else
            {
                addPointToNode(node->right, point, depth + 1);
            }
        }
    }

    void addPoint(int point_id, Eigen::Vector3d new_point)
    {
        Point point;
        point.pointID = point_id;
        point.position = new_point;

        // skip if point already exist
        if (std::find(point_list.begin(), point_list.end(), point) != point_list.end()) return;

        // store data
        point_list.push_back(point);

        // conditional rebuild
        if (point_list.size() >= size_at_last_rebuild * rebuild_threshold)
        {
            rebuild();
            size_at_last_rebuild = point_list.size();
        }
        else
        {
            addPointToNode(root, point, 0);
        }
        
        // // always rebuild
        // rebuild();
    }

    void radiusSearch(std::shared_ptr<Node> node, const Eigen::Vector3d& target, double radius, double radius_squared, int depth, std::set<int>& result)
    {
        // null node
        if (node == nullptr)
        {
            return;
        }

        // current node
        if ((node->point.position - target).squaredNorm() <= radius_squared)
        {
            if (!node->deleted) result.insert(node->point.pointID);
        }

        // children nodes   
        int axis = depth % 3;
        double diff = target[axis] - node->point.position[axis];

        if (diff < 0)
        {
            radiusSearch(node->left, target, radius, radius_squared, depth + 1, result);
            if (fabs(diff) <= radius)
            {
                radiusSearch(node->right, target, radius, radius_squared, depth + 1, result);
            }
        }
        else
        {
            radiusSearch(node->right, target, radius, radius_squared, depth + 1, result);
            if (fabs(diff) <= radius)
            {
                radiusSearch(node->left, target, radius, radius_squared, depth + 1, result);
            }
        }

        // radiusSearch(node->left, target, radius, radius_squared, depth + 1, result);
        // radiusSearch(node->right, target, radius, radius_squared, depth + 1, result);
    }

    std::set<int> radiusSearch(Eigen::Vector3d target, double radius)
    {
        // initialize
        std::set<int> result;

        // search
        radiusSearch(root, target, radius, radius*radius, 0, result);

        // return
        return result;
    }

    std::shared_ptr<Node> search(std::shared_ptr<Node> node, Point point, int depth) 
    {
        // null node
        if (node == nullptr) 
        {
            return nullptr;
        }

        // current node
        if (node->point == point) 
        {
            return node;
        }

        // children nodes
        int axis = depth % 3;
        if (point.position[axis] < node->point.position[axis])
        {
            return search(node->left, point, depth + 1);
        } 
        else
        {
            return search(node->right, point, depth + 1);
        }
    }

    void deletePoint(int pointID, Eigen::Vector3d new_point)
    {
        Point point;
        point.pointID = pointID;
        point.position = new_point;

        // delete from point_list
        point_list.erase(std::remove(point_list.begin(), point_list.end(), point), point_list.end());

        // delete from tree
        auto node = search(root, point, 0);
        if (node != nullptr)
        {
            node->deleted = true;
        }
    }

    void print_tree()
    {
        std::cout << "Printing tree" << std::endl;
        print_node(root);
    }

    void print_node(std::shared_ptr<Node> node)
    {
        if (node == nullptr)
        {
            return;
        }

        std::cout << "Node: " << node->point.pointID << std::endl;
        print_node(node->left);
        print_node(node->right);
    }

    void print_size()
    {
        std::cout << "Size: " << point_list.size() << std::endl;
    }

    std::vector<Point> point_list;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::shared_ptr<Node> root;
    std::size_t size_at_last_rebuild;
    int rebuild_threshold;
};
