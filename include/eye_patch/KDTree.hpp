#pragma once

#include <set>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <memory>


struct Node {
    Node(){}
    int pointID;
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

    void sort_points_in_axis(std::vector<int>& point_list, int axis, int start, int mid, int end)
    {
        std::nth_element(point_list.begin() + start, point_list.begin() + mid, point_list.begin() + end, 
            [&](const int& point_a, const int& point_b) 
            {
                Eigen::Vector3d vector_a = point_to_vector3d_map.at(point_a);
                Eigen::Vector3d vector_b = point_to_vector3d_map.at(point_b);
                return vector_a[axis] < vector_b[axis];
            });
    }

    std::shared_ptr<Node> build_node(std::vector<int> point_list, int depth)
    {
        auto node = std::make_shared<Node>();
        
        int size = point_list.size();
        if (size == 0)
        {
            return nullptr;
        }
        else if (size == 1)
        {
            node->pointID = point_list[0];
        }
        else
        {
            int axis = depth % 3;
            int start = 0;
            int end = point_list.size();
            int mid = end / 2;
            sort_points_in_axis(point_list, axis, start, mid, end);

            node->pointID = point_list[mid];
            node->left = build_node(std::vector<int>(point_list.begin(), point_list.begin() + mid), depth + 1);
            node->right = build_node(std::vector<int>(point_list.begin() + mid + 1, point_list.end()), depth + 1);
        }

        return node;
    }

    void addPointToNode(std::shared_ptr<Node>& node, int point_id, int depth)
    {
        // null node
        if (node == nullptr)
        {
            node = std::make_shared<Node>();
            node->pointID = point_id;
            return;
        }
        // branch/leaf node
        else if (node->pointID == point_id)
        {
            node->deleted = false;
            return;
        }
        // children nodes
        else
        {
            int axis = depth % 3;
            if (point_to_vector3d_map.at(point_id)[axis] < point_to_vector3d_map.at(node->pointID)[axis])
            {
                addPointToNode(node->left, point_id, depth + 1);
            }
            else
            {
                addPointToNode(node->right, point_id, depth + 1);
            }
        }
    }

    void addPoint(Eigen::Vector3d new_point, int point_id)
    {
        // skip if point already exist
        if (std::find(point_list.begin(), point_list.end(), point_id) != point_list.end()) return;

        // store data
        point_list.push_back(point_id);
        point_to_vector3d_map[point_id] = new_point;

        // conditional rebuild
        if (point_list.size() >= size_at_last_rebuild * rebuild_threshold)
        {
            rebuild();
            size_at_last_rebuild = point_list.size();
        }
        else
        {
            addPointToNode(root, point_id, 0);
        }
        
        // // always rebuild
        // rebuild();
    }

    void radiusSearch(std::shared_ptr<Node> node, Eigen::Vector3d target, double radius, int depth, std::set<int>& result)
    {
        // null node
        if (node == nullptr)
        {
            return;
        }

        // current node
        if ((point_to_vector3d_map.at(node->pointID) - target).norm() <= radius)
        {
            if (!node->deleted) result.insert(node->pointID);
        }

        // children nodes   
        int axis = depth % 3;
        double diff = target[axis] - point_to_vector3d_map.at(node->pointID)[axis];

        if (diff < 0)
        {
            radiusSearch(node->left, target, radius, depth + 1, result);
            if (fabs(diff) <= radius)
            {
                radiusSearch(node->right, target, radius, depth + 1, result);
            }
        }
        else
        {
            radiusSearch(node->right, target, radius, depth + 1, result);
            if (fabs(diff) <= radius)
            {
                radiusSearch(node->left, target, radius, depth + 1, result);
            }
        }

        // radiusSearch(node->left, target, radius, depth + 1, result);
        // radiusSearch(node->right, target, radius, depth + 1, result);
    }

    std::set<int> radiusSearch(Eigen::Vector3d target, double radius)
    {
        // initialize
        std::set<int> result;

        // search
        radiusSearch(root, target, radius, 0, result);

        // return
        return result;
    }

    std::shared_ptr<Node> search(std::shared_ptr<Node> node, int pointID, int depth) 
    {
        // null node
        if (node == nullptr) 
        {
            return nullptr;
        }

        // current node
        if (node->pointID == pointID) 
        {
            return node;
        }

        // children nodes
        int axis = depth % 3;
        if (point_to_vector3d_map.at(pointID)[axis] < point_to_vector3d_map.at(node->pointID)[axis])
        {
            return search(node->left, pointID, depth + 1);
        } 
        else
        {
            return search(node->right, pointID, depth + 1);
        }
    }

    void deletePoint(int pointID)
    {
        // error check
        if (point_to_vector3d_map.find(pointID) == point_to_vector3d_map.end())
        {
            throw std::invalid_argument("PointID given does not have corresponding vector3d data.");
        }

        // delete from point_list
        point_list.erase(std::remove(point_list.begin(), point_list.end(), pointID), point_list.end());

        // delete from tree
        auto node = search(root, pointID, 0);
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

        std::cout << "Node: " << node->pointID << std::endl;
        print_node(node->left);
        print_node(node->right);
    }

    std::vector<int> point_list;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::shared_ptr<Node> root;
    std::size_t size_at_last_rebuild;
    int rebuild_threshold;
};
