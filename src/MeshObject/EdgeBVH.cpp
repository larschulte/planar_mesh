#include "MeshObject/EdgeBVH.hpp"
#include <iostream>

#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Surface.hpp"


EdgeBVH::BoundingBox::BoundingBox()
    : min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
      max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

void EdgeBVH::BoundingBox::expand(const Eigen::Vector3d& point) 
{
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);
}

bool EdgeBVH::BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const 
{
    for (int i = 0; i < 3; ++i) 
    {
        double invD = 1.0 / dir[i];
        double t0 = (min[i] - orig[i]) * invD;
        double t1 = (max[i] - orig[i]) * invD;
        if (invD < 0.0) std::swap(t0, t1);
        tMin = std::max(tMin, t0);
        tMax = std::min(tMax, t1);
        if (tMax < tMin) return false;
    }
    return true;
}

bool EdgeBVH::BoundingBox::intersect(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const 
{
    // count as intersect if both points are inside
    bool a_inside = (a[0] >= min[0] && a[0] <= max[0] &&
                     a[1] >= min[1] && a[1] <= max[1] &&
                     a[2] >= min[2] && a[2] <= max[2]);
    bool b_inside = (b[0] >= min[0] && b[0] <= max[0] &&
                     b[1] >= min[1] && b[1] <= max[1] &&
                     b[2] >= min[2] && b[2] <= max[2]);
    if (a_inside && b_inside) return true;

    // Check for intersection if not both points are inside
    double tMin = 0.0;
    double tMax = 1.0;
    Eigen::Vector3d dir = b - a;
    return intersect(a, dir, tMin, tMax) && tMax >= 0.0 && tMin <= 1.0;
}

int EdgeBVH::BoundingBox::get_longest_axis()
{
    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

bool EdgeBVH::Node::isLeaf() const 
{
    return !left && !right;
}


double EdgeBVH::sort_edge_list_in_axis(std::vector<std::shared_ptr<Edge>>& edge_list, int axis, int start, int mid, int end)
{
    std::sort(edge_list.begin() + start, edge_list.begin() + end, 
        [&](const std::shared_ptr<Edge>& edge_a, const std::shared_ptr<Edge>& edge_b) 
        {
            return edge_a->get_center()[axis] < edge_b->get_center()[axis];
        });
    return edge_list[mid]->get_center()[axis];
}

void EdgeBVH::expand_node_box(const std::shared_ptr<Node>& node, const std::shared_ptr<Edge>& edge)
{
    node->box.expand(edge->get_max());
    node->box.expand(edge->get_min());
}

bool EdgeBVH::node_intersect_edge(const std::shared_ptr<Node>& node, const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
{
    bool intersected = node->box.intersect(vertex0->get_position(), vertex1->get_position());
    if (!intersected) return false;
    
    if (!node->isLeaf())
    {
        if (node_intersect_edge(node->left, vertex0, vertex1)) return true;
        if (node_intersect_edge(node->right, vertex0, vertex1)) return true;   
    }
    else
    {
        for (const std::shared_ptr<Edge>& edge : node->edges)
        {
            if (edge->intersects_edge(surface_, vertex0, vertex1)) return true;
        }
    }

    return false;
}

void EdgeBVH::convert_leaf_to_branch(const std::shared_ptr<Node>& node)
{
    int start = 0;
    int end = node->edges.size();
    int mid = (start + end) / 2;
    int axis = node->box.get_longest_axis();
    double split_value = sort_edge_list_in_axis(node->edges, axis, start, mid, end);
    
    node->split_axis = axis;
    node->split_value = split_value;
    node->left = build_node(node->edges, start, mid);
    node->right = build_node(node->edges, mid, end);

    node->edges.clear();
}

std::shared_ptr<EdgeBVH::Node> EdgeBVH::build_node(const std::vector<std::shared_ptr<Edge>>& edge_list, const int& start, const int& end)
{
    auto node = std::make_shared<Node>();

    // expand box
    for (int i = start; i < end; ++i)
    {
        expand_node_box(node, edge_list[i]);
    }

    // store vector
    node->edges = std::vector<std::shared_ptr<Edge>>(edge_list.begin() + start, edge_list.begin() + end);

    // convert to branch if necessary
    if (node->edges.size() > 4) convert_leaf_to_branch(node);

    return node;
}

void EdgeBVH::node_add_edge(const std::shared_ptr<Node>& node, const std::shared_ptr<Edge>& edge)
{
    expand_node_box(node, edge);

    if (!node->isLeaf())
    {   
        if (edge->get_center()[node->split_axis] < node->split_value)
        {
            node_add_edge(node->left, edge);
        }
        else 
        {
            node_add_edge(node->right, edge);
        }
    }
    else
    {
        node->edges.push_back(edge);
        if (node->edges.size() > 4) convert_leaf_to_branch(node);
    }
}

bool EdgeBVH::node_delete_edge(const std::shared_ptr<Node>& node, const std::shared_ptr<Edge>& edge)
{
    if (!node->isLeaf())
    {
        if (edge->get_center()[node->split_axis] < node->split_value)
        {
            return node_delete_edge(node->left, edge);
        }
        else if (edge->get_center()[node->split_axis] > node->split_value)
        {
            return node_delete_edge(node->right, edge);
        }
        else
        {
            return node_delete_edge(node->left, edge) || node_delete_edge(node->right, edge);
        }
    }
    else
    {
        auto it = std::remove(node->edges.begin(), node->edges.end(), edge);
        if (it != node->edges.end())
        {
            node->edges.erase(it, node->edges.end());
            return true;
        }
        else
        {
            return false;
        }
    }
}

void EdgeBVH::node_print(const std::shared_ptr<Node>& node, int level) const
{
    if (!node->isLeaf())
    {
        node_print(node->left, level+1);
        node_print(node->right, level+1);
    }
    else
    {
        for (const std::shared_ptr<Edge>& edge : node->edges)
        {
            std::cout << "Level: " <<  level << " | ID: " << edge->get_id() << " | Center: " << edge->get_center().transpose() << std::endl;
        }
        std::cout << std::endl;
    }
}

void EdgeBVH::node_flatten(const std::shared_ptr<EdgeBVH::Node>& node, std::vector<std::shared_ptr<Edge>>& edge_list) const
{
    if (!node->isLeaf())
    {
        node_flatten(node->left, edge_list);
        node_flatten(node->right, edge_list);
    }
    else
    {
        edge_list.insert(edge_list.end(), node->edges.begin(), node->edges.end());
    }
}

std::vector<std::shared_ptr<Edge>> EdgeBVH::get_edge_list() const
{
    std::vector<std::shared_ptr<Edge>> edge_list;
    node_flatten(root, edge_list);
    return edge_list;
}

void EdgeBVH::tree_delete_edge(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to delete expired edge.");

    // decrement size
    edge_size--;
    
    // delete from BVH
    if (!node_delete_edge(root, edge)) throw std::runtime_error("Edge not found in BVH.");
}

EdgeBVH::EdgeBVH()
    : rebuild_threshold(2),
      size_at_last_rebuild(0), 
      edge_size(0)
{
    rebuild();
}

void EdgeBVH::set_surface(const std::shared_ptr<Surface>& surface)
{
    surface_ = surface;
}

void EdgeBVH::rebuild()
{
    if (edge_size == 0)
    {
        root = build_node(std::vector<std::shared_ptr<Edge>>(), 0, 0);
    }
    else
    {
        std::vector<std::shared_ptr<Edge>> edge_list = get_edge_list();
        root = build_node(edge_list, 0, edge_list.size());
    }
}

void EdgeBVH::tree_add_edge(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to add expired edge.");

    // increment count
    edge_size++;

    node_add_edge(root, edge);

    if (edge_size > size_at_last_rebuild * rebuild_threshold)
    {    
        rebuild();
        size_at_last_rebuild = edge_size;
    }
}

bool EdgeBVH::tree_intersect_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
{
    // check input
    if (vertex0->is_expired() || vertex1->is_expired()) throw std::runtime_error("Attempts to intersect with expired vertex.");

    return node_intersect_edge(root, vertex0, vertex1);
}

void EdgeBVH::tree_print() const
{
    node_print(root, 0);
}
