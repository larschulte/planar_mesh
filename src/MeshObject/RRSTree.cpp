#include "MeshObject/RRSTree.hpp"
#include "MeshObject/Vertex.hpp"

RRSTree::BoundingBox::BoundingBox() : 
    min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
    max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

void RRSTree::BoundingBox::expand(const Eigen::Vector3d& point)
{
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);
}

bool RRSTree::BoundingBox::contains(const Eigen::Vector3d& point)
{
    return (point[0] >= min[0] && point[0] <= max[0] &&
            point[1] >= min[1] && point[1] <= max[1] &&
            point[2] >= min[2] && point[2] <= max[2]);
}

int RRSTree::BoundingBox::get_longest_axis()
{
    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

bool RRSTree::Node::isLeaf() const
{
    return !left && !right;
}

double RRSTree::sort_boundary_vertex_list_in_axis(std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, int axis, int start, int mid, int end)
{
    std::sort(boundary_vertex_list.begin() + start, boundary_vertex_list.begin() + end, 
        [&](const std::shared_ptr<Vertex>& boundary_vertex_a, const std::shared_ptr<Vertex>& boundary_vertex_b) 
        {
            return boundary_vertex_a->get_position()[axis] < boundary_vertex_b->get_position()[axis];
        });
    return boundary_vertex_list[mid]->get_position()[axis];
}

void RRSTree::expand_node_box(const std::shared_ptr<Node>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    node->box.expand(boundary_vertex->get_min());
    node->box.expand(boundary_vertex->get_max());
}

void RRSTree::convert_leaf_to_branch(const std::shared_ptr<Node>& node)
{
    int start = 0;
    int end = node->boundary_vertices.size();
    int mid = (start + end) / 2;
    int axis = node->box.get_longest_axis();
    double split_value = sort_boundary_vertex_list_in_axis(node->boundary_vertices, axis, start, mid, end);
    
    node->split_axis = axis;
    node->split_value = split_value;
    node->left = build_node(node->boundary_vertices, start, mid);
    node->right = build_node(node->boundary_vertices, mid, end);
    node->boundary_vertices.clear();
}

std::shared_ptr<RRSTree::Node> RRSTree::build_node(const std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, const int& start, const int& end)
{
    auto node = std::make_shared<Node>();

    // expand box
    for (int i = start; i < end; i++)
    {
        expand_node_box(node, boundary_vertex_list[i]);
    }

    // store vertices
    node->boundary_vertices = std::vector<std::shared_ptr<Vertex>>(boundary_vertex_list.begin() + start, boundary_vertex_list.begin() + end);

    // convert to branch
    if (node->boundary_vertices.size() > leaf_size) convert_leaf_to_branch(node);

    return node;
}

void RRSTree::node_add_vertex(const std::shared_ptr<Node>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    expand_node_box(node, boundary_vertex);

    if (!node->isLeaf())
    {    
        if (boundary_vertex->get_position()[node->split_axis] < node->split_value)
        {
            node_add_vertex(node->left, boundary_vertex);
        }
        else 
        {
            node_add_vertex(node->right, boundary_vertex);
        }
    }
    else
    {
        node->boundary_vertices.push_back(boundary_vertex);
        if (node->boundary_vertices.size() > leaf_size) convert_leaf_to_branch(node);
    }
}

void RRSTree::node_increase_radius(const std::shared_ptr<Node>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    expand_node_box(node, boundary_vertex);

    if (!node->isLeaf())
    {
        if (boundary_vertex->get_position()[node->split_axis] < node->split_value)
        {
            node_increase_radius(node->left, boundary_vertex);
        }
        else 
        {
            node_increase_radius(node->right, boundary_vertex);
        }
    }
}

bool RRSTree::node_delete_vertex(const std::shared_ptr<Node>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    if (!node->isLeaf())
    {
        if (boundary_vertex->get_position()[node->split_axis] < node->split_value)
        {
            return node_delete_vertex(node->left, boundary_vertex);
        }
        else if (boundary_vertex->get_position()[node->split_axis] > node->split_value)
        {
            return node_delete_vertex(node->right, boundary_vertex);
        }
        else
        {
            return node_delete_vertex(node->left, boundary_vertex) || node_delete_vertex(node->right, boundary_vertex);
        }
    }
    else
    {
        auto it = std::remove(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex);
        if (it != node->boundary_vertices.end())
        {
            node->boundary_vertices.erase(it, node->boundary_vertices.end());
            return true;
        }
        else
        {
            return false;
        }
    }
}

void RRSTree::node_reverse_radius_search(const std::shared_ptr<Node>& node, const Eigen::Vector3d& point, std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& search_results)
{
    bool contained = node->box.contains(point);
    if (!contained) return;
    
    if (!node->isLeaf())
    {
        if (node->left->box.contains(point))
        {
            node_reverse_radius_search(node->left, point, search_results);
        }
        if (node->right->box.contains(point))
        {
            node_reverse_radius_search(node->right, point, search_results);
        }
    }
    else
    {
        for (const std::shared_ptr<Vertex>& boundary_vertex : node->boundary_vertices)
        {
            if (boundary_vertex->approx_contains(point)) 
            {
                search_results.insert(boundary_vertex);
            }
        }
    }
}

void RRSTree::node_flattern(const std::shared_ptr<Node>& node, std::vector<std::shared_ptr<Vertex>>& flatten_list)
{
    if (!node->isLeaf())
    {
        node_flattern(node->left, flatten_list);
        node_flattern(node->right, flatten_list);
    }
    else
    {
        flatten_list.insert(flatten_list.end(), node->boundary_vertices.begin(), node->boundary_vertices.end());
    }
}

std::vector<std::shared_ptr<Vertex>> RRSTree::compute_vertices_list()
{
    std::vector<std::shared_ptr<Vertex>> flatten_list;
    node_flattern(root, flatten_list);
    return flatten_list;
}

void RRSTree::node_print(const std::shared_ptr<Node>& node, int level) const
{
    if (!node->isLeaf())
    {
        node_print(node->left, level+1);
        node_print(node->right, level+1);
    }
    else
    {
        for (const std::shared_ptr<Vertex>& boundary_vertex : node->boundary_vertices)
        {
            std::cout << "Level: " <<  level << " | ID: " << boundary_vertex->get_id() << " | Position: " << boundary_vertex->get_position().transpose() << " | Radius: " << boundary_vertex->get_radius() << std::endl;
        }
        std::cout << std::endl;
    }
}

RRSTree::RRSTree() : rebuild_threshold(5), size_at_last_rebuild(0), tree_size(0), leaf_size(64)
{
    rebuild();
}

void RRSTree::rebuild()
{
    if (tree_size == 0)
    {
        root = build_node(std::vector<std::shared_ptr<Vertex>>(), 0, 0);
    }
    else
    {
        std::vector<std::shared_ptr<Vertex>> boundary_vertex_list = compute_vertices_list();
        root = build_node(boundary_vertex_list, 0, boundary_vertex_list.size());
    }
}

bool RRSTree::can_reverse_radius_search()
{
    return (tree_size > 0);
}

void RRSTree::tree_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{   
    // check input
    if (boundary_vertex->is_expired()) throw std::invalid_argument("Invalid vertex in boundary_vertex_list");

    // increase size
    tree_size++;

    node_add_vertex(root, boundary_vertex);

    // add to tree
    if (tree_size > size_at_last_rebuild * rebuild_threshold)
    {    
        rebuild();
        size_at_last_rebuild = tree_size;
    }
}

void RRSTree::tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // check input
    if (boundary_vertex->is_expired()) throw std::invalid_argument("Invalid vertex in boundary_vertex_list");

    // decrease size
    tree_size--;

    // delete from BVH
    if (!node_delete_vertex(root, boundary_vertex)) throw std::invalid_argument("Vertex not found in BVH.");
}

void RRSTree::tree_reverse_radius_search(const Eigen::Vector3d& point, std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& search_results)
{
    node_reverse_radius_search(root, point, search_results);
}

void RRSTree::tree_increase_radius(std::shared_ptr<Vertex> boundary_vertex)
{   
    // check input
    if (boundary_vertex->is_expired()) throw std::invalid_argument("Invalid vertex in boundary_vertex_list");

    node_increase_radius(root, boundary_vertex);
}

void RRSTree::tree_print() const
{
    node_print(root, 0);
}

void RRSTree::print_size()
{
    std::cout << "Size: " << tree_size << std::endl;
}