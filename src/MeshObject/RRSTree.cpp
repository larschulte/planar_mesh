#include "MeshObject/RRSTree.hpp"
#include "MeshObject/Vertex.hpp"

#include "MeshObject/Surface.hpp"

std::ostream& operator<<(std::ostream& os, const RRSReturnType& type)
{
    switch (type)
    {
        case RRSReturnType::INTERSECTED:
            os << "_";
            break;
        case RRSReturnType::SKIP:
            os << "-";
            break;
        case RRSReturnType::ABORT:
            os << "X";
            break;
        default:
            os << "?";
            break;
    }
    return os;
}

RRSBoundingBox::RRSBoundingBox() : 
    min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
    max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

RRSBoundingBox::RRSBoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max) : min(min), max(max) {}

bool RRSBoundingBox::expand(const Eigen::Vector3d& point)
{
    Eigen::Vector3d oldMin = min;
    Eigen::Vector3d oldMax = max;

    // Update min and max to include the new point
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);

    // Check if min or max changed
    bool changed = (min != oldMin || max != oldMax);
    
    return changed;
}

bool RRSBoundingBox::contains(const Eigen::Vector3d& point)
{
    return (point[0] >= min[0] && point[0] <= max[0] &&
            point[1] >= min[1] && point[1] <= max[1] &&
            point[2] >= min[2] && point[2] <= max[2]);
}

int RRSBoundingBox::get_longest_axis()
{
    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

void RRSNode::recursive_unlock()
{
    omp_unset_nest_lock(&omp_lock);
    if (locked_children)
    {
        left->recursive_unlock();
        right->recursive_unlock();
        locked_children = false;
    }
}

void RRSNode::recursive_expand_parent_box()
{
    if (parent)
    {
        // expanded
        const bool expanded = parent->box.expand(box.min) || parent->box.expand(box.max);

        // recursive update
        if (expanded)
        {
            parent->recursive_expand_parent_box();
        }
    }
}

void RRSNode::recursive_shrink_parent_box()
{
    if (parent)
    {
        // old parent box
        RRSBoundingBox old_parent_box = parent->box;

        // new parent box
        RRSBoundingBox new_parent_box = RRSBoundingBox();
        if (!std::isnan(parent->left->box.min[0])) new_parent_box.expand(parent->left->box.min);
        if (!std::isnan(parent->left->box.max[0])) new_parent_box.expand(parent->left->box.max);
        if (!std::isnan(parent->right->box.min[0])) new_parent_box.expand(parent->right->box.min);
        if (!std::isnan(parent->right->box.max[0])) new_parent_box.expand(parent->right->box.max);
        
        // shrunk
        const bool shrunk = new_parent_box.min[0] > old_parent_box.min[0] &&
                            new_parent_box.min[1] > old_parent_box.min[1] &&
                            new_parent_box.min[2] > old_parent_box.min[2] &&
                            new_parent_box.max[0] < old_parent_box.max[0] &&
                            new_parent_box.max[1] < old_parent_box.max[1] &&
                            new_parent_box.max[2] < old_parent_box.max[2];

        // recursive update
        if (shrunk) 
        {
            parent->box = new_parent_box;
            parent->recursive_shrink_parent_box();
        }
    }
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

void RRSTree::expand_node_box(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    node->box.expand(boundary_vertex->get_min());
    node->box.expand(boundary_vertex->get_max());
}

void RRSTree::convert_leaf_to_branch(const std::shared_ptr<RRSNode>& node)
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
    node->left->parent = node;
    node->right->parent = node;
    node->left->sibling = node->right;
    node->right->sibling = node->left;
    node->boundary_vertices.clear();

    // set lock before setting isLeaf to false to prevent other threads from accessing the children node
    omp_set_nest_lock(&node->left->omp_lock);
    omp_set_nest_lock(&node->right->omp_lock);
    node->locked_children = true;

    // update isLeaf after locking the children
    node->isLeaf = false;
}

std::shared_ptr<RRSNode> RRSTree::build_node(const std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, const int& start, const int& end)
{
    auto node = std::make_shared<RRSNode>();

    // expand box
    for (int i = start; i < end; i++)
    {
        expand_node_box(node, boundary_vertex_list[i]);
    }

    // store vertices
    node->boundary_vertices = std::vector<std::shared_ptr<Vertex>>(boundary_vertex_list.begin() + start, boundary_vertex_list.begin() + end);

    // store node pointer in vertices
    for (const std::shared_ptr<Vertex>& vertex : node->boundary_vertices)
    {
        vertex->node = node;
    }

    // convert to branch
    if (node->boundary_vertices.size() > leaf_size) convert_leaf_to_branch(node);

    return node;
}

void RRSTree::node_add_vertex(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    expand_node_box(node, boundary_vertex);

    if (!node->isLeaf)
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

        // add node pointer to vertex
        boundary_vertex->node = node;

        if (node->boundary_vertices.size() > leaf_size) convert_leaf_to_branch(node);
    }
}

void RRSTree::node_increase_radius(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    expand_node_box(node, boundary_vertex);

    if (!node->isLeaf)
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

bool RRSTree::node_delete_vertex(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    if (!node->isLeaf)
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
            boundary_vertex->node = nullptr;
            return true;
        }
        else
        {
            return false;
        }
    }
}

RRSReturnType RRSTree::node_reverse_radius_search(const std::shared_ptr<RRSNode>& node, const Eigen::Vector3d& point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    // skip if not contained
    if (!node->box.contains(point))
    {
        return RRSReturnType::SKIP;
    }

    // branch if not leaf
    if (!node->isLeaf)
    {
        // search left and right
        RRSReturnType left_return = node_reverse_radius_search(node->left, point, search_results);
        if (left_return == RRSReturnType::ABORT) return RRSReturnType::ABORT;
        RRSReturnType right_return = node_reverse_radius_search(node->right, point, search_results);
        if (right_return == RRSReturnType::ABORT) return RRSReturnType::ABORT;

        // skip if both is skip
        if (left_return == RRSReturnType::SKIP && right_return == RRSReturnType::SKIP) return RRSReturnType::SKIP;
        // intersected if any is intersected
        return RRSReturnType::INTERSECTED;
    }
    else
    {
        // abort if can't lock node
        if (!omp_test_nest_lock(&node->omp_lock))
        {
            // std::cout << "_ _ _ _ X _ _ _" << std::endl;
            return RRSReturnType::ABORT;
        }
        
        // skip if no vertices
        if (node->boundary_vertices.size() == 0)
        {
            omp_unset_nest_lock(&node->omp_lock);
            return RRSReturnType::SKIP;
        }

        // skip if boundary vertex is not searchable
        if (!node->boundary_vertices[0]->is_searchable())
        {
            omp_unset_nest_lock(&node->omp_lock);
            return RRSReturnType::SKIP;
        }

        // skip if not contained
        if (!node->boundary_vertices[0]->approx_contains(point))
        {
            omp_unset_nest_lock(&node->omp_lock);
            return RRSReturnType::SKIP;
        }

        // abort if can't lock vertex's surface
        const std::shared_ptr<Vertex>& boundary_vertex = node->boundary_vertices[0];
        const std::shared_ptr<Surface>& surface = boundary_vertex->get_surface_check();
        if (!omp_test_nest_lock(&surface->lock)) 
        {
            omp_unset_nest_lock(&node->omp_lock);
            // std::cout << "_ _ _ _ _ _ X _" << std::endl;
            return RRSReturnType::ABORT;
        }

        // return
        search_results.push_back(node->boundary_vertices[0]);
        return RRSReturnType::INTERSECTED;
    }
}

RRSReturnType RRSTree::node_find_leaf_node(const std::shared_ptr<RRSNode>& node, const Eigen::Vector3d& point, std::shared_ptr<RRSNode>& return_node)
{    
    // branch if not leaf
    if (!node->isLeaf)
    {
        if (point[node->split_axis] < node->split_value)
        {
            return node_find_leaf_node(node->left, point, return_node);
        }
        else
        {
            return node_find_leaf_node(node->right, point, return_node);
        }
    }
    else
    {
        // abort if can't lock leaf node
        if (!omp_test_nest_lock(&node->omp_lock))
        {
            // std::cout << "_ _ _ _ _ X _ _" << std::endl;
            return RRSReturnType::ABORT;
        }

        // return the leaf node locked if no vertex
        const bool no_vertex = node->boundary_vertices.size() == 0;
        if (no_vertex)
        {
            return_node = node;
            return RRSReturnType::INTERSECTED;
        }
        
        // if there is vertex, return new empty leaf node
        
        // temperary vertex
        std::shared_ptr<Vertex> temp_vertex = std::make_shared<Vertex>(); 
        temp_vertex->temp_initialize(point, 0);

        //  add and branch
        node->boundary_vertices.push_back(temp_vertex);
        convert_leaf_to_branch(node); // this may cause error in Application when locking surface's node

        // delete
        return_node = point[node->split_axis] < node->split_value ? node->left : node->right;
        return_node->boundary_vertices.pop_back();

        // locks
        omp_set_nest_lock(&return_node->omp_lock);
        node->recursive_unlock();

        // return
        return RRSReturnType::INTERSECTED;
    }
}

void RRSTree::node_flattern(const std::shared_ptr<RRSNode>& node, std::vector<std::shared_ptr<Vertex>>& flatten_list)
{
    if (!node->isLeaf)
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

void RRSTree::node_print(const std::shared_ptr<RRSNode>& node, int level) const
{
    if (!node->isLeaf)
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

RRSTree::RRSTree() : rebuild_threshold(5), size_at_last_rebuild(0), tree_size(0), leaf_size(1)
{
    rebuild();
}

void RRSTree::check_rebuild()
{
    if (tree_size > size_at_last_rebuild * rebuild_threshold)
    {
        std::cout << "Rebuilding RRS tree ...." << std::endl;
        rebuild();
        std::cout << "Rebuilding RRS tree done" << std::endl;
        size_at_last_rebuild = tree_size;
    }

    // release lock
    if (root->left) root->left->recursive_unlock();
    if (root->right) root->right->recursive_unlock();
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
}

void RRSTree::tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // check input
    if (boundary_vertex->is_expired()) throw std::invalid_argument("Invalid vertex in boundary_vertex_list");

    // decrease size
    tree_size--;

    // get vertex's node reference
    const std::shared_ptr<RRSNode>& node = boundary_vertex->node;
    if (node == nullptr) throw std::invalid_argument("Vertex not found in BVH.");

    // delete from node
    node->boundary_vertices.erase(std::remove(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex), node->boundary_vertices.end());
    boundary_vertex->node = nullptr;
}

RRSReturnType RRSTree::tree_reverse_radius_search(const Eigen::Vector3d& point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    return node_reverse_radius_search(root, point, search_results);
}

RRSReturnType RRSTree::tree_find_leaf_node(const Eigen::Vector3d& point, std::shared_ptr<RRSNode>& return_node)
{
    return node_find_leaf_node(root, point, return_node);
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