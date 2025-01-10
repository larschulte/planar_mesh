#include "MeshObject/EdgeBVH.hpp"
#include <iostream>

#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Surface.hpp"


EdgeBVH::BoundingBox::BoundingBox()
    : min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
      max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

EdgeBVH::BoundingBox::BoundingBox(const std::shared_ptr<Edge>& edge) 
    : min(edge->get_min()), max(edge->get_max()) {}

void EdgeBVH::BoundingBox::expand(const Eigen::Vector3d& point) 
{
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);
}

void EdgeBVH::BoundingBox::expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max) 
{
    if (input_min[0] < min[0]) min[0] = input_min[0];
    if (input_min[1] < min[1]) min[1] = input_min[1];
    if (input_min[2] < min[2]) min[2] = input_min[2];

    if (input_max[0] > max[0]) max[0] = input_max[0];
    if (input_max[1] > max[1]) max[1] = input_max[1];
    if (input_max[2] > max[2]) max[2] = input_max[2];
}

void EdgeBVH::BoundingBox::expand_box_no_return(const std::shared_ptr<Edge>& edge) 
{
    expand(edge->get_min());
    expand(edge->get_max());
}

void EdgeBVH::BoundingBox::expand_box_no_return(const EdgeBVH::BoundingBox& box) 
{
    expand_box_no_return(box.min, box.max);
}

bool EdgeBVH::BoundingBox::expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max) 
{
    Eigen::Vector3d oldMin = min;
    Eigen::Vector3d oldMax = max;

    expand_box_no_return(input_min, input_max);

    return (min != oldMin || max != oldMax);
}

bool EdgeBVH::BoundingBox::expand_box(const EdgeBVH::BoundingBox& box) 
{
    return expand_box(box.min, box.max);
}

bool EdgeBVH::BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir_inv, double& tMin, double& tMax) const 
{
    for (int i = 0; i < 3; ++i) 
    {
        const double& invD = dir_inv[i];
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
    // count as intersect if either point is inside
    const bool a_inside = (a[0] >= min[0] && a[0] <= max[0] &&
                     a[1] >= min[1] && a[1] <= max[1] &&
                     a[2] >= min[2] && a[2] <= max[2]);
    const bool b_inside = (b[0] >= min[0] && b[0] <= max[0] &&
                     b[1] >= min[1] && b[1] <= max[1] &&
                     b[2] >= min[2] && b[2] <= max[2]);
    if (a_inside || b_inside) return true;

    // Check for intersection if no point is inside
    double tMin = 0.0;
    double tMax = 1.0;
    const Eigen::Vector3d dir = b - a;
    const Eigen::Vector3d dir_inv = dir.cwiseInverse();
    return intersect(a, dir_inv, tMin, tMax) && tMax >= 0.0 && tMin <= 1.0;
}

int EdgeBVH::BoundingBox::get_longest_axis()
{
    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

double EdgeBVH::BoundingBox::compute_surface_area() const
{
    Eigen::Vector3d dimensions = max - min;
    double area = 2.0 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[0]);
    if (std::isnan(area)) throw std::runtime_error("BoundingBox::compute_surface_area() returned nan."); // throw if nan
    return area;
}

const double& EdgeBVH::BoundingBox::get_surface_area()
{
    // if min and max are not updated, return the stored value
    if (min == min_used_for_surface_area && max == max_used_for_surface_area) return surface_area;

    // else update and return
    min_used_for_surface_area = min;
    max_used_for_surface_area = max;
    surface_area = compute_surface_area();
    return surface_area;
}

bool EdgeBVH::Node::isLeaf() const 
{
    return !left && !right;
}

void EdgeBVH::Node::recursive_expand_parent_box()
{
    if (parent)
    {
        // expanded
        const bool expanded = parent->box.expand_box(box);

        // recursive update
        if (expanded)
        {
            parent->recursive_expand_parent_box();
        }
    }
}

void EdgeBVH::Node::recursive_shrink_parent_box()
{
    if (parent)
    {
        // old parent box
        EdgeBVH::BoundingBox old_parent_box = parent->box;

        // new parent box
        EdgeBVH::BoundingBox new_parent_box = EdgeBVH::BoundingBox();
        new_parent_box.expand_box_no_return(parent->left->box);
        new_parent_box.expand_box_no_return(parent->right->box);
        
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

void EdgeBVH::sort_edge_list_in_axis(std::vector<std::shared_ptr<Edge>>& edge_list, int axis, int start, int end)
{
    std::sort(edge_list.begin() + start, edge_list.begin() + end, 
        [&](const std::shared_ptr<Edge>& edge_a, const std::shared_ptr<Edge>& edge_b) 
        {
            return edge_a->get_center()[axis] < edge_b->get_center()[axis];
        });
}

void EdgeBVH::expand_node_box(const std::shared_ptr<EdgeBVH::Node>& node, const std::shared_ptr<Edge>& edge)
{
    node->box.expand(edge->get_max());
    node->box.expand(edge->get_min());
}

bool EdgeBVH::node_intersect_edge(const std::shared_ptr<EdgeBVH::Node>& node, const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
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

std::shared_ptr<EdgeBVH::Node> EdgeBVH::find_best_node(const std::shared_ptr<EdgeBVH::Node>& node, const std::shared_ptr<Edge>& edge)
{
    // Store a queue of nodes to process using raw pointers
    std::queue<std::pair<EdgeBVH::Node*, double>> queue;
    queue.push(std::make_pair(node.get(), 0));

    // initialize
    double best_cost = std::numeric_limits<double>::infinity();
    EdgeBVH::Node* best_node = nullptr;

    // while queue is not empty
    while (!queue.empty())
    {
        // get the first element
        std::pair<EdgeBVH::Node*, double> current = queue.front();
        queue.pop();

        // get the node and inherited cost
        EdgeBVH::Node* current_node = current.first;
        double inherited_cost = current.second;

        // cost to branch from current node
        double cost;
        {
            // cost of creating a new branch node
            EdgeBVH::BoundingBox new_branch_box = current_node->box;
            new_branch_box.expand_box_no_return(edge);
            double new_branch_node_cost = new_branch_box.get_surface_area();

            // total cost
            cost = inherited_cost + new_branch_node_cost;
        }

        if (cost < best_cost)
        {
            best_cost = cost;
            best_node = current_node;
        }

        // check if it is worth it to go to the children
        if (!current_node->isLeaf())
        {
            // compute change to inherited cost
            EdgeBVH::BoundingBox expanded_box = current_node->box;
            expanded_box.expand_box_no_return(edge);
            double change_to_inherited = expanded_box.get_surface_area() - current_node->box.get_surface_area();

            // compute lower bound cost to add branch node
            EdgeBVH::BoundingBox smallest_branch_box(edge);
            double lower_bound_cost = smallest_branch_box.get_surface_area();
            
            if (inherited_cost + change_to_inherited + lower_bound_cost > best_cost)
            {
                // it is not worth it to go to the children
                continue;
            }
            else
            {
                // we should go into the children
                queue.push(std::make_pair(current_node->left.get(), inherited_cost + change_to_inherited));
                queue.push(std::make_pair(current_node->right.get(), inherited_cost + change_to_inherited));
            }
        }
    }

    // Return the best node as a shared_ptr
    if (best_node)
    {
        // Since best_node is part of the existing tree, we can create a shared_ptr using std::shared_ptr aliasing constructor
        // This shares ownership with the original node's shared_ptr
        return std::shared_ptr<EdgeBVH::Node>(node, best_node);
    }
    else
    {
        return nullptr;
    }
}


void EdgeBVH::node_add_edge(const std::shared_ptr<EdgeBVH::Node>& node, const std::shared_ptr<Edge>& edge)
{    
    // create new node
    std::shared_ptr<EdgeBVH::Node> new_node = std::make_shared<EdgeBVH::Node>();
    {
        new_node->box.expand_box_no_return(edge);
        new_node->edges.push_back(edge);
        edge->node = new_node;
    }
    
    // create duplicate node of the current node added to
    std::shared_ptr<EdgeBVH::Node> duplicate_node = std::make_shared<EdgeBVH::Node>();
    {
        duplicate_node->box = node->box;
        duplicate_node->split_value = node->split_value;
        duplicate_node->split_axis = node->split_axis;
        
        // if node is leaf, copy boundary vertices
        if (node->isLeaf())
        {
            duplicate_node->edges = node->edges;
            for (const std::shared_ptr<Edge>& edge : duplicate_node->edges)
            {
                edge->node = duplicate_node;
            }
        }
        else
        {
            // else, copy children
            duplicate_node->left = node->left;
            duplicate_node->right = node->right;

            node->left->parent = duplicate_node;
            node->right->parent = duplicate_node;
        }
    }
    
    // make new node and duplicate node children of the current node
    {
        // expand current node box
        node->box.expand_box_no_return(new_node->box);
        node->recursive_expand_parent_box();

        // get split axis
        node->split_axis = node->box.get_longest_axis();

        // get split value  
        node->split_value = edge->get_center()[node->split_axis];

        // put into children
        node->left = duplicate_node;
        node->right = new_node;

        duplicate_node->parent = node;
        new_node->parent = node;
    }    
}

void EdgeBVH::node_print(const std::shared_ptr<EdgeBVH::Node>& node, int level) const
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

    // get node
    std::shared_ptr<Node> node = edge->node;
    if (!node) throw std::runtime_error("Edge does not belong to any node.");
    
    // delete from BVH
    auto it = std::remove(node->edges.begin(), node->edges.end(), edge);
    node->edges.erase(it, node->edges.end());
    edge->node = nullptr;
    
    // recompute box
    node->box = BoundingBox();
    for (const std::shared_ptr<Edge>& edge : node->edges)
    {
        expand_node_box(node, edge);
    }

    // shrink parent box
    node->recursive_shrink_parent_box();
}

EdgeBVH::EdgeBVH()
    : rebuild_threshold(2),
      size_at_last_rebuild(0), 
      edge_size(0)
{
    root = std::make_shared<EdgeBVH::Node>();
}

void EdgeBVH::set_surface(const std::shared_ptr<Surface>& surface)
{
    surface_ = surface;
}

void EdgeBVH::tree_add_edge(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to add expired edge.");

    // increment count
    edge_size++;

    // find best node to add
    std::shared_ptr<Node> best_node = find_best_node(root, edge);

    // add to best node
    node_add_edge(best_node, edge);
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
