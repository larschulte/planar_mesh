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
    if (std::isnan(area))
    {
        // this means the box is deleted, thus should have zero area
        return 0;
    }
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

void EdgeBVH::Node::recursive_expand_parent_box()
{
    if (parent_)
    {
        // expanded
        const bool expanded = parent_->box_.expand_box(box_);

        // recursive update
        if (expanded)
        {
            parent_->recursive_expand_parent_box();
        }
    }
}

void EdgeBVH::Node::recursive_shrink_parent_box()
{
    if (parent_)
    {
        // old parent box
        EdgeBVH::BoundingBox old_parent_box = parent_->box_;

        // new parent box
        EdgeBVH::BoundingBox new_parent_box = EdgeBVH::BoundingBox();

        // copy boxes of children of parent
        EdgeBVH::BoundingBox parent_left_box;
        EdgeBVH::BoundingBox parent_right_box;
        {
            // read lock
            std::shared_lock<std::shared_mutex> lock_left(parent_->left_->rwlock_box_, std::defer_lock);
            std::shared_lock<std::shared_mutex> lock_right(parent_->right_->rwlock_box_, std::defer_lock);
            std::lock(lock_left, lock_right);

            // copy boxes
            parent_left_box = parent_->left_->box_;
            parent_right_box = parent_->right_->box_;
        }

        // expand new parent box
        new_parent_box.expand_box_no_return(parent_left_box);
        new_parent_box.expand_box_no_return(parent_right_box);
        
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
            parent_->box_ = new_parent_box;
            parent_->recursive_shrink_parent_box();
        }
    }
}

EdgeBVH::EdgeBVHReturnType EdgeBVH::Node::node_intersect_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1, std::vector<std::shared_ptr<Edge>>& edges_encountered)
{
    // skip if not intersected
    if (!box_.intersect(vertex0->get_position(), vertex1->get_position())) return EdgeBVHReturnType::SKIP;
        
    if (!isLeaf_)
    {
        // search left and right
        EdgeBVHReturnType left_return = left_->node_intersect_edge(vertex0, vertex1, edges_encountered);
        EdgeBVHReturnType right_return = right_->node_intersect_edge(vertex0, vertex1, edges_encountered);

        // skip if both is skip
        if (left_return == EdgeBVHReturnType::SKIP && right_return == EdgeBVHReturnType::SKIP) return EdgeBVHReturnType::SKIP;
        return EdgeBVHReturnType::INTERSECTED;
    }
    else
    {
        // store
        edges_encountered.push_back(edge_);

        // return
        return EdgeBVHReturnType::INTERSECTED;
    }
}

std::shared_ptr<EdgeBVH::Node> EdgeBVH::find_best_node(const std::shared_ptr<EdgeBVH::Node>& node, const std::shared_ptr<Edge>& edge)
{
    // Store a queue of nodes to process using raw pointers
    std::queue<std::pair<EdgeBVH::Node*, double>> queue;
    queue.push(std::make_pair(node.get(), 0));

    // initialize
    double best_cost = std::numeric_limits<double>::infinity();
    EdgeBVH::Node* best_node = nullptr;

    // compute lower bound cost to add branch node
    EdgeBVH::BoundingBox smallest_branch_box(edge);
    const double lower_bound_cost = smallest_branch_box.get_surface_area();

    // while queue is not empty
    while (!queue.empty())
    {
        // get the first element
        std::pair<EdgeBVH::Node*, double> current = queue.front();
        queue.pop();

        // get the node and inherited cost
        EdgeBVH::Node* current_node = current.first;
        double inherited_cost = current.second;

        // skip if inherited cost is already greater than best cost
        if (inherited_cost > best_cost) continue;

        // cost to add to the node directly if the node is leaf and does not have edge
        if (current_node->isLeaf_ && current_node->edge_ == nullptr)
        {
            // total cost
            double cost = inherited_cost;

            // update best cost and best node
            if (cost < best_cost)
            {
                best_cost = cost;
                best_node = current_node;
            }
        }

        // skip if inherited cost is already greater than best cost
        if (inherited_cost + lower_bound_cost > best_cost) continue;

        // cost to add a branch that contains current node and leaf node
        {
            // cost of creating a new branch node
            EdgeBVH::BoundingBox new_branch_box = current_node->box_;
            new_branch_box.expand_box_no_return(edge);
            double new_branch_node_cost = new_branch_box.get_surface_area();

            // total cost
            double cost = inherited_cost + new_branch_node_cost;

            // update best cost and best node
            if (cost < best_cost)
            {
                best_cost = cost;
                best_node = current_node;
            }
        }

        // check if it is worth it to go to the children
        if (!current_node->isLeaf_)
        {
            // compute change to inherited cost
            EdgeBVH::BoundingBox expanded_box = current_node->box_;
            expanded_box.expand_box_no_return(edge);
            double change_to_inherited = expanded_box.get_surface_area() - current_node->box_.get_surface_area();

            if (inherited_cost + change_to_inherited + lower_bound_cost > best_cost)
            {
                // it is not worth it to go to the children
                continue;
            }
            else
            {
                // we should go into the children
                queue.push(std::make_pair(current_node->left_.get(), inherited_cost + change_to_inherited));
                queue.push(std::make_pair(current_node->right_.get(), inherited_cost + change_to_inherited));
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


void EdgeBVH::Node::node_add_edge(const std::shared_ptr<Edge>& edge)
{    
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_node_);

    // if node is leaf, and do not have edge, simply add edge to node
    if (isLeaf_ && edge_ == nullptr)
    {
        edge_ = edge;
        edge_->node = shared_from_this();
        box_.expand_box_no_return(edge->get_min(), edge->get_max());
        recursive_expand_parent_box();
        return;
    }

    // create new node
    std::shared_ptr<EdgeBVH::Node> new_leaf_node = std::make_shared<EdgeBVH::Node>();
    {
        new_leaf_node->box_.expand_box_no_return(edge);
        new_leaf_node->edge_ = edge;
        new_leaf_node->edge_->node = new_leaf_node;
    }
    
    // create duplicate node of the current node added to
    std::shared_ptr<EdgeBVH::Node> duplicate_node = std::make_shared<EdgeBVH::Node>();
    {
        duplicate_node->box_ = box_;
        
        // if node is leaf, copy boundary vertices
        if (isLeaf_)
        {
            duplicate_node->edge_ = edge_;
            if (duplicate_node->edge_) duplicate_node->edge_->node = duplicate_node;
        }
        else
        {
            // else, copy children
            duplicate_node->left_ = left_;
            duplicate_node->right_ = right_;

            left_->parent_ = duplicate_node;
            right_->parent_ = duplicate_node;

            duplicate_node->isLeaf_.store(isLeaf_.load());
        }
    }
    
    // make new node and duplicate node children of the current node
    {
        // expand current node box
        box_.expand_box_no_return(new_leaf_node->box_);
        recursive_expand_parent_box();

        // put into children
        left_ = duplicate_node;
        right_ = new_leaf_node;

        duplicate_node->parent_ = shared_from_this();
        new_leaf_node->parent_ = shared_from_this();

        // change to branch
        isLeaf_.store(false);
    }    
}

void EdgeBVH::Node::node_delete_edge(const std::shared_ptr<Edge>& edge)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_node_);

    // remove node from edge
    edge->node = nullptr;

    // edge from node
    edge_ = nullptr;
    
    // keep in parent's left or right

    // reset box
    {
        // write lock
        std::unique_lock<std::shared_mutex> lock_box(rwlock_box_);

        // reset box
        box_ = BoundingBox();
    }
    
    // shrink parent box
    recursive_shrink_parent_box();
}

void EdgeBVH::Node::node_print(int level) const
{
    if (!isLeaf_)
    {
        left_->node_print(level+1);
        right_->node_print(level+1);
    }
    else
    {
        std::cout << "Level: " <<  level << " | ID: " << edge_->get_id() << " | Center: " << edge_->get_center().transpose() << std::endl;

        std::cout << std::endl;
    }
}

void EdgeBVH::Node::node_flatten(std::vector<std::shared_ptr<Edge>>& edge_list) const
{
    if (!isLeaf_)
    {
        left_->node_flatten(edge_list);
        right_->node_flatten(edge_list);
    }
    else
    {
        edge_list.push_back(edge_);
    }
}

std::vector<std::shared_ptr<Edge>> EdgeBVH::get_edge_list() const
{
    std::vector<std::shared_ptr<Edge>> edge_list;
    root_->node_flatten(edge_list);
    return edge_list;
}

void EdgeBVH::tree_delete_edge(const std::shared_ptr<Edge>& edge)
{
    // get node
    std::shared_ptr<Node> node = edge->node;

    // skip if edge not in tree
    if (!node) return;
    
    // remove from node
    node->node_delete_edge(edge);

    // decrement size
    edge_size_--;
}

EdgeBVH::EdgeBVH() : edge_size_(0)
{
    root_ = std::make_shared<EdgeBVH::Node>();
}

void EdgeBVH::tree_add_edge(const std::shared_ptr<Edge>& edge)
{
    // get node
    std::shared_ptr<Node> node = edge->node;

    // skip if edge already in tree
    if (node) return;

    // find best node to add
    std::shared_ptr<Node> best_node = find_best_node(root_, edge);

    // add to best node
    best_node->node_add_edge(edge);

    // increment count
    edge_size_++;
}

EdgeBVH::EdgeBVHReturnType EdgeBVH::tree_intersect_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1, std::vector<std::shared_ptr<Edge>>& edges_encountered)
{
    // check input
    if (vertex0->is_expired() || vertex1->is_expired()) throw std::runtime_error("Attempts to intersect with expired vertex.");

    return root_->node_intersect_edge(vertex0, vertex1, edges_encountered);
}

void EdgeBVH::tree_print() const
{
    root_->node_print(0);
}
