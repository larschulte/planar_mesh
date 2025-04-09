#include "MeshObject/RRSTree.hpp"
#include "MeshObject/Vertex.hpp"

#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include <queue>
#include <thread>

Settings RRSTree::settings_;

std::ostream& operator<<(std::ostream& os, const RRSReturnType& type)
{
    switch (type)
    {
        case RRSReturnType::INTERSECTED:
            os << "I _ _";
            break;
        case RRSReturnType::SKIP:
            os << "_ S _";
            break;
        default:
            os << "? ? ?";
            break;
    }
    return os;
}

RRSBoundingBox& RRSBoundingBox::operator=(const RRSBoundingBox& other)
{
    min = other.min;
    max = other.max;
    min_used_for_surface_area = other.min_used_for_surface_area;
    max_used_for_surface_area = other.max_used_for_surface_area;
    surface_area = other.surface_area;
    
    return *this;
}

RRSBoundingBox::RRSBoundingBox(const RRSBoundingBox& other)
{
    // std::shared_lock lock(other.mutex_); // Acquire shared (read) lock

    min = other.min;
    max = other.max;
    min_used_for_surface_area = other.min_used_for_surface_area;
    max_used_for_surface_area = other.max_used_for_surface_area;
    surface_area = other.surface_area;
} 

RRSBoundingBox::RRSBoundingBox() : 
    min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
    max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

RRSBoundingBox::RRSBoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max) : min(min), max(max) {}

bool RRSBoundingBox::expand(const Eigen::Vector3d& point)
{
    // std::unique_lock lock(mutex_); // Acquire unique (write) lock

    Eigen::Vector3d oldMin = min;
    Eigen::Vector3d oldMax = max;

    // Update min and max to include the new point
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);

    // Check if min or max changed
    bool changed = (min != oldMin || max != oldMax);
    
    return changed;
}

void RRSBoundingBox::expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max)
{
    // std::unique_lock lock(mutex_); // Acquire unique (write) lock

    // Component-wise min and max without calling Eigen's functions
    if (input_min[0] < min[0]) min[0] = input_min[0];
    if (input_min[1] < min[1]) min[1] = input_min[1];
    if (input_min[2] < min[2]) min[2] = input_min[2];

    if (input_max[0] > max[0]) max[0] = input_max[0];
    if (input_max[1] > max[1]) max[1] = input_max[1];
    if (input_max[2] > max[2]) max[2] = input_max[2];
}

void RRSBoundingBox::expand_box_no_return(const RRSBoundingBox& box)
{
    // std::shared_lock lock(box.mutex_); // Acquire shared (read) lock
    expand_box_no_return(box.min, box.max);
}

bool RRSBoundingBox::expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max)
{
    // make copy of old min and max
    Eigen::Vector3d oldMin = min;
    Eigen::Vector3d oldMax = max;

    // Update min and max to include the new box
    expand_box_no_return(input_min, input_max);

    // Check if min or max changed
    bool changed = (min != oldMin || max != oldMax);

    return changed;
}

bool RRSBoundingBox::expand_box(const RRSBoundingBox& box)
{
    // std::shared_lock lock(box.mutex_); // Acquire shared (read) lock
    return expand_box(box.min, box.max);
} 

bool RRSBoundingBox::contains(const Eigen::Vector3d& point)
{
    // std::shared_lock lock(mutex_); // Acquire shared (read) lock

    return (point[0] >= min[0] && point[0] <= max[0] &&
            point[1] >= min[1] && point[1] <= max[1] &&
            point[2] >= min[2] && point[2] <= max[2]);
}

int RRSBoundingBox::get_longest_axis()
{
    // std::shared_lock lock(mutex_); // Acquire shared (read) lock

    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

double RRSBoundingBox::get_surface_area()
{
    // std::unique_lock lock(mutex_); // Lock for read/write

    if (min != min_used_for_surface_area || max != max_used_for_surface_area) 
    {
        Eigen::Vector3d dimensions = max - min;
        double area = 2.0 * (dimensions[0] * dimensions[1] + 
                             dimensions[1] * dimensions[2] + 
                             dimensions[2] * dimensions[0]);

        surface_area = std::isnan(area) ? 0.0 : area;
        min_used_for_surface_area = min;
        max_used_for_surface_area = max;
    }

    return surface_area; // Return by value
}

std::atomic<unsigned int> RRSNode::counter_ = 0;

RRSNode::RRSNode()
{
    id_ = counter_.fetch_add(1);
}

void RRSNode::recursive_expand_parent_box()
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

void RRSNode::recursive_shrink_parent_box()
{
    if (parent_)
    {
        // old parent box
        RRSBoundingBox old_parent_box = parent_->box_;

        // new parent box
        RRSBoundingBox new_parent_box = RRSBoundingBox();

        // expand new parent box
        {
            // // read lock
            // std::shared_lock<std::shared_mutex> lock(parent_->rwlock_node_);

            // expand left box
            if (parent_->left_)
            {
                // // read lock
                // std::shared_lock<std::shared_mutex> lock_left(parent_->left_->rwlock_node_);

                new_parent_box.expand_box_no_return(parent_->left_->box_);
            }

            // expand right box
            if (parent_->right_)
            {
                // // read lock
                // std::shared_lock<std::shared_mutex> lock_right(parent_->right_->rwlock_node_);

                new_parent_box.expand_box_no_return(parent_->right_->box_);
            }
        }
        
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

std::shared_ptr<RRSNode> RRSTree::find_best_node(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    // Store a queue of nodes to process using raw pointers
    std::queue<std::pair<RRSNode*, double>> queue;
    queue.push(std::make_pair(node.get(), 0));

    // initialize
    double best_cost = std::numeric_limits<double>::infinity();
    RRSNode* best_node = nullptr;

    // compute lower bound cost to add branch node
    RRSBoundingBox smallest_branch_box(boundary_vertex->get_min(), boundary_vertex->get_max());
    const double lower_bound_cost = smallest_branch_box.get_surface_area();

    // while queue is not empty
    while (!queue.empty())
    {
        // get the first element
        std::pair<RRSNode*, double> current = queue.front();
        queue.pop();

        // get the node and inherited cost
        RRSNode* current_node = current.first;
        double inherited_cost = current.second;

        // skip if inherited cost is already greater than best cost
        if (inherited_cost > best_cost) continue;

        // cost to add to the node directly if the node is leaf and does not have boundary vertex
        if (current_node->isLeaf_ && current_node->boundary_vertex_ == nullptr)
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
            RRSBoundingBox new_branch_box = current_node->box_;
            new_branch_box.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());
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
            RRSBoundingBox expanded_box = current_node->box_;
            expanded_box.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());
            double change_to_inherited = expanded_box.get_surface_area() - current_node->box_.get_surface_area();

            // lower-bound cost
            double cost = inherited_cost + change_to_inherited + lower_bound_cost;
            
            // update best cost and best node
            if (cost < best_cost)
            {
                // we should check the children
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
        return std::shared_ptr<RRSNode>(node, best_node);
    }
    else
    {
        return nullptr;
    }
}


//
// PROBLEM
//
// if we do not have find leaf node function.
// we need to lock the node here which will cause deadlock.

// 
// FACTS
// 
// within one thread, as all searches are performed before we process point
// thus addition of point to search tree can be done after we process points
// 
// thus we add the point into the search tree because we want future searches to find this point, not the current process point function        

//
// SOLUTION
// 
// only add vertex to tree or faces to tree after all surface lock and node lock are released. 
// when adding vertex from storage, add to a queue. that is processed after all locks are released
void RRSNode::node_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_node_);

    // if node is leaf, and do not have boundary vertex, simply add boundary_vertex to node
    if (isLeaf_ && boundary_vertex_ == nullptr)
    {
        boundary_vertex_ = boundary_vertex;
        boundary_vertex_->node = shared_from_this();
        box_.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());
        recursive_expand_parent_box();
        return;
    }
    
    // create new node
    std::shared_ptr<RRSNode> new_leaf_node = std::make_shared<RRSNode>();
    {
        new_leaf_node->box_.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());
        new_leaf_node->boundary_vertex_ = boundary_vertex;
        new_leaf_node->boundary_vertex_->node = new_leaf_node;
    }
    
    // create duplicate node
    std::shared_ptr<RRSNode> duplicate_node = std::make_shared<RRSNode>();
    {
        duplicate_node->box_ = box_;
        
        // if node is leaf, copy boundary vertices
        if (isLeaf_)
        {
            duplicate_node->boundary_vertex_ = boundary_vertex_;
            if (duplicate_node->boundary_vertex_) duplicate_node->boundary_vertex_->node = duplicate_node;
        }
        else
        {
            // else, copy children
            duplicate_node->left_ = left_;
            left_->parent_ = duplicate_node;

            duplicate_node->right_ = right_;
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
        duplicate_node->parent_ = shared_from_this();
        right_ = new_leaf_node;
        new_leaf_node->parent_ = shared_from_this();

        // // assign sibling
        // duplicate_node->sibling_ = sibling_;
        // new_leaf_node->sibling_ = duplicate_node;

        // change to branch
        isLeaf_ = false;
    }    
}

void RRSNode::node_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_node_);

    // this node could be locked by other thread and add new vertex before we obtain the lock
    if (!isLeaf_)
    {
        // obtain new node and delete vertex from it
        std::shared_ptr<RRSNode> new_node = boundary_vertex->node;

        // remove from node
        new_node->node_delete_vertex(boundary_vertex);

        return;
    }

    // throw if have no parent
    if (!parent_)
    {
        throw std::runtime_error("node_delete_vertex, Node has no parent.");
    }

    // actual deletion
    {
        // remove node from vertices
        boundary_vertex->node = nullptr;

        // remove vertex from node
        boundary_vertex_ = nullptr;

        // reset box
        box_ = RRSBoundingBox(); // here

        // shrink parent box
        recursive_shrink_parent_box();
    }
}

void RRSNode::node_delete_self()
{
    // skip if contains vertex (after the vertex is deleted, the node could be given another vertex)
    if (boundary_vertex_) return;

    // skip if is branch now (after storing antoher vertex, the node could become a branch node)
    if (!isLeaf_) return;

    // throw if have no parent (unless this is truely root node, in which case i should implement something)
    if (!parent_) throw std::runtime_error("node_delete_self, Node has no parent.");

    // remove from parent and simplify tree
    parent_->node_delete_child(shared_from_this());
}

void RRSNode::node_delete_child(const std::shared_ptr<RRSNode> child)
{
    // replace this parent with the other child, instead of copying it, this allow the other child, if also removed vertex, to be deleted as well
    
    // throw if have no parent (unless this is truely root node, in which case i should implement something)
    if (!parent_) throw std::runtime_error("node_delete_child, Node has no parent.");
    
    // throw if any left or right is null (since this function is called by children)
    if (!left_ || !right_) throw std::runtime_error("One or more Child is null.");

    // get grandparent
    std::shared_ptr<RRSNode> grandparent = parent_;

    // get other child
    std::shared_ptr<RRSNode> other_child;
    if (left_ == child)
    {
        other_child = right_;
    }
    else if (right_ == child)
    {
        other_child = left_;
    }
    else
    {
        throw std::runtime_error("Input child is not a child of its parent.");
    }

    // disconnect grandparent
    {
        // grandparent -> parent
        if (grandparent->left_ == shared_from_this())
        {
            grandparent->left_ = nullptr;
        }
        else if (grandparent->right_ == shared_from_this())
        {
            grandparent->right_ = nullptr;
        }
        else
        {
            throw std::runtime_error("This node is not a child of the grandparent.");
        }

        // parent -> grandparent
        parent_ = nullptr;
    }

    // disconnect child and other child
    {
        // parent -> children
        left_ = nullptr;
        right_ = nullptr;

        // children -> parent
        child->parent_ = nullptr;
        other_child->parent_ = nullptr;
    }

    // connect other child with grandparent
    {
        // grandparent -> other child
        if (grandparent->left_ == nullptr)
        {
            grandparent->left_ = other_child;
        }
        else if (grandparent->right_ == nullptr)
        {
            grandparent->right_ = other_child;
        }
        else
        {
            throw std::runtime_error("Grandparent already has two children.");
        }

        // other child -> grandparent
        other_child->parent_ = grandparent;
    }
}

void RRSNode::node_update_vertex_box(const std::shared_ptr<Vertex>& boundary_vertex)
{
    double old_box_size;
    double new_box_size;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_node_);

        // current box size
        old_box_size = box_.max[0] - box_.min[0];
        new_box_size = boundary_vertex->get_max()[0] - boundary_vertex->get_min()[0];

        // skip if no change
        if (old_box_size == new_box_size) return;

        // store new box size
        box_.min = boundary_vertex->get_min();
        box_.max = boundary_vertex->get_max();
    }

    // recursive update
    if (new_box_size > old_box_size)
    {
        recursive_expand_parent_box();
    }
    else
    {
        recursive_shrink_parent_box();
    }
}

RRSReturnType RRSNode::node_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    // skip if not contained
    if (!box_.contains(generic_point->get_position()))
    {
        return RRSReturnType::SKIP;
    }

    // branch if not leaf
    if (!isLeaf_)
    {
        // search left and right
        RRSReturnType left_return = left_->node_reverse_radius_search(generic_point, search_results);
        RRSReturnType right_return = right_->node_reverse_radius_search(generic_point, search_results);

        // skip if both is skip
        if (left_return == RRSReturnType::SKIP && right_return == RRSReturnType::SKIP) return RRSReturnType::SKIP;
        // intersected if any is intersected
        return RRSReturnType::INTERSECTED;
    }
    else
    {        
        // store
        search_results.push_back(boundary_vertex_);

        // return
        return RRSReturnType::INTERSECTED;
    }
}

void RRSNode::node_flattern(std::vector<std::shared_ptr<Vertex>>& flatten_list)
{
    if (!isLeaf_)
    {
        left_->node_flattern(flatten_list);
        right_->node_flattern(flatten_list);
    }
    else
    {
        flatten_list.push_back(boundary_vertex_);
    }
}

void RRSNode::node_count(unsigned int& count) const
{
    count++;

    if (!isLeaf_)
    {
        left_->node_count(count);
        right_->node_count(count);
    }
}

std::vector<std::shared_ptr<Vertex>> RRSTree::compute_vertices_list()
{
    std::vector<std::shared_ptr<Vertex>> flatten_list;
    root->node_flattern(flatten_list);
    return flatten_list;
}

void RRSNode::node_print(int level) const
{
    if (!isLeaf_)
    {
        left_->node_print(level+1);
        right_->node_print(level+1);
    }
    else
    {
        std::cout << "Level: " <<  level << " | ID: " << boundary_vertex_->get_id() << " | Position: " << boundary_vertex_->get_position().transpose() << " | Radius: " << boundary_vertex_->get_radius() << std::endl;
        std::cout << std::endl;
    }
}

RRSTree::RRSTree() : tree_size(0), leaf_size(1)
{
    root = std::make_shared<RRSNode>();
}

void RRSTree::tree_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{   
    // get vertex's node
    std::shared_ptr<RRSNode> node = boundary_vertex->node;

    // skip if vertex already in tree
    if (node != nullptr) return;

    // find best node to add
    std::shared_ptr<RRSNode> best_node = find_best_node(root, boundary_vertex);

    // add to best node
    best_node->node_add_vertex(boundary_vertex);

    // increase size
    tree_size++;
}

void RRSTree::tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // get vertex's node
    std::shared_ptr<RRSNode> node = boundary_vertex->node;

    // skip if vertex not in tree
    if (node == nullptr) return;

    // remove from node
    node->node_delete_vertex(boundary_vertex);

    // decrease size
    tree_size--;
}

void RRSTree::tree_update_vertex_box(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // get vertex's node
    std::shared_ptr<RRSNode> node = boundary_vertex->node;

    // skip if vertex not in tree
    if (node == nullptr) return;

    // update node box
    node->node_update_vertex_box(boundary_vertex);
}

RRSReturnType RRSTree::tree_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    return root->node_reverse_radius_search(generic_point, search_results);
}

void RRSTree::tree_print() const
{
    root->node_print(0);
}

void RRSTree::print_size()
{
    std::cout << "Size: " << tree_size << std::endl;
}

unsigned int RRSTree::get_size() const
{
    return get_node_size();
}

unsigned int RRSTree::get_node_size() const
{
    unsigned int count = 0;
    root->node_count(count);

    return count;
}