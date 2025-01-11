#include "MeshObject/RRSTree.hpp"
#include "MeshObject/Vertex.hpp"

#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"

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

void RRSBoundingBox::expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max)
{
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
    return expand_box(box.min, box.max);
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

double RRSBoundingBox::compute_surface_area() const
{
    Eigen::Vector3d dimensions = max - min;
    double area = 2.0 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[0]);
    if (std::isnan(area)) 
    {
        // this means the box is being deleted, thus should have zero area.
        return 0;
    }
    return area;
}

const double& RRSBoundingBox::get_surface_area()
{
    // if min and max are not updated, return the stored value
    if (min == min_used_for_surface_area && max == max_used_for_surface_area) return surface_area;

    // else update and return
    min_used_for_surface_area = min;
    max_used_for_surface_area = max;
    surface_area = compute_surface_area();
    return surface_area;
}

void RRSNode::recursive_expand_parent_box()
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

void RRSNode::recursive_shrink_parent_box()
{
    if (parent)
    {
        // old parent box
        RRSBoundingBox old_parent_box = parent->box;

        // new parent box
        RRSBoundingBox new_parent_box = RRSBoundingBox();
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
        if (inherited_cost + lower_bound_cost > best_cost) continue;

        // cost to add a branch that contains current node and leaf node
        {
            // cost of creating a new branch node
            RRSBoundingBox new_branch_box = current_node->box;
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
        if (!current_node->isLeaf)
        {
            // compute change to inherited cost
            RRSBoundingBox expanded_box = current_node->box;
            expanded_box.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());
            double change_to_inherited = expanded_box.get_surface_area() - current_node->box.get_surface_area();

            // lower-bound cost
            double cost = inherited_cost + change_to_inherited + lower_bound_cost;
            
            // update best cost and best node
            if (cost < best_cost)
            {
                // we should check the children
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
    // after locking the node
    
    // create new node
    std::shared_ptr<RRSNode> new_node = std::make_shared<RRSNode>();
    {
        new_node->box.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());
        new_node->boundary_vertices.push_back(boundary_vertex);
        boundary_vertex->node = new_node;
    }
    
    // create duplicate node
    std::shared_ptr<RRSNode> duplicate_node = std::make_shared<RRSNode>();
    {
        duplicate_node->box = box;
        duplicate_node->split_value = split_value;
        duplicate_node->split_axis = split_axis;
        
        // if node is leaf, copy boundary vertices
        if (isLeaf)
        {
            duplicate_node->boundary_vertices = boundary_vertices;
            for (const std::shared_ptr<Vertex>& vertex : duplicate_node->boundary_vertices)
            {
                vertex->node = duplicate_node;
            }
        }
        else
        {
            // else, copy children
            duplicate_node->left = left;
            left->parent = duplicate_node;

            duplicate_node->right = right;
            right->parent = duplicate_node;

            duplicate_node->isLeaf.store(isLeaf.load());
        }
    }
    
    // make new node and duplicate node children of the current node
    {
        // expand current node box
        box.expand_box_no_return(new_node->box);
        recursive_expand_parent_box();

        // get split axis
        split_axis = box.get_longest_axis();

        // get split value  
        split_value = boundary_vertex->get_position()[split_axis];

        // put into children
        left = duplicate_node;
        duplicate_node->parent = shared_from_this();
        right = new_node;
        new_node->parent = shared_from_this();

        // assign sibling
        duplicate_node->sibling = sibling;
        new_node->sibling = duplicate_node;

        // change to branch
        isLeaf = false;
    }    
}

bool RRSNode::node_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    if (!isLeaf)
    {
        if (boundary_vertex->get_position()[split_axis] < split_value)
        {
            return left->node_delete_vertex(boundary_vertex);
        }
        else if (boundary_vertex->get_position()[split_axis] > split_value)
        {
            return right->node_delete_vertex(boundary_vertex);
        }
        else
        {
            return left->node_delete_vertex(boundary_vertex) || right->node_delete_vertex(boundary_vertex);
        }
    }
    else
    {
        auto it = std::remove(boundary_vertices.begin(), boundary_vertices.end(), boundary_vertex);
        if (it != boundary_vertices.end())
        {
            boundary_vertices.erase(it, boundary_vertices.end());
            boundary_vertex->node = nullptr;
            return true;
        }
        else
        {
            return false;
        }
    }
}

RRSReturnType RRSNode::node_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    // skip if not contained
    if (!box.contains(generic_point->get_position()))
    {
        return RRSReturnType::SKIP;
    }

    // branch if not leaf
    if (!isLeaf)
    {
        // search left and right
        RRSReturnType left_return = left->node_reverse_radius_search(generic_point, search_results);
        RRSReturnType right_return = right->node_reverse_radius_search(generic_point, search_results);

        // skip if both is skip
        if (left_return == RRSReturnType::SKIP && right_return == RRSReturnType::SKIP) return RRSReturnType::SKIP;
        // intersected if any is intersected
        return RRSReturnType::INTERSECTED;
    }
    else
    {
        // skip if no vertices
        if (boundary_vertices.size() == 0)
        {
            return RRSReturnType::SKIP;
        }
        const std::shared_ptr<Vertex>& boundary_vertex = boundary_vertices[0];

        // Double-Checked Locking
        if (boundary_vertex->is_expired() || !boundary_vertex->approx_contains(generic_point->get_position()))
        {
            return RRSReturnType::SKIP;
        }
        
        // we lock the vertex during vertex->delete_() call, and then release it before locking it again when removing it from the rrs tree
        // thus this thread may have lock this vertex during the gap
        // and see an expired vertex in the search tree
        // thus we need to check if the vertex is expired and skip if it is
        if (boundary_vertex->is_expired() || !boundary_vertex->approx_contains(generic_point->get_position()))
        {
            return RRSReturnType::SKIP;
        }

        const std::shared_ptr<Surface>& surface = boundary_vertex->get_surface_check();
        generic_point->intersected_surfaces.insert(surface);

        // return
        search_results.push_back(boundary_vertex);

        return RRSReturnType::INTERSECTED;
    }
}

void RRSNode::node_flattern(std::vector<std::shared_ptr<Vertex>>& flatten_list)
{
    if (!isLeaf)
    {
        left->node_flattern(flatten_list);
        right->node_flattern(flatten_list);
    }
    else
    {
        flatten_list.insert(flatten_list.end(), boundary_vertices.begin(), boundary_vertices.end());
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
    if (!isLeaf)
    {
        left->node_print(level+1);
        right->node_print(level+1);
    }
    else
    {
        for (const std::shared_ptr<Vertex>& boundary_vertex : boundary_vertices)
        {
            std::cout << "Level: " <<  level << " | ID: " << boundary_vertex->get_id() << " | Position: " << boundary_vertex->get_position().transpose() << " | Radius: " << boundary_vertex->get_radius() << std::endl;
        }
        std::cout << std::endl;
    }
}

RRSTree::RRSTree() : tree_size(0), leaf_size(1)
{
    root = std::make_shared<RRSNode>();
}

void RRSTree::tree_add_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{   
    // check input
    if (boundary_vertex->is_expired()) throw std::invalid_argument("Invalid vertex in boundary_vertex_list");

    // increase size
    tree_size++;

    // find best node to add
    std::shared_ptr<RRSNode> best_node = find_best_node(root, boundary_vertex);

    // add to best node
    best_node->node_add_vertex(boundary_vertex);
}

void RRSTree::tree_delete_vertex(const std::shared_ptr<Vertex>& boundary_vertex)
{
    // check input
    // if (boundary_vertex->is_expired()) throw std::invalid_argument("Invalid vertex in boundary_vertex_list");

    // decrease size
    tree_size--;

    // get vertex's node reference
    const std::shared_ptr<RRSNode>& node = boundary_vertex->node;
    if (node == nullptr) throw std::invalid_argument("node is null.");

    // make copy for later release
    const std::shared_ptr<RRSNode> locked_node = node;

    // throw if not found in node->boundary_vertices
    const bool found = std::find(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex) != node->boundary_vertices.end();
    if (!found) throw std::invalid_argument("Vertex not found in node's boundary_vertices.");

    // shrink bounding box
    node->box = RRSBoundingBox(); // here
    node->recursive_shrink_parent_box();

    // delete from node
    node->boundary_vertices.erase(std::remove(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex), node->boundary_vertices.end());
    boundary_vertex->node = nullptr;
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
    return tree_size;
}