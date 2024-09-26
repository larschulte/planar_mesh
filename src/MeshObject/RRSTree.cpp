#include "MeshObject/RRSTree.hpp"
#include "MeshObject/Vertex.hpp"

#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"

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
        case RRSReturnType::ABORT:
            os << "_ _ A";
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
    if (std::isnan(area)) throw std::runtime_error("RRSBoundingBox::compute_surface_area() returned nan."); // throw if nan
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

void RRSNode::recursive_unlock()
{
    omp_unset_nest_lock(&omp_lock);
    if (locked_children)
    {
        if (!left || !right) throw std::runtime_error("RRSNode has locked children but one of them is null.");
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

double RRSTree::sort_boundary_vertex_list_in_axis(std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, int axis, int start, int mid, int end)
{
    std::sort(boundary_vertex_list.begin() + start, boundary_vertex_list.begin() + end, 
        [&](const std::shared_ptr<Vertex>& boundary_vertex_a, const std::shared_ptr<Vertex>& boundary_vertex_b) 
        {
            return boundary_vertex_a->get_position()[axis] < boundary_vertex_b->get_position()[axis];
        });
    return boundary_vertex_list[mid]->get_position()[axis];
}

void RRSTree::sort_boundary_vertex_list_in_axis(std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, int axis, int start, int end)
{
    std::sort(boundary_vertex_list.begin() + start, boundary_vertex_list.begin() + end, 
        [&](const std::shared_ptr<Vertex>& boundary_vertex_a, const std::shared_ptr<Vertex>& boundary_vertex_b) 
        {
            return boundary_vertex_a->get_position()[axis] < boundary_vertex_b->get_position()[axis];
        });
}

void RRSTree::convert_leaf_to_branch(const std::shared_ptr<RRSNode>& node)
{
    int start = 0;
    int end = node->boundary_vertices.size();
    int split_axis;     // value to be computed
    double split_value; // value to be computed
    int split_index;    // value to be computed

    // fill in value to be computed
    const bool use_sah = node->boundary_vertices.size() > 2;
    if (!use_sah)
    {
        // use simple median split
        split_index = (start + end) / 2;
        split_axis = node->box.get_longest_axis();
        split_value = sort_boundary_vertex_list_in_axis(node->boundary_vertices, split_axis, start, split_index, end);
    }
    else
    {
        // use SAH
        
        // loop through all axes to find the optimal split axis and value
        double min_sah_cost = std::numeric_limits<double>::infinity();
        int best_split_axis = -1;
        double best_split_value = 0;
        int best_split_index = -1;
        for (int axis = 0; axis < 3; axis++)  // Assuming 3 axes (x, y, z)
        {
            // Sort boundary vertices along this axis
            sort_boundary_vertex_list_in_axis(node->boundary_vertices, axis, start, end);

            // Recompute the min_suffix and max_suffix after sorting
            std::vector<Eigen::Vector3d> min_suffix(end - start);
            std::vector<Eigen::Vector3d> max_suffix(end - start);

            // Initialize suffix arrays with the last element
            min_suffix[end - start - 1] = node->boundary_vertices[end - 1]->get_min();
            max_suffix[end - start - 1] = node->boundary_vertices[end - 1]->get_max();

            // Compute suffix bounding boxes for the right side
            for (int j = end - 2; j >= start; --j)
            {
                min_suffix[j - start] = min_suffix[j - start + 1].cwiseMin(node->boundary_vertices[j]->get_min());
                max_suffix[j - start] = max_suffix[j - start + 1].cwiseMax(node->boundary_vertices[j]->get_max());
            }

            // Initialize the left bounding box
            RRSBoundingBox left_box;
            Eigen::Vector3d min_left = node->boundary_vertices[start]->get_min();
            Eigen::Vector3d max_left = node->boundary_vertices[start]->get_max();

            // Iterate through potential split points and evaluate SAH cost
            for (int i = start + 1; i < end; i++)
            {
                // Update left bounding box incrementally
                min_left = min_left.cwiseMin(node->boundary_vertices[i - 1]->get_min());
                max_left = max_left.cwiseMax(node->boundary_vertices[i - 1]->get_max());
                left_box.expand_box_no_return(min_left, max_left);

                // Use precomputed right bounding box from suffix arrays
                RRSBoundingBox right_box;
                right_box.expand_box_no_return(min_suffix[i - start], max_suffix[i - start]);

                // Calculate SAH cost for this split
                double sah_cost = calculate_sah(node->box, left_box, right_box, i - start, end - i);

                if (sah_cost < min_sah_cost)
                {
                    min_sah_cost = sah_cost;
                    best_split_axis = axis;
                    best_split_value = node->boundary_vertices[i]->get_position()[axis];
                    best_split_index = i;
                }
            }
        }
        sort_boundary_vertex_list_in_axis(node->boundary_vertices, best_split_axis, start, end);

        // fill in computed value
        split_axis = best_split_axis;
        split_value = best_split_value;
        split_index = best_split_index;
    }

    // Create left and right child nodes
    node->split_axis = split_axis;
    node->split_value = split_value;
    node->left = build_node(node->boundary_vertices, start, split_index);
    node->right = build_node(node->boundary_vertices, split_index, end);
    node->left->parent = node;
    node->right->parent = node;
    node->left->sibling = node->right;
    node->right->sibling = node->left;
    node->boundary_vertices.clear();

    // Locking as before
    node->locked_children = true;
    node->isLeaf = false;

    // return
    return;
}

double RRSTree::calculate_sah(RRSBoundingBox& parent_box, RRSBoundingBox& left_box, RRSBoundingBox& right_box, int left_count, int right_count)
{
    double S_parent = parent_box.get_surface_area();  // Surface area of parent node
    double S_left = left_box.get_surface_area();     // Surface area of left child
    double S_right = right_box.get_surface_area();   // Surface area of right child

    // The cost of traversing the node, typically set to a constant, e.g., 1.0
    const double traversal_cost = 1.0;

    // SAH cost formula
    return traversal_cost + (S_left / S_parent) * left_count + (S_right / S_parent) * right_count;
}

std::shared_ptr<RRSNode> RRSTree::build_node(const std::vector<std::shared_ptr<Vertex>>& boundary_vertex_list, const int& start, const int& end)
{
    auto node = std::make_shared<RRSNode>();

    // lock before creating link between vertex and node
    while (!omp_test_nest_lock(&node->omp_lock))
    {
        std::cout << "RRS lock node inside build node waiting ..." << std::endl;
    };

    // expand box
    for (int i = start; i < end; i++)
    {
        node->box.expand_box_no_return(boundary_vertex_list[i]->get_min(), boundary_vertex_list[i]->get_max());
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
void RRSTree::node_add_vertex(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    node->box.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());

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
        // lock before adding vertex
        while (!omp_test_nest_lock(&node->omp_lock))
        {
            std::cout << "RRS lock node inside add vertex waiting ... " << std::endl;
        };

        // node could become branch while waiting (when other thread add to the same node), hence need a new node_add_vertex call
        if (!node->isLeaf)
        {
            omp_unset_nest_lock(&node->omp_lock);
            node_add_vertex(node, boundary_vertex);
            return;
        }

        // connect to this node
        node->boundary_vertices.push_back(boundary_vertex);
        boundary_vertex->node = node;

        // convert to branch
        if (node->boundary_vertices.size() > leaf_size) convert_leaf_to_branch(node);

        // unlock
        node->recursive_unlock();        
    }
}

void RRSTree::node_increase_radius(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<Vertex>& boundary_vertex)
{
    node->box.expand_box_no_return(boundary_vertex->get_min(), boundary_vertex->get_max());

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

RRSReturnType RRSTree::node_reverse_radius_search(const std::shared_ptr<RRSNode>& node, const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    // skip if not contained
    if (!node->box.contains(generic_point->get_position()))
    {
        return RRSReturnType::SKIP;
    }

    // branch if not leaf
    if (!node->isLeaf)
    {
        // search left and right
        RRSReturnType left_return = node_reverse_radius_search(node->left, generic_point, search_results);
        if (left_return == RRSReturnType::ABORT) return RRSReturnType::ABORT;
        RRSReturnType right_return = node_reverse_radius_search(node->right, generic_point, search_results);
        if (right_return == RRSReturnType::ABORT) return RRSReturnType::ABORT;

        // skip if both is skip
        if (left_return == RRSReturnType::SKIP && right_return == RRSReturnType::SKIP) return RRSReturnType::SKIP;
        // intersected if any is intersected
        return RRSReturnType::INTERSECTED;
    }
    else
    {
        // abort if can't lock node
        // node lock prevents vertex from being added or deleted from the node
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
        const std::shared_ptr<Vertex>& boundary_vertex = node->boundary_vertices[0];

        // abort if can't lock vertex
        // vertex lock prevent vertex's properties from being changed when we are accessing
        if (!omp_test_nest_lock(&boundary_vertex->vertex_lock))
        {
            omp_unset_nest_lock(&node->omp_lock);
            return RRSReturnType::ABORT;
        }
        
        // skip if boundary vertex is not searchable ( this mean the vertex should be removed from node but has not been removed yet)
        if (!node->boundary_vertices[0]->is_searchable())
        {
            omp_unset_nest_lock(&boundary_vertex->vertex_lock);
            omp_unset_nest_lock(&node->omp_lock);
            return RRSReturnType::SKIP;
        }

        // skip if not contained
        if (!node->boundary_vertices[0]->approx_contains(generic_point->get_position()))
        {
            omp_unset_nest_lock(&boundary_vertex->vertex_lock);
            omp_unset_nest_lock(&node->omp_lock);
            return RRSReturnType::SKIP;
        }

        // abort if can't lock vertex's surface
        const std::shared_ptr<Surface>& surface = boundary_vertex->get_surface_check();
        generic_point->intersected_surfaces.insert(surface);
        if (!omp_test_nest_lock(&surface->lock)) 
        {
            generic_point->contented_surfaces[surface]++;
            omp_unset_nest_lock(&boundary_vertex->vertex_lock);
            omp_unset_nest_lock(&node->omp_lock);
            // std::cout << "_ _ _ _ _ _ X _" << std::endl;
            return RRSReturnType::ABORT;
        }

        // return
        search_results.push_back(node->boundary_vertices[0]);

        omp_unset_nest_lock(&boundary_vertex->vertex_lock);
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
        node->box.expand_box_no_return(temp_vertex->get_min(), temp_vertex->get_max());
        node->boundary_vertices.push_back(temp_vertex);
        convert_leaf_to_branch(node); // this may cause error in Application when locking surface's node

        // delete
        return_node = point[node->split_axis] < node->split_value ? node->left : node->right;

        // throw if return_node is not the same as temp_vertex->node
        if (return_node != temp_vertex->node) throw std::runtime_error("RRSTree::node_find_leaf_node() returned wrong node.");

        return_node->boundary_vertices.pop_back();

        // locks
        while (!omp_test_nest_lock(&return_node->omp_lock))
        {
            std::cout << "RRS lock return_node waiting ..." << std::endl;
        };
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

    // release lock
    root->recursive_unlock();
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
    if (node == nullptr) throw std::invalid_argument("node is null.");

    // lock the node
    while (!omp_test_nest_lock(&node->omp_lock)) // don't copy node, as node can change when other thread is adding point to the node
    {
        std::cout << "delete vertex waiting ..." << std::endl;
    }

    // make copy for later release
    const std::shared_ptr<RRSNode> locked_node = node;

    // throw if not found in node->boundary_vertices
    const bool found = std::find(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex) != node->boundary_vertices.end();
    if (!found) throw std::invalid_argument("Vertex not found in node's boundary_vertices.");

    // shrink bounding box
    node->box = RRSBoundingBox();
    node->recursive_shrink_parent_box();

    // delete from node
    node->boundary_vertices.erase(std::remove(node->boundary_vertices.begin(), node->boundary_vertices.end(), boundary_vertex), node->boundary_vertices.end());
    boundary_vertex->node = nullptr;

    // release
    omp_unset_nest_lock(&locked_node->omp_lock);
}

RRSReturnType RRSTree::tree_reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& search_results)
{
    return node_reverse_radius_search(root, generic_point, search_results);
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