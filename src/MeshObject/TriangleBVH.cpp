#include "MeshObject/TriangleBVH.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Surface.hpp" // Include the header file for the 'Surface' class
#include "MeshObject/GenericPoint.hpp" // Include the header file for the 'GenericPoint' class
#include <iostream>

Settings TriangleBVH::settings_;

std::ostream& operator<<(std::ostream& os, const BVHReturnType& type)
{
    switch (type)
    {
        case BVHReturnType::INTERSECTED:
            os << "I _ _";
            break;
        case BVHReturnType::SKIP:
            os << "_ S _";
            break;
        default:
            os << "? ? ?";
            break;
    }
    return os;
}

bool ray_triangle_intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
    const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
    Eigen::Vector3d& outIntersection)
{
    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    Eigen::Vector3d pvec = dir.cross(edge2);
    double det = edge1.dot(pvec);
    if (std::fabs(det) < EPSILON) {
        return false;
    }

    double invDet = 1.0 / det;

    Eigen::Vector3d tvec = orig - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) return false;
    
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = dir.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return false;

    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) return false;

    outIntersection = orig + dir * t;
    return true;
}

BoundingBox& BoundingBox::operator=(const BoundingBox& other)
{
    min = other.min;
    max = other.max;
    min_used_for_surface_area = other.min_used_for_surface_area;
    max_used_for_surface_area = other.max_used_for_surface_area;
    surface_area = other.surface_area;
    
    return *this;
}

BoundingBox::BoundingBox(const BoundingBox& other)
{
    // std::shared_lock lock_other(other.mutex_); // Acquire shared (read) lock for other object

    min = other.min;
    max = other.max;
    min_used_for_surface_area = other.min_used_for_surface_area;
    max_used_for_surface_area = other.max_used_for_surface_area;
    surface_area = other.surface_area;
}

BoundingBox::BoundingBox()
    : min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
      max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

BoundingBox::BoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max) : min(min), max(max) {}

bool BoundingBox::expand(const Eigen::Vector3d& point) 
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

void BoundingBox::expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max)
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

void BoundingBox::expand_box_no_return(const BoundingBox& box)
{
    // std::shared_lock lock(box.mutex_); // Acquire shared (read) lock
    expand_box_no_return(box.min, box.max);
}

bool BoundingBox::expand_box(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max)
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

bool BoundingBox::expand_box(const BoundingBox& box)
{
    // std::shared_lock lock(box.mutex_); // Acquire shared (read) lock
    return expand_box(box.min, box.max);
} 

bool BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& invDir, double& tMin, double& tMax) const 
{
    // std::shared_lock lock(mutex_); // Acquire shared (read) lock

    for (int i = 0; i < 3; ++i) 
    {
        double t0 = std::min((min[i] - orig[i]) * invDir[i], (max[i] - orig[i]) * invDir[i]);
        double t1 = std::max((min[i] - orig[i]) * invDir[i], (max[i] - orig[i]) * invDir[i]);
        if (t1 < tMin || t0 > tMax) return false;  // Early exit if miss
        tMin = std::max(tMin, t0);
        tMax = std::min(tMax, t1);
    }
    return true;
}

bool BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) const 
{
    double tMin = -std::numeric_limits<double>::infinity();
    double tMax = std::numeric_limits<double>::infinity();
    return intersect(orig, dir, tMin, tMax);
}

int BoundingBox::get_longest_axis()
{
    // std::shared_lock lock(mutex_); // Acquire shared (read) lock

    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

double BoundingBox::get_surface_area()
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

void Node::recursive_expand_parent_box()
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

void Node::recursive_shrink_parent_box()
{
    if (parent_)
    {
        // old parent box
        BoundingBox old_parent_box = parent_->box_;

        // new parent box
        BoundingBox new_parent_box = BoundingBox();

        // copy boxes of children of parent
        BoundingBox parent_left_box;
        BoundingBox parent_right_box;
        {
            std::shared_lock lock_left(parent_->left_->rwlock_box_, std::defer_lock);
            std::shared_lock lock_right(parent_->right_->rwlock_box_, std::defer_lock);
            std::lock(lock_left, lock_right);

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

BVHReturnType Node::node_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const
{       
    // skip if not intersected
    if (!box_.intersect(generic_point->get_position(), generic_point->get_inv_direction()))
    {
        return BVHReturnType::SKIP;
    }

    // branch if not leaf
    if (!isLeaf_)
    {
        // search left and right
        BVHReturnType left_return = left_->node_intersection_search(generic_point, faces_intersected);
        BVHReturnType right_return = right_->node_intersection_search(generic_point, faces_intersected);

        // skip if both is skip
        if (left_return == BVHReturnType::SKIP && right_return == BVHReturnType::SKIP) return BVHReturnType::SKIP;
        // intersected if any is intersected
        return BVHReturnType::INTERSECTED;
    }
    else
    {        
        // store
        faces_intersected.push_back(face_);

        // return
        return BVHReturnType::INTERSECTED;
    }
}

std::shared_ptr<Node> TriangleBVH::find_best_node(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face)
{
    // store a list of queue to process
    std::queue<std::pair<Node*, double>> queue;
    queue.push(std::make_pair(node.get(), 0));

    // initialize
    double best_cost = std::numeric_limits<double>::infinity();
    Node* best_node = nullptr;

    // compute lower bound cost to add branch node
    RRSBoundingBox smallest_branch_box(face->get_min(), face->get_max());
    const double lower_bound_cost = smallest_branch_box.get_surface_area();

    // while queue is not empty
    while (!queue.empty())
    {
        // get the first element
        std::pair<Node*, double> current = queue.front();
        queue.pop();

        // get the node and inherited cost
        Node* current_node = current.first;
        double inherited_cost = current.second;

        // skip if inherited cost is already greater than best cost
        if (inherited_cost > best_cost) continue;

        // cost to add to the node directly if the node is leaf and does not have face
        if (current_node->isLeaf_ && current_node->face_ == nullptr)
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
            BoundingBox new_branch_box = current_node->box_;
            new_branch_box.expand_box_no_return(face->get_min(), face->get_max());
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
            BoundingBox expanded_box = current_node->box_;
            expanded_box.expand_box_no_return(face->get_min(), face->get_max());
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
        return std::shared_ptr<Node>(node, best_node);
    }
    else
    {
        return nullptr;
    }
}

void Node::node_add_face(const std::shared_ptr<Face>& face)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_node_);

    // if node is leaf, and do not have face, simply add face to node
    if (isLeaf_ && face_ == nullptr)
    {
        face_ = face;
        face_->node = shared_from_this();
        box_.expand_box_no_return(face->get_min(), face->get_max());
        recursive_expand_parent_box();
        return;
    }
    
    // create new node
    std::shared_ptr<Node> new_leaf_node = std::make_shared<Node>();
    {
        new_leaf_node->box_.expand_box_no_return(face->get_min(), face->get_max());
        new_leaf_node->face_ = face;
        new_leaf_node->face_->node = new_leaf_node;
    }
    
    // create duplicate node
    std::shared_ptr<Node> duplicate_node = std::make_shared<Node>();
    {
        duplicate_node->box_ = box_;
        
        // if node is leaf, copy boundary vertices
        if (isLeaf_)
        {
            duplicate_node->face_ = face_;
            if (duplicate_node->face_) duplicate_node->face_->node = duplicate_node;
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

        // assign sibling
        duplicate_node->sibling_ = sibling_;
        new_leaf_node->sibling_ = duplicate_node;

        // change to branch
        isLeaf_ = false;
    }
}

void Node::node_delete_face(const std::shared_ptr<Face>& face)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_node_);

    // face is only stored in leaf node

    // remove node from face
    face->node = nullptr;

    // remove face from node
    face_ = nullptr;

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

void Node::node_print(int level) const
{
    if (!isLeaf_)
    {
        left_->node_print(level+1);
        right_->node_print(level+1);
    }
    else
    {
        std::cout << "Level: " <<  level << " | ID: " << face_->get_id() << " | Center: " << face_->get_center().transpose() << std::endl;
        std::cout << std::endl;
    }
}

void Node::node_flatten(std::vector<std::shared_ptr<Face>>& face_list) const
{
    if (!isLeaf_)
    {
        left_->node_flatten(face_list);
        right_->node_flatten(face_list);
    }
    else
    {
        face_list.push_back(face_);
    }
}

std::vector<std::shared_ptr<Face>> TriangleBVH::get_face_list() const
{
    std::vector<std::shared_ptr<Face>> face_list;
    root->node_flatten(face_list);
    return face_list;
}

void TriangleBVH::tree_delete_face(std::shared_ptr<Face> face)
{    
    // get face's node reference
    std::shared_ptr<Node> node = face->node;

    // skip if face not in tree
    if (node == nullptr) return; 

    // delete face from node
    node->node_delete_face(face);

    // decrement face size
    face_size--;
}

TriangleBVH::TriangleBVH() :
        face_size(0),
        leaf_size(1)
{
    root = std::make_shared<Node>();
}

void TriangleBVH::tree_add_face(std::shared_ptr<Face> face)
{
    // get face's node reference
    std::shared_ptr<Node> node = face->node;

    // skip if face already in tree
    if (node != nullptr) return;

    // find best node to add
    std::shared_ptr<Node> best_node = find_best_node(root, face);

    // add to best node
    best_node->node_add_face(face);

    // increment face size
    face_size++;
}

BVHReturnType TriangleBVH::tree_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const
{
    return root->node_intersection_search(generic_point, faces_intersected);
}

void TriangleBVH::tree_print() const
{
    root->node_print(0);
}

unsigned int TriangleBVH::get_size() const
{
    return face_size;
}