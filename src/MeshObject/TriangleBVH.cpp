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
        case BVHReturnType::ABORT:
            os << "_ _ A";
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

BoundingBox::BoundingBox()
    : min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
      max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

BoundingBox::BoundingBox(const Eigen::Vector3d& min, const Eigen::Vector3d& max) : min(min), max(max) {}

bool BoundingBox::expand(const Eigen::Vector3d& point) 
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

void BoundingBox::expand_box_no_return(const Eigen::Vector3d& input_min, const Eigen::Vector3d& input_max)
{
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
    return expand_box(box.min, box.max);
} 

bool BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& invDir, double& tMin, double& tMax) const 
{
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
    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

double BoundingBox::compute_surface_area() const
{
    Eigen::Vector3d dimensions = max - min;
    double area = 2.0 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[0]);
    if (std::isnan(area)) throw std::runtime_error("BVHBoundingBox::compute_surface_area() returned nan."); // throw if nan
    return area;
}

const double& BoundingBox::get_surface_area()
{
    // if min and max are not updated, return the stored value
    if (min == min_used_for_surface_area && max == max_used_for_surface_area) return surface_area;

    // else update and return
    min_used_for_surface_area = min;
    max_used_for_surface_area = max;
    surface_area = compute_surface_area();
    return surface_area;
}

void Node::recursive_unlock()
{
    omp_unset_nest_lock(&omp_lock);
    if (locked_children)
    {
        if (!left || !right) throw std::runtime_error("Node has locked children but one of them is null.");
        left->recursive_unlock();
        right->recursive_unlock();
        locked_children = false;
    }
}

void Node::recursive_expand_parent_box()
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

void Node::recursive_shrink_parent_box()
{
    if (parent)
    {
        // old parent box
        BoundingBox old_parent_box = parent->box;

        // new parent box
        BoundingBox new_parent_box = BoundingBox();
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

double TriangleBVH::sort_face_list_in_axis(std::vector<std::shared_ptr<Face>>& face_list, int axis, int start, int mid, int end)
{
    std::sort(face_list.begin() + start, face_list.begin() + end, 
        [&](const std::shared_ptr<Face>& triangle_a, const std::shared_ptr<Face>& triangle_b) 
        {
            return triangle_a->get_first_vertex()->get_position()[axis] < triangle_b->get_first_vertex()->get_position()[axis];
        });
    return face_list[mid]->get_first_vertex()->get_position()[axis];
}

void TriangleBVH::sort_face_list_in_axis(std::vector<std::shared_ptr<Face>>& face_list, int axis, int start, int end)
{
    std::sort(face_list.begin() + start, face_list.begin() + end, 
        [&](const std::shared_ptr<Face>& triangle_a, const std::shared_ptr<Face>& triangle_b) 
        {
            return triangle_a->get_first_vertex()->get_position()[axis] < triangle_b->get_first_vertex()->get_position()[axis];
        });
}

BVHReturnType TriangleBVH::node_intersection_search(const std::shared_ptr<Node>& node, const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const
{    
    // lock the node for exclusive read/write (abort if can't lock node -> someone else is reading/writing this node)
    if (!omp_test_nest_lock(&node->omp_lock))
    {
        // std::cout << "_ X waiting ..." << std::endl;
        return BVHReturnType::ABORT;
    }
    
    // skip if not intersected
    if (!node->box.intersect(generic_point->get_position(), generic_point->get_inv_direction()))
    {
        omp_unset_nest_lock(&node->omp_lock);
        return BVHReturnType::SKIP;
    }
    
    // branch if not leaf
    if (!node->isLeaf)
    {
        // get copy of left and right node
        std::shared_ptr<Node> left = node->left;
        std::shared_ptr<Node> right = node->right;
        // unlock node
        omp_unset_nest_lock(&node->omp_lock);

        // search left and right
        BVHReturnType left_return = node_intersection_search(left, generic_point, faces_intersected);
        if (left_return == BVHReturnType::ABORT) return BVHReturnType::ABORT;
        BVHReturnType right_return = node_intersection_search(right, generic_point, faces_intersected);
        if (right_return == BVHReturnType::ABORT) return BVHReturnType::ABORT;

        // skip if both is skip
        if (left_return == BVHReturnType::SKIP && right_return == BVHReturnType::SKIP) return BVHReturnType::SKIP;
        // intersected if any is intersected
        return BVHReturnType::INTERSECTED;
    }
    else
    {
        // skip if no faces
        if (node->faces.size() == 0)
        {
            omp_unset_nest_lock(&node->omp_lock);
            return BVHReturnType::SKIP;
        }

        // abort if can't lock face lock
        const std::shared_ptr<Face>& face = node->faces[0];
        if (!omp_test_nest_lock(&face->face_lock))
        {
            omp_unset_nest_lock(&node->omp_lock);
            return BVHReturnType::ABORT;
        }

        // skip if face is not searchable
        if (!node->faces[0]->is_searchable())
        {
            omp_unset_nest_lock(&face->face_lock);
            omp_unset_nest_lock(&node->omp_lock);
            return BVHReturnType::SKIP;
        }

        // skip if not intersected
        if (!node->faces[0]->intersects_point(generic_point->get_origin(), generic_point->get_direction()))
        {
            omp_unset_nest_lock(&face->face_lock);
            omp_unset_nest_lock(&node->omp_lock);
            return BVHReturnType::SKIP;
        }

        // abort if can't lock face's surface
        const std::shared_ptr<Surface>& surface = face->get_surface();
        generic_point->intersected_surfaces.insert(surface);
        if (!omp_test_nest_lock(&surface->lock)) // nest lock here since a ray could intersect two faces of the same surface in two nodes
        {
            generic_point->contented_surfaces[surface]++;
            omp_unset_nest_lock(&face->face_lock);
            omp_unset_nest_lock(&node->omp_lock);
            // std::cout << "_ _ X _ _ _ _ _" << std::endl;
            return BVHReturnType::ABORT;
        }

        // return
        faces_intersected.push_back(node->faces[0]);

        omp_unset_nest_lock(&face->face_lock);
        return BVHReturnType::INTERSECTED;
    }
}

void TriangleBVH::convert_leaf_to_branch(const std::shared_ptr<Node>& node)
{
    int start = 0;
    int end = node->faces.size();
    int split_axis;     // value to be computed
    double split_value; // value to be computed
    int split_index;    // value to be computed

    // fill in value to be computed
    const bool use_sah = (node->faces.size() > 2) && !settings_.turn_off_sah;
    if (!use_sah)
    {
        // use simple median split
        split_index = (start + end) / 2;
        split_axis = node->box.get_longest_axis();
        split_value = sort_face_list_in_axis(node->faces, split_axis, start, split_index, end);
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
            sort_face_list_in_axis(node->faces, axis, start, end);

            // Recompute the min_suffix and max_suffix after sorting
            std::vector<Eigen::Vector3d> min_suffix(end - start);
            std::vector<Eigen::Vector3d> max_suffix(end - start);

            // Initialize suffix arrays with the last element
            min_suffix[end - start - 1] = node->faces[end - 1]->get_min();
            max_suffix[end - start - 1] = node->faces[end - 1]->get_max();

            // Compute suffix bounding boxes for the right side
            for (int j = end - 2; j >= start; --j)
            {
                min_suffix[j - start] = min_suffix[j - start + 1].cwiseMin(node->faces[j]->get_min());
                max_suffix[j - start] = max_suffix[j - start + 1].cwiseMax(node->faces[j]->get_max());
            }

            // Initialize the left bounding box
            BoundingBox left_box;
            Eigen::Vector3d min_left = node->faces[start]->get_min();
            Eigen::Vector3d max_left = node->faces[start]->get_max();

            // Iterate through potential split points and evaluate SAH cost
            for (int i = start + 1; i < end; i++)
            {
                // Update left bounding box incrementally
                min_left = min_left.cwiseMin(node->faces[i - 1]->get_min());
                max_left = max_left.cwiseMax(node->faces[i - 1]->get_max());
                left_box.expand_box_no_return(min_left, max_left);

                // Use precomputed right bounding box from suffix arrays
                BoundingBox right_box;
                right_box.expand_box_no_return(min_suffix[i - start], max_suffix[i - start]);

                // Calculate SAH cost for this split
                double sah_cost = calculate_sah(node->box, left_box, right_box, i - start, end - i);

                if (sah_cost < min_sah_cost)
                {
                    min_sah_cost = sah_cost;
                    best_split_axis = axis;
                    best_split_value = node->faces[i]->get_first_vertex()->get_position()[axis];
                    best_split_index = i;
                }
            }
        }
        sort_face_list_in_axis(node->faces, best_split_axis, start, end);

        // fill in computed value
        split_axis = best_split_axis;
        split_value = best_split_value;
        split_index = best_split_index;
    }

    // Create left and right child nodes
    node->split_axis = split_axis;
    node->split_value = split_value;
    node->left = build_node(node->faces, start, split_index);
    node->right = build_node(node->faces, split_index, end);
    node->left->parent = node;
    node->right->parent = node;
    node->left->sibling = node->right;
    node->right->sibling = node->left;
    node->faces.clear();

    // Locking as before
    node->locked_children = true;
    node->isLeaf = false;

    // return
    return;
}

std::shared_ptr<Node> TriangleBVH::build_node(const std::vector<std::shared_ptr<Face>>& face_list, const int& start, const int& end)
{
    auto node = std::make_shared<Node>();
    
    // lock before creating link between face and node
    while (!omp_test_nest_lock(&node->omp_lock))
    {
        std::cout << "BVH lock node inside build node waiting ..." << std::endl;
    };

    // expand box
    for (int i = start; i < end; i++)
    {
        node->box.expand_box_no_return(face_list[i]->get_min(), face_list[i]->get_max());
    }

    // store faces    
    node->faces = std::vector<std::shared_ptr<Face>>(face_list.begin() + start, face_list.begin() + end);

    // store node pointer in faces
    for (const std::shared_ptr<Face>& face : node->faces)
    {
        face->node = node;
    }

    // convert to branch
    if (node->faces.size() > leaf_size) convert_leaf_to_branch(node);

    return node;
}

std::shared_ptr<Node> TriangleBVH::find_best_node(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face)
{
    if (!node->isLeaf)
    {
        if (face->get_first_vertex()->get_position()[node->split_axis] < node->split_value)
        {
            return find_best_node(node->left, face);
        }
        else
        {
            return find_best_node(node->right, face);
        }
    }
    else
    {
        return node;
    }
}

void TriangleBVH::node_add_face(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face)
{
    // lock before adding vertex
    while (!omp_test_nest_lock(&node->omp_lock))
    {
        std::cout << "BVH lock node inside add face waiting ... " << std::endl;
    };

    // after locking the node
    
    // create new node
    std::shared_ptr<Node> new_node = std::make_shared<Node>();
    {
        new_node->box.expand_box_no_return(face->get_min(), face->get_max());
        new_node->faces.push_back(face);
        face->node = new_node;
    }
    
    // create duplicate node
    std::shared_ptr<Node> duplicate_node = std::make_shared<Node>();
    {
        duplicate_node->box = node->box;
        duplicate_node->split_value = node->split_value;
        duplicate_node->split_axis = node->split_axis;
        
        // if node is leaf, copy boundary vertices
        if (node->isLeaf)
        {
            duplicate_node->faces = node->faces;
            for (const std::shared_ptr<Face>& sub_face : duplicate_node->faces)
            {
                sub_face->node = duplicate_node;
            }
        }
        else
        {
            // else, copy children
            duplicate_node->left = node->left;
            node->left->parent = duplicate_node;

            duplicate_node->right = node->right;
            node->right->parent = duplicate_node;

            duplicate_node->isLeaf.store(node->isLeaf.load());
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
        node->split_value = face->get_first_vertex()->get_position()[node->split_axis];

        // put into children
        node->left = duplicate_node;
        duplicate_node->parent = node;
        node->right = new_node;
        new_node->parent = node;

        // assign sibling
        duplicate_node->sibling = node->sibling;
        new_node->sibling = duplicate_node;

        // change to branch
        node->isLeaf = false;
    }
    
    // unlock
    node->recursive_unlock();
}

double TriangleBVH::calculate_sah(BoundingBox& parent_box, BoundingBox& left_box, BoundingBox& right_box, int left_count, int right_count)
{
    double S_parent = parent_box.get_surface_area();  // Surface area of parent node
    double S_left = left_box.get_surface_area();     // Surface area of left child
    double S_right = right_box.get_surface_area();   // Surface area of right child

    // The cost of traversing the node, typically set to a constant, e.g., 1.0
    const double traversal_cost = 1.0;

    // SAH cost formula
    return traversal_cost + (S_left / S_parent) * left_count + (S_right / S_parent) * right_count;
}

bool TriangleBVH::node_delete_face(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face)
{
    if (!node->isLeaf)
    {
        if (face->get_first_vertex()->get_position()[node->split_axis] < node->split_value)
        {
            return node_delete_face(node->left, face);
        }
        else if (face->get_first_vertex()->get_position()[node->split_axis] > node->split_value)
        {
            return node_delete_face(node->right, face);
        }
        else
        {
            return node_delete_face(node->left, face) || node_delete_face(node->right, face);
        }
    }
    else
    {
        auto it = std::remove(node->faces.begin(), node->faces.end(), face);
        if (it != node->faces.end())
        {
            node->faces.erase(it, node->faces.end());
            face->node = nullptr;
            return true;
        }
        else
        {
            return false;
        }
    }
}

void TriangleBVH::node_print(const std::shared_ptr<Node> &node, int level) const
{
    if (!node->isLeaf)
    {
        node_print(node->left, level+1);
        node_print(node->right, level+1);
    }
    else
    {
        for (const std::shared_ptr<Face>& face : node->faces)
        {
            std::cout << "Level: " <<  level << " | ID: " << face->get_id() << " | Center: " << face->get_center().transpose() << std::endl;
        }
        std::cout << std::endl;
    }
}

void TriangleBVH::node_flatten(const std::shared_ptr<Node>& node, std::vector<std::shared_ptr<Face>>& face_list) const
{
    if (!node->isLeaf)
    {
        node_flatten(node->left, face_list);
        node_flatten(node->right, face_list);
    }
    else
    {
        face_list.insert(face_list.end(), node->faces.begin(), node->faces.end());
    }
}

std::vector<std::shared_ptr<Face>> TriangleBVH::get_face_list() const
{
    std::vector<std::shared_ptr<Face>> face_list;
    node_flatten(root, face_list);
    return face_list;
}

void TriangleBVH::tree_delete_face(std::shared_ptr<Face> face)
{
    // check input
    // if (face->is_expired()) throw std::runtime_error("Attempts to delete expired face.");

    // decrement face size
    face_size--;
    
    // get face's node reference
    const std::shared_ptr<Node>& node = face->node;
    if (node == nullptr) throw std::invalid_argument("Vertex not found in BVH.");

    // lock node
    while (!omp_test_nest_lock(&node->omp_lock)) 
    {
        std::cout << "delete face waiting for " << std::endl;
    }

    // make copy for later release
    const std::shared_ptr<Node> locked_node = node;

    // throw if not found in node->faces
    const bool found = std::find(node->faces.begin(), node->faces.end(), face) != node->faces.end();
    if (!found) throw std::invalid_argument("Face not found in BVH.");

    node->box = BoundingBox();
    node->recursive_shrink_parent_box();

    // delete from node
    node->faces.erase(std::remove(node->faces.begin(), node->faces.end(), face), node->faces.end());
    face->node = nullptr;

    // release lock
    omp_unset_nest_lock(&locked_node->omp_lock);
}

TriangleBVH::TriangleBVH()
        :rebuild_threshold(2),
        size_at_last_rebuild(0),
        face_size(0),
        leaf_size(1)
{
    rebuild();
}

void TriangleBVH::check_rebuild()
{
    if (face_size > size_at_last_rebuild * rebuild_threshold)
    {
        std::cout << "Rebuilding BVH ...." << std::endl;
        rebuild();
        std::cout << "Rebuilding BVH done" << std::endl;
        size_at_last_rebuild = face_size;
    }
}

void TriangleBVH::rebuild()
{
    if (face_size == 0)
    {
        root = build_node(std::vector<std::shared_ptr<Face>>(), 0, 0);
    }
    else
    {
        std::vector<std::shared_ptr<Face>> face_list = get_face_list();
        root = build_node(face_list, 0, face_list.size());
    }

    // release lock
    root->recursive_unlock();
}

void TriangleBVH::tree_add_face(std::shared_ptr<Face> face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to add expired face.");

    // increment face size
    face_size++;

    // find best node to add
    std::shared_ptr<Node> best_node = find_best_node(root, face);

    // add to best node
    node_add_face(best_node, face);
}

BVHReturnType TriangleBVH::tree_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& faces_intersected) const
{
    return node_intersection_search(root, generic_point, faces_intersected);
}

void TriangleBVH::tree_print() const
{
    node_print(root, 0);
}

unsigned int TriangleBVH::get_size() const
{
    return face_size;
}