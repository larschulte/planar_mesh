#include "MeshObject/TriangleBVH.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Surface.hpp" // Include the header file for the 'Surface' class
#include <iostream>

std::ostream& operator<<(std::ostream& os, const BVHReturnType& type)
{
    switch (type)
    {
        case BVHReturnType::INTERSECTED:
            os << "_";
            break;
        case BVHReturnType::SKIP:
            os << "-";
            break;
        case BVHReturnType::ABORT:
            os << "X";
            break;
        default:
            os << "?";
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

bool BoundingBox::expand_box(const BoundingBox& box) 
{
    Eigen::Vector3d oldMin = min;
    Eigen::Vector3d oldMax = max;

    // Update min and max to include the new box
    min = min.cwiseMin(box.min);
    max = max.cwiseMax(box.max);

    // Check if min or max changed
    bool changed = (min != oldMin || max != oldMax);

    return changed;
}

bool BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const 
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

void Node::recursive_unlock()
{
    omp_unset_nest_lock(&omp_lock);
    if (locked_children)
    {
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
        new_parent_box.expand_box(parent->left->box);
        new_parent_box.expand_box(parent->right->box);
                
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

void TriangleBVH::expand_node_box(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face)
{
    node->box.expand(face->get_vertex(0)->get_position());
    node->box.expand(face->get_vertex(1)->get_position());
    node->box.expand(face->get_vertex(2)->get_position());
}

BVHReturnType TriangleBVH::node_intersection_search(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, std::vector<std::shared_ptr<Face>>& faces_intersected, std::set<std::shared_ptr<Surface>, MeshObjectCompare>& intersected_surfaces) const
{    
    // skip if not intersected
    if (!node->box.intersect(orig, dir))
    {
        return BVHReturnType::SKIP;
    }
    
    // branch if not leaf
    if (!node->isLeaf)
    {
        // search left and right
        BVHReturnType left_return = node_intersection_search(node->left, orig, dir, faces_intersected, intersected_surfaces);
        if (left_return == BVHReturnType::ABORT) return BVHReturnType::ABORT;
        BVHReturnType right_return = node_intersection_search(node->right, orig, dir, faces_intersected, intersected_surfaces);
        if (right_return == BVHReturnType::ABORT) return BVHReturnType::ABORT;

        // skip if both is skip
        if (left_return == BVHReturnType::SKIP && right_return == BVHReturnType::SKIP) return BVHReturnType::SKIP;
        // intersected if any is intersected
        return BVHReturnType::INTERSECTED;
    }
    else
    {
        // abort if can't lock node
        if (!omp_test_nest_lock(&node->omp_lock))
        {
            return BVHReturnType::ABORT;
        }

        // skip if no faces
        if (node->faces.size() == 0)
        {
            omp_unset_nest_lock(&node->omp_lock);
            return BVHReturnType::SKIP;
        }

        // skip if not intersected
        if (!node->faces[0]->intersects_point(orig, dir))
        {
            omp_unset_nest_lock(&node->omp_lock);
            return BVHReturnType::SKIP;
        }

        // abort if can't lock face's surface
        const std::shared_ptr<Face>& face = node->faces[0];
        const std::shared_ptr<Surface>& surface = face->get_surface();
        intersected_surfaces.insert(surface);
        if (!omp_test_nest_lock(&surface->lock)) // nest lock here since a ray could intersect two faces of the same surface in two nodes
        {
            omp_unset_nest_lock(&node->omp_lock);
            // std::cout << "_ _ X _ _ _ _ _" << std::endl;
            return BVHReturnType::ABORT;
        }

        // return
        faces_intersected.push_back(node->faces[0]);
        return BVHReturnType::INTERSECTED;
    }
}

BVHReturnType TriangleBVH::node_find_leaf_node(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& endPoint, std::shared_ptr<Node>& return_node)
{
    // branch if not leaf
    if (!node->isLeaf)
    {   
        // release lock if not leaf
        if (endPoint[node->split_axis] < node->split_value)
        {
            return node_find_leaf_node(node->left, orig, endPoint, return_node);
        }
        else 
        {
            return node_find_leaf_node(node->right, orig, endPoint, return_node);
        }
    }
    else
    {
        // abort if can't lock node
        if (!omp_test_nest_lock(&node->omp_lock))
        {
            // std::cout << "X _ _ _ _ _ _ _" << std::endl;
            return BVHReturnType::ABORT;
        }
        
        // return the leaf node locked if no face
        const bool no_face = node->faces.size() == 0;
        if (no_face)
        {
            return_node = node;
            return BVHReturnType::INTERSECTED;
        }
        
        // if there is face, return new empty leaf node

        // temperary face
        std::shared_ptr<Face> temp_face = std::make_shared<Face>(); 
        temp_face->temp_initialize(endPoint);

        //  add and branch
        node->faces.push_back(temp_face);
        convert_leaf_to_branch(node); // this may cause error in Application when locking surface's node

        // delete
        return_node = endPoint[node->split_axis] < node->split_value ? node->left : node->right;
        return_node->faces.pop_back();

        // locks
        while (!omp_test_nest_lock(&return_node->omp_lock))
        {
            std::cout << "BVH lock return_node waiting ..." << std::endl;
        };
        node->recursive_unlock();

        // return
        return BVHReturnType::INTERSECTED;
    }
}

void TriangleBVH::convert_leaf_to_branch(const std::shared_ptr<Node>& node)
{
    int start = 0;
    int end = node->faces.size();
    int mid = (start + end) / 2;
    int axis = node->box.get_longest_axis();
    double split_value = sort_face_list_in_axis(node->faces, axis, start, mid, end);

    node->split_axis = axis;
    node->split_value = split_value;
    node->left = build_node(node->faces, start, mid);
    node->right = build_node(node->faces, mid, end);
    node->left->parent = node;
    node->right->parent = node;
    node->left->sibling = node->right;
    node->right->sibling = node->left;
    node->faces.clear();

    // set lock before setting isLeaf to false to prevent other threads from accessing the children node
    while (!omp_test_nest_lock(&node->left->omp_lock))
    {
        std::cout << "BVH lock left children waiting ..." << std::endl;
    };
    while (!omp_test_nest_lock(&node->right->omp_lock))
    {
        std::cout << "BVH lock right children waiting ..." << std::endl;
    };
    node->locked_children = true;

    // update isLeaf after locking the children
    node->isLeaf = false;
}

std::shared_ptr<Node> TriangleBVH::build_node(const std::vector<std::shared_ptr<Face>>& face_list, const int& start, const int& end)
{
    auto node = std::make_shared<Node>();

    // expand box
    for (int i = start; i < end; i++)
    {
        expand_node_box(node, face_list[i]);
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

void TriangleBVH::node_add_face(const std::shared_ptr<Node>& node, const std::shared_ptr<Face>& face)
{
    expand_node_box(node, face);

    if (!node->isLeaf)
    {    
        if (face->get_first_vertex()->get_position()[node->split_axis] < node->split_value)
        {
            node_add_face(node->left, face);
        }
        else 
        {
            node_add_face(node->right, face);
        }
    }
    else
    {
        node->faces.push_back(face);

        // store node pointer in face
        face->node = node;

        if (node->faces.size() > leaf_size) convert_leaf_to_branch(node);
    }
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
    if (face->is_expired()) throw std::runtime_error("Attempts to delete expired face.");

    // decrement face size
    face_size--;
    
    // get face's node reference
    const std::shared_ptr<Node>& node = face->node;
    if (node == nullptr) throw std::invalid_argument("Vertex not found in BVH.");

    // delete from node
    node->faces.erase(std::remove(node->faces.begin(), node->faces.end(), face), node->faces.end());
    face->node = nullptr;
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
    if (root->left) root->left->recursive_unlock();
    if (root->right) root->right->recursive_unlock();
}

void TriangleBVH::tree_add_face(std::shared_ptr<Face> face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to add expired face.");

    // increment face size
    face_size++;

    node_add_face(root, face);
}

BVHReturnType TriangleBVH::tree_intersection_search(Eigen::Vector3d origin, Eigen::Vector3d endPoint, std::vector<std::shared_ptr<Face>>& faces_intersected, std::set<std::shared_ptr<Surface>, MeshObjectCompare>& intersected_surfaces) const
{
    Eigen::Vector3d dir = (endPoint - origin).normalized();
    return node_intersection_search(root, origin, dir, faces_intersected, intersected_surfaces);
}

BVHReturnType TriangleBVH::tree_find_leaf_node(const Eigen::Vector3d& origin, const Eigen::Vector3d& endPoint, std::shared_ptr<Node>& return_node)
{
    return node_find_leaf_node(root, origin, endPoint, return_node);
}

void TriangleBVH::tree_print() const
{
    node_print(root, 0);
}

unsigned int TriangleBVH::get_size() const
{
    return face_size;
}