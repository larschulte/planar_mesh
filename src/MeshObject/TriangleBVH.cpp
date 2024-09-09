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

void BoundingBox::expand(const Eigen::Vector3d& point) 
{
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);
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

bool Node::isLeaf() const 
{
    return !left && !right;
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

BVHReturnType TriangleBVH::node_intersection_search(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, std::vector<std::shared_ptr<Face>>& faces_intersected) const
{    
    // abort if can't lock node
    if (!omp_test_nested_lock_with_log(node->lock, "node lock")) 
    {
        // abort message
        std::cout << "Can't lock node" << std::endl;
        return BVHReturnType::ABORT;
    }

    // skip if not intersected
    if (!node->box.intersect(orig, dir))
    {
        omp_unset_nested_lock_with_log(node->lock, "node lock (node not intersected)");
        return BVHReturnType::SKIP;
    }
    
    // branch if not leaf
    if (!node->isLeaf())
    {
        // release lock if not leaf
        omp_unset_nested_lock_with_log(node->lock, "node lock (node not leaf)");

        // search left and right
        BVHReturnType left_return = node_intersection_search(node->left, orig, dir, faces_intersected);
        BVHReturnType right_return = node_intersection_search(node->right, orig, dir, faces_intersected);

        // abort if any is abort
        if (left_return == BVHReturnType::ABORT || right_return == BVHReturnType::ABORT) return BVHReturnType::ABORT;
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
            omp_unset_nested_lock_with_log(node->lock, "node lock (no faces)");
            return BVHReturnType::SKIP;
        }

        // skip if not intersected
        if (!node->faces[0]->intersects_point(orig, dir))
        {
            omp_unset_nested_lock_with_log(node->lock, "node lock (face not intersected)");
            return BVHReturnType::SKIP;
        }

        // abort if can't lock face's surface
        if (!omp_test_nested_lock_with_log(node->faces[0]->get_surface()->lock, "face's surface lock")) // nest lock here since a ray could intersect two faces of the same surface in two nodes
        {
            omp_unset_nested_lock_with_log(node->lock, "node lock (can't lock face's surface)");

            // abort message
            std::cout << "Can't lock face's surface" << std::endl;
            return BVHReturnType::ABORT;
        }
        
        // hold both the node lock and the face's surface lock
        faces_intersected.push_back(node->faces[0]);
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
    node->faces.clear();
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

    if (!node->isLeaf())
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
    if (!node->isLeaf())
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
    if (!node->isLeaf())
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
    if (!node->isLeaf())
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
    
    // delete from BVH
    if (!node_delete_face(root, face)) throw std::runtime_error("Face not found in BVH.");
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
}

void TriangleBVH::tree_add_face(std::shared_ptr<Face> face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to add expired face.");

    // increment face size
    face_size++;

    node_add_face(root, face);
}

BVHReturnType TriangleBVH::tree_intersection_search(Eigen::Vector3d origin, Eigen::Vector3d endPoint, std::vector<std::shared_ptr<Face>>& faces_intersected) const
{
    Eigen::Vector3d dir = (endPoint - origin).normalized();
    return node_intersection_search(root, origin, dir, faces_intersected);
}

void TriangleBVH::tree_print() const
{
    node_print(root, 0);
}
