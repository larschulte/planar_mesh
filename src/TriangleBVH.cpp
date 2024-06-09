#include "eye_patch/TriangleBVH.hpp"
#include <iostream>

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

TriangleBVH::BoundingBox::BoundingBox()
    : min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
      max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) {}

void TriangleBVH::BoundingBox::expand(const Eigen::Vector3d& point) 
{
    min = min.cwiseMin(point);
    max = max.cwiseMax(point);
}

bool TriangleBVH::BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const 
{
    for (int i = 0; i < 3; ++i) 
    {
        double invD = 1.0 / dir[i];
        double t0 = (min[i] - orig[i]) * invD;
        double t1 = (max[i] - orig[i]) * invD;
        if (invD < 0.0) std::swap(t0, t1);
        tMin = std::max(tMin, t0);
        tMax = std::min(tMax, t1);
        if (tMax <= tMin) return false;
    }
    return true;
}

bool TriangleBVH::BoundingBox::intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) const 
{
    double tMin = -std::numeric_limits<double>::infinity();
    double tMax = std::numeric_limits<double>::infinity();
    return intersect(orig, dir, tMin, tMax);
}

int TriangleBVH::BoundingBox::get_longest_axis()
{
    Eigen::Vector3d diagonal_line = max - min;
    int axis = 0;
    if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
    if (diagonal_line[2] > diagonal_line[axis]) axis = 2;
    return axis;
}

bool TriangleBVH::Node::isLeaf() const 
{
    return !left && !right;
}


double TriangleBVH::sort_face_list_in_axis(std::vector<std::weak_ptr<Face>>& face_list, int axis, int start, int mid, int end)
{
    std::nth_element(face_list.begin() + start, face_list.begin() + mid, face_list.begin() + end, 
        [&](const std::weak_ptr<Face>& triangle_a, const std::weak_ptr<Face>& triangle_b) 
        {
            return triangle_a.lock()->get_center()[axis] < triangle_b.lock()->get_center()[axis];
        });
    return face_list[mid].lock()->get_center()[axis];
}

void TriangleBVH::expand_node_box(std::shared_ptr<Node> node, std::weak_ptr<Face> face)
{
    for (std::weak_ptr<Vertex> vertex : face.lock()->get_vertices())
    {
        node->box.expand(vertex.lock()->get_position());
    }
}

std::set<std::weak_ptr<Face>> TriangleBVH::intersectHierarchy(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) 
{
    bool intersected = node->box.intersect(orig, dir);
    if (!intersected) return std::set<std::weak_ptr<Face>>();
    
    if (node->isLeaf())
    {
        std::set<std::weak_ptr<Face>> intersected_faces;
        for (std::weak_ptr<Face> face : node->faces)
        {
            if (face.lock()->intersects_point(orig, dir)) intersected_faces.insert(face);
        }
        return intersected_faces;
    }
    else
    {
        std::set<std::weak_ptr<Face>> faces_left = intersectHierarchy(node->left, orig, dir);
        std::set<std::weak_ptr<Face>> faces_right = intersectHierarchy(node->right, orig, dir);
        faces_left.insert(faces_right.begin(), faces_right.end());
        return faces_left;
    }
}

void TriangleBVH::convert_leaf_to_branch(std::shared_ptr<Node> node)
{
    std::vector<std::weak_ptr<Face>> face_list(node->faces.begin(), node->faces.end());
    int start = 0;
    int end = face_list.size();
    int mid = (start + end) / 2;
    int axis = node->box.get_longest_axis();
    double split_value = sort_face_list_in_axis(face_list, axis, start, mid, end);
    
    node->faces.clear();
    node->split_axis = axis;
    node->split_value = split_value;
    node->left = build_node(std::vector<std::weak_ptr<Face>>(face_list.begin(), face_list.begin() + mid));
    node->right = build_node(std::vector<std::weak_ptr<Face>>(face_list.begin() + mid, face_list.end()));
}

std::shared_ptr<TriangleBVH::Node> TriangleBVH::build_node(std::vector<std::weak_ptr<Face>> face_list)
{
    auto node = std::make_shared<Node>();
    for (std::weak_ptr<Face> face : face_list) 
    {
        expand_node_box(node, face);
        node->faces.insert(face);
    }

    if (node->faces.size() > 4) convert_leaf_to_branch(node);

    return node;
}

void TriangleBVH::add_face_to_node(std::shared_ptr<Node> node, std::weak_ptr<Face> face)
{
    if (node->isLeaf())
    {
        expand_node_box(node, face);
        node->faces.insert(face);
        if (node->faces.size() > 4) convert_leaf_to_branch(node);
    }
    else
    {
        expand_node_box(node, face);

        if (face.lock()->get_center()[node->split_axis] < node->split_value)
        {
            add_face_to_node(node->left, face);
        }
        else 
        {
            add_face_to_node(node->right, face);
        }
    }
}

void TriangleBVH::delete_face_from_node(std::shared_ptr<Node> node, std::weak_ptr<Face> face)
{
    if (node->isLeaf())
    {
        node->faces.erase(face);
    }
    else
    {
        if (face.lock()->get_center()[node->split_axis] < node->split_value)
        {
            delete_face_from_node(node->left, face);
        }
        else
        {
            delete_face_from_node(node->right, face);
        }
    }
}

void TriangleBVH::print_node(const std::shared_ptr<Node> &node, int level) const
{
    if (node->isLeaf())
    {
        for (std::weak_ptr<Face> face : node->faces)
        {
            std::cout << "Level: " <<  level << " | ID: " << face.lock()->get_id() << " | Center: " << face.lock()->get_center().transpose() << std::endl;
        }
        std::cout << std::endl;
    }
    else
    {
        print_node(node->left, level+1);
        print_node(node->right, level+1);
    }
}

void TriangleBVH::delete_face(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to delete expired face.");

    // check if face exists
    if (face_set.find(face) == face_set.end()) return;

    // delete from face set
    face_set.erase(face);

    // delete from face list using erase-remove idiom
    face_list.erase(std::remove(face_list.begin(), face_list.end(), face), face_list.end());
    
    // delete from BVH
    delete_face_from_node(root, face);
}

TriangleBVH::TriangleBVH()
    : rebuild_threshold(2),
      size_at_last_rebuild(0)
{
    rebuild();
}

void TriangleBVH::rebuild()
{
    root = build_node(face_list);
}

void TriangleBVH::add_face(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to add expired face.");

    // check if face already exists
    if (face_set.find(face) != face_set.end()) return;

    // add to face set
    face_set.insert(face);

    // add to face list
    face_list.push_back(face);

    if (face_list.size() > size_at_last_rebuild * rebuild_threshold)
    {
        rebuild();
        size_at_last_rebuild = face_list.size();
    }
    else
    {
        add_face_to_node(root, face);
    }
}

std::set<std::weak_ptr<Face>> TriangleBVH::intersectionSearch(Eigen::Vector3d origin, Eigen::Vector3d endPoint)
{
    Eigen::Vector3d dir = (endPoint - origin).normalized();
    return intersectHierarchy(root, origin, dir);
}

void TriangleBVH::print() const
{
    print_node(root, 0);
}
