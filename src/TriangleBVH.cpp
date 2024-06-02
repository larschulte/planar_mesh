#include "eye_patch/TriangleBVH.hpp"

bool rayTriangleIntersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
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

Eigen::Vector3d TriangleBVH::compute_triangle_center(int triangleID)
{
    int point0_id = triangle_to_indices_map.at(triangleID)[0];
    int point1_id = triangle_to_indices_map.at(triangleID)[1];
    int point2_id = triangle_to_indices_map.at(triangleID)[2];
    Eigen::Vector3d v0 = point_to_vector3d_map.at(point0_id);
    Eigen::Vector3d v1 = point_to_vector3d_map.at(point1_id);
    Eigen::Vector3d v2 = point_to_vector3d_map.at(point2_id);
    return (v0 + v1 + v2) / 3.0;
}

double TriangleBVH::sort_triangle_list_in_axis(std::vector<Triangle>& triangle_list, int axis, int start, int mid, int end)
{
    std::nth_element(triangle_list.begin() + start, triangle_list.begin() + mid, triangle_list.begin() + end, 
        [&](const Triangle& triangle_a, const Triangle& triangle_b) 
        {
            return triangle_a.center[axis] < triangle_b.center[axis];
        });
    return triangle_list[mid].center[axis];
}

void TriangleBVH::expand_node_box(std::shared_ptr<Node> node, Triangle triangle)
{
    node->box.expand(triangle.v0);
    node->box.expand(triangle.v1);
    node->box.expand(triangle.v2);
}

std::set<int> TriangleBVH::intersectHierarchy(const std::shared_ptr<Node>& node, const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) 
{
    bool intersected = node->box.intersect(orig, dir);
    if (!intersected) return std::set<int>();
    
    if (node->isLeaf())
    {
        std::set<int> intersected_triangleIDs;
        for (Triangle triangle : node->triangles)
        {
            Eigen::Vector3d v0 = triangle.v0;
            Eigen::Vector3d v1 = triangle.v1;
            Eigen::Vector3d v2 = triangle.v2;
            Eigen::Vector3d intersection;
            if (rayTriangleIntersect(orig, dir, v0, v1, v2, intersection)) intersected_triangleIDs.insert(triangle.triangleID);
        }
        return intersected_triangleIDs;
    }
    else
    {
        std::set<int> triangles_left = intersectHierarchy(node->left, orig, dir);
        std::set<int> triangles_right = intersectHierarchy(node->right, orig, dir);
        triangles_left.insert(triangles_right.begin(), triangles_right.end());
        return triangles_left;
    }
}

void TriangleBVH::convert_leaf_to_branch(std::shared_ptr<Node> node)
{
    std::vector<Triangle> triangle_list(node->triangles.begin(), node->triangles.end());
    int start = 0;
    int end = triangle_list.size();
    int mid = (start + end) / 2;
    int axis = node->box.get_longest_axis();
    double split_value = sort_triangle_list_in_axis(triangle_list, axis, start, mid, end);
    
    node->triangles.clear();
    node->split_axis = axis;
    node->split_value = split_value;
    node->left = build_node(std::vector<Triangle>(triangle_list.begin(), triangle_list.begin() + mid));
    node->right = build_node(std::vector<Triangle>(triangle_list.begin() + mid, triangle_list.end()));
}

std::shared_ptr<TriangleBVH::Node> TriangleBVH::build_node(std::vector<Triangle> triangle_list)
{
    auto node = std::make_shared<Node>();
    for (Triangle triangle : triangle_list) 
    {
        expand_node_box(node, triangle);
        node->triangles.insert(triangle);
    }

    if (node->triangles.size() > 4) convert_leaf_to_branch(node);

    return node;
}

void TriangleBVH::addTriangleToNode(std::shared_ptr<Node> node, Triangle triangle)
{
    if (node->isLeaf())
    {
        expand_node_box(node, triangle);
        node->triangles.insert(triangle);
        if (node->triangles.size() > 4) convert_leaf_to_branch(node);
    }
    else
    {
        expand_node_box(node, triangle);

        if (triangle.center[node->split_axis] < node->split_value)
        {
            addTriangleToNode(node->left, triangle);
        }
        else 
        {
            addTriangleToNode(node->right, triangle);
        }
    }
}

void TriangleBVH::deleteTriangleFromNode(std::shared_ptr<Node> node, Triangle triangle)
{
    if (node->isLeaf())
    {
        node->triangles.erase(triangle);
    }
    else
    {
        if (triangle.center[node->split_axis] < node->split_value)
        {
            deleteTriangleFromNode(node->left, triangle);
        }
        else
        {
            deleteTriangleFromNode(node->right, triangle);
        }
    }
}

void TriangleBVH::deleteTriangle(int triangleID, std::array<int, 3> indices, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
{
    Triangle triangle;
    triangle.triangleID = triangleID;
    triangle.v0 = v0;
    triangle.v1 = v1;
    triangle.v2 = v2;
    triangle.center = (v0 + v1 + v2) / 3.0;

    // delete from triangle list using erase-remove idiom
    triangle_list.erase(std::remove(triangle_list.begin(), triangle_list.end(), triangle), triangle_list.end());
    
    // delete from BVH
    deleteTriangleFromNode(root, triangle);
}

TriangleBVH::TriangleBVH()
    : rebuild_threshold(2),
      size_at_last_rebuild(0)
{
    rebuild();
}

void TriangleBVH::addData(std::vector<int> _triangle_list, std::map<int, std::array<int, 3>> _triangle_to_indices_map, std::map<int, Eigen::Vector3d> _point_to_vector3d_map)
{
    for (int triangleID : _triangle_list)
    {
        Triangle triangle;
        triangle.triangleID = triangleID;
        triangle.v0 = _point_to_vector3d_map.at(_triangle_to_indices_map.at(triangleID)[0]);
        triangle.v1 = _point_to_vector3d_map.at(_triangle_to_indices_map.at(triangleID)[1]);
        triangle.v2 = _point_to_vector3d_map.at(_triangle_to_indices_map.at(triangleID)[2]);
        triangle.center = (triangle.v0 + triangle.v1 + triangle.v2) / 3.0;
        triangle_list.push_back(triangle);
    }
}

void TriangleBVH::rebuild()
{
    root = build_node(triangle_list);
}

void TriangleBVH::addTriangle(int triangleID, std::array<int, 3> indices, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
{
    Triangle triangle;
    triangle.triangleID = triangleID;
    triangle.v0 = v0;
    triangle.v1 = v1;
    triangle.v2 = v2;
    triangle.center = (v0 + v1 + v2) / 3.0;

    triangle_list.push_back(triangle);

    if (triangle_list.size() > size_at_last_rebuild * rebuild_threshold)
    {
        rebuild();
        size_at_last_rebuild = triangle_list.size();
    }
    else
    {
        addTriangleToNode(root, triangle);
    }
}

std::set<int> TriangleBVH::intersectionSearch(Eigen::Vector3d origin, Eigen::Vector3d endPoint)
{
    Eigen::Vector3d dir = (endPoint - origin).normalized();
    return intersectHierarchy(root, origin, dir);
}
