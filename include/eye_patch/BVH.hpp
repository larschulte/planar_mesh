#include <vector>
#include <memory>
#include <limits>
#include <iostream>
#include <Eigen/Dense>

// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
bool rayTriangleIntersect(
    const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
    const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
    Eigen::Vector3d& outIntersection)
{
    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    // compute determinant
    Eigen::Vector3d pvec = dir.cross(edge2);
    double det = edge1.dot(pvec);
    if (det > -EPSILON && det < EPSILON) {
        return false; // This ray is parallel to this triangle.
    }

    // compute inverse determinant
    double invDet = 1.0 / det;

    // compute u
    Eigen::Vector3d tvec = orig - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) return false;
    
    // compute v
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = dir.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return false;

    // compute t
    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) return false; // This means that there is a line intersection but not a ray intersection.

    // comptue intersection
    outIntersection = orig + dir * t;
    return true;
}

// https://chatgpt.com/share/96c43118-6cb4-4549-9478-2725dee3b44d
struct AABB 
{
    Eigen::Vector3d min;
    Eigen::Vector3d max;

    AABB() : 
        min(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity())),
        max(Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity())) 
    {}

    void expand(const Eigen::Vector3d& point) 
    {
        min = min.cwiseMin(point);
        max = max.cwiseMax(point);
    }

    bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const 
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
};

struct BVHNode 
{
    AABB box;
    std::shared_ptr<BVHNode> left;
    std::shared_ptr<BVHNode> right;
    std::vector<int> triangleIndices;

    bool isLeaf() const 
    {
        return !left && !right;
    }
};

// the start and end are indices of the triangles list
std::shared_ptr<BVHNode> buildBVH(std::map<int, Eigen::Vector3d>& vertex_to_vector3d_map, std::map<int, std::array<int, 3>>& triangle_to_indices_map, std::vector<int>& triangle_id_list, int start, int end) 
{
    // create node
    auto node = std::make_shared<BVHNode>();
    
    // compute bounding box
    for (int i = start; i < end; ++i) 
    {
        std::array<int, 3>& indices = triangle_to_indices_map[triangle_id_list[i]];
        node->box.expand(vertex_to_vector3d_map[indices[0]]);
        node->box.expand(vertex_to_vector3d_map[indices[1]]);
        node->box.expand(vertex_to_vector3d_map[indices[2]]);

    }

    // create leaf node when we have less than 4 triangles
    if (end - start <= 4) 
    { 
        // reserve space
        node->triangleIndices.reserve(end - start);

        // store triangle id
        for (int i = start; i < end; ++i) 
        {
            node->triangleIndices.push_back(triangle_id_list[i]);
        }

        // return node
        return node;
    }

    // find current box longest axis
    Eigen::Vector3d extent = node->box.max - node->box.min;
    int axis = 0;
    if (extent[1] > extent[0]) axis = 1;
    if (extent[2] > extent[axis]) axis = 2;

    // sort triangle by its center along the longest axis
    // nth_element is a partial sort, it will put the mid element in the mid position, everything smaller on the left but not sorted, and everything bigger on the right but not sorted
    int mid = (start + end) / 2;
    std::nth_element(triangle_id_list.begin() + start, triangle_id_list.begin() + mid, triangle_id_list.begin() + end, 
        [&](const int& a, const int& b) 
        {
            double centerA = (vertex_to_vector3d_map[triangle_to_indices_map[a][0]][axis] + vertex_to_vector3d_map[triangle_to_indices_map[a][1]][axis] + vertex_to_vector3d_map[triangle_to_indices_map[a][2]][axis]) / 3.0;
            double centerB = (vertex_to_vector3d_map[triangle_to_indices_map[b][0]][axis] + vertex_to_vector3d_map[triangle_to_indices_map[b][1]][axis] + vertex_to_vector3d_map[triangle_to_indices_map[b][2]][axis]) / 3.0;
            return centerA < centerB;
        });

    // recursive build children nodes
    node->left = buildBVH(vertex_to_vector3d_map, triangle_to_indices_map, triangle_id_list, start, mid);
    node->right = buildBVH(vertex_to_vector3d_map, triangle_to_indices_map, triangle_id_list, mid, end);

    // return node
    return node;
}

// no const for map since map[i] may create new element if not present
bool intersectBVH(
    const std::shared_ptr<BVHNode>& node, 
    const Eigen::Vector3d& orig, 
    const Eigen::Vector3d& dir, 
    std::map<int, Eigen::Vector3d>& vertex_to_vector3d_map, 
    std::map<int, std::array<int, 3>>& triangle_to_indices_map, 
    std::vector<int>& outTriangleIdList,
    std::vector<Eigen::Vector3d>& outIntersectionList) 
{
    // skip if ray doesn't intersect bounding box
    double tMin = -std::numeric_limits<double>::infinity();
    double tMax = std::numeric_limits<double>::infinity();
    if (!node->box.intersect(orig, dir, tMin, tMax)) return false;

    // if not leaf, recursively check children
    if (!node->isLeaf())
    {
        std::vector<int> outTriangleIdListLeft;
        std::vector<int> outTriangleIdListRight;
        std::vector<Eigen::Vector3d> outIntersectionListLeft;
        std::vector<Eigen::Vector3d> outIntersectionListRight;

        bool hitLeft = intersectBVH(node->left, orig, dir, vertex_to_vector3d_map, triangle_to_indices_map, outTriangleIdListLeft, outIntersectionListLeft);
        bool hitRight = intersectBVH(node->right, orig, dir, vertex_to_vector3d_map, triangle_to_indices_map, outTriangleIdListRight, outIntersectionListRight);

        outTriangleIdList.insert(outTriangleIdList.end(), outTriangleIdListLeft.begin(), outTriangleIdListLeft.end());
        outTriangleIdList.insert(outTriangleIdList.end(), outTriangleIdListRight.begin(), outTriangleIdListRight.end());
        outIntersectionList.insert(outIntersectionList.end(), outIntersectionListLeft.begin(), outIntersectionListLeft.end());
        outIntersectionList.insert(outIntersectionList.end(), outIntersectionListRight.begin(), outIntersectionListRight.end());

        return hitLeft || hitRight;
    }

    // when leaf, check each triangle in node
    bool hasIntersection = false;
    for (int node_triangle : node->triangleIndices) 
    {
        // each triangle
        const std::array<int, 3> indices = triangle_to_indices_map[node_triangle];
        Eigen::Vector3d v0 = vertex_to_vector3d_map[indices[0]];
        Eigen::Vector3d v1 = vertex_to_vector3d_map[indices[1]];
        Eigen::Vector3d v2 = vertex_to_vector3d_map[indices[2]];

        // check intersection
        Eigen::Vector3d intersectPoint;
        bool intersected = rayTriangleIntersect(orig, dir, v0, v1, v2, intersectPoint);

        // return if intersected
        if (intersected) 
        {
            outTriangleIdList.push_back(node_triangle);
            outIntersectionList.push_back(intersectPoint);
            hasIntersection = true;
        }
    }
    return hasIntersection;
    
}

