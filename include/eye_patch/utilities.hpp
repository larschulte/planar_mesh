#include <tuple>
#include <Eigen/Dense>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>


std::tuple<int, int, int> valueToJet(float value) 
{
    // Ensure value is within [0, 1]
    if (value < 0.0f) value = 0.0f;
    if (value > 1.0f) value = 1.0f;

    float r = 0, g = 0, b = 0;

    if (value < 0.25f) 
    {
        r = 0;
        g = 4 * value;
        b = 1;
    } 
    else if (value < 0.5f) 
    {
        r = 0;
        g = 1;
        b = 1 - 4 * (value - 0.25f);
    } 
    else if (value < 0.75f) 
    {
        r = 4 * (value - 0.5f);
        g = 1;
        b = 0;
    } 
    else 
    {
        r = 1;
        g = 1 - 4 * (value - 0.75f);
        b = 0;
    }

    return std::make_tuple(static_cast<int>(r * 255), static_cast<int>(g * 255), static_cast<int>(b * 255));
}

Eigen::Vector3d merge_means(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2) 
{
    return (size1 * mean1 + size2 * mean2) / (size1 + size2);
}

Eigen::Matrix3d merge_covariances(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                int size1, int size2) 
{
    Eigen::Vector3d combined_mean = merge_means(mean1, mean2, size1, size2);
    Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
    Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
    Eigen::Matrix3d combined_covariance = (size1 * cov1 + size2 * cov2 + size1 * mean_diff1 + size2 * mean_diff2) / (size1 + size2);

    return combined_covariance;
}

std::array<int, 3> sortThreeInts(int a, int b, int c) {
    std::array<int, 3> sorted = {a, b, c};
    if (sorted[0] > sorted[1]) std::swap(sorted[0], sorted[1]);
    if (sorted[0] > sorted[2]) std::swap(sorted[0], sorted[2]);
    if (sorted[1] > sorted[2]) std::swap(sorted[1], sorted[2]);
    return sorted;
}

std::array<int, 2> sortTwoInts(int a, int b) {
    std::array<int, 2> sorted = {a, b};
    if (sorted[0] > sorted[1]) std::swap(sorted[0], sorted[1]);
    return sorted;
}


// Function to find the orientation of the ordered triplet (p, q, r).
int orientation(const Eigen::Vector2d &p, const Eigen::Vector2d &q, const Eigen::Vector2d &r) {
    double val = (q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y());
    if (val == 0) return 0;           // collinear
    return (val > 0) ? 1 : 2;         // clock or counterclockwise
}

// Function to check if point q lies on segment pr excluding endpoints.
bool onSegment(const Eigen::Vector2d &p, const Eigen::Vector2d &q, const Eigen::Vector2d &r) {
    if (q.x() < std::max(p.x(), r.x()) && q.x() > std::min(p.x(), r.x()) &&
        q.y() < std::max(p.y(), r.y()) && q.y() > std::min(p.y(), r.y()))
        return true;
    return false;
}

// Function to check if two segments (p1, q1) and (p2, q2) intersect.
bool doIntersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &q1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q2) {
    // Find the four orientations needed for the general and special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return false;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return false;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return false;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return false;

    return false; // Doesn't fall in any of the above cases
}

// ray plane intersection
Eigen::Vector3d ray_plane_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayEndPoint, Eigen::Vector3d planeMean, Eigen::Vector3d planeNormal)
{   
    // if perpendicular, return NaN
    Eigen::Vector3d rayDirection = (rayEndPoint - rayOrigin).normalized();
    if (planeNormal.dot(rayDirection) == 0)
    {
        // terminate programe by throwing an error
        throw std::invalid_argument("Ray and plane are perpendicular");
    }

    // compute intersection
    double distance = (planeMean - rayOrigin).dot(planeNormal) / rayDirection.dot(planeNormal);
    Eigen::Vector3d intersection = rayOrigin + rayDirection * distance;

    // return
    return intersection;
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
{
    // transform point cloud
    typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
    pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
    return transformed_cloud;
}

// point_in_triangle - generated by copilot
bool point_in_triangle(Eigen::Vector2d point, Eigen::Vector2d v0, Eigen::Vector2d v1, Eigen::Vector2d v2) 
{
    Eigen::Vector2d v0v1 = v1 - v0;
    Eigen::Vector2d v0v2 = v2 - v0;
    Eigen::Vector2d v0p = point - v0;

    double dot00 = v0v2.dot(v0v2);
    double dot01 = v0v2.dot(v0v1);
    double dot02 = v0v2.dot(v0p);
    double dot11 = v0v1.dot(v0v1);
    double dot12 = v0v1.dot(v0p);

    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Return true if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v <= 1);
}
