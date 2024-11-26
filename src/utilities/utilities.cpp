#include "utilities/utilities.hpp"

std::tuple<int, int, int> valueToJet(double value) 
{
    // Ensure value is within [0, 1]
    if (value < 0.0) value = 0.0;
    if (value > 1.0) value = 1.0;

    double r = 0.0, g = 0.0, b = 0.0;

    if (value < 0.25) 
    {
        r = 0.0;
        g = 4.0 * value;
        b = 1.0;
    } 
    else if (value < 0.5) 
    {
        r = 0.0;
        g = 1.0;
        b = 1.0 - 4.0 * (value - 0.25);
    } 
    else if (value < 0.75) 
    {
        r = 4.0 * (value - 0.5);
        g = 1.0;
        b = 0.0;
    } 
    else 
    {
        r = 1.0;
        g = 1.0 - 4.0 * (value - 0.75);
        b = 0.0;
    }

    return std::make_tuple(static_cast<int>(r * 255), static_cast<int>(g * 255), static_cast<int>(b * 255));
}

// is_point_in_triangle(surface_coordinate0, surface_coordinate1, surface_coordinate2, surface_coordinate)
bool is_point_in_triangle(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c, const Eigen::Vector2d& p)
{
    // Compute vectors        
    Eigen::Vector2d v0 = c - a;
    Eigen::Vector2d v1 = b - a;
    Eigen::Vector2d v2 = p - a;

    // Compute dot products
    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    // Compute barycentric coordinates
    double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u > 0) && (v > 0) && (u + v < 1);
}

Eigen::Matrix3d generate_orthogonal_basis(const Eigen::Vector3d& unit_vector)
{
    // Ensure the input vector is normalized
    Eigen::Vector3d z = unit_vector.normalized();

    // Create an arbitrary vector that is not parallel to z
    Eigen::Vector3d arbitrary_vector;
    if (std::fabs(z.x()) > std::fabs(z.y()))
        arbitrary_vector = Eigen::Vector3d(-z.z(), 0, z.x());
    else
        arbitrary_vector = Eigen::Vector3d(0, z.z(), -z.y());

    // Compute the first basis vector (x) by normalizing the cross product of z and the arbitrary vector
    Eigen::Vector3d x = z.cross(arbitrary_vector).normalized();

    // Compute the second basis vector (y) as the cross product of z and x
    Eigen::Vector3d y = z.cross(x);

    // Construct the orthogonal basis matrix
    Eigen::Matrix3d basis;
    basis.col(0) = x;
    basis.col(1) = y;
    basis.col(2) = z;

    return basis;
}

// generate unit vector covariance
Eigen::Matrix3d generate_unit_vector_covariance(const Eigen::Vector3d& unit_vector, double angular_uncertainty_in_radian, double epsilon)
{
    // unrotated covariance matrix
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    cov(0, 0) = std::pow(angular_uncertainty_in_radian, 2);
    cov(1, 1) = std::pow(angular_uncertainty_in_radian, 2);
    cov(2, 2) = std::pow(epsilon, 2);

    // generate orthogonal basis
    Eigen::Matrix3d basis = generate_orthogonal_basis(unit_vector);

    // rotated covariance matrix
    return basis * cov * basis.transpose();
}