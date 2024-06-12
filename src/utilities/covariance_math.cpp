#include "utilities/covariance_math.hpp"

// add mean2 to mean1
Eigen::Vector3d merge_mean(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2) 
{
    return (size1 * mean1 + size2 * mean2) / (size1 + size2);
}

// add cov2 to cov1
Eigen::Matrix3d merge_covariance(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                int size1, int size2) 
{
    // Handle the edge case where one of the sizes is zero
    if (size1 == 0 && size2 == 0) throw std::invalid_argument("Both sizes are zero");
    if (size1 == 0) return cov2;
    if (size2 == 0) return cov1;

    Eigen::Vector3d combined_mean = merge_mean(mean1, mean2, size1, size2);
    Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
    Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
    Eigen::Matrix3d combined_covariance = (size1 * cov1 + size2 * cov2 + size1 * mean_diff1 + size2 * mean_diff2) / (size1 + size2);

    return combined_covariance;
}

// remove mean2 from combined_mean
Eigen::Vector3d remove_mean(const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, int combined_size, int size2) 
{
    if (combined_size < size2) throw std::invalid_argument("Size of the combined_mean must be greater than size of mean2");
    int size1 = combined_size - size2;
    if (size1 == 0) return Eigen::Vector3d::Zero();
    
    return (combined_size * combined_mean - size2 * mean2) / size1;
}

// remove cov2 from combined_covariance
Eigen::Matrix3d remove_covariance(const Eigen::Matrix3d& combined_covariance, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, 
                                int combined_size, int size2) 
{
    if (combined_size < size2) throw std::invalid_argument("Size of the combined_cov must be greater than size of the cov2");
    int size1 = combined_size - size2;
    if (size1 == 0) return Eigen::Matrix3d::Zero();

    Eigen::Vector3d mean1 = remove_mean(combined_mean, mean2, combined_size, size2);

    Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
    Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
    Eigen::Matrix3d cov1 = (combined_covariance * combined_size - size2 * cov2 - size1 * mean_diff1 - size2 * mean_diff2) / size1;

    return cov1;
}