#include "utilities/covariance_math.hpp"
#include <numeric>

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


// weighted add mean2 to mean1
Eigen::Vector3d weighted_merge_mean(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, double weight1, double weight2) 
{
    return (weight1 * mean1 + weight2 * mean2) / (weight1 + weight2);
}

// weighted add cov2 to cov1
Eigen::Matrix3d weighted_merge_covariance(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                double weight1, double weight2) 
{
    // Handle the edge case where one of the sizes is zero
    if (weight1 == 0 && weight2 == 0) throw std::invalid_argument("Both sizes are zero");
    if (weight1 == 0) return cov2;
    if (weight2 == 0) return cov1;

    Eigen::Vector3d combined_mean = merge_mean(mean1, mean2, weight1, weight2);
    Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
    Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
    Eigen::Matrix3d combined_covariance = (weight1 * cov1 + weight2 * cov2 + weight1 * mean_diff1 + weight2 * mean_diff2) / (weight1 + weight2);

    return combined_covariance;
}

// weighted remove mean2 from combined_mean
Eigen::Vector3d weighted_remove_mean(const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, double combined_weight, double weight2) 
{
    double epsilon = 1e-6;

    if (combined_weight - weight2 < -epsilon) throw std::invalid_argument("Size of the combined_weight " + std::to_string(combined_weight) + " must be greater than size of weight2 " + std::to_string(weight2));
    double weight1 = combined_weight - weight2;
    if (std::abs(weight1) < epsilon) return Eigen::Vector3d::Zero();
    
    return (combined_weight * combined_mean - weight2 * mean2) / weight1;
}

// weighted remove cov2 from combined_covariance
Eigen::Matrix3d weighted_remove_covariance(const Eigen::Matrix3d& combined_covariance, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, 
                                double combined_weight, double weight2) 
{
    double epsilon = 1e-6;

    if (combined_weight - weight2 < -epsilon) throw std::invalid_argument("Size of the combined_weight " + std::to_string(combined_weight) + " must be greater than size of weight2 " + std::to_string(weight2));
    double weight1 = combined_weight - weight2;
    if (std::abs(weight1) < epsilon) return Eigen::Matrix3d::Zero();

    Eigen::Vector3d mean1 = remove_mean(combined_mean, mean2, combined_weight, weight2);

    Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
    Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
    Eigen::Matrix3d cov1 = (combined_covariance * combined_weight - weight2 * cov2 - weight1 * mean_diff1 - weight2 * mean_diff2) / weight1;

    return cov1;
}

double merge_information_weighted_mean(const double& mean1, const double& mean2, const double& information1, const double& information2) 
{
    return (information1 * mean1 + information2 * mean2) / (information1 + information2);
}

double merge_information(const double& information1, const double& information2) 
{
    return information1 + information2;
}

double remove_information_weighted_mean(const double& combined_mean, const double& mean2, const double& combined_information, const double& information2) 
{
    double information1 = combined_information - information2;
    return (combined_mean * combined_information - information2 * mean2) / information1;
}

double remove_information(const double& combined_information, const double& information2) 
{
    return combined_information - information2;
}

// Function to compute the mean of a vector
double compute_mean(const std::vector<double>& data) 
{
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Function to compute the standard deviation of a vector
double compute_std(const std::vector<double>& data)
{
    double mean = compute_mean(data);
    double sq_sum = std::accumulate(data.begin(), data.end(), 0.0, 
        [mean](double acc, double val) 
        {
            return acc + (val - mean) * (val - mean);
        });
    return std::sqrt(sq_sum / data.size());
}