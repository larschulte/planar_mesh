#include <Eigen/Dense>

Eigen::Vector3d merge_mean(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2);
Eigen::Matrix3d merge_covariance(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                int size1, int size2);

Eigen::Vector3d remove_mean(const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, int combined_size, int size2);
Eigen::Matrix3d remove_covariance(const Eigen::Matrix3d& combined_covariance, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, 
                                int combined_size, int size2);

Eigen::Vector3d weighted_merge_mean(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, double weight1, double weight2);
Eigen::Matrix3d weighted_merge_covariance(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                double weight1, double weight2);

Eigen::Vector3d weighted_remove_mean(const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, double combined_weight, double weight2);
Eigen::Matrix3d weighted_remove_covariance(const Eigen::Matrix3d& combined_covariance, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& combined_mean, const Eigen::Vector3d& mean2, 
                                double combined_weight, double weight2);

double merge_information_weighted_mean(const double& mean1, const double& mean2, const double& information1, const double& information2);
double merge_information(const double& information1, const double& information2);

double remove_information_weighted_mean(const double& combined_mean, const double& mean2, const double& combined_information, const double& information2);
double remove_information(const double& combined_information, const double& information2);

double compute_mean(const std::vector<double>& data);
double compute_std(const std::vector<double>& data);