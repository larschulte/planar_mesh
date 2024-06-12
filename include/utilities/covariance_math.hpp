#include <Eigen/Dense>

Eigen::Vector3d merge_mean(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2);
Eigen::Matrix3d merge_covariance(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                int size1, int size2);

Eigen::Vector3d remove_mean(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2);
Eigen::Matrix3d remove_covariance(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                int size1, int size2);