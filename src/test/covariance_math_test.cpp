#include "utilities/covariance_math.hpp"
#include <iostream>

void test1()
{
    std::cout << "=============================== start test 1 ===============================" << std::endl;

    // generate random points
    std::vector<Eigen::Vector3d> points_set1;
    for (std::size_t i = 0; i < 10; i++)
    {
        Eigen::Vector3d point;
        point << rand() % 100, rand() % 100, rand() % 100;
        points_set1.push_back(point);
    }

    // compute mean and covariance
    Eigen::Vector3d mean1(0, 0, 0);
    Eigen::Matrix3d covariance1 = Eigen::Matrix3d::Zero();
    for (std::size_t i = 0; i < points_set1.size(); i++)
    {
        mean1 += points_set1[i];
    }
    mean1 /= points_set1.size();
    for (std::size_t i = 0; i < points_set1.size(); i++)
    {
        Eigen::Vector3d diff = points_set1[i] - mean1;
        covariance1 += diff * diff.transpose();
    }
    covariance1 /= points_set1.size();

    // generate random points
    std::vector<Eigen::Vector3d> points_set2;
    for (std::size_t i = 0; i < 10; i++)
    {
        Eigen::Vector3d point;
        point << rand() % 100, rand() % 100, rand() % 100;
        points_set2.push_back(point);
    }

    // compute mean and covariance
    Eigen::Vector3d mean2(0, 0, 0);
    Eigen::Matrix3d covariance2 = Eigen::Matrix3d::Zero();
    for (std::size_t i = 0; i < points_set2.size(); i++)
    {
        mean2 += points_set2[i];
    }
    mean2 /= points_set2.size();
    for (std::size_t i = 0; i < points_set2.size(); i++)
    {
        Eigen::Vector3d diff = points_set2[i] - mean2;
        covariance2 += diff * diff.transpose();
    }
    covariance2 /= points_set2.size();


    // compute merged mean and cov
    // get combined set
    std::vector<Eigen::Vector3d> combined_points;
    for (std::size_t i = 0; i < points_set1.size(); i++)
    {
        combined_points.push_back(points_set1[i]);
    }
    for (std::size_t i = 0; i < points_set2.size(); i++)
    {
        combined_points.push_back(points_set2[i]);
    }

    // compute mean and covariance
    Eigen::Vector3d mean_combined(0, 0, 0);
    Eigen::Matrix3d covariance_combined = Eigen::Matrix3d::Zero();
    for (std::size_t i = 0; i < combined_points.size(); i++)
    {
        mean_combined += combined_points[i];
    }
    mean_combined /= combined_points.size();
    for (std::size_t i = 0; i < combined_points.size(); i++)
    {
        Eigen::Vector3d diff = combined_points[i] - mean_combined;
        covariance_combined += diff * diff.transpose();
    }
    covariance_combined /= combined_points.size();

    // merge mean and covariance
    Eigen::Vector3d merged_mean = merge_mean(mean1, mean2, points_set1.size(), points_set2.size());
    Eigen::Matrix3d merged_cov = merge_covariance(covariance1, covariance2, mean1, mean2, points_set1.size(), points_set2.size());
    int merged_size = points_set1.size() + points_set2.size();

    // print
    std::cout << "merge_means: " << std::endl;
    std::cout << merged_mean.transpose() << std::endl;
    std::cout << "combined_mean: " << std::endl;
    std::cout << mean_combined.transpose() << std::endl;

    std::cout << "merge_covariances: " << std::endl;
    std::cout << merged_cov << std::endl;
    std::cout << "combined_covariance: " << std::endl;
    std::cout << covariance_combined << std::endl;


    // remove_mean
    Eigen::Vector3d removed_mean = remove_mean(merged_mean, mean2, merged_size, points_set2.size());
    Eigen::Matrix3d removed_cov = remove_covariance(merged_cov, covariance2, merged_mean, mean2, merged_size, points_set2.size());

    // print
    std::cout << "remove_mean: " << std::endl;
    std::cout << removed_mean.transpose() << std::endl;
    std::cout << "mean1: " << std::endl;
    std::cout << mean1.transpose() << std::endl;
    std::cout << "remove_covariance: " << std::endl;
    std::cout << removed_cov << std::endl;
    std::cout << "covariance1: " << std::endl;
    std::cout << covariance1 << std::endl;
}

void test2()
{
    std::cout << "=============================== start test 2 ===============================" << std::endl;

    // surface
    Eigen::Vector3d mean1(3, 3, 3);
    Eigen::Matrix3d cov1;
    cov1 << 1, 2, 3,
            4, 5, 6,
            7, 8, 9;
    int size1 = 10;

    // point
    Eigen::Vector3d mean2(2, 2, 2);
    Eigen::Matrix3d cov2;
    cov2 << 0, 0, 0,
            0, 0, 0,
            0, 0, 0;
    int size2 = 1;

    Eigen::Vector3d merged_mean = merge_mean(mean1, mean2, size1, size2);
    Eigen::Matrix3d merged_cov = merge_covariance(cov1, cov2, mean1, mean2, size1, size2);
    int merged_size = size1 + size2;

    // print
    std::cout << "merge_means: " << std::endl;
    std::cout << merged_mean.transpose() << std::endl;
    std::cout << "merged_covariance: " << std::endl;
    std::cout << merged_cov << std::endl;
    std::cout << "merged_size: " << std::endl;
    std::cout << merged_size << std::endl;


    // recover mean1, cov1 and size1
    Eigen::Vector3d removed_mean = remove_mean(merged_mean, mean2, merged_size, size2);
    Eigen::Matrix3d removed_cov = remove_covariance(merged_cov, cov2, merged_mean, mean2, merged_size, size2);
    int removed_size = merged_size - size2;

    // print
    std::cout << "remove_mean: " << std::endl;
    std::cout << removed_mean.transpose() << std::endl;
    std::cout << "removed_covariance: " << std::endl;
    std::cout << removed_cov << std::endl;
    std::cout << "removed_size: " << std::endl;
    std::cout << removed_size << std::endl;
}

int main()
{
    test1();
    test2();
    return 0;
}