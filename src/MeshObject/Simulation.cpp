#include "MeshObject/Simulation.hpp"
#include <iostream> 
#include <random>

Simulation::Simulation()
{
}

void Simulation::set_object(int id)
{
    id_ = id;
}

void Simulation::set_noise(double noise_std)
{
    noise_std_ = noise_std;
}

void Simulation::get_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position)
{
    if (id_ == 0)
    {
        get_cube_data_pair(origin, position);
    }
    if (id_ == 1)
    {
        get_plane_data_pair(origin, position);
    }
    if (id_ == 2)
    {
        get_gap_data_pair(origin, position);
    }
    else
    {
        std::cout << "Error: invalid object id" << std::endl;
    }

    // add noise

    // gaussian noise with mean and variance
    double mean = 0.0; // mean of the Gaussian distribution
    double std = noise_std_; // variance of the Gaussian distribution

    // generate random noise value from Gaussian distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(mean, std);
    double noise_value = distribution(gen);

    // convert to noise vector and add to position
    Eigen::Vector3d noise_vector = (position - origin).normalized() * noise_value;
    position += noise_vector;
}

void Simulation::get_cube_data_pair(Eigen::Vector3d& out_origin, Eigen::Vector3d& out_position)
{
    // settings on the cube
    double size = 1.0;
    Eigen::Vector3d origin = Eigen::Vector3d(size, size, size);

    // compute position
    // xy, yz, xz
    
    // random double between 0 and size
    double u = (double)rand() / RAND_MAX * size;
    double v = (double)rand() / RAND_MAX * size;

    // set position
    // random generate an integer between 0 and 2
    int axis = rand() % 3;
    Eigen::Vector3d position = Eigen::Vector3d(0, 0, 0);
    if (axis == 0)
    {
        position = Eigen::Vector3d(u, v, 0);
    }
    else if (axis == 1)
    {
        position = Eigen::Vector3d(0, u, v);
    }
    else if (axis == 2)
    {
        position = Eigen::Vector3d(u, 0, v);
    }

    // return
    out_origin = origin;
    out_position = position;

    // log
    std::cout << "origin: " << origin.transpose() << ", position: " << position.transpose() << std::endl;
}

void Simulation::get_plane_data_pair(Eigen::Vector3d& out_origin, Eigen::Vector3d& out_position)
{
    // settings on the cube
    double size = 2.0;
    Eigen::Vector3d origin = Eigen::Vector3d(1, 1, 1);

    // compute position
    // xy, yz, xz
    
    // random double between 0 and size
    double u = (double)rand() / RAND_MAX * size;
    double v = (double)rand() / RAND_MAX * size;

    // set position
    // random generate an integer between 0 and 2
    Eigen::Vector3d position = Eigen::Vector3d(u, v, 0);

    // return
    out_origin = origin;
    out_position = position;

    // log
    std::cout << "origin: " << origin.transpose() << ", position: " << position.transpose() << std::endl;
}

void Simulation::get_gap_data_pair(Eigen::Vector3d& out_origin, Eigen::Vector3d& out_position)
{
    // settings on the cube
    double size = 2.0;
    Eigen::Vector3d origin = Eigen::Vector3d(1, 1, 1);

    // compute position
    // xy, yz, xz
    
    // random double between 0 and size
    double u = (double)rand() / RAND_MAX * size;
    double v = (double)rand() / RAND_MAX * size;

    double z;
    if (u < 0.6 || u > 1.4)
    {
        z = 0.0;
    }
    else
    {
        z = -0.1;
    }

    // set position
    Eigen::Vector3d position = Eigen::Vector3d(u, v, z);

    // return
    out_origin = origin;
    out_position = position;

    // log
    std::cout << "origin: " << origin.transpose() << ", position: " << position.transpose() << std::endl;
}