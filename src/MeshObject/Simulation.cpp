#include "MeshObject/Simulation.hpp"
#include <iostream> 

Simulation::Simulation()
{
}

void Simulation::set_object(int id)
{
    id_ = id;
}

void Simulation::get_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position)
{
    if (id_ == 0)
    {
        return get_cube_data_pair(origin, position);
    }
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