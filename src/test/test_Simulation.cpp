#include "MeshObject/Simulation.hpp"

int main()
{
    Simulation sim;
    sim.set_object(0);
    Eigen::Vector3d origin, position;
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    sim.get_data_pair(origin, position);
    return 0;
}