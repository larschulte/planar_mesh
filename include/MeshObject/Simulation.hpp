#include <Eigen/Dense>

class Simulation
{
public:
    Simulation();

    void set_object(int id);
    void get_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_cube_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_plane_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);

private:
    int id_;
};
