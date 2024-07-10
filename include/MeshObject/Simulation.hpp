#include <Eigen/Dense>

class Simulation
{
public:
    Simulation();

    void set_object(int id);
    void set_noise(double noise_std);
    void get_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_cube_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_plane_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);
    void get_gap_data_pair(Eigen::Vector3d& origin, Eigen::Vector3d& position);

private:
    int id_;
    double noise_std_;
};
