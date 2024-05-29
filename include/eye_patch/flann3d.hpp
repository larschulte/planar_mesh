#ifndef FLANN3D_H
#define FLANN3D_H

#include <set>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <flann/flann.hpp>
#include <flann/util/matrix.h>

class flann3d
{
public:
    flann3d();
    
    void set_input(std::set<int> point_set, std::map<int, Eigen::Vector3d> point_to_vector3d_map);
    void addPoint(Eigen::Vector3d new_point, int point_id);
    std::set<int> radiusSearch(Eigen::Vector3d searchPoint, double radius);

    int flann_last_id;

private:
    std::vector<double> flann_data_storage;
    flann::Index<flann::L2_Simple<double>> flann_tree;
    std::vector<int> index_to_pointID; // index to point id correspondence
};

#endif // FLANN3D_H
