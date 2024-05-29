#include "eye_patch/flann3d.hpp"

flann3d::flann3d() : flann_tree(flann::KDTreeIndexParams(1))
{
    // placeholder
    std::set<int> point_set = {-1};
    std::map<int, Eigen::Vector3d> point_to_vector3d_map = {{-1, Eigen::Vector3d(0, 0, 0)}};
    set_input(point_set, point_to_vector3d_map);
}

void flann3d::set_input(std::set<int> point_set, std::map<int, Eigen::Vector3d> point_to_vector3d_map)
{
    // reserve
    flann_data_storage.reserve(point_set.size() * 3);
    index_to_pointID.reserve(point_set.size());

    // process
    for (int point_id : point_set)
    {
        // add to data storage
        flann_data_storage.push_back(point_to_vector3d_map.at(point_id)[0]);
        flann_data_storage.push_back(point_to_vector3d_map.at(point_id)[1]);
        flann_data_storage.push_back(point_to_vector3d_map.at(point_id)[2]);

        // add to index
        index_to_pointID.push_back(point_id);
    }

    // add to flann
    flann_tree.buildIndex(flann::Matrix<double>(flann_data_storage.data(), point_set.size(), 3));

    // update id
    flann_last_id = point_set.size() - 1;
}

void flann3d::addPoint(Eigen::Vector3d new_point, int point_id)
{
    // add to data storage
    flann_data_storage.push_back(new_point[0]);
    flann_data_storage.push_back(new_point[1]);
    flann_data_storage.push_back(new_point[2]);

    // add to index
    index_to_pointID.push_back(point_id);

    // add to flann
    flann_tree.addPoints(flann::Matrix<double>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3));

    // update id
    flann_last_id++;
}

std::set<int> flann3d::radiusSearch(Eigen::Vector3d searchPoint, double radius)
{
    // initialize
    std::vector<int> searched_indices;
    std::vector<double> searched_dists;

    // convert to vector
    std::vector<double> query_point = {searchPoint[0], searchPoint[1], searchPoint[2]}; 

    // initialize
    std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>());
    std::vector<std::vector<double>> list_of_search_dists(1, std::vector<double>());

    // search
    flann_tree.radiusSearch(flann::Matrix<double>(query_point.data(), 1, 3), list_of_search_indices, list_of_search_dists, radius * radius, flann::SearchParams(-1, 0));

    // extract
    searched_indices = list_of_search_indices[0];
    searched_dists = list_of_search_dists[0];

    // convert to id
    std::set<int> searched_ids;
    for (const int& index : searched_indices)
    {
        searched_ids.insert(index_to_pointID[index]);
    }

    // remove placeholder
    searched_ids.erase(-1);
    
    // return
    return searched_ids;
}
