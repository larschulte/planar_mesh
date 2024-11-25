#include "utilities/DataLoader.hpp"
#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"

std::vector<std::string> read_under_folder(std::string pcd_file_folder)
{
    // initialize
    std::vector<std::string> pcd_file_list;

    // read
    DIR *dirstream = opendir(pcd_file_folder.c_str());
    struct dirent *entry;
    while ((entry = readdir(dirstream)) != NULL) 
    {
        // obtain file name
        std::string file_name = entry->d_name;
        if (file_name.find(".pcd") == std::string::npos) continue; // skip if not pcd file

        // pushback
        pcd_file_list.push_back(pcd_file_folder + file_name);
    }
    closedir(dirstream);
    
    // sort
    std::sort(pcd_file_list.begin(), pcd_file_list.end());

    // return 
    return pcd_file_list;
}

Eigen::Affine3d find_pose(std::string pcd_file, std::string pose_file)
{
    // get sec and nsec from pcd_file, using split
    std::string pcd_file_name = pcd_file.substr(pcd_file.find_last_of("/\\") + 1);
    std::string sec_str = pcd_file_name.substr(6, 10);
    std::string nsec_str = pcd_file_name.substr(17, 9);

    // find pose 
    std::ifstream pose_stream(pose_file);
    std::string line;
    double pose_x, pose_y, pose_z, pose_qx, pose_qy, pose_qz, pose_qw;
    bool pose_found = false;
    while (std::getline(pose_stream, line))
    {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "VERTEX_SE3:QUAT_TIME")
        {
            int vertex_id;
            double x, y, z, qx, qy, qz, qw;
            int timestamp_sec, timestamp_nsec;
            iss >> vertex_id >> x >> y >> z >> qx >> qy >> qz >> qw >> timestamp_sec >> timestamp_nsec;

            if (timestamp_sec == std::stoi(sec_str) && timestamp_nsec == std::stoi(nsec_str))
            {
                pose_x = x;
                pose_y = y;
                pose_z = z;
                pose_qx = qx;
                pose_qy = qy;
                pose_qz = qz;
                pose_qw = qw;
                pose_found = true;
                break;
            }
        }
    }
    if (!pose_found)
    {
        std::cout << "pose not found for time " << sec_str << " " << nsec_str << std::endl;
        return Eigen::Isometry3d::Identity();
    }

    // convert pose to eigen::isometry3d
    Eigen::Affine3d pose_eigen = Eigen::Isometry3d::Identity();
    pose_eigen.translation() << pose_x, pose_y, pose_z;
    Eigen::Quaterniond q(pose_qw, pose_qx, pose_qy, pose_qz);
    pose_eigen.rotate(q);

    return pose_eigen;
}

std::map<std::string, Eigen::Affine3d> create_file_to_pose_map(std::vector<std::string> pcd_file_list, std::string pose_file_path)
{
    // initialize
    std::map<std::string, Eigen::Affine3d> file_to_pose_map;

    // create
    for (std::string pcd_file : pcd_file_list)
    {
        Eigen::Affine3d pose = find_pose(pcd_file, pose_file_path);
        file_to_pose_map[pcd_file] = pose;
    }

    // return 
    return file_to_pose_map;
}

template<typename PointT> 
typename pcl::PointCloud<PointT>::Ptr load_pointcloud(std::string pcd_file)
{
    // load the pcd file
    typename pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);
    if (pcl::io::loadPCDFile<PointT> (pcd_file, *cloud) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
    }
    return cloud;
}

template <typename PointT>
DataLoader<PointT>::DataLoader(){}

template <typename PointT>
void DataLoader<PointT>::load_dataset(std::string pcd_file_folder, std::string pose_file_path)
{
    pcd_file_list_ = read_under_folder(pcd_file_folder);
    file_to_pose_map_ = create_file_to_pose_map(pcd_file_list_, pose_file_path);
}

template <typename PointT>
DataLoader<PointT>::DataLoader(std::string pcd_file_folder, std::string pose_file_path)
{
    pcd_file_list_ = read_under_folder(pcd_file_folder);
    file_to_pose_map_ = create_file_to_pose_map(pcd_file_list_, pose_file_path);
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr DataLoader<PointT>::get_cloud(int i)
{
    return load_pointcloud<PointT>(pcd_file_list_[i]);
}

template <typename PointT>
Eigen::Affine3d DataLoader<PointT>::get_pose(int i)
{
    return file_to_pose_map_[pcd_file_list_[i]];
}

template <typename PointT>
int DataLoader<PointT>::size()
{
    return pcd_file_list_.size();
}

// Explicit template instantiation
template class DataLoader<VilensPointT>;
template class DataLoader<BagPointT>;