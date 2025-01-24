#include "utilities/DataLoader.hpp"
#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"
#include <pcl/filters/passthrough.h>

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

using Pose = std::tuple<double, double, double, double, double, double, double>;

Pose parse_g2o_line(const std::string& line, const std::string& sec_str, const std::string& nsec_str, bool& pose_found) 
{
    std::istringstream iss(line);
    std::string type;
    iss >> type;

    if (type == "VERTEX_SE3:QUAT_TIME") {
        int vertex_id;
        double x, y, z, qx, qy, qz, qw;
        int timestamp_sec, timestamp_nsec;
        iss >> vertex_id >> x >> y >> z >> qx >> qy >> qz >> qw >> timestamp_sec >> timestamp_nsec;

        if (timestamp_sec == std::stoi(sec_str) && timestamp_nsec == std::stoi(nsec_str)) {
            pose_found = true;
            return {x, y, z, qx, qy, qz, qw};
        }
    }

    return {0, 0, 0, 0, 0, 0, 0}; // Default pose
}

Pose parse_csv_line(const std::string& line, const std::string& sec_str, const std::string& nsec_str, bool& pose_found) 
{
    if (line[0] == '#') 
    {
        return {0, 0, 0, 0, 0, 0, 0}; // Default pose
    }

    std::istringstream iss(line);
    std::string token;
    std::vector<std::string> tokens;

    while (std::getline(iss, token, ',')) 
    {
        tokens.push_back(token);
    }

    if (tokens.size() >= 10) 
    {
        int timestamp_sec = std::stoi(tokens[1]);
        int timestamp_nsec = std::stoi(tokens[2]);

        if (timestamp_sec == std::stoi(sec_str) && timestamp_nsec == std::stoi(nsec_str)) 
        {
            double x = std::stod(tokens[3]);
            double y = std::stod(tokens[4]);
            double z = std::stod(tokens[5]);
            double qx = std::stod(tokens[6]);
            double qy = std::stod(tokens[7]);
            double qz = std::stod(tokens[8]);
            double qw = std::stod(tokens[9]);
            pose_found = true;
            return {x, y, z, qx, qy, qz, qw};
        }
    }

    return {0, 0, 0, 0, 0, 0, 0}; // Default pose
}

Pose parse_tum_line(const std::string& line, const std::string& sec_str, const std::string& nsec_str, bool& pose_found) 
{
    // get line
    std::istringstream iss(line);

    // get individual values
    std::string timestamp, sec, nsec;
    double x, y, z, qx, qy, qz, qw;
    iss >> timestamp >> x >> y >> z >> qx >> qy >> qz >> qw;
    sec = timestamp.substr(0, timestamp.find('.'));
    nsec = timestamp.substr(timestamp.find('.') + 1);

    // return
    if (sec == sec_str && nsec == nsec_str) 
    {
        pose_found = true;
        return {x, y, z, qx, qy, qz, qw};
    }
    else
    {
        return {0, 0, 0, 0, 0, 0, 0};
    }
}

Eigen::Affine3d find_pose(const std::string& pcd_file, const std::string& pose_file) 
{
    // Extract sec and nsec from the pcd_file name
    std::string pcd_file_name = pcd_file.substr(pcd_file.find_last_of("/\\") + 1);
    std::string sec_str = pcd_file_name.substr(6, 10);
    std::string nsec_str = pcd_file_name.substr(17, 9);

    // Determine file type and set parser function
    std::ifstream pose_stream(pose_file);
    std::string extension = pose_file.substr(pose_file.find_last_of(".") + 1);

    std::function<Pose(const std::string&, const std::string&, const std::string&, bool&)> parser;
    if (extension == "g2o" || extension == "slam") 
    {
        parser = parse_g2o_line;
    } 
    else if (extension == "csv") 
    {
        parser = parse_csv_line;
    } 
    else if (extension == "txt")
    {
        parser = parse_tum_line;
    }
    else 
    {
        std::cerr << "Unsupported file format: " << extension << std::endl;
        return Eigen::Isometry3d::Identity();
    }

    // Parse the file and find the pose
    std::string line;
    bool pose_found = false;
    Pose pose;

    while (std::getline(pose_stream, line)) 
    {
        pose = parser(line, sec_str, nsec_str, pose_found);
        if (pose_found) 
        {
            break;
        }
    }

    if (!pose_found) 
    {
        std::cout << "Pose not found for time " << sec_str << " " << nsec_str << std::endl;
        return Eigen::Isometry3d::Identity();
    }

    // Convert the pose to Eigen::Affine3d
    Eigen::Affine3d pose_eigen = Eigen::Isometry3d::Identity();
    pose_eigen.translation() << std::get<0>(pose), std::get<1>(pose), std::get<2>(pose);
    Eigen::Quaterniond q(std::get<6>(pose), std::get<3>(pose), std::get<4>(pose), std::get<5>(pose));
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
void DataLoader<PointT>::load_dataset(DataLoader_Settings settings)
{
    // set settings
    azimuth_resolution_ = settings.azimuth_resolution;
    altitude_resolution_ = settings.altitude_resolution;
    remove_double_return_flag_ = settings.remove_double_return_flag;
    filter_low_intensity_flag_ = settings.filter_low_intensity_flag;

    // load dataset
    pcd_file_list_ = read_under_folder(settings.pcd_file_folder);
    file_to_pose_map_ = create_file_to_pose_map(pcd_file_list_, settings.pose_file_path);
}

template <typename PointT>
DataLoader<PointT>::DataLoader(std::string pcd_file_folder, std::string pose_file_path)
{
    pcd_file_list_ = read_under_folder(pcd_file_folder);
    file_to_pose_map_ = create_file_to_pose_map(pcd_file_list_, pose_file_path);
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr DataLoader<PointT>::remove_double_return(typename pcl::PointCloud<PointT>::Ptr input_pointcloud)
{
    // initialize mapping grid
    std::map<std::pair<double, double>, std::pair<double, int>> point_map; // key: azimuth and altitude, value: range and index

    // for each point, compute azimuth and altitude and store in the map if it is closer
    for (std::size_t i = 0; i < input_pointcloud->size(); i++)
    {
        // get point
        PointT point = input_pointcloud->points[i];

        // compute azimuth and altitude
        double range = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
        double azimuth = atan2(point.y, point.x) / M_PI * 180;
        double altitude = asin(point.z / range) / M_PI * 180;

        // round to the nearest resolution
        azimuth = round(azimuth / azimuth_resolution_) * azimuth_resolution_;
        altitude = round(altitude / altitude_resolution_) * altitude_resolution_;

        // if same azimuth and altitude not already in the map
        if (point_map.find(std::make_pair(azimuth, altitude)) == point_map.end())
        {
            // store azimuth and altitude as key, range and the index of the point as value
            point_map[std::make_pair(azimuth, altitude)] = std::make_pair(range, i);
        }
        else
        {
            // compare range, keep the closer one
            if (range < point_map[std::make_pair(azimuth, altitude)].first)
            {
                point_map[std::make_pair(azimuth, altitude)] = std::make_pair(range, i);
            }
        }        
    }

    // from the map, create a new pointcloud
    typename pcl::PointCloud<PointT>::Ptr output_pointcloud(new pcl::PointCloud<PointT>);
    for (const auto& point : point_map)
    {
        output_pointcloud->push_back(input_pointcloud->points[point.second.second]);
    }

    // return 
    return output_pointcloud;
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr DataLoader<PointT>::filter_low_intensity(typename pcl::PointCloud<PointT>::Ptr input_pointcloud)
{
    // get intensity field name
    std::string intensity_field_name;
    if (std::is_same<PointT, VilensPointT>::value)
    {
        intensity_field_name = "curvature";
    }
    else if (std::is_same<PointT, BagPointT>::value)
    {
        intensity_field_name = "intensity";
    }

    // filter point cloud by intensity
    pcl::PassThrough<PointT> filter_low_intensity_points;
    filter_low_intensity_points.setInputCloud(input_pointcloud);
    filter_low_intensity_points.setFilterFieldName (intensity_field_name);
    filter_low_intensity_points.setFilterLimits(150, 255);
    filter_low_intensity_points.filter(*input_pointcloud);

    return input_pointcloud;
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr DataLoader<PointT>::get_cloud(int i)
{
    typename pcl::PointCloud<PointT>::Ptr loaded_pointcloud = load_pointcloud<PointT>(pcd_file_list_[i]);

    if (remove_double_return_flag_) loaded_pointcloud = remove_double_return(loaded_pointcloud);
    if (filter_low_intensity_flag_) loaded_pointcloud = filter_low_intensity(loaded_pointcloud);

    return loaded_pointcloud;
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