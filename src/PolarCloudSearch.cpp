#include "eye_patch/PolarCloudSearch.hpp"

PolarCloudSearch::PolarCloudSearch()
{
}

void PolarCloudSearch::setCloudPointer(pcl::PointCloud<PointT>::Ptr old_pointcloud)
{
    old_pointcloud_ = old_pointcloud;
}

void PolarCloudSearch::setResolution(float res_azimuth, float res_altitude)
{
    // input check
    if (res_azimuth <= 0 || res_altitude <= 0) {
        throw std::invalid_argument("azimuth and altitude resolutions must be > 0");
    }

    // compute
    res_azimuth_ = res_azimuth;
    res_altitude_ = res_altitude;
    number_of_azimuth_ = std::ceil(360.f/res_azimuth_);
    number_of_altitude_ = std::ceil(180.f/res_altitude_);

    // resize
    grid_.resize(number_of_azimuth_, std::vector<pcl::Indices>(number_of_altitude_)); 
}

void PolarCloudSearch::prepareSearch(const Eigen::Isometry3f &T_world2lidar)
{
    clearGrid();
    
    // Note: This can't be parallelized since the push_back() command means
    // the loop iterations are order-dependent.
    // #pragma omp parallel for
    for (int i = 0; i < old_pointcloud_->size(); i++)
    {
        GridIndexBounds bounds = calculateGridIndexBounds(old_pointcloud_->points[i], T_world2lidar);
        addIndexToBounds(bounds, i);
    }
}

void PolarCloudSearch::clearGrid()
{
    for (auto &row : grid_)
    {
        for (auto &column : row)
        {
            column.clear();
        }
    }
}

void PolarCloudSearch::addPoint(const PointT& lidar_point, const Eigen::Isometry3f &T_world2lidar)
{
    int i = old_pointcloud_->size();
    old_pointcloud_->push_back(lidar_point);

    GridIndexBounds bounds = calculateGridIndexBounds(lidar_point, T_world2lidar);
    addIndexToBounds(bounds, i);
}

pcl::Indices PolarCloudSearch::extractNearbyIndices(const PointT& lidar_point,
                                                    const Eigen::Isometry3f &T_world2lidar)
{
    GridIndexBounds bounds = calculateGridIndexBounds(lidar_point, T_world2lidar);
    return extractIndicesFromBounds(bounds);
}

GridIndexBounds PolarCloudSearch::calculateGridIndexBounds(const PointT& surfel, const Eigen::Isometry3f& T_world2lidar) 
{
    Eigen::Vector3f normal = surfel.getNormalVector3fMap();
    Eigen::Vector3f position = surfel.getArray3fMap();
    float radius = surfel.radius;

    Eigen::Vector3f viewpoint = T_world2lidar.translation();

    Eigen::Vector3f view_to_surfel = position - viewpoint;

    Eigen::Vector3f major_axis = normal.cross(view_to_surfel).normalized();
    Eigen::Vector3f minor_axis = normal.cross(major_axis).normalized();

    Eigen::Vector3f point1 = position + radius * (  major_axis + minor_axis);
    Eigen::Vector3f point2 = position + radius * (  major_axis - minor_axis);
    Eigen::Vector3f point3 = position + radius * (- major_axis + minor_axis);
    Eigen::Vector3f point4 = position + radius * (- major_axis - minor_axis);

    Eigen::Vector3f point1_L = T_world2lidar * point1;
    Eigen::Vector3f point2_L = T_world2lidar * point2; 
    Eigen::Vector3f point3_L = T_world2lidar * point3;
    Eigen::Vector3f point4_L = T_world2lidar * point4;

    Eigen::Vector3f spherical_1 = computeSphericalCoordinate(point1_L);
    Eigen::Vector3f spherical_2 = computeSphericalCoordinate(point2_L);
    Eigen::Vector3f spherical_3 = computeSphericalCoordinate(point3_L);
    Eigen::Vector3f spherical_4 = computeSphericalCoordinate(point4_L);
    
    GridIndexBounds bounds;
    bounds.start_azimuth  = index_of_azimuth(std::min({spherical_1[1], spherical_2[1], spherical_3[1], spherical_4[1]}));
    bounds.end_azimuth    = index_of_azimuth(std::max({spherical_1[1], spherical_2[1], spherical_3[1], spherical_4[1]}));
    bounds.start_altitude = index_of_altitude(std::min({spherical_1[2], spherical_2[2], spherical_3[2], spherical_4[2]}));
    bounds.end_altitude   = index_of_altitude(std::max({spherical_1[2], spherical_2[2], spherical_3[2], spherical_4[2]}));

    return bounds;
}

void PolarCloudSearch::addIndexToBounds(const GridIndexBounds& bounds, int point_index) 
{
    // this doesn't consider the case when azimuth wraps around!!!
    for (int i_azimuth = bounds.start_azimuth; i_azimuth <= bounds.end_azimuth; i_azimuth ++) 
    {
        for (int i_altitude = bounds.start_altitude; i_altitude <= bounds.end_altitude; i_altitude ++) 
        {
            grid_[i_azimuth][i_altitude].push_back(point_index);
        }
    }
}

pcl::Indices PolarCloudSearch::extractIndicesFromBounds(const GridIndexBounds& bounds) 
{
    // this doesn't consider the case when azimuth wraps around!!!
    pcl::Indices indices_within_bounds;
    for (int i_azimuth = bounds.start_azimuth; i_azimuth <= bounds.end_azimuth; i_azimuth ++) 
    {
        for (int i_altitude = bounds.start_altitude; i_altitude <= bounds.end_altitude; i_altitude ++)
        {
            indices_within_bounds.insert(indices_within_bounds.end(), grid_[i_azimuth][i_altitude].begin(), grid_[i_azimuth][i_altitude].end());
        }
    }
    return indices_within_bounds;
}

Eigen::Vector3f PolarCloudSearch::computeSphericalCoordinate(const Eigen::Vector3f &XL)
{
    const float x = XL[0];
    const float y = XL[1];
    const float z = XL[2];

    // spherical
    static constexpr float rad2deg = 180 / M_PI;
    const float range = XL.norm();
    const float azimuth = std::atan2(y, x) * rad2deg;
    const float altitude = std::asin(z / range) * rad2deg;
    
    Eigen::Vector3f spherical(range, azimuth, altitude);

    return spherical;
}

int PolarCloudSearch::index_of_azimuth(float azimuth_angle)
{
    // compute shperical returns -180 and 180, adding 180 makes all positive. 
    int index = (azimuth_angle + 180)/res_azimuth_; 
    // this is to wrap around azimuth when 180 +- 0.1 degree is passed in
    if (index < 0) index += number_of_azimuth_; 
    if (index >= number_of_azimuth_) index -= number_of_azimuth_;

    return index;
}

int PolarCloudSearch::index_of_altitude(float altitude_angle)
{
    // compute shperical returns -90 and 90, adding 90 makes all positive. 
    int index = (altitude_angle + 90)/res_altitude_; 
    // no wrap around for altitude, but instead caps
    if (index < 0) index = 0; 
    if (index >= number_of_altitude_) index = number_of_altitude_;

    return index;
}
