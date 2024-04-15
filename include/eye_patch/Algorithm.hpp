#pragma once
#define PCL_NO_PRECOMPILE

#include <delaunator.hpp>
#include <pcl/common/transforms.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/extract_indices.h>


// point to triangle map for d2
std::map<int, std::vector<int>> 
point_to_triangle_map(delaunator::Delaunator d)
{
    std::map<int, std::vector<int>> pt_map;
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        // indices
        int i1 = d.triangles[i];
        int i2 = d.triangles[i + 1];
        int i3 = d.triangles[i + 2];

        // store
        pt_map[i1].push_back(i);
        pt_map[i2].push_back(i);
        pt_map[i3].push_back(i);
    }

    return pt_map;
}

// convert to 2d polar cloud
template<typename PointT>
pcl::PointCloud<pcl::PointXY>::Ptr 
obtain_2d_polar_cloud(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    pcl::PointCloud<pcl::PointXY>::Ptr cloud_polar (new pcl::PointCloud<pcl::PointXY>);
    cloud_polar->resize(cloud->size());
    for (std::size_t i = 0; i < cloud->size(); i+=1)
    {
        pcl::PointXY point;

        // compute azimuth and altitude
        float x = cloud->points[i].x;
        float y = cloud->points[i].y;
        float z = cloud->points[i].z;
        double r = sqrt(x * x + y * y + z * z);
        double azimuth = atan2(y, x) * 180 / M_PI;
        double altitude = asin(z / r) * 180 / M_PI;

        point.x = azimuth;
        point.y = altitude;
        cloud_polar->points[i] = point;
    }

    return cloud_polar;
}

// check if within triangle
bool 
is_inside_triangle(Eigen::Vector3f p0, Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f point)
{
    Eigen::Vector3f v0 = p1 - p0;
    Eigen::Vector3f v1 = p2 - p0;
    Eigen::Vector3f v2 = point - p0;
    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);
    float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return (u >= 0) && (v >= 0) && (u + v < 1);
}

// compute ray to triangle intersection
Eigen::Vector3f 
ray_triangle_intersection(Eigen::Vector3f ray_origin, Eigen::Vector3f ray_direction, Eigen::Vector3f p0, Eigen::Vector3f p1, Eigen::Vector3f p2)
{
    // compute triangle normal
    Eigen::Vector3f normal = (p1 - p0).cross(p2 - p0);
    
    // if parallel, return NaN
    if (normal.dot(ray_direction) == 0)
    {
        return Eigen::Vector3f(NAN, NAN, NAN);
    }

    // compute intersection
    double distance = (p0 - ray_origin).dot(normal) / ray_direction.dot(normal);
    Eigen::Vector3f intersection = ray_origin + ray_direction * distance;

    // return
    return intersection;
}

// compute direction vector
template<typename PointT>
std::vector<Eigen::Vector3f> 
compute_point_directions(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    // initialize
    std::vector<Eigen::Vector3f> directions;
    directions.resize(cloud->size());

    // compute
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        Eigen::Vector3f direction = cloud->points[i].getVector3fMap().normalized();
        directions[i] = direction;
    }

    // return
    return directions;
}


// get k nearest neighbor in tree
template<typename PointT>
pcl::Indices 
get_k_nearest_neighbor(pcl::KdTreeFLANN<PointT> kdtree, PointT point, int k)
{
    pcl::Indices nearest_neighbor(k);
    std::vector<float> nearest_neighbor_distance(k);
    kdtree.nearestKSearch(point, k, nearest_neighbor, nearest_neighbor_distance);
    return nearest_neighbor;
}

// compute triangle center
template<typename PointT>
typename pcl::PointCloud<pcl::PointXYZ>::Ptr 
computer_triangle_center(typename pcl::PointCloud<PointT>::Ptr vertex_cloud, delaunator::Delaunator d)
{
    // initialize
    typename pcl::PointCloud<pcl::PointXYZ>::Ptr center_cloud (new typename pcl::PointCloud<pcl::PointXYZ>);
    center_cloud->resize(d.triangles.size() / 3);

    // compute triangle centers
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        // vertcies index
        int v1_index = d.triangles[i];
        int v2_index = d.triangles[i + 1];
        int v3_index = d.triangles[i + 2];

        // vertices
        Eigen::Vector3f v1 = vertex_cloud->points[v1_index].getVector3fMap();
        Eigen::Vector3f v2 = vertex_cloud->points[v2_index].getVector3fMap();
        Eigen::Vector3f v3 = vertex_cloud->points[v3_index].getVector3fMap();

        // center
        Eigen::Vector3f center = (v1 + v2 + v3) / 3;

        // store center
        center_cloud->points[i / 3].x = center(0);
        center_cloud->points[i / 3].y = center(1);
        center_cloud->points[i / 3].z = center(2);
    }

    // return
    return center_cloud;
}

// compute triangle center to vertices index map
// map[triangle_index in d] = [v1_index, v2_index, v3_index in vertex_cloud]
std::map<int, std::vector<int>> 
compute_center_to_vertices_index_map(delaunator::Delaunator d)
{
    // initialize
    std::map<int, std::vector<int>> map;
    
    // compute
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        // vertcies index
        int v1_index = d.triangles[i];
        int v2_index = d.triangles[i + 1];
        int v3_index = d.triangles[i + 2];

        // store center to vertex indices map
        map[i / 3].push_back(v1_index);
        map[i / 3].push_back(v2_index);
        map[i / 3].push_back(v3_index);
    }

    // return
    return map;
}

// triangulation
// triangles list stores [trig1-vertex1, trig1-vertex2, trig1-vertex3, trig2-vertex1, trig2-vertex2, trig2-vertex3, ...]
template<typename PointT>
delaunator::Delaunator 
obtain_triangulation(typename pcl::PointCloud<PointT>::Ptr cloud)
{

    // compute azimuth and altitude coordinates
    std::vector<double>* coords = new std::vector<double>();
    coords->reserve(cloud->size() * 2);
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        float x = cloud->points[i].x;
        float y = cloud->points[i].y;
        float z = cloud->points[i].z;
        
        double r = sqrt(x * x + y * y + z * z);
        double azimuth = atan2(y, x) * 180 / M_PI;
        double altitude = asin(z / r) * 180 / M_PI;

        coords->push_back(azimuth);
        coords->push_back(altitude);
    }

    // create delaunay triangles given x and y
    delaunator::Delaunator d(*coords);
    // std::cout << "Number of triangles: " << d.triangles.size() << std::endl;

    return d;
}


// https://stackoverflow.com/questions/11719168/how-do-i-find-the-sphere-center-from-3-points-and-radius
bool compute_sphere_centers(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3, float r, Eigen::Vector3f &sphere_center1, Eigen::Vector3f &sphere_center2)
{
    // edge and normal
    Eigen::Vector3f e1 = p2 - p1;
    Eigen::Vector3f e2 = p3 - p1;
    Eigen::Vector3f normal = e1.cross(e2);

    // circle center
    Eigen::Vector3f c = (e1.squaredNorm() * e2 - e2.squaredNorm() * e1).cross(normal) / ( 2 * normal.squaredNorm() ) + p1;

    // compute t
    Eigen::Vector3f cp1 = c - p1;
    float t_squared = (r*r - cp1.squaredNorm()) / normal.squaredNorm();
    if (t_squared < 0) return false;
    float t = std::sqrt(t_squared);

    // compute sphere center
    sphere_center1 = c + normal * t;
    sphere_center2 = c - normal * t;
    return true;
}


// https://stackoverflow.com/questions/5883169/intersection-between-a-line-and-a-sphere
bool line_sphere_intersection(Eigen::Vector3f point, Eigen::Vector3f direction, Eigen::Vector3f sphere_center, double sphere_radius, Eigen::Vector3f& solution1, Eigen::Vector3f& solution2)
{
    // http://www.codeproject.com/Articles/19799/Simple-Ray-Tracing-in-C-Part-II-Triangles-Intersec

    // compute A B C
    double cx = sphere_center[0];
    double cy = sphere_center[1];
    double cz = sphere_center[2];
    double r = sphere_radius;
    double px = point[0];
    double py = point[1];
    double pz = point[2];
    double dx = direction[0];
    double dy = direction[1];
    double dz = direction[2];
    double A = dx * dx + dy * dy + dz * dz;
    double B = 2.0 * (px * dx + py * dy + pz * dz - dx * cx - dy * cy - dz * cz);
    double C = px * px - 2.0 * px * cx + cx * cx + py * py - 2.0 * py * cy + cy * cy + pz * pz - 2.0 * pz * cz + cz * cz - r * r;

    // discriminant
    double det = B * B - 4.0 * A * C;
    if ( det < 0 ) return false; // no solution

    // compute solutions
    double t1 = 1/(2.0*A) * (- B - std::sqrt(det));
    double t2 = 1/(2.0*A) * (- B + std::sqrt(det));
    solution1 = point + t1 * direction;    
    solution2 = point + t2 * direction;

    return true;
}


bool eye_patch_intersection(Eigen::Vector3f current_point,
                            Eigen::Vector3f current_point_direction,
                            Eigen::Vector3f v1, 
                            Eigen::Vector3f v2,
                            Eigen::Vector3f v3,
                            double vertex_variance,
                            double sphere_radius,
                            Eigen::Vector3f& likelihood_point,
                            double& likelihood_variance
                            )
{
    
    Eigen::Vector3f sphere_center1, sphere_center2;
    bool sphere_exists = compute_sphere_centers(v1, v2, v3, sphere_radius, sphere_center1, sphere_center2);
    if (!sphere_exists)
    {
        return false;
    }

    // compute line and sphere intersection
    Eigen::Vector3f intersection11, intersection12;
    Eigen::Vector3f intersection21, intersection22;
    bool intersection_exists1 = line_sphere_intersection(current_point, current_point_direction, sphere_center1, sphere_radius, intersection11, intersection12);
    bool intersection_exists2 = line_sphere_intersection(current_point, current_point_direction, sphere_center2, sphere_radius, intersection21, intersection22);
    if (!intersection_exists1 || !intersection_exists2)
    {
        return false;
    }

    // find the two inner intersection
    Eigen::Vector3f inner_intersection1, inner_intersection2;
    double distance1121 = (intersection11 - intersection21).norm();
    double distance1122 = (intersection11 - intersection22).norm();
    double distance1221 = (intersection12 - intersection21).norm();
    double distance1222 = (intersection12 - intersection22).norm();
    if (distance1121 < distance1122 && distance1121 < distance1221 && distance1121 < distance1222)
    {
        inner_intersection1 = intersection11;
        inner_intersection2 = intersection21;
    }
    else if (distance1122 < distance1121 && distance1122 < distance1221 && distance1122 < distance1222)
    {
        inner_intersection1 = intersection11;
        inner_intersection2 = intersection22;
    }
    else if (distance1221 < distance1121 && distance1221 < distance1122 && distance1221 < distance1222)
    {
        inner_intersection1 = intersection12;
        inner_intersection2 = intersection21;
    }
    else
    {
        inner_intersection1 = intersection12;
        inner_intersection2 = intersection22;
    }

    // compute variance
    double uniform_variance = (inner_intersection1 - inner_intersection2).squaredNorm() / 12;
    double noise_variance = vertex_variance;

    // assignment
    likelihood_point = (inner_intersection1 + inner_intersection2) / 2;
    likelihood_variance = uniform_variance + noise_variance;

    return true;
}


// compute updated pointcloud
template <typename PointT>
void 
update_pointcloud(
    typename pcl::PointCloud<PointT>::Ptr old_cloud, 
    const std::vector<Eigen::Vector3f> old_cloud_direction, 
    std::vector<float>& old_cloud_variance, 
    const typename pcl::PointCloud<PointT>::Ptr new_cloud, 
    std::vector<int>& used_points,
    double range_std
    )
{
    // triangulate new cloud
    delaunator::Delaunator d = obtain_triangulation<PointT>(new_cloud);
    
    // compute triangle centers 
    pcl::PointCloud<pcl::PointXYZ>::Ptr center_cloud = computer_triangle_center<PointT>(new_cloud, d);

    // compute triangle center to vertices index map
    std::map<int, std::vector<int>> center_to_vertices_index_map = compute_center_to_vertices_index_map(d);
    
    // prepare kd tree search
    pcl::KdTreeFLANN<pcl::PointXY> new_cloud_triangle_center_polar_kdtree;
    new_cloud_triangle_center_polar_kdtree.setInputCloud(obtain_2d_polar_cloud<pcl::PointXYZ>(center_cloud));
    
    // iterate through polar of old cloud
    std::vector<int> used_triangles;
    pcl::PointCloud<pcl::PointXY>::Ptr old_cloud_polar = obtain_2d_polar_cloud<PointT>(old_cloud);
    for (std::size_t i = 0; i < old_cloud_polar->size(); i++)
    {
        // current point
        Eigen::Vector3f current_point = old_cloud->points[i].getVector3fMap();
        pcl::PointXY current_point_polar = old_cloud_polar->points[i];
        Eigen::Vector3f current_point_direction = old_cloud_direction[i];
        
        // 1. search closest triangle center
        int K = 4;
        pcl::Indices center_indices_searched = get_k_nearest_neighbor(new_cloud_triangle_center_polar_kdtree, current_point_polar, K);
        
        // 2. find intersections to those triangles
        std::vector<int> all_intersected_triangle_indices;
        std::vector<Eigen::Vector3f> all_intersections;
        std::vector<double> all_likelihood_variance;
        for (int center_index : center_indices_searched)
        {
            // get triangle vertices xyz
            Eigen::Vector3f v1 = new_cloud->points[center_to_vertices_index_map[center_index][0]].getVector3fMap();
            Eigen::Vector3f v2 = new_cloud->points[center_to_vertices_index_map[center_index][1]].getVector3fMap();
            Eigen::Vector3f v3 = new_cloud->points[center_to_vertices_index_map[center_index][2]].getVector3fMap();

            // compute intersection
            Eigen::Vector3f likelihood_point;
            double likelihood_variance;
            double sphere_radius = 0.05;
            bool intersection_exists = eye_patch_intersection(current_point, current_point_direction, v1, v2, v3, std::pow(range_std, 2), sphere_radius, likelihood_point, likelihood_variance);
            if (!intersection_exists) continue;

            // check if intersection is inside the triangle
            bool inside = is_inside_triangle(v1, v2, v3, likelihood_point);
            if (inside)
            {
                all_intersected_triangle_indices.push_back(center_index);
                all_intersections.push_back(likelihood_point);
                all_likelihood_variance.push_back(likelihood_variance);
            }
        }
        bool no_intersections = all_intersections.size() == 0;
        if (no_intersections) continue;

        // 3. find cloest intersection
        auto iterator = std::min_element(all_intersections.begin(), all_intersections.end(),
            [current_point](const Eigen::Vector3f& a, const Eigen::Vector3f& b) {return (a - current_point).norm() < (b - current_point).norm();});
        int closest_intersected_triangle_index = all_intersected_triangle_indices[std::distance(all_intersections.begin(), iterator)];
        double closest_intersected_triangle_likelihood_variance = all_likelihood_variance[std::distance(all_intersections.begin(), iterator)];
        Eigen::Vector3f closest_likelihood_point = *iterator;
        
        // 4. compute posterior
        Eigen::Vector3f prior_point = current_point;
        float prior_variance = old_cloud_variance[i];

        Eigen::Vector3f likelihood_point = closest_likelihood_point;
        float likelihood_variance = closest_intersected_triangle_likelihood_variance;

        Eigen::Vector3f posterior_point = (prior_point / prior_variance + likelihood_point / likelihood_variance) / (1 / prior_variance + 1 / likelihood_variance);
        float posterior_variance = 1 / (1 / prior_variance + 1 / likelihood_variance);

        float update_threshold_distance = 3 * std::pow(prior_variance, 0.5) + 3 * std::pow(likelihood_variance, 0.5);
        bool too_far = (prior_point - likelihood_point).norm() > update_threshold_distance;
        if (too_far) continue;


        // 5. update old cloud
        old_cloud->points[i].x = posterior_point(0);
        old_cloud->points[i].y = posterior_point(1);
        old_cloud->points[i].z = posterior_point(2);
        old_cloud_variance[i] = posterior_variance;
        used_triangles.push_back(closest_intersected_triangle_index); // record used triangles
    }

    // record used_points
    for (int used_triangle : used_triangles)
    {
        for (int vertex_index : center_to_vertices_index_map[used_triangle])
        {
            used_points.push_back(vertex_index);
        }
    }
}

// unique a vector
template<typename T>
void 
unique_vector(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}


// remove used points
template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
remove_used_points(typename pcl::PointCloud<PointT>::Ptr cloud, std::vector<int> used_points)
{
    // convert into indices
    pcl::PointIndices::Ptr used_points_indices (new pcl::PointIndices);
    used_points_indices->indices = used_points;

    // remove used points
    pcl::ExtractIndices<PointT> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(used_points_indices);
    extract.setNegative(true);
    typename pcl::PointCloud<PointT>::Ptr cloud_copied (new typename pcl::PointCloud<PointT>);
    extract.filter(*cloud_copied);
    
    return cloud_copied;
}

// remove used direction
std::vector<Eigen::Vector3f> 
remove_used_directions(std::vector<Eigen::Vector3f> directions, std::vector<int> used_points)
{
    std::vector<Eigen::Vector3f> directions_copied;
    for (std::size_t i = 0; i < directions.size(); i++)
    {
        // if i not in used points
        if (std::find(used_points.begin(), used_points.end(), i) == used_points.end())
        {
            directions_copied.push_back(directions[i]);
        }
    }
    return directions_copied;
}

// remove used variance
std::vector<float> 
remove_used_variance(std::vector<float> variances, std::vector<int> used_points)
{
    std::vector<float> variances_copied;
    for (std::size_t i = 0; i < variances.size(); i++)
    {
        // if i not in used points
        if (std::find(used_points.begin(), used_points.end(), i) == used_points.end())
        {
            variances_copied.push_back(variances[i]);
        }
    }
    return variances_copied;
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
{
    // transform point cloud
    typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
    pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
    return transformed_cloud;
}

// convert to target frame using pose
template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_to_frame(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d current_pose, Eigen::Affine3d target_pose)
{
    Eigen::Affine3d pose_current2target = target_pose.inverse() * current_pose; // the cloud will first be transformed to current pose, then to target pose
    return transform_to_global<PointT>(cloud, pose_current2target);
}


template <typename PointT>
class Algorithm
{
public:
    // constructor, destructor
    Algorithm(double range_std) : sensor_range_std_(range_std), old_cloud_(new typename pcl::PointCloud<PointT>){};

    // input
    void add_pointcloud_and_pose(typename pcl::PointCloud<PointT>::Ptr new_cloud, Eigen::Affine3d new_pose)
    {
        // compute near cloud and far for new cloud
        std::vector<Eigen::Vector3f> new_cloud_direction = compute_point_directions<PointT>(new_cloud);
        std::vector<float> new_cloud_variance(new_cloud->size(), std::pow(sensor_range_std_, 2));

        // transform old cloud to local frame
        typename pcl::PointCloud<PointT>::Ptr old_cloud_local = transform_to_frame<PointT>(old_cloud_, Eigen::Isometry3d::Identity(), new_pose);

        // collect used point, update old cloud
        std::vector<int> used_points;
        update_pointcloud<PointT>(old_cloud_local, old_cloud_direction_, old_cloud_variance_, new_cloud, used_points, sensor_range_std_);
        unique_vector<int>(used_points);
        
        // add ununsed points / directions to old cloud
        typename pcl::PointCloud<PointT>::Ptr filtered_new_cloud = remove_used_points<PointT>(new_cloud, used_points);
        *old_cloud_local += *filtered_new_cloud;

        std::vector<Eigen::Vector3f> filtered_directions = remove_used_directions(new_cloud_direction, used_points);
        old_cloud_direction_.insert(old_cloud_direction_.end(), filtered_directions.begin(), filtered_directions.end());

        std::vector<float> filtered_varaince = remove_used_variance(new_cloud_variance, used_points);
        old_cloud_variance_.insert(old_cloud_variance_.end(), filtered_varaince.begin(), filtered_varaince.end());

        // transform old cloud to global frame
        typename pcl::PointCloud<PointT>::Ptr old_cloud_global = transform_to_frame<PointT>(old_cloud_local, new_pose, Eigen::Isometry3d::Identity());

        // update old cloud
        *old_cloud_ = *old_cloud_global;
    }

    // output
    typename pcl::PointCloud<PointT>::Ptr get_old_cloud()
    {
        return old_cloud_;
    }

private:
    // sensor parameters
    double sensor_range_std_;

    // algorithm data storage
    typename pcl::PointCloud<PointT>::Ptr old_cloud_;
    std::vector<Eigen::Vector3f> old_cloud_direction_;
    std::vector<float> old_cloud_variance_;
};