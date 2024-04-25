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

void compute_azimuth_and_altitude(const Eigen::Vector3f& point, float& azimuth, float& altitude)
{
    // x y z
    float x = point[0];
    float y = point[1];
    float z = point[2];

    // r azimuth altitude
    double r = sqrt(x * x + y * y + z * z);
    azimuth = atan2(y, x) * 180 / M_PI;
    altitude = asin(z / r) * 180 / M_PI;
}

// convert to 2d polar cloud
template<typename PointT>
pcl::PointCloud<pcl::PointXY>::Ptr 
compute_2d_polar_cloud(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    // initialize
    pcl::PointCloud<pcl::PointXY>::Ptr cloud_polar (new pcl::PointCloud<pcl::PointXY>);
    cloud_polar->resize(cloud->size());

    // process
    for (std::size_t i = 0; i < cloud->size(); i+=1)
    {
        // initialize
        float azimuth, altitude;

        // process
        compute_azimuth_and_altitude(cloud->points[i].getVector3fMap(), azimuth, altitude);

        // store
        cloud_polar->points[i].x = azimuth;
        cloud_polar->points[i].y = altitude;
    }

    // return
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
d_to_triangle_map(delaunator::Delaunator d)
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

// compute triangle center
template<typename PointT>
typename pcl::PointCloud<pcl::PointXYZ>::Ptr 
compute_triangle_center_cloud(typename pcl::PointCloud<PointT>::Ptr cloud, std::map<int, std::vector<int>> triangle_map)
{
    // initialize
    typename pcl::PointCloud<pcl::PointXYZ>::Ptr center_cloud (new typename pcl::PointCloud<pcl::PointXYZ>);

    // compute triangle centers
    for (const auto& entry : triangle_map)
    {
        // vertices index
        int v1_index = entry.second[0];
        int v2_index = entry.second[1];
        int v3_index = entry.second[2];
        
        // vertices vector
        Eigen::Vector3f v1 = cloud->points[v1_index].getVector3fMap();
        Eigen::Vector3f v2 = cloud->points[v2_index].getVector3fMap();
        Eigen::Vector3f v3 = cloud->points[v3_index].getVector3fMap();

        // center vector
        Eigen::Vector3f center = (v1 + v2 + v3) / 3;

        // center point
        pcl::PointXYZ center_point(center(0), center(1), center(2));

        // store
        center_cloud->push_back(center_point);
    }

    // return
    return center_cloud;
}

template<typename PointT>
typename pcl::PointCloud<pcl::PointXYZ>::Ptr 
compute_triangle_center_cloud(typename pcl::PointCloud<PointT>::Ptr cloud, delaunator::Delaunator d)
{
    return compute_triangle_center_cloud<PointT>(cloud, d_to_triangle_map(d));
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
        // std::cout << "Sphere does not exist" << std::endl;
        return false;
    }

    // compute line and sphere intersection
    Eigen::Vector3f intersection11, intersection12;
    Eigen::Vector3f intersection21, intersection22;
    bool intersection_exists1 = line_sphere_intersection(current_point, current_point_direction, sphere_center1, sphere_radius, intersection11, intersection12);
    bool intersection_exists2 = line_sphere_intersection(current_point, current_point_direction, sphere_center2, sphere_radius, intersection21, intersection22);
    if (!intersection_exists1 || !intersection_exists2)
    {
        // std::cout << "Ray does not intersect with sphere" << std::endl;
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


// compute vertex normal using discrete vector area
template <typename PointT>
std::map<int, Eigen::Vector3f> compute_vertex_to_normal_map(typename pcl::PointCloud<PointT>::Ptr cloud, delaunator::Delaunator d)
{   
    // initialize
    std::map<int, Eigen::Vector3f> vertex_to_normal_map;

    // compute normal
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        // vertex indices
        int v1_index = d.triangles[i];
        int v2_index = d.triangles[i + 1];
        int v3_index = d.triangles[i + 2];
        // vertex vectors
        Eigen::Vector3f v1 = cloud->points[v1_index].getVector3fMap();
        Eigen::Vector3f v2 = cloud->points[v2_index].getVector3fMap();
        Eigen::Vector3f v3 = cloud->points[v3_index].getVector3fMap();

        // find the order and push cross products
        Eigen::Vector3f v21 = v2 - v1;
        Eigen::Vector3f v31 = v3 - v1;
        // order
        if (v21.cross(v31).dot(v1) < 0)
        {
            vertex_to_normal_map[v1_index] += v2.cross(v3);
            vertex_to_normal_map[v2_index] += v3.cross(v1);
            vertex_to_normal_map[v3_index] += v1.cross(v2);
        }
        else
        {
            vertex_to_normal_map[v1_index] += v3.cross(v2);
            vertex_to_normal_map[v2_index] += v1.cross(v3);
            vertex_to_normal_map[v3_index] += v2.cross(v1);
        }
    }
    for (auto& [vertex_index, normal] : vertex_to_normal_map)
    {
        normal.normalize();
    }
    
    // return
    return vertex_to_normal_map;
}

float compute_edge_curvature(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f n1, Eigen::Vector3f n2)
{
    // compute curvature
    Eigen::Vector3f p = p2 - p1;
    Eigen::Vector3f n = n2 - n1;
    float curvature = (n).dot(p) / p.squaredNorm();
    return curvature;
}


// compute updated pointcloud
template <typename PointT>
void 
update_pointcloud(
    typename pcl::PointCloud<PointT>::Ptr old_cloud, 
    std::vector<Eigen::Vector3f>& old_cloud_direction, 
    std::vector<float>& old_cloud_variance, 
    const typename pcl::PointCloud<PointT>::Ptr new_cloud, 
    const std::vector<Eigen::Vector3f> new_cloud_direction, 
    const std::vector<float> new_cloud_variance, 
    double range_std
    )
{
    // triangulate new cloud
    delaunator::Delaunator d = obtain_triangulation<PointT>(new_cloud);

    // compute vertex to normal map
    std::map<int, Eigen::Vector3f> vertex_to_normal_map = compute_vertex_to_normal_map<PointT>(new_cloud, d);
    
    // compute triangle centers 
    pcl::PointCloud<pcl::PointXYZ>::Ptr center_cloud = compute_triangle_center_cloud<PointT>(new_cloud, d);

    // compute triangle center to vertices index map
    std::map<int, std::vector<int>> center_to_vertices_index_map = d_to_triangle_map(d);
    
    // new point to old point map
    std::map<int, std::vector<int>> new_to_updated_old_point_map;

    // prepare kd tree search
    pcl::KdTreeFLANN<pcl::PointXY> new_cloud_triangle_center_polar_kdtree;
    new_cloud_triangle_center_polar_kdtree.setInputCloud(compute_2d_polar_cloud<pcl::PointXYZ>(center_cloud));
    
    // update old cloud
    pcl::PointCloud<pcl::PointXY>::Ptr old_cloud_polar = compute_2d_polar_cloud<PointT>(old_cloud);
    for (std::size_t i = 0; i < old_cloud_polar->size(); i++)
    {
        // current point
        Eigen::Vector3f current_point = old_cloud->points[i].getVector3fMap();
        pcl::PointXY current_point_polar = old_cloud_polar->points[i];
        Eigen::Vector3f current_point_direction = old_cloud_direction[i];
        
        // 1. search closest triangle center
        int K = 4;
        pcl::Indices center_indices_searched = get_k_nearest_neighbor(new_cloud_triangle_center_polar_kdtree, current_point_polar, K); // the returned indices is ordered by closest
        
        // 2. find intersections to those triangles
        std::vector<int> all_intersected_triangle_indices;
        std::vector<Eigen::Vector3f> all_likelihood_points;
        std::vector<double> all_likelihood_variance;
        for (int center_index : center_indices_searched)
        {
            // get triangle vertices index
            int v1_index = center_to_vertices_index_map[center_index][0];
            int v2_index = center_to_vertices_index_map[center_index][1];
            int v3_index = center_to_vertices_index_map[center_index][2];

            // get triangle vertices xyz
            Eigen::Vector3f v1 = new_cloud->points[v1_index].getVector3fMap();
            Eigen::Vector3f v2 = new_cloud->points[v2_index].getVector3fMap();
            Eigen::Vector3f v3 = new_cloud->points[v3_index].getVector3fMap();

            // get triangle curvature
            Eigen::Vector3f v1_normal = vertex_to_normal_map[v1_index];
            Eigen::Vector3f v2_normal = vertex_to_normal_map[v2_index];
            Eigen::Vector3f v3_normal = vertex_to_normal_map[v3_index];

            float curvature_12 = compute_edge_curvature(v1, v2, v1_normal, v2_normal);
            float curvature_23 = compute_edge_curvature(v2, v3, v2_normal, v3_normal);
            float curvature_31 = compute_edge_curvature(v3, v1, v3_normal, v1_normal);

            float mean_curvature = (curvature_12 + curvature_23 + curvature_31) / 3;

            // get sphere radius
            float sphere_radius = 1 / mean_curvature;

            // compute intersection
            Eigen::Vector3f likelihood_point;
            double likelihood_variance;
            bool intersection_exists = eye_patch_intersection(current_point, current_point_direction, v1, v2, v3, std::pow(range_std, 2), sphere_radius, likelihood_point, likelihood_variance);
            if (!intersection_exists) continue;

            // check if intersection is inside the triangle
            bool inside = is_inside_triangle(v1, v2, v3, likelihood_point);
            if (!inside) continue;
            
            // store intersection
            all_intersected_triangle_indices.push_back(center_index);
            all_likelihood_points.push_back(likelihood_point);
            all_likelihood_variance.push_back(likelihood_variance);
        }
        if (all_likelihood_points.size() != 1) 
        {
            continue; 
            // only proceed if there is one and only one intersectoin (which is most of the cases)
            // with more than one intersection, there is ambiguity and is better to wait for more aligned scan to update the point
        }

        // 3. compute posterior
        Eigen::Vector3f prior_point = current_point;
        float prior_variance = old_cloud_variance[i];

        Eigen::Vector3f likelihood_point = all_likelihood_points[0];
        float likelihood_variance = all_likelihood_variance[0];

        Eigen::Vector3f posterior_point = (prior_point / prior_variance + likelihood_point / likelihood_variance) / (1 / prior_variance + 1 / likelihood_variance);
        float posterior_variance = 1 / (1 / prior_variance + 1 / likelihood_variance);

        float update_threshold_distance = 3 * std::pow(prior_variance, 0.5) + 3 * std::pow(likelihood_variance, 0.5);
        bool too_far = (prior_point - likelihood_point).norm() > update_threshold_distance;
        if (too_far) continue;


        // 4. update old cloud
        old_cloud->points[i].x = posterior_point(0);
        old_cloud->points[i].y = posterior_point(1);
        old_cloud->points[i].z = posterior_point(2);
        old_cloud_variance[i] = posterior_variance;

        int intersected_triangle_index = all_intersected_triangle_indices[0];
        new_to_updated_old_point_map[center_to_vertices_index_map[intersected_triangle_index][0]].push_back(i);
        new_to_updated_old_point_map[center_to_vertices_index_map[intersected_triangle_index][1]].push_back(i);
        new_to_updated_old_point_map[center_to_vertices_index_map[intersected_triangle_index][2]].push_back(i);
    }

    // add new cloud
    for (std::size_t i = 0; i < new_cloud->size(); i++)
    {
        // get updated old points
        std::vector<int> updated_old_points_indices = new_to_updated_old_point_map[i];

        // if no old point is updated by this new point
        if (updated_old_points_indices.size() == 0) 
        {
            // add this new point to the old cloud
            old_cloud->push_back(new_cloud->points[i]);
            old_cloud_direction.push_back(new_cloud_direction[i]);
            old_cloud_variance.push_back(new_cloud_variance[i]);
            continue;
        }

        // // if there are old points updated by this new point, but the new point is far away from the old point
        // bool far_away = true;
        // for (int updated_old_point_index : updated_old_points_indices)
        // {
        //     Eigen::Vector3f updated_old_point = old_cloud->points[updated_old_point_index].getVector3fMap();
        //     double updated_old_point_variance = old_cloud_variance[updated_old_point_index];
        //     double updated_old_point_std = std::pow(updated_old_point_variance, 0.5);

        //     Eigen::Vector3f new_point = new_cloud->points[i].getVector3fMap();
        //     float distance = (updated_old_point - new_point).norm();
        //     if (distance < 3 * updated_old_point_std)
        //     {
        //         far_away = false;
        //         break;
        //     }
        // }
        // if (far_away)
        // {
        //     // add this new point to the old cloud
        //     old_cloud->push_back(new_cloud->points[i]);
        //     old_cloud_direction.push_back(new_cloud_direction[i]);
        //     old_cloud_variance.push_back(new_cloud_variance[i]);
        //     continue;
        // }
    }
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
{
    // transform point cloud
    typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
    pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
    return transformed_cloud;
}

// transform direction to global
std::vector<Eigen::Vector3f> transform_direction_to_global(const std::vector<Eigen::Vector3f>& direction, const Eigen::Affine3d& pose_eigen)
{
    // initialize
    std::vector<Eigen::Vector3f> direction_transformed = direction;

    // process
    for (std::size_t i = 0; i < direction.size(); i++)
    {
        // directional vector doesn't require translation, only rotation
        direction_transformed[i] = pose_eigen.cast<float>().rotation() * direction[i]; 
    }

    // return
    return direction_transformed;
}

// convert to target frame using pose
template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr 
transform_cloud_to_frame(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d current_pose, Eigen::Affine3d target_pose)
{
    Eigen::Affine3d pose_current2target = target_pose.inverse() * current_pose; // the cloud will first be transformed to current pose, then to target pose
    return transform_cloud_to_global<PointT>(cloud, pose_current2target);
}

// convert direction to target frame using pose
std::vector<Eigen::Vector3f>
transform_direction_to_frame(const std::vector<Eigen::Vector3f>& direction, Eigen::Affine3d current_pose, Eigen::Affine3d target_pose)
{
    Eigen::Affine3d pose_current2target = target_pose.inverse() * current_pose; // the cloud will first be transformed to current pose, then to target pose
    return transform_direction_to_global(direction, pose_current2target);
}

template <typename PointT>
class Algorithm
{
public:
    // constructor, destructor
    Algorithm(double range_std) 
        : 
        sensor_range_std_(range_std), 
        old_cloud_(new typename pcl::PointCloud<PointT>),
        control_cloud_(new typename pcl::PointCloud<PointT>)
        {};

    // input
    void input(typename pcl::PointCloud<PointT>::Ptr new_cloud, Eigen::Affine3d new_pose)
    {
        // generate new cloud direction and variance
        std::vector<Eigen::Vector3f> new_cloud_direction = compute_point_directions<PointT>(new_cloud);
        std::vector<float> new_cloud_variance(new_cloud->size(), std::pow(sensor_range_std_, 2));

        // global to local
        typename pcl::PointCloud<PointT>::Ptr old_cloud_local = transform_cloud_to_frame<PointT>(old_cloud_, Eigen::Isometry3d::Identity(), new_pose);
        std::vector<Eigen::Vector3f> old_cloud_direction_local = transform_direction_to_frame(old_cloud_direction_, Eigen::Isometry3d::Identity(), new_pose);
        
        // process
        update_pointcloud<PointT>(old_cloud_local, old_cloud_direction_local, old_cloud_variance_, new_cloud, new_cloud_direction, new_cloud_variance, sensor_range_std_);

        // local to global
        typename pcl::PointCloud<PointT>::Ptr old_cloud_global = transform_cloud_to_frame<PointT>(old_cloud_local, new_pose, Eigen::Isometry3d::Identity());
        std::vector<Eigen::Vector3f> old_cloud_direction_global = transform_direction_to_frame(old_cloud_direction_local, new_pose, Eigen::Isometry3d::Identity());

        // store
        *old_cloud_ = *old_cloud_global;
        old_cloud_direction_ = old_cloud_direction_global;
        *control_cloud_ += *transform_cloud_to_global<PointT>(new_cloud, new_pose);
    }

    // output
    typename pcl::PointCloud<PointT>::Ptr get_old_cloud()
    {
        return old_cloud_;
    }

    typename pcl::PointCloud<PointT>::Ptr get_control_cloud()
    {
        return control_cloud_;
    }

    typename pcl::PointCloud<PointT>::Ptr get_near_cloud()
    {
        typename pcl::PointCloud<PointT>::Ptr near_cloud (new pcl::PointCloud<PointT>);
        for (std::size_t i = 0; i < old_cloud_->size(); i++)
        {
            Eigen::Vector3f point = old_cloud_->points[i].getVector3fMap();
            Eigen::Vector3f direction = old_cloud_direction_[i];
            float variance = old_cloud_variance_[i];
            float shift_distance = 3 * std::pow(variance, 0.5);
            Eigen::Vector3f near_point = point - direction * shift_distance;

            PointT near_point_pcl;
            near_point_pcl.x = near_point(0);
            near_point_pcl.y = near_point(1);
            near_point_pcl.z = near_point(2);

            near_cloud->push_back(near_point_pcl);
        }

        return near_cloud;
    }

    typename pcl::PointCloud<PointT>::Ptr get_far_cloud()
    {
        typename pcl::PointCloud<PointT>::Ptr far_cloud (new pcl::PointCloud<PointT>);
        for (std::size_t i = 0; i < old_cloud_->size(); i++)
        {
            Eigen::Vector3f point = old_cloud_->points[i].getVector3fMap();
            Eigen::Vector3f direction = old_cloud_direction_[i];
            float variance = old_cloud_variance_[i];
            float shift_distance = 3 * std::pow(variance, 0.5);
            Eigen::Vector3f far_point = point + direction * shift_distance;

            PointT far_point_pcl;
            far_point_pcl.x = far_point(0);
            far_point_pcl.y = far_point(1);
            far_point_pcl.z = far_point(2);

            far_cloud->push_back(far_point_pcl);
        }

        return far_cloud;
    }

private:
    // sensor parameters
    double sensor_range_std_;

    // algorithm data storage
    typename pcl::PointCloud<PointT>::Ptr old_cloud_;
    std::vector<Eigen::Vector3f> old_cloud_direction_;
    std::vector<float> old_cloud_variance_;

    // control cloud
    typename pcl::PointCloud<PointT>::Ptr control_cloud_;
};