#include "MeshObject/Application.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Storage.hpp"
#include "utilities/utilities.hpp"
#include "utilities/covariance_math.hpp"

#include <iostream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include "MeshObject/MeshObject.hpp"

template class Application<VilensPointT>;

template <typename PointT>
Application<PointT>::Application() 
{
    // Initialization code
    std::map<std::string, std::pair<std::string, std::string>> dataset_map;
    dataset_map["room"] = std::make_pair(
        "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/",
        "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam"
    );
    dataset_map["osney"] = std::make_pair(
        "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/",
        "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam"
    );
    dataset_map["blenheim"] = std::make_pair(
        "/home/jiahao/datasets/2024-03-14-09-09-02-lenord-walk-for-lintong/individual_clouds/",
        "/home/jiahao/datasets/2024-03-14-09-09-02-lenord-walk-for-lintong/slam_pose_graph.g2o"
    );

    distance_threshold = 0.05;
    fit_plane_threshold = 10;
    merged_eigenvalue_threshold = 15e-5;
    dataset = "room";
    ith_cloud = 50;
    ith_point = 0;
    shuffle_pointcloud = false;
    pointcloud_fraction = 1;
    distance_to_radius_ratio = tan(4 * M_PI / 180);

    storage_ = std::make_shared<Storage>();
    data_loader.load_dataset(dataset_map.at(dataset).first, dataset_map.at(dataset).second);
    load_point_cloud();
}

template <typename PointT>
Eigen::Matrix3d Application<PointT>::merge_covariances_of_surfaces(std::shared_ptr<Surface> surface1, std::shared_ptr<Surface> surface2) 
{
    const Eigen::Matrix3d& cov1 = surface1->get_covariance();
    const Eigen::Matrix3d& cov2 = surface2->get_covariance();
    const Eigen::Vector3d& mean1 = surface1->get_mean();
    const Eigen::Vector3d& mean2 = surface2->get_mean();
    int size1 = surface1->get_total_point_size();
    int size2 = surface2->get_total_point_size();
    return merge_covariance(cov1, cov2, mean1, mean2, size1, size2);
}

// compute eigen value of merged surfaces
template <typename PointT>
double Application<PointT>::compute_eigenvalue_of_merged_surfaces(std::shared_ptr<Surface> surface1, std::shared_ptr<Surface> surface2)
{
    Eigen::Matrix3d covariance_matrix = merge_covariances_of_surfaces(surface1, surface2);
    return Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance_matrix).eigenvalues()[0];
}

template <typename PointT>
void Application<PointT>::try_merge_surfaces(std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& surfaces_to_merge)
{
    while (true) 
    {
        std::set<std::pair<std::shared_ptr<Surface>, std::shared_ptr<Surface>>> surface_pairs;
        for (std::shared_ptr<Surface> surface1 : surfaces_to_merge) 
        {
            for (std::shared_ptr<Surface> surface2 : surfaces_to_merge) 
            {
                if (surface1 >= surface2) continue;
                surface_pairs.insert(std::make_pair(surface1, surface2));
            }
        }
        
        bool again = false;
        for (const auto& pairs : surface_pairs) 
        {
            Eigen::Matrix3d covariance_matrix = merge_covariances_of_surfaces(pairs.first, pairs.second);
            double eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance_matrix).eigenvalues()[0];
            if (eigenvalue > merged_eigenvalue_threshold) continue;

            surfaces_to_merge.erase(pairs.second);
            pairs.first->merge_surface(pairs.second);

            again = true;
            break;
        }
        if (!again) break;
    }
}

template <typename PointT>
void Application<PointT>::add_point_by_radius_search(const std::shared_ptr<GenericPoint>& generic_point)
{
    // the reverse radius search provides a set of vertices that think the new point should be in
    // the new point to surface fit is then performed -> represented by point to surface projective distance, point uncertainty and plane uncertainty
    // by adding a point to a surface, we are essentially expanding the surface's boundary

    // each vertex can have multiple surfaces, but each edge and face and interiror point can only have one surface

    // when a surface have very few points, the surface have high self uncertainty, and the point is likely to be added to the surface with few points
    // if the new point can be added into multiple surfaces, but not all surfaces, added to the surfaces, reduce the search radius of the neighboring vertices from different surfaces
    // if the new point can be added into one surface, add it to that surface, reduce the search radius of the neighboring vertices from different surfaces
    // if the new point can not be added into any surface, create a new surface and add it to that surface, reduce the search radius of the neighboring vertices from different surfaces

    // what about surface merging
    // after adding the new point to surfaces, if the new point have multiple surfaces, try merge them

    // create new vertex
    std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(generic_point);

    // when can not search
    if (!storage_->can_reverse_radius_search())
    {
        std::shared_ptr<Surface> new_surface = storage_->add_surface();
        new_surface->connect(new_vertex);
        return;
    }

    // get neighboring vertices
    std::map<int, double> point_to_radius_map;
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices = storage_->reverse_radius_search(new_vertex);

    // when no search results
    if (neighboring_vertices.size() == 0)
    {
        std::shared_ptr<Surface> new_surface = storage_->add_surface();
        new_surface->connect(new_vertex);
        return;
    }

    // get neighboring surfaces
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> neighboring_surfaces; 
    for (std::shared_ptr<Vertex> vertex : neighboring_vertices)
    {
        neighboring_surfaces.insert(vertex->get_surfaces().begin(), vertex->get_surfaces().end());
    }

    // split into surfaces to add to and surfaces to not add to
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_to_add_to;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_to_not_add_to;
    for (std::shared_ptr<Surface> surface : neighboring_surfaces) 
    {
        // add if surface is small in size
        if (surface->get_total_point_size() < fit_plane_threshold) 
        {
            surfaces_to_add_to.insert(surface);
            continue;
        }
        
        // now left with surfaces with large size

        // skip if large distance
        double distance = std::fabs(surface->compute_point_to_surface_distance(new_vertex));
        if (distance >= distance_threshold) 
        {
            surfaces_to_not_add_to.insert(surface);
            continue;
        }

        // skip if not in the same direction
        Eigen::Vector3d normal = surface->get_normal();
        Eigen::Vector3d direction = new_vertex->get_origin() - new_vertex->get_position();
        if (normal.dot(direction) < 0) 
        {
            surfaces_to_not_add_to.insert(surface);
            continue;
        }

        // now left with surfaces with large size, small distance, and in the same direction
        surfaces_to_add_to.insert(surface);
        continue;
    }

    // for surfaces not to add to
    for (std::shared_ptr<Surface> surface : surfaces_to_not_add_to)
    {
        // find neighboring vertices from the surface
        for (std::shared_ptr<Vertex> vertex : neighboring_vertices)
        {
            // skip if vertex is not in the surface
            if (vertex->get_surfaces().find(surface) == vertex->get_surfaces().end()) continue;

            // reduce the search radius
            vertex->reduce_reverse_radius_search_radius((vertex->get_position() - new_vertex->get_position()).norm());
        }
    }

    // for surfaces to add to
    for (std::shared_ptr<Surface> surface : surfaces_to_add_to)
    {
        // add the new point to the surface
        surface->connect_by_edges_and_faces(new_vertex, neighboring_vertices);
    }

    // if new_vertex is in multiple surfaces, try merge them
    if (new_vertex->get_surfaces().size() > 1)
    {
        new_vertex->try_merge_surfaces();
    }

    // if new_vertex is not in any surface, create a new surface and add it to that surface
    if (new_vertex->get_surfaces().size() == 0)
    {
        std::shared_ptr<Surface> new_surface = storage_->add_surface();
        new_surface->connect(new_vertex);
    }

    // // try refine the candidate surface this point added to
    // if (candidate_surfaces.size() > 0) 
    // {
    //     if (new_vertex->get_surface()->get_eigenvalues()[0] > merged_eigenvalue_threshold) 
    //     {
    //         new_vertex->get_surface()->refine_surface();
    //     }
    // }
}

template <typename PointT>
void Application<PointT>::load_point_cloud()
{
    if (ith_cloud < 0)
    {
        std::cout << "reached the first pointcloud" << std::endl;
        ith_cloud = 0;
    }
    if (ith_cloud >= data_loader.size())
    {
        std::cout << "reached the last pointcloud" << std::endl;
        ith_cloud = data_loader.size() - 1;
    }
    
    typename pcl::PointCloud<PointT>::Ptr pointcloud_local = data_loader.get_cloud(ith_cloud);
    Eigen::Affine3d pose = data_loader.get_pose(ith_cloud);
    pointcloud = transform_cloud_to_global<PointT>(pointcloud_local, pose);
    origin = pose.translation();
    ith_size = pointcloud->size() * pointcloud_fraction;

    std::cout << "loaded pointcloud " << ith_cloud << " with " << pointcloud->size() << " points" << std::endl;

    if (shuffle_pointcloud) 
    {
        std::random_shuffle(pointcloud->points.begin(), pointcloud->points.end());
    }
}

template <typename PointT>
void Application<PointT>::process_point(const std::shared_ptr<GenericPoint>& generic_point)
{
    // the face search provides a set of faces that are pointed by the new point

    // if the new face is penetrated by the new point, the penetrated face should be deleted
    // if the new face is behind the new point, nothing should be done
    // if the new face can contain the new point, add the new point to the face

    // if the new point can be added into multiple faces, if these faces are from the same surface
    // if the new point can be added into multiple faces, if these faces are from different surface

    // by adding a point to a surface, we are essentially refining the surface's plane estimate

    // searched faces
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> searched_faces = storage_->face_intersection_search(generic_point);

    // searched surfaces
    std::map<std::shared_ptr<Surface>, std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>> searched_surface_to_searched_faces;
    for (const std::shared_ptr<Face>& face : searched_faces)
    {
        searched_surface_to_searched_faces[face->get_surface()].insert(face);
    }

    bool point_added = false;
    for (const auto& pair : searched_surface_to_searched_faces)
    {
        const std::shared_ptr<Surface>& surface = pair.first;
        const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& searched_faces = pair.second;

        double distance = surface->compute_point_to_surface_distance(generic_point);
        bool points_before_surface = distance > distance_threshold;
        bool points_behind_surface = distance < -distance_threshold;
        bool points_within_surface = !points_before_surface && !points_behind_surface;
        
        if (points_behind_surface)
        {
            storage_->set_penetrating_point(generic_point);
            for (const std::shared_ptr<Face>& face : searched_faces) storage_->delete_face(face);
            storage_->clear_penetrating_point();
        }
        else if (points_within_surface)
        {
            if (!point_added)
            {
                storage_->add_interior_point(*searched_faces.begin(), generic_point);
                point_added = true;
            }
            else
            {
                std::cout << "point within multiple surface" << std::endl;
            }
        }
        else if (points_before_surface)
        {
            continue;
        }
    }
    if (!point_added) add_point_by_radius_search(generic_point);

    if (ith_point == ith_size) 
    {   
        ith_cloud += 1;
        ith_point = 0;
        load_point_cloud();
    }
}

template <typename PointT>
void Application<PointT>::step()
{        
    Eigen::Vector3d thisPointVEC = pointcloud->points[ith_point].getVector3fMap().cast<double>();
    Eigen::Vector3d thisPointOriginVEC = origin;
    ith_point++;

    // log
    std::cout << "==================================================================== Processing point " << ith_point << " of cloud " << ith_cloud << std::endl;

    const std::shared_ptr<GenericPoint>& generic_point = storage_->add_generic_point(thisPointVEC, thisPointOriginVEC);
    process_point(generic_point);
    storage_->delete_generic_point(generic_point);
}

template <typename PointT>
void Application<PointT>::add_back_generic_points()
{
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> copy_of_generic_points = storage_->get_generic_points();
    for (const std::shared_ptr<GenericPoint>& generic_point : copy_of_generic_points)
    {
        process_point(generic_point);
        storage_->delete_generic_point(generic_point);
    }
}

template <typename PointT>
void Application<PointT>::loop()
{
    step();
    while (ith_point != 0)
    {
        step();
    }
}

template <typename PointT>
void Application<PointT>::restart()
{
    storage_ = std::make_shared<Storage>();
    ith_cloud = 50;
    ith_point = 0;
    load_point_cloud();
}

template <typename PointT>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr Application<PointT>::compute_generic_point_pointcloud()
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    for (const std::shared_ptr<GenericPoint>& generic_point : storage_->get_generic_points())
    {
        pcl::PointXYZRGB point;
        point.x = generic_point->get_position()[0];
        point.y = generic_point->get_position()[1];
        point.z = generic_point->get_position()[2];
        point.r = 255;
        point.g = 0;
        point.b = 0;
        cloud->push_back(point);
    }
    return cloud;
}

template <typename PointT>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr Application<PointT>::compute_interior_point_pointcloud(bool show_projected_point, bool show_error_color)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    for (const std::shared_ptr<InteriorPoint>& interior_point : storage_->get_interior_points())
    {
        pcl::PointXYZRGB point;
        if (show_projected_point)
        {
            point.x = interior_point->get_projected_position()[0];
            point.y = interior_point->get_projected_position()[1];
            point.z = interior_point->get_projected_position()[2];
        }
        else
        {
            point.x = interior_point->get_position()[0];
            point.y = interior_point->get_position()[1];
            point.z = interior_point->get_position()[2];
        }
        if (show_error_color)
        {
            double distance = std::fabs(interior_point->get_projected_distance() / 0.05);
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else
        {
            const std::tuple<int, int, int>& color = interior_point->get_surface()->get_color();
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        cloud->push_back(point);
    }
    return cloud;
}

template <typename PointT>
std::map<std::shared_ptr<Vertex>, int> Application<PointT>::get_vertex_to_cloud_indices_map()
{
    return vertex_to_cloud_indices_map;
} 

template <typename PointT>
const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Application<PointT>::get_faces() 
{
    return storage_->get_faces();
}

template <typename PointT>
const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Application<PointT>::get_edges() 
{
    return storage_->get_edges();
}

template <typename PointT>
std::vector<std::shared_ptr<Vertex>> Application<PointT>::get_rrs_vertices() 
{
    return storage_->get_rrs_vertices();
}

template <typename PointT>
std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> Application<PointT>::get_boundary_edges() 
{
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> boundary_edges;
    for (const std::shared_ptr<Edge>& edge : storage_->get_edges())
    {
        if (edge->is_boundary()) 
        {
            boundary_edges.insert(edge);
        }
    }
    return boundary_edges;
}

template <typename PointT>
void Application<PointT>::refine_surfaces()
{
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> copy_of_surfaces = storage_->get_surfaces();
    for (const std::shared_ptr<Surface>& surface : copy_of_surfaces)
    {
        surface->refine_surface();
    }
    std::cout << "number of generic points after refine: " << storage_->get_generic_points().size() << std::endl;
}

template <typename PointT>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr Application<PointT>::compute_vertex_point_pointcloud(bool show_projected_point, bool show_error_color)
{
    vertex_to_cloud_indices_map.clear();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    for (const std::shared_ptr<Vertex>& vertex : storage_->get_vertices())
    {
        pcl::PointXYZRGB point;
        if (show_projected_point)
        {
            point.x = vertex->get_projected_position()[0];
            point.y = vertex->get_projected_position()[1];
            point.z = vertex->get_projected_position()[2];
        }
        else
        {
            point.x = vertex->get_position()[0];
            point.y = vertex->get_position()[1];
            point.z = vertex->get_position()[2];
        }
        if (show_error_color)
        {
            double distance = std::abs(vertex->get_projected_distance() / 0.05);
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else
        {
            const std::tuple<int, int, int>& color = vertex->get_surface()->get_color();
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        cloud->push_back(point);
        vertex_to_cloud_indices_map[vertex] = cloud->size() - 1;
    }
    return cloud;
}

template <typename PointT>
void Application<PointT>::change_color()
{
    for (const std::shared_ptr<Surface>& surface : storage_->get_surfaces())
    {
        surface->set_random_color();
    }
}