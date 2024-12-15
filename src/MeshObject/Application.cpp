#define PCL_NO_PRECOMPILE

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
#include "MeshObject/Simulation.hpp"

#include "utilities/omp_utilities.hpp"
#include "MeshObject/RRSTree.hpp"

#include "MeshObject/TriangleBVH.hpp"
#include <pcl/filters/passthrough.h>

#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"

template class Application<VilensPointT>;
template class Application<BagPointT>;


template <typename PointT>
Application<PointT>::Application() 
{
    restart();
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

// template <typename PointT>
// void Application<PointT>::try_merge_surfaces(std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& surfaces_to_merge)
// {
    // while (true) 
    // {
    //     std::set<std::pair<std::shared_ptr<Surface>, std::shared_ptr<Surface>>> surface_pairs;
    //     for (std::shared_ptr<Surface> surface1 : surfaces_to_merge) 
    //     {
    //         for (std::shared_ptr<Surface> surface2 : surfaces_to_merge) 
    //         {
    //             if (surface1 >= surface2) continue;
    //             surface_pairs.insert(std::make_pair(surface1, surface2));
    //         }
    //     }
        
    //     bool again = false;
    //     for (const auto& pairs : surface_pairs) 
    //     {
    //         Eigen::Matrix3d covariance_matrix = merge_covariances_of_surfaces(pairs.first, pairs.second);
    //         double eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance_matrix).eigenvalues()[0];
    //         if (eigenvalue > settings_.merged_eigenvalue_threshold) continue;

    //         surfaces_to_merge.erase(pairs.second);
    //         pairs.first->merge_surface(pairs.second);

    //         again = true;
    //         break;
    //     }
    //     if (!again) break;
    // }
// }


template <typename PointT>
void Application<PointT>::load_point_cloud()
{
    // fix ith_cloud
    if (ith_cloud < 0)
    {
        if (settings_.log.load_point_cloud) std::cout << "reached the first pointcloud" << std::endl;
        ith_cloud = 0;
    }
    if (ith_cloud >= data_loader.size())
    {
        if (settings_.log.load_point_cloud) std::cout << "reached the last pointcloud" << std::endl;
        ith_cloud = data_loader.size() - 1;
    }
    
    // get cloud and pose
    typename pcl::PointCloud<PointT>::Ptr pointcloud_local = data_loader.get_cloud(ith_cloud);
    Eigen::Affine3d pose = data_loader.get_pose(ith_cloud);

    // get number of points in the current cloud
    ith_size = pointcloud_local->size();

    // transform cloud to global
    pointcloud = transform_cloud_to_global<PointT>(pointcloud_local, pose);

    // update origin and distance traveled
    Eigen::Vector3d previous_origin = origin;
    origin = pose.translation();
    distance_travelled_ += (origin - previous_origin).norm();
    std::cout << "distance traveled: " << distance_travelled_ << std::endl;

    // shuffle pointcloud
    if (settings_.shuffle_pointcloud) 
    {
        std::random_shuffle(pointcloud->points.begin(), pointcloud->points.end());
    }

    // log
    if (settings_.log.load_point_cloud) std::cout << "loaded pointcloud " << ith_cloud << " with " << pointcloud->size() << " points" << std::endl;
}

template <typename PointT>
void Application<PointT>::process_point(const std::shared_ptr<GenericPoint>& generic_point)
{
    //
    // find new storage leaf node
    //
    std::vector<std::shared_ptr<Node>> locked_bvh_nodes;
    std::vector<std::shared_ptr<RRSNode>> locked_rrs_nodes;
    std::vector<std::shared_ptr<Surface>> locked_surfaces;
    std::set<std::shared_ptr<Surface>, MeshObjectCompare> prelocked_surfaces;

    // get candidate prelock surfaces
    std::vector<std::shared_ptr<Surface>> prelock_surface_candidates;
    for (const auto& surface : generic_point->intersected_surfaces)
    {
        prelock_surface_candidates.push_back(surface);
    }
    
    // lock
    for (const std::shared_ptr<Surface>& surface : prelock_surface_candidates)
    {
        const bool can_lock = omp_test_nest_lock(&surface->lock);
        if (can_lock)
        {
            prelocked_surfaces.insert(surface);
        }
        else
        {
            // std::cout << "X _ _ _" << std::endl;

            // increment contention count
            generic_point->contented_surfaces[surface]++;

            for (const std::shared_ptr<Surface>& surface : prelocked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");
            if (generic_point->get_contention_count() > settings_.retry_threshold)
            {
                storage_->add_to_abort_queue(generic_point);
            }
            else
            {
                storage_->add_to_queue(generic_point);
            }
            return;
        }
    }

    std::vector<std::shared_ptr<Face>> bvh_results;
    BVHReturnType BVH_return = storage_->face_intersection_search(generic_point, bvh_results);
    for (const std::shared_ptr<Face>& face : bvh_results) locked_bvh_nodes.emplace_back(face->node); // store the locked nodes
    for (const std::shared_ptr<Face>& face : bvh_results) locked_surfaces.emplace_back(face->get_surface()); // store the surface

    if (BVH_return == BVHReturnType::ABORT)
    {
        // std::cout << "_ X _ _" << std::endl;
        for (const std::shared_ptr<Surface>& surface : prelocked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");
        for (const std::shared_ptr<Surface>& surface : locked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");
        for (const std::shared_ptr<Node>& node : locked_bvh_nodes) omp_unset_nest_lock(&node->omp_lock);
        for (const std::shared_ptr<RRSNode>& node : locked_rrs_nodes) omp_unset_nest_lock(&node->omp_lock);
        
        if (generic_point->get_contention_count() > settings_.retry_threshold)
        {
            storage_->add_to_abort_queue(generic_point);
        }
        else
        {
            storage_->add_to_queue(generic_point);
        }
        return;
    }

    std::vector<std::shared_ptr<Vertex>> rrs_results;
    RRSReturnType RRS_return = storage_->reverse_radius_search(generic_point, rrs_results);    
    for (const std::shared_ptr<Vertex>& vertex : rrs_results) locked_rrs_nodes.emplace_back(vertex->node);
    for (const std::shared_ptr<Vertex>& vertex : rrs_results) locked_surfaces.emplace_back(vertex->get_surface());

    if (RRS_return == RRSReturnType::ABORT)
    {
        // std::cout << "_ _ _ X" << std::endl;
        for (const std::shared_ptr<Surface>& surface : prelocked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");
        for (const std::shared_ptr<Surface>& surface : locked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");
        for (const std::shared_ptr<Node>& node : locked_bvh_nodes) omp_unset_nest_lock(&node->omp_lock);
        for (const std::shared_ptr<RRSNode>& node : locked_rrs_nodes) omp_unset_nest_lock(&node->omp_lock);

        if (generic_point->get_contention_count() > settings_.retry_threshold)
        {
            storage_->add_to_abort_queue(generic_point);
        }
        else
        {
            storage_->add_to_queue(generic_point);
        }
        return;
    }

    // unlock prelocked surfaces
    for (const std::shared_ptr<Surface>& surface : prelocked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");

    // convert surfaces vector to set to lock nodes
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> locked_surfaces_set(locked_surfaces.begin(), locked_surfaces.end());

    // update stats
    num_of_concurrent_processes++;
    accumulated_points++;
    std::chrono::time_point<std::chrono::high_resolution_clock> t_now = std::chrono::high_resolution_clock::now();

    // compute duration
    std::chrono::duration<double> duration_accumulated = t_now - t_init;
    std::chrono::duration<double> duration_instantaneous;
    {
        // lock
        std::unique_lock<std::mutex> lock(t_last_mutex);

        // compute duration    
        duration_instantaneous = t_now - t_last;

        // update timestamp
        t_last = t_now;
    }

    // compute speed
    double speed_accumulated = accumulated_points / duration_accumulated.count();
    double speed_instantaneous = 1.0 / duration_instantaneous.count();
    
    // // output number of concurrent threads that reaches this point
    if (settings_.log.num_of_concurrent_processes)
    {
        std::stringstream ss;
        ss  << " | number of concurrent processes = " << num_of_concurrent_processes
            << " | processing point " << generic_point->get_id()
            << " | by thread " << omp_get_thread_num()
            << " | BVH tree size = " << storage_->get_bvh_size()
            << std::endl;
        std::cout << ss.str();
    }

    if (settings_.log.total_processed_points)
    {
        std::stringstream ss;
        ss  << " | total processed points = " << accumulated_points
            << " | speed accumulated = " << speed_accumulated
            << " | speed instantaneous = " << speed_instantaneous
            << std::endl;
        std::cout << ss.str();
    }
    
    // std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // process
    bool added_to_new_surface = false;
    double radius = settings_.radius_value;
    std::shared_ptr<Surface> added_surface;
    if (add_point_by_intersection_search(generic_point, radius, bvh_results, added_surface))
    {
        // if point added, go to end to unlock all locks
    }
    else if (add_point_by_radius_search(generic_point, radius, rrs_results, added_surface))
    {
        // if point added, go to end to unlock all locks
    }
    else
    {
        add_point_by_new_surface(generic_point, radius, added_surface);
        added_to_new_surface = true;
        // if point added, go to end to unlock all locks
    }

    // throw if added_surface is nullptr
    if (added_surface == nullptr) throw std::runtime_error("added surface is nullptr");

    // don't reduce radius if the point is added to a new surface
    if (!added_to_new_surface)
    {
        // rrs search radius reduction
        for (const std::shared_ptr<Vertex>& vertex : rrs_results)
        {
            // skip if expired
            if (vertex->is_expired()) continue;

            // skip if the same surface
            if (vertex->get_surface() == added_surface) continue;

            // reduce using shortest distance to ray
            const double distance = shortest_distance_to_line_segment(generic_point->get_origin(), generic_point->get_position(), vertex->get_position());
            vertex->reduce_reverse_radius_search_radius(distance);
        }

        // bvh search radius reduction
        for (const std::shared_ptr<Face>& face : bvh_results)
        {
            // skip if expired
            if (face->is_expired()) continue;

            // skip if the same surface
            if (face->get_surface() == added_surface) continue;

            // don't reduce radius of face if counter is not zero
            if (face->get_reduce_radius_counter() > 0)
            {
                face->decrement_reduce_radius_counter();
                continue;
            }

            // get copy of vertices and interior points (as they might be deleted)
            std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices = face->get_vertices();
            std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points = face->get_interior_points();

            // reduce using shortest distance to ray
            for (const std::shared_ptr<Vertex>& vertex : vertices)
            {
                // skip if expired
                if (vertex->is_expired()) continue;

                // reduce using shortest distance to ray
                const double distance = shortest_distance_to_line_segment(generic_point->get_origin(), generic_point->get_position(), vertex->get_position());
                vertex->reduce_reverse_radius_search_radius(distance);
            }
            for (const std::shared_ptr<InteriorPoint>& interior_point : interior_points)
            {
                // skip if expired
                if (interior_point->is_expired()) continue;

                // reduce using shortest distance to ray
                const double distance = shortest_distance_to_line_segment(generic_point->get_origin(), generic_point->get_position(), interior_point->get_position());
                interior_point->reduce_reverse_radius_search_radius(distance);
            }
        }
    }

    num_of_concurrent_processes--;
    
    // after unlocking all locks, add the point in queue to the search tree
    storage_->add_points_in_affected_vertices_set();
    storage_->add_faces_in_affected_faces_set();

    //
    // they all lead to return, thus unlock all locks
    //

    // unlock surface
    for (const std::shared_ptr<Surface>& surface : locked_surfaces) omp_unset_nested_lock_with_log(surface->lock, "unlock surface");

    // unlock bvh nodes
    for (const std::shared_ptr<Node>& node : locked_bvh_nodes) node->recursive_unlock();

    // unlock rrs nodes
    for (const std::shared_ptr<RRSNode>& node : locked_rrs_nodes) node->recursive_unlock();

    // lock.unlock();
}

template <typename PointT>
bool Application<PointT>::add_point_by_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, double& new_point_radius, std::vector<std::shared_ptr<Face>>& bvh_results, std::shared_ptr<Surface>& surface_to_add_to)
{
    // don't reduce radius if the point is added as new surface, thus need to move the reduce radius part outside this function.

    // // add point by radius search or add to new surface (when no search results)
    if (bvh_results.size() == 0) return false;

    // all surfaces
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> all_surfaces;
    for (const std::shared_ptr<Face>& face : bvh_results)
    {
        all_surfaces.insert(face->get_surface());
    }

    // check relative position with the searched surfaces
    // update the smallest distance to surface not within
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_seed;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_in_front;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_within;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_behind;
    for (const std::shared_ptr<Surface>& surface : all_surfaces)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // seed surface
        if (surface->get_total_point_size() < settings_.fit_plane_threshold)
        {
            surfaces_seed.insert(surface);
            continue;
        }

        // relative position
        RelativePosition relative_position = surface->check_relative_position(generic_point);
        if (relative_position == RelativePosition::WITHIN)
        {
            surfaces_within.insert(surface);
        }
        else if (relative_position == RelativePosition::BEHIND)
        {
            surfaces_behind.insert(surface);
        }
        else if (relative_position == RelativePosition::IN_FRONT)
        {
            surfaces_in_front.insert(surface);
        }
    }

    // update new point radius
    for (const std::shared_ptr<Surface>& surface : surfaces_in_front)
    {
        // update radius
        double distance = surface->compute_point_to_plane_distance(generic_point->get_position());
        if (distance < new_point_radius) new_point_radius = distance;
    }
    for (const std::shared_ptr<Surface>& surface : surfaces_behind)
    {
        // update radius
        double distance = surface->compute_point_to_plane_distance(generic_point->get_position());
        if (distance < new_point_radius) new_point_radius = distance;
    }
    
    // log
    if (settings_.log.process_point) std::cout << ">> found " << bvh_results.size() << " searched faces grouped into " << all_surfaces.size() << " searched surfaces" << std::endl;

    // decide surface_to_add_to
    if (surfaces_within.size() > 0)
    {
        // add to largest surface
        unsigned int largest_surface_size = 0;
        for (std::shared_ptr<Surface> surface : surfaces_within)
        {
            if (surface->get_total_point_size() > largest_surface_size)
            {
                largest_surface_size = surface->get_total_point_size();
                surface_to_add_to = surface;
            }
        }
    }
    else if (surfaces_seed.size() > 0)
    {
        // add to the closest surface measured by point to point distance
        double smallest_distance = std::numeric_limits<double>::max();

        // for each face
        for (std::shared_ptr<Face> face : bvh_results)
        {
            // skip if expired
            if (face->is_expired()) continue;

            // this surface
            std::shared_ptr<Surface> this_surface = face->get_surface();
            
            // skip if the face is not in surfaces_seed
            if (surfaces_seed.find(this_surface) == surfaces_seed.end()) continue;

            // for each vertex of the face
            for (std::shared_ptr<Vertex> vertex : face->get_vertices())
            {
                // compute distance
                double distance = (vertex->get_position() - generic_point->get_position()).norm();
                if (distance < smallest_distance)
                {
                    smallest_distance = distance;
                    surface_to_add_to = this_surface;
                }
            }
        }
    }
    else
    {
        // add point by radius search or new surface (when no surfaces within nor seed)
        return false;
    }

    // log
    if (settings_.log.process_point) std::cout << "========================== within surface " << surface_to_add_to->get_id() << std::endl;

    // add to surface_to_add_to
    const std::shared_ptr<InteriorPoint>& new_interior_point = storage_->add_interior_point(generic_point);
    new_interior_point->reduce_reverse_radius_search_radius(new_point_radius);
    new_interior_point->reduce_previous_radius(new_point_radius);
    new_interior_point->connect(surface_to_add_to);
    for (const std::shared_ptr<Face>& face : bvh_results)
    {
        // skip if the face is not from the surface
        if (face->get_surface() != surface_to_add_to) continue;

        // skip if the face is expired
        if (face->is_expired()) continue;

        // connect
        new_interior_point->connect(face);

        // increment reduce radius counter
        face->increment_reduce_radius_counter();
        
        break;
    }    
    return true;
}

template <typename PointT>
bool Application<PointT>::add_point_by_radius_search(const std::shared_ptr<GenericPoint>& generic_point, double& new_point_radius, std::vector<std::shared_ptr<Vertex>>& rrs_results, std::shared_ptr<Surface>& surface_to_add_to)
{
    // log
    if (settings_.log.add_point_by_radius_search) std::cout << ">> found " << rrs_results.size() << " neighboring vertices" << std::endl;

    // convert to set
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> all_vertices(rrs_results.begin(), rrs_results.end());

    // // add to new surface (when no search results)
    if (all_vertices.size() == 0) return false;
    
    // all_surfaces
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> all_surfaces; 
    for (std::shared_ptr<Vertex> vertex : all_vertices)
    {
        // only add to neighboring surface if the vertex is boundary in that surface
        if (vertex->is_boundary())
        {
            all_surfaces.insert(vertex->get_surface());
        }
    }
    if (settings_.log.add_point_by_radius_search) std::cout << ">> grouped into " << all_surfaces.size() << " neighboring surfaces" << std::endl;

    // split surfaces into seed, within, and outside
    std::set<std::shared_ptr<Surface>> surfaces_seed;
    std::set<std::shared_ptr<Surface>> surfaces_within;
    std::set<std::shared_ptr<Surface>> surfaces_outside;
    for (std::shared_ptr<Surface> surface : all_surfaces)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // seed surface
        if (surface->get_total_point_size() < settings_.fit_plane_threshold)
        {
            surfaces_seed.insert(surface);
            continue;
        }
        
        // check relative position
        RelativePosition relative_position = surface->check_relative_position(generic_point);
        if (relative_position == RelativePosition::WITHIN)
        { 
            surfaces_within.insert(surface);
        }
        else
        {
            surfaces_outside.insert(surface);
        }
    }

    // update radius
    for (std::shared_ptr<Vertex> vertex : all_vertices)
    {
        // skip if the vertex is expired
        if (vertex->is_expired()) continue;

        // skip if the vertex is not in surfaces_with_point_not_within
        if (surfaces_outside.find(vertex->get_surface()) == surfaces_outside.end()) continue;

        // compute distance
        double distance = (vertex->get_position() - generic_point->get_position()).norm();

        // reduce new_point_radius
        if (distance < new_point_radius) new_point_radius = distance;
    }

    // decide surface_to_add_to
    if (surfaces_within.size() > 0)
    {
        // add to largest surface
        unsigned int largest_surface_size = 0;
        for (std::shared_ptr<Surface> surface : surfaces_within)
        {
            if (surface->get_total_point_size() > largest_surface_size)
            {
                largest_surface_size = surface->get_total_point_size();
                surface_to_add_to = surface;
            }
        }
    }
    else if (surfaces_seed.size() > 0)
    {
        // add to the closest surface measured by point to point distance
        double smallest_distance = std::numeric_limits<double>::max();
        
        // for each vertex
        for (std::shared_ptr<Vertex> vertex : all_vertices)
        {
            // skip if expired
            if (vertex->is_expired()) continue;

            // this surface
            std::shared_ptr<Surface> this_surface = vertex->get_surface();
            
            // skip if the face is not in surfaces_seed
            if (surfaces_seed.find(this_surface) == surfaces_seed.end()) continue;

            // compute distance
            double distance = (vertex->get_position() - generic_point->get_position()).norm();
            if (distance < smallest_distance)
            {
                smallest_distance = distance;
                surface_to_add_to = this_surface;
            }
        }
    }
    else
    {
        // add to new surface (when no surfaces within nor seed)
        return false;
    }

    // add to surface_to_add_to
    std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(surface_to_add_to, generic_point);
    new_vertex->reduce_reverse_radius_search_radius(new_point_radius);
    new_vertex->reduce_previous_radius(new_point_radius);
    surface_to_add_to->connect_by_edges_and_faces(new_vertex, all_vertices);
    return true;
}

template <typename PointT>
void Application<PointT>::add_point_by_new_surface(const std::shared_ptr<GenericPoint>& generic_point, double& radius, std::shared_ptr<Surface>& added_surface)
{
    std::shared_ptr<Surface> new_surface = storage_->add_surface();
    std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(new_surface, generic_point);
    new_vertex->reduce_reverse_radius_search_radius(radius);
    new_vertex->reduce_previous_radius(radius);
    
    added_surface = new_surface;
}

template <typename PointT>
void Application<PointT>::get_lidar_data(Eigen::Vector3d& origin, Eigen::Vector3d& position)
{
    origin = this->origin;
    position = pointcloud->points[ith_point].getVector3fMap().template cast<double>();
    ith_point += settings_.process_every_n_points;

    if (ith_point >= ith_size)
    {
        ith_cloud += 1;
        ith_point = 0;
        load_point_cloud();
    }
}

template <typename PointT>
void Application<PointT>::get_sim_data(Eigen::Vector3d& origin, Eigen::Vector3d& position)
{
    Simulation sim;
    sim.set_object(settings_.sim_object);
    sim.set_noise(settings_.range_precision, settings_.range_accuracy);
    sim.get_data_pair(origin, position);
}

template <typename PointT>
void Application<PointT>::step()
{        
    Eigen::Vector3d thisPointVEC;
    Eigen::Vector3d thisPointOriginVEC;
    if (settings_.use_sim_data)
    {
        get_sim_data(thisPointOriginVEC, thisPointVEC);
    }
    else
    {
        get_lidar_data(thisPointOriginVEC, thisPointVEC);
    }

    // log
    if (settings_.log.step) std::cout << "==================================================================== Processing point " << ith_point << " of cloud " << ith_cloud << std::endl;

    std::shared_ptr<GenericPoint> generic_point = std::make_shared<GenericPoint>();
    generic_point->initialize_(storage_, thisPointVEC, thisPointOriginVEC, distance_travelled_);
    process_point(generic_point);
}

template <typename PointT>
void Application<PointT>::loop()
{
    // add all points from the cloud to the storage's processing queue
    for (std::size_t i = 0; i < ith_size; i++)
    {
        Eigen::Vector3d thisPointVEC = pointcloud->points[i].getVector3fMap().template cast<double>();
        Eigen::Vector3d thisPointOriginVEC = this->origin;

        storage_->add_to_main_queue(thisPointVEC, thisPointOriginVEC, distance_travelled_);
    }

    // add points from repeated queue to queue
    storage_->add_points_in_smaller_repeated_queues_to_main_queue();

    // split the queue into smaller queues
    storage_->split_main_queue_into_smaller_queues();

    // reset stats
    t_init = std::chrono::high_resolution_clock::now();
    accumulated_points = 0;

    // process all points in the queue
    unsigned int num_iteration = 0;
    while (true)
    {
        #pragma omp parallel num_threads(settings_.num_threads)
        {
            while (true)
            {
                std::shared_ptr<GenericPoint> generic_point = nullptr;
                unsigned int total_points_in_queue = 0;

                total_points_in_queue = storage_->get_queue_size();
                if (total_points_in_queue > 0)
                {
                    generic_point = storage_->pop_from_queue();
                    if (settings_.log.step) 
                    {
                        std::stringstream ss;
                        ss << " | remaining point " << total_points_in_queue << " | Processing point " << generic_point->get_id() << " | by thread " << omp_get_thread_num() << std::endl;
                        std::cout << ss.str();
                    }
                }

                // If the queue is empty, break the loop
                if (!generic_point) break;

                // Process the point (outside the critical section)
                process_point(generic_point);
            }
        }
        num_iteration ++;

        // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double> duration = t_end - t_init;
        // std::cout << "==================================================================== Processed " << accumulated_points << " points in " << duration.count() << " s" << std::endl;

        // add all points from the abort queue to the main queue
        storage_->add_points_in_smaller_abort_queues_to_main_queue();

        storage_->print_main_queue_stats();

        // break loop if no more points
        if (storage_->get_main_queue_size() == 0)
        {
            // clear queue before loading next point cloud
            storage_->clear_queues();
            break;
        }

        // print aborted point stats
        storage_->split_main_queue_into_smaller_queues_by_contention();
    }

    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = t_end - t_init;
    std::cout << "==================================================================== Processed " << accumulated_points << " points in " << duration.count() << " s" << std::endl;

    // store time and ith_cloud into file
    if (settings_.output_time)
    {
        std::ofstream file;
        file.open(settings_.output_file_name, std::ios_base::app);
        file << ith_cloud << " " << duration.count() << std::endl;
        file.close();
    }

    // print repeated queue size
    if (settings_.log.step) std::cout << "==================================================================== repeated queue size: " << storage_->get_repeated_queue_size() << std::endl;

    // // check tree rebuild
    // storage_->check_tree_rebuild();

    // load next cloud
    ith_cloud += 1;
    load_point_cloud();
}

template <typename PointT>
void Application<PointT>::process_the_rest()
{
    for (int i = ith_cloud; i < data_loader.size(); i++)
    {
        loop();
    }
}

template <typename PointT>
void Application<PointT>::restart()
{
    // reset objects
    storage_ = std::make_shared<Storage>();
    data_loader.load_dataset(settings_.data_loader_settings);

    // reset state
    ith_cloud = settings_.data_loader_settings.start_cloud;
    ith_point = 0;
    ith_size = 0;
    origin = Eigen::Vector3d(0, 0, 0);
    distance_travelled_ = 0;
    
    // load point cloud
    load_point_cloud();
}

template <typename PointT>
void Application<PointT>::rebuild_tree()
{
    storage_->rebuild_tree();
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
pcl::PointCloud<pcl::PointXYZRGB>::Ptr Application<PointT>::compute_interior_point_pointcloud(const Settings& settings)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    for (const std::shared_ptr<InteriorPoint>& interior_point : storage_->get_interior_points())
    {
        // skip if not confirmed
        if (settings.show_confirmed_only && !interior_point->is_confirmed()) continue;
        
        pcl::PointXYZRGB point;
        if (settings.show_projected_point)
        {
            point.x = interior_point->buffer_compute_projected_position()[0];
            point.y = interior_point->buffer_compute_projected_position()[1];
            point.z = interior_point->buffer_compute_projected_position()[2];
        }
        else
        {
            point.x = interior_point->get_position()[0];
            point.y = interior_point->get_position()[1];
            point.z = interior_point->get_position()[2];
        }
        if (settings.color_mode == ColorMode::ID)
        {
            const std::tuple<int, int, int>& color = interior_point->get_surface()->get_color();
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::POSITIONAL_UNCERTAINTY)
        {
            double distance = std::fabs(interior_point->buffer_compute_projected_distance() / 0.05);
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::POSITIONAL_UNCERTAINTY_NORMALIZED)
        {
            // print
            std::cout << "uncertainty: " << interior_point->get_projected_uncertainty() << std::endl;

            double distance = std::abs(interior_point->buffer_compute_projected_distance()) / interior_point->get_projected_uncertainty() / 3.f;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::PROJECTED_UNCERTAINTY)
        {
            double distance = interior_point->get_projected_uncertainty() / 0.1f;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::RADIUS)
        {
            double ratio = interior_point->get_radius() / settings.radius_denominator;
            std::tuple<int, int, int> color = valueToJet(1.0-ratio);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::SURFACE_UNCERTAINTY)
        {
            double distance = interior_point->get_surface()->get_surface_position_std_in_normal_direction() / settings.positional_uncertainty_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::CONTENTION)
        {
            double distance = surface_to_contention_count[interior_point->get_surface()] / settings.contention_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::DISTANCE_TRAVELLED)
        {
            double distance = interior_point->get_distance_travelled() / 20.f;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == ColorMode::WEIGHT)
        {
            double distance = interior_point->weight_ / 10000.f;
            std::tuple<int, int, int> color = valueToJet(distance);
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
        surface->remove_unmatched_points();
        surface->remove_singular_components();
        surface->split_surface_by_connected_components();
    }
    if (settings_.log.refine_surfaces) std::cout << "number of generic points after refine: " << storage_->get_generic_points().size() << std::endl;
}

template <typename PointT>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr Application<PointT>::compute_vertex_point_pointcloud(const Settings& setting)
{
    vertex_to_cloud_indices_map.clear();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    for (const std::shared_ptr<Vertex>& vertex : storage_->get_vertices())
    {
        // skip if not confirmed
        if (setting.show_confirmed_only && !vertex->is_confirmed()) continue;

        // skip if singular
        if (!setting.show_singular_vertex && vertex->is_singular()) continue;

        pcl::PointXYZRGB point;
        if (setting.show_projected_point)
        {
            point.x = vertex->buffer_compute_projected_position()[0];
            point.y = vertex->buffer_compute_projected_position()[1];
            point.z = vertex->buffer_compute_projected_position()[2];
        }
        else
        {
            point.x = vertex->get_position()[0];
            point.y = vertex->get_position()[1];
            point.z = vertex->get_position()[2];
        }
        if (setting.color_mode == ColorMode::ID)
        {
            const std::tuple<int, int, int>& color = vertex->get_surface()->get_color();
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::POSITIONAL_UNCERTAINTY)
        {
            double distance = std::abs(vertex->buffer_compute_projected_distance() / 0.05);
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::POSITIONAL_UNCERTAINTY_NORMALIZED)
        {
            double distance = std::abs(vertex->buffer_compute_projected_distance()) / vertex->get_projected_uncertainty() / 3.f;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::PROJECTED_UNCERTAINTY)
        {
            double distance = vertex->get_projected_uncertainty() / 0.1f;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::RADIUS)
        {
            double ratio = vertex->get_radius() / setting.radius_denominator;
            std::tuple<int, int, int> color = valueToJet(1.0-ratio);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::SURFACE_UNCERTAINTY)
        {
            double distance = vertex->get_surface()->get_surface_position_std_in_normal_direction() / setting.positional_uncertainty_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::CONTENTION)
        {
            double distance = surface_to_contention_count[vertex->get_surface()] / setting.contention_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::DISTANCE_TRAVELLED)
        {
            double distance = vertex->get_distance_travelled() / 20.f;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == ColorMode::WEIGHT)
        {
            double distance = vertex->weight_ / 10000.f;
            std::tuple<int, int, int> color = valueToJet(distance);
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