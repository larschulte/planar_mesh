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
    add_point_to_map(generic_point, bvh_results, rrs_results);

    num_of_concurrent_processes--;
    
    // after unlocking all locks, add the point in queue to the search tree
    storage_->add_or_remove_vertices_from_rrs_tree();
    storage_->add_or_remove_faces_from_bvh_tree();

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
void Application<PointT>::add_point_to_map(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& bvh_results, std::vector<std::shared_ptr<Vertex>>& rrs_results)
{
    // from bvh results and rrs results, get surfaces
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> bvh_surfaces;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> rrs_surfaces;
    for (const std::shared_ptr<Face>& face : bvh_results)
    {
        // skip if expired
        if (face->is_expired()) continue;

        bvh_surfaces.insert(face->get_surface());
    }
    for (const std::shared_ptr<Vertex>& vertex : rrs_results)
    {
        // skip if expired
        if (vertex->is_expired()) continue;
        
        rrs_surfaces.insert(vertex->get_surface());
    }

    // split surfaces into bvh and rrs seed, within, behind, in front
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_bvh_seed;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_bvh_within;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_bvh_behind;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_bvh_in_front;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_rrs_seed;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_rrs_within;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_rrs_behind;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_rrs_in_front;
    for (const std::shared_ptr<Surface>& surface : bvh_surfaces)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // check relative position
        RelativePosition relative_position = surface->check_relative_position(generic_point);

        // skip if no relative position
        if (relative_position == RelativePosition::NO_RELATIVE_POSITION) surfaces_bvh_seed.insert(surface);

        // add to within
        if (relative_position == RelativePosition::WITHIN) surfaces_bvh_within.insert(surface);

        // add to behind
        if (relative_position == RelativePosition::BEHIND) surfaces_bvh_behind.insert(surface);

        // add to in front
        if (relative_position == RelativePosition::IN_FRONT) surfaces_bvh_in_front.insert(surface);
    }
    for (const std::shared_ptr<Surface>& surface : rrs_surfaces)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // check relative position
        RelativePosition relative_position = surface->check_relative_position(generic_point);

        // skip if no relative position
        if (relative_position == RelativePosition::NO_RELATIVE_POSITION) surfaces_rrs_seed.insert(surface);

        // add to within
        if (relative_position == RelativePosition::WITHIN) surfaces_rrs_within.insert(surface);

        // add to behind
        if (relative_position == RelativePosition::BEHIND) surfaces_rrs_behind.insert(surface);

        // add to in front
        if (relative_position == RelativePosition::IN_FRONT) surfaces_rrs_in_front.insert(surface);
    }

    // construct a list of surfaces to add to
    std::vector<std::shared_ptr<Surface>> list_of_surfaces_to_add_to;

    // 1. surfaces_bvh_within
    if (surfaces_bvh_within.size() > 0)
    {
        // sort by surface size
        std::vector<std::shared_ptr<Surface>> surfaces_bvh_within_sorted(surfaces_bvh_within.begin(), surfaces_bvh_within.end());
        std::sort(surfaces_bvh_within_sorted.begin(), surfaces_bvh_within_sorted.end(), 
            [](const std::shared_ptr<Surface>& a, const std::shared_ptr<Surface>& b) { return a->get_total_point_size() > b->get_total_point_size(); });

        // add to list
        list_of_surfaces_to_add_to.insert(list_of_surfaces_to_add_to.end(), surfaces_bvh_within_sorted.begin(), surfaces_bvh_within_sorted.end());
    }

    // 2. surfaces_rrs_within
    if (surfaces_rrs_within.size() > 0)
    {
        // sort by surface size
        std::vector<std::shared_ptr<Surface>> surfaces_rrs_within_sorted(surfaces_rrs_within.begin(), surfaces_rrs_within.end());
        std::sort(surfaces_rrs_within_sorted.begin(), surfaces_rrs_within_sorted.end(), 
            [](const std::shared_ptr<Surface>& a, const std::shared_ptr<Surface>& b) { return a->get_total_point_size() > b->get_total_point_size(); });

        // add to list
        list_of_surfaces_to_add_to.insert(list_of_surfaces_to_add_to.end(), surfaces_rrs_within_sorted.begin(), surfaces_rrs_within_sorted.end());
    }

    // 3. surfaces_rrs_seed
    if (surfaces_rrs_seed.size() > 0)
    {        
        // sort by surface closeness
        std::unordered_map<std::shared_ptr<Surface>, double, MeshObjectHash> surfaces_rrs_seed_distances;
        for (std::shared_ptr<Surface> surface : surfaces_rrs_seed)
        {
            for (std::shared_ptr<Vertex> vertex : rrs_results)
            {
                // skip if the vertex is expired
                if (vertex->is_expired()) continue;

                // skip if the vertex is not in surfaces_rrs_seed
                if (vertex->get_surface() != surface) continue;

                // compute distance
                double distance_vertex = (vertex->get_position() - generic_point->get_position()).norm();
                if (distance_vertex < surfaces_rrs_seed_distances[surface])
                {
                    surfaces_rrs_seed_distances[surface] = distance_vertex;
                }
            }
        }
        std::vector<std::shared_ptr<Surface>> surfaces_rrs_seed_sorted(surfaces_rrs_seed.begin(), surfaces_rrs_seed.end());
        std::sort(surfaces_rrs_seed_sorted.begin(), surfaces_rrs_seed_sorted.end(), 
            [&surfaces_rrs_seed_distances](const std::shared_ptr<Surface>& a, const std::shared_ptr<Surface>& b) { return surfaces_rrs_seed_distances[a] < surfaces_rrs_seed_distances[b]; });

        // add to list
        list_of_surfaces_to_add_to.insert(list_of_surfaces_to_add_to.end(), surfaces_rrs_seed_sorted.begin(), surfaces_rrs_seed_sorted.end());
    }

    // 4. new surface
    std::shared_ptr<Surface> new_surface = storage_->add_surface();
    list_of_surfaces_to_add_to.push_back(new_surface);

    // add to surface according to the list
    std::shared_ptr<Surface> surface_to_add_to;
    std::shared_ptr<InteriorPoint> new_interior_point = nullptr;
    std::shared_ptr<Vertex> new_vertex = nullptr;
    for (std::shared_ptr<Surface> surface : list_of_surfaces_to_add_to)
    {
        // store current surface to add to
        surface_to_add_to = surface;

        // if added as interior point
        if (bvh_surfaces.find(surface_to_add_to) != bvh_surfaces.end())
        {
            // add as interior point
            new_interior_point = storage_->add_interior_point(generic_point);

            // get the first face
            std::shared_ptr<Face> face_to_add_to;
            for (const std::shared_ptr<Face>& face : bvh_results)
            {
                // skip if the face is expired
                if (face->is_expired()) continue;

                // skip if the face is not from the surface
                if (face->get_surface() != surface_to_add_to) continue;

                // set face_to_add_to
                face_to_add_to = face;
                break;
            }

            // add
            new_interior_point->connect(surface_to_add_to);
            new_interior_point->connect(face_to_add_to);

            // delete new surface
            storage_->delete_surface(new_surface);

            // return
            break;
        }
        
        // if added as vertex
        new_vertex = storage_->add_vertex(surface_to_add_to, generic_point);

        // reduce radius of the new vertex
        for (std::shared_ptr<Vertex> vertex : rrs_results)
        {
            // skip if the vertex is expired
            if (vertex->is_expired()) continue;

            // skip if the vertex is the same as the new vertex
            if (vertex->get_surface() == surface_to_add_to) continue;

            // add neighboring vertex
            new_vertex->add_vertex_point_distance_publisher(vertex);
        }

        // try update
        new_vertex->try_update_radius();

        // if new surface
        if (surface_to_add_to == new_surface)
        {
            break;
        }

        // if old surface
        std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> rrs_results_set = std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>(rrs_results.begin(), rrs_results.end());
        const bool connected = surface_to_add_to->connect_by_edges_and_faces(new_vertex, rrs_results_set);
        if (connected)
        {
            // delete new surface
            storage_->delete_surface(new_surface);
            break;
        }
        else
        {
            // need to prevent this vertex from generating a new generic point
            new_vertex->can_create_generic_point(false);
            storage_->delete_vertex(new_vertex);
            continue;
        }
    }

    // delete within and penetrated that are from differnt surface
    for (const std::shared_ptr<Face>& face : bvh_results)
    {
        // skip if expired
        if (face->is_expired()) continue;

        // skip if same surface
        if (face->get_surface() == surface_to_add_to) continue;

        // skip if no relative position
        if (surfaces_bvh_seed.find(face->get_surface()) != surfaces_bvh_seed.end()) continue;

        // skip if in front 
        if (surfaces_bvh_in_front.find(face->get_surface()) != surfaces_bvh_in_front.end()) continue;

        // delete penetrated face
        storage_->delete_face(face);
    }
    
    // reduce radius of nearby vertices as ray
    for (const std::shared_ptr<Vertex>& vertex : rrs_results)
    {
        // skip if expired
        if (vertex->is_expired()) continue;

        // skip if same surface
        if (vertex->get_surface() == surface_to_add_to) continue;

        if (new_interior_point)
        {
            // skip if no relative position
            if (surfaces_rrs_seed.find(vertex->get_surface()) != surfaces_rrs_seed.end()) continue;

            // skip if in front
            if (surfaces_rrs_in_front.find(vertex->get_surface()) != surfaces_rrs_in_front.end()) continue;

            // add as subscriber
            new_interior_point->add_interior_ray_distance_subscriber(vertex);
        }

        if (new_vertex)
        {
            // add neighboring vertex
            new_vertex->add_vertex_point_distance_subscriber(vertex);
        }
    }    
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
                if (generic_point)
                {
                    process_point(generic_point);
                }
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
    total_duration += duration.count();
    total_loops += 1;
    std::cout << "==================================================================== Processed " << accumulated_points << " points in " << duration.count() << " s, " << "average duration: " << total_duration / total_loops << std::endl;

    // print size of rrs, vertices, and boundary vertices
    std::cout << "rrs size: " << storage_->get_rrs_size() << " | b-vertices size: " << storage_->get_boundary_vertices().size() << " | vertices size: " << storage_->get_vertices().size() << " | total point size: " << storage_->get_vertices().size() + storage_->get_interior_points().size() << std::endl;
    // print size of bvh, faces
    std::cout << "bvh size: " << storage_->get_bvh_size() << " | faces size: " << storage_->get_faces().size() << std::endl;

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

    total_duration = 0;
    total_loops = 0;
    
    // load point cloud
    load_point_cloud();
}

template <typename PointT>
void Application<PointT>::rebuild_tree()
{
    storage_->rebuild_tree();
}

template <typename PointT>
std::shared_ptr<Storage> Application<PointT>::get_storage()
{
    return storage_;
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
        if (settings.point_mode == PointMode::PROJECTED)
        {
            point.x = interior_point->buffer_compute_projected_position()[0];
            point.y = interior_point->buffer_compute_projected_position()[1];
            point.z = interior_point->buffer_compute_projected_position()[2];
        }
        else if (settings.point_mode == PointMode::ORIGINAL)
        {
            point.x = interior_point->get_original_position()[0];
            point.y = interior_point->get_original_position()[1];
            point.z = interior_point->get_original_position()[2];
        }
        else if (settings.point_mode == PointMode::USED)
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
            // the lowest uncertainty is 0.02 which is the range precision
            // 1/(0.02*0.02) = 2500
            // the best weight is 2500
            // all weight should be lower than this number
            const double best_weight = 1.f / (settings.range_precision * settings.range_precision);
            double distance = interior_point->weight_ / best_weight; 
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
        // skip if internal
        if (!setting.show_internal_vertices && !vertex->is_boundary()) continue;

        // skip if not confirmed
        if (setting.show_confirmed_only && !vertex->is_confirmed()) continue;

        // skip if singular
        if (!setting.show_singular_vertex && vertex->is_singular()) continue;

        pcl::PointXYZRGB point;
        if (setting.point_mode == PointMode::PROJECTED)
        {
            point.x = vertex->buffer_compute_projected_position()[0];
            point.y = vertex->buffer_compute_projected_position()[1];
            point.z = vertex->buffer_compute_projected_position()[2];
        }
        else if (setting.point_mode == PointMode::ORIGINAL)
        {
            point.x = vertex->get_original_position()[0];
            point.y = vertex->get_original_position()[1];
            point.z = vertex->get_original_position()[2];
        }
        else if (setting.point_mode == PointMode::USED)
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
            const double best_weight = 1.f / (settings_.range_precision * settings_.range_precision);
            double distance = vertex->weight_ / best_weight;
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