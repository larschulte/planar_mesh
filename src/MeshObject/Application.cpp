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

template class Application<VilensPointT>;

template <typename PointT>
Application<PointT>::Application() 
{
    ith_cloud = settings_.start_cloud;
    ith_point = settings_.start_point;

    storage_ = std::make_shared<Storage>();
    data_loader.load_dataset(settings_.cloud_path, settings_.pose_path);
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
    
    typename pcl::PointCloud<PointT>::Ptr pointcloud_local = data_loader.get_cloud(ith_cloud);
    Eigen::Affine3d pose = data_loader.get_pose(ith_cloud);
    pointcloud = transform_cloud_to_global<PointT>(pointcloud_local, pose);
    origin = pose.translation();
    ith_size = pointcloud->size() * settings_.pointcloud_fraction;

    if (settings_.log.load_point_cloud) std::cout << "loaded pointcloud " << ith_cloud << " with " << pointcloud->size() << " points" << std::endl;

    if (settings_.shuffle_pointcloud) 
    {
        std::random_shuffle(pointcloud->points.begin(), pointcloud->points.end());
    }
}

template <typename PointT>
void Application<PointT>::process_point(const std::shared_ptr<GenericPoint>& generic_point)
{
    double radius = settings_.radius_value;
    if (add_point_by_intersection_search(generic_point, radius)) 
    {
        return;
    }
    else if (add_point_by_radius_search(generic_point, radius))
    {
        return;
    }
    else
    {
        add_point_by_new_surface(generic_point, radius);
        return;
    }
}

template <typename PointT>
bool Application<PointT>::add_point_by_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, double& radius)
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

    // skip if too close to any vertices of any faces
    for (const std::shared_ptr<Face>& face : searched_faces)
    {
        for (const std::shared_ptr<Vertex>& vertex : face->get_vertices())
        {
            double distance = (vertex->get_position() - generic_point->get_position()).norm();
            if (distance < settings_.duplicated_point_distance_threshold)
            {
                if (settings_.log.duplicated_point) std::cout << ">> point too close to existing vertex of face, not adding" << std::endl;
                return true;
            }
        }
    }

    // searched surfaces
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> searched_surfaces;
    for (const std::shared_ptr<Face>& face : searched_faces)
    {
        searched_surfaces.insert(face->get_surface());
    }

    // check relative position with the searched surfaces
    // update the smallest distance to surface not within
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_with_point_within;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_with_point_behind;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_with_point_infront;
    for (const std::shared_ptr<Surface>& surface : searched_surfaces)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // if within surface
        if (surface->check_relative_position(generic_point) == RelativePosition::WITHIN)
        {
            surfaces_with_point_within.insert(surface);
        }
        // if behind surface
        else if (surface->check_relative_position(generic_point) == RelativePosition::BEHIND)
        {
            surfaces_with_point_behind.insert(surface);
            double projective_distance = surface->compute_point_projective_distance(generic_point);
            if (projective_distance < radius)
            {
                radius = projective_distance;
            }
        }
        else
        {
            surfaces_with_point_infront.insert(surface);
            double projective_distance = surface->compute_point_projective_distance(generic_point);
            if (projective_distance < radius)
            {
                radius = projective_distance;
            }
        }
    }

    // // delete abnormal surfaces
    // bool delete_abnormal_surfaces = false;
    // for (const auto& pair : searched_surface_to_searched_faces)
    // {
    //     const std::shared_ptr<Surface>& surface = pair.first;

    //     if (surface->is_abnormal())
    //     {
    //         // delete
    //         std::cout << ">> removing abnormal surface during intersection search" << surface->get_id() << std::endl; // log
    //         storage_->delete_surface(surface);
    //         delete_abnormal_surfaces = true;
    //     }
    // }

    // if (delete_abnormal_surfaces)
    // {
    //     // recompute searched faces
    //     std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> searched_faces_copy = searched_faces;
    //     searched_faces.clear();
    //     for (const std::shared_ptr<Face>& face : searched_faces_copy)
    //     {
    //         // skip if the face is expired
    //         if (face->is_expired()) continue;

    //         searched_faces.insert(face);
    //     }

    //     // recompute searched surfaces to searched faces map
    //     searched_surface_to_searched_faces.clear();
    //     for (const std::shared_ptr<Face>& face : searched_faces)
    //     {
    //         searched_surface_to_searched_faces[face->get_surface()].insert(face);   
    //     }
    // }
    
    // log
    if (settings_.log.process_point) std::cout << ">> found " << searched_faces.size() << " searched faces grouped into " << searched_surfaces.size() << " searched surfaces" << std::endl;

    // process point behind surface
    for (const std::shared_ptr<Surface>& surface : surfaces_with_point_behind)
    {
        // log
        if (settings_.log.process_point) std::cout << "========================== behind surface " << surface->get_id() << std::endl;

        // delete penetrated faces
        for (const std::shared_ptr<Face>& face : searched_faces)
        {
            // skip if face is not from the surface
            if (face->get_surface() != surface) continue;

            // skip if face is expired from previous delete face operation
            if (face->is_expired()) continue;

            // log
            if (settings_.log.process_point) std::cout << ">> disconnect penetrated face " << face->get_id() << " from surface " << surface->get_id() << std::endl;
            
            // update radius of points in the face
            face->update_radius(generic_point);

            // delete face
            storage_->delete_face(face);
        }
    }

    // process points within surface
    // find the surface with the largest size
    bool added_as_interior_point = false;
    std::shared_ptr<Surface> largest_surface;
    std::size_t largest_surface_size = 0;
    for (const std::shared_ptr<Surface>& surface : surfaces_with_point_within)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // update largest surface
        std::size_t surface_size = surface->get_total_point_size();
        if (surface_size > largest_surface_size)
        {
            largest_surface = surface;
            largest_surface_size = surface_size;
        }
    }

    // add to the largest surface
    if (largest_surface != nullptr)
    {
        // log
        if (settings_.log.process_point) std::cout << "========================== within surface " << largest_surface->get_id() << std::endl;

        // add as interior point
        const std::shared_ptr<InteriorPoint>& temp_interior_point = storage_->add_interior_point(generic_point);
        temp_interior_point->connect(largest_surface);
        temp_interior_point->reduce_reverse_radius_search_radius(radius);

        // get faces with points within
        std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_with_points_within;
        for (const std::shared_ptr<Face>& face : searched_faces)
        {
            // skip if the face is not from the surface
            if (face->get_surface() != largest_surface) continue;

            // skip if the face is expired
            if (face->is_expired()) continue;

            // insert
            faces_with_points_within.insert(face);
        }

        // connect to the first face
        for (const std::shared_ptr<Face>& face : faces_with_points_within)
        {
            // log
            if (settings_.log.process_point) std::cout << ">> adding interior point to face " << face->get_id() << std::endl;

            // connect
            temp_interior_point->connect(face);

            // update flag
            added_as_interior_point = true;

            // break so only added to the first surface
            break;
        }
    }

    // if can't be added as interior point, add as vertex
    if (!added_as_interior_point)
    {
        return false;
    }
    else
    {
        return true;
    }
}

template <typename PointT>
bool Application<PointT>::add_point_by_radius_search(const std::shared_ptr<GenericPoint>& generic_point, double& radius)
{
    /*

    contribution
    - without
        - voxelation
        - outlier
        - denoising
        - radius
        - neighborhood information
    - generalize 

    background:
        - the reverse radius search provides a set of vertices that think the new point should be in
        - the new point to surface fit is then performed -> represented by point to surface projective distance, point uncertainty and plane uncertainty
        - by adding a point to a surface, we are essentially expanding the surface's boundary

    settings:
        - each vertex can have multiple surfaces, but each edge and face and interiror point can only have one surface

    steps:
        - for surface with low confidence, add new point to surface
        - for surface with high confidence and match, add new point to surface (if edge intersects, reduce radius of the other edge vertex)
        - for surface with high confidence and mismatch, reduce the search radius of the neighboring vertices from the surface

    merge neighboring surface
        - after adding the new point to surfaces, if the new point have multiple surfaces, try merge them

    average projected distance
        - for surface with only three points, the average projective distance is zero

    the surface the new point connected to may be duplicate

    what does duplicate means
        - if a vertex/edge/face/point is connected to more than one surface, it is duplicated

    measurement points from planear structure don't have duplicate
        - since edge and vertex of real surface have zero area, they can not be measured

    for curved surface, duplicate may exist


    a diff method between surfaces

    new surface are created when a new point found no existing surface to add to


    remove duplicate surface
        - if a low confidence surface is also penetrated

    
    computing projective distance std when the surface is accurate, would in fact give the std of lidar range sensor
        - perhaps use a metric to determine if the distribution of std is gaussian like, which can be used check if the surface is single modal
        - graph the distribution of projective distance! perhaps we can see multiple peaks, which can be used to classify the surface

    */

   // 1. when we review a point, we should treat it as a new point, then decide which surface to keep / remove
    // this is to ensure indepdenence of point adding order - if a reviewed point is disconnected from some surface, a new duplicate point 
    // at the same location should not connect to the disconnected surfaces

    // 2. confidence flag
    // indicate whether the normal of a surface is reliable
    
    // 3. we add a point to the largest confidence surface
    // (largest number of addtioanl points that agree with the surface)

    // 4. if no confidence surface, we add it to all non confidence surface nearby 
    // (this is the same as grouping a set of nearby points and try to fit a plane to them)

    // 5. if a point is not added to a confidence surface, it will start a new seed

    // 4 and 5 combined gives us a slightly more comprehensive method of grouping a set of nearby points
    // (the best should be N choose K the set of nearby points)

    // 6. abnormal test
    // the projective distance of a sampled point to the sampled surface should follow the distribution of lidar sensor noise (fact)
    // thus if the projective distance of a set of points to a surface have a distribution unlike the lidar sensor noise, unlike means
        // gaussian but shifted origin
        // gaussian but shifted origin and changed std
        // changed std
        // non gaussian
    // that means the point to plane projective distance variable is due to more than just sensor noise, possible source includes
        // shifted surface position
        // shifted surface position and changed surface normal
        // changed surface normal
        // surface is not planar
    // (non confidence surface should not attend the abnormal test)

    // multi surface idea
        // when at large n point threhsold threshold, it would be difficult for n points to be on the same surface, 
        // which is the only way the proposed surface won't be treated as abnormal

        // the proposed adding new point to multiple surface does not work in this case, 
        // because each point only start a new seed from that point on
        // which mean the next n consecutive points added to that point need to be on the same surface
        // for the final surface to be treated as normal

        // the adding to multiple surface idea is only possible if we solve the combination N choose P problem 

    // only update surface close to new observation, such that the update time will be consistant
    // otherwise the update time will be proportional to number of surface added
    

    // when can not search
    if (!storage_->can_reverse_radius_search())
    {
        if (settings_.log.add_point_by_radius_search) std::cout << ">> no points to search, adding new surface" << std::endl;
        return false;
    }

    // get neighboring vertices
    std::map<int, double> point_to_radius_map;
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices = storage_->reverse_radius_search(generic_point);
    if (settings_.log.add_point_by_radius_search) std::cout << ">> found " << neighboring_vertices.size() << " neighboring vertices" << std::endl;

    // when no search results
    if (neighboring_vertices.size() == 0)
    {
        if (settings_.log.add_point_by_radius_search) std::cout << ">> no search results, adding new surface" << std::endl;
        return false;
    }

    // if there is search result, check if the current point is too close to existing point
    // if too close, then the new point will not add any information to the surface
    // thus we should not add the point
    for (const std::shared_ptr<Vertex>& vertex : neighboring_vertices)
    {
        // skip if the vertex is expired
        if (vertex->is_expired()) continue;

        // compute distance
        double distance = (vertex->get_position() - generic_point->get_position()).norm();

        // update closest distance
        if (distance < settings_.duplicated_point_distance_threshold)
        {
            // don't add the point
            if (settings_.log.duplicated_point) std::cout << ">> point too close to existing point, not adding" << std::endl;
            return true;
        }
    }

    // add point logic
        // add to matched surfaces
        // if no matched surfaces, add to all low confidence surfaces
    
    // review point logic
        // if low confidence surface changed into high confidence surface, check if match, if mismatch then disconnect
        // if have matched surfaces, remove from low confidence surfaces
        // remain in matched surfaces

    // bool review_point = false;
    // if (review_point)
    // {
    //     // review point
    //     // make a copy of the neighboring vertices as we will be modifying it
    //     std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices_copy = neighboring_vertices;
    //     for (const std::shared_ptr<Vertex>& vertex : neighboring_vertices_copy)
    //     {
    //         // some sibling vertex may be expired during previous review 
    //         if (vertex->is_expired()) continue;

    //         vertex->review_surfaces();
    //     }

    //     // recompute neighboring vertices
    //     neighboring_vertices.clear();
    //     for (std::shared_ptr<Vertex> vertex : neighboring_vertices_copy)
    //     {
    //         // skip if the vertex is expired
    //         if (vertex->is_expired()) continue;

    //         neighboring_vertices.insert(vertex);
    //     }
    // }

    // get neighboring surfaces
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> neighboring_surfaces; 
    for (std::shared_ptr<Vertex> vertex : neighboring_vertices)
    {
        // only add to neighboring surface if the vertex is boundary in that surface
        if (vertex->is_boundary())
        {
            neighboring_surfaces.insert(vertex->get_surface());
        }
    }
    if (settings_.log.add_point_by_radius_search) std::cout << ">> grouped into " << neighboring_surfaces.size() << " neighboring surfaces" << std::endl;

    // // delete abnormal surfaces
    // bool delete_abnormal_surfaces = false;
    // for (const std::shared_ptr<Surface>& surface : neighboring_surfaces)
    // {
    //     if (surface->is_abnormal())
    //     {
    //         // delete
    //         if (settings_.log.add_point_by_radius_search) std::cout << ">> removing abnormal surface during intersection search" << surface->get_id() << std::endl; // log
    //         storage_->delete_surface(surface);
    //         delete_abnormal_surfaces = true;
    //     }
    // }

    // if (delete_abnormal_surfaces)
    // {
    //     // recompute neighboring vertices
    //     std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices_copy = neighboring_vertices;
    //     neighboring_vertices.clear();
    //     for (std::shared_ptr<Vertex> vertex : neighboring_vertices_copy)
    //     {
    //         // skip if the vertex is expired
    //         if (vertex->is_expired()) continue;

    //         neighboring_vertices.insert(vertex);
    //     }
        
    //     // recompute neighboring surfaces
    //     neighboring_surfaces.clear();
    //     for (std::shared_ptr<Vertex> vertex : neighboring_vertices)
    //     {
    //         // only add to neighboring surface if the vertex is boundary in that surface
    //         if (vertex->is_boundary())
    //         {
    //             neighboring_surfaces.insert(vertex->get_surface());
    //         }
    //     }
    //     if (settings_.log.add_point_by_radius_search) std::cout << ">> grouped into " << neighboring_surfaces.size() << " neighboring surfaces" << std::endl;
    // }

    // remove unmatched points in the surfaces
    // bool removed_unmatched_points = false;
    // for (const std::shared_ptr<Surface>& surface : neighboring_surfaces)
    // {
    //     // the remove unmatched point creates generic points that are not added back
    //     if (surface->remove_unmatched_points())
    //     {
    //         removed_unmatched_points = true;
    //         // remove singular components
    //         bool remove_singular_components = true;
    //         if (remove_singular_components)
    //         {
    //             // remove singular components
    //             surface->remove_singular_components();
    //         }
    //         surface->split_surface_by_connected_components();
    //     }
    // }

    // if (removed_unmatched_points)
    // {
    //     // recompute neighboring vertices
    //     std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices_copy = neighboring_vertices;
    //     neighboring_vertices.clear();
    //     for (std::shared_ptr<Vertex> vertex : neighboring_vertices_copy)
    //     {
    //         // skip if the vertex is expired
    //         if (vertex->is_expired()) continue;

    //         neighboring_vertices.insert(vertex);
    //     }
        
    //     // recompute neighboring surfaces
    //     neighboring_surfaces.clear();
    //     for (std::shared_ptr<Vertex> vertex : neighboring_vertices)
    //     {
    //         // only add to neighboring surface if the vertex is boundary in that surface
    //         if (vertex->is_boundary())
    //         {
    //             neighboring_surfaces.insert(vertex->get_surface());
    //         }
    //     }
    //     if (settings_.log.add_point_by_radius_search) std::cout << ">> grouped into " << neighboring_surfaces.size() << " neighboring surfaces" << std::endl;
    // }

    // for each surface, check if confidence surface, then check if new point is within
    std::set<std::shared_ptr<Surface>> surfaces_with_low_confidence;
    std::set<std::shared_ptr<Surface>> surfaces_with_point_within;
    std::set<std::shared_ptr<Surface>> surfaces_with_point_not_within;
    for (std::shared_ptr<Surface> surface : neighboring_surfaces)
    {
        // skip if the surface is expired
        if (surface->is_expired()) continue;

        // if low confidence surface
        bool low_confidence = surface->get_total_point_size() < settings_.fit_plane_threshold;
        if (low_confidence)
        {
            surfaces_with_low_confidence.insert(surface);
        }
        // else high confidence surface
        else
        {
            // if within high confidence surface
            bool within = surface->check_relative_position(generic_point) == RelativePosition::WITHIN;
            if (within)
            { 
                surfaces_with_point_within.insert(surface);
            }
            else
            {
                surfaces_with_point_not_within.insert(surface);
            }
        }
    }

    // for points belong to surface_with_point_not_within, update reverse search radius
    for (std::shared_ptr<Vertex> vertex : neighboring_vertices)
    {
        // skip if the vertex is expired
        if (vertex->is_expired()) continue;

        // skip if the vertex is not in surfaces_with_point_not_within
        if (surfaces_with_point_not_within.find(vertex->get_surface()) == surfaces_with_point_not_within.end()) continue;

        // compute distance
        double distance = (vertex->get_position() - generic_point->get_position()).norm();

        // reduce reverse search radius
        vertex->reduce_reverse_radius_search_radius(distance);

        // update smallest distance
        if (distance < radius)
        {
            radius = distance;
        }
    }
    
    // reset search radius if surface changes
    // the search radius is a function of, 
        // 1. the point location
        // 2. the surface normal we are basing of
        // 3. the point that do not belong to the surface we are basing of
    // thus, the search radius should reset if any of the above is changed
        // the point location is unchanged
        // the surface normal we are basing of could change if that point location is no longer considered in the surface
        // the point that are not in the surface we are basing of is irrelevant

    // the reverse search radius should be a function of location and surface at that location
    // if the deleted point is added to a different surface than the original one, reset the search radius of the vertex
    // however, as the new point is also not within some surfaces regardless of the new surface added, the new radius should reduce as well


    // if there is surface with points within, add to these surfaces then do a review
    // if there is no surface with points within, add to all low confidence surfaces
    // if there is no surface with points within and no low confidence surfaces, add new surface
    // [todo] if a point is within a surface, but can't connect to it dues to intersecting edge, ignore it for now
    if (surfaces_with_point_within.size() > 0)
    {
        // sort by number of points ([todo] replace by better metric later)
        std::vector<std::shared_ptr<Surface>> sorted_surfaces_with_point_within;
        sorted_surfaces_with_point_within.insert(sorted_surfaces_with_point_within.end(), surfaces_with_point_within.begin(), surfaces_with_point_within.end());
        std::sort(sorted_surfaces_with_point_within.begin(), sorted_surfaces_with_point_within.end(), 
            [](const std::shared_ptr<Surface>& a, const std::shared_ptr<Surface>& b) -> bool
            {
                return a->get_total_point_size() > b->get_total_point_size();
            });
        
        // add to the smallest uncertainty surface as current surface
        std::shared_ptr<Surface> current_surface = sorted_surfaces_with_point_within[0];
        std::shared_ptr<Vertex> current_vertex = storage_->add_vertex(generic_point);

        current_vertex->reduce_reverse_radius_search_radius(radius);
        current_vertex->reduce_previous_radius(radius);
        current_surface->connect_by_edges_and_faces(current_vertex, neighboring_vertices);
        current_surface->connect(current_vertex);
        
        // // connect and merge to next surface if possible
        // for (std::size_t i = 1; i < sorted_surfaces_with_point_within.size(); i++)
        // {
        //     std::shared_ptr<Surface> next_surface = sorted_surfaces_with_point_within[i];
            
        //     // check if mergable
        //     // merge happened between overlapped surfaces causes error!!!
        //     bool can_merge = current_surface->can_merge(next_surface);
        //     if (can_merge)
        //     {
        //         // add to surface
        //         std::shared_ptr<Vertex> next_vertex = storage_->add_vertex(generic_point);
        //         next_vertex->reduce_reverse_radius_search_radius(radius);
        //         next_surface->connect_by_edges_and_faces(next_vertex, neighboring_vertices);
        //         next_surface->connect(next_vertex);
        //         // merge surfaces
        //         current_vertex->absorbs(next_vertex);
        //         DeletedPointStorage original_name = storage_->get_deleted_points_storage_name();
        //         storage_->set_deleted_points_storage_name(DeletedPointStorage::NONE);
        //         storage_->delete_vertex(next_vertex);
        //         storage_->set_deleted_points_storage_name(original_name);
        //     }
        // }
    }
    else if (surfaces_with_low_confidence.size() > 0)
    {
        // add to the smallest surface with low confidence

        // sort by number of points ([todo] replace by better metric later)
        std::vector<std::shared_ptr<Surface>> sorted_surfaces_with_low_confidence;
        sorted_surfaces_with_low_confidence.insert(sorted_surfaces_with_low_confidence.end(), surfaces_with_low_confidence.begin(), surfaces_with_low_confidence.end());
        std::sort(sorted_surfaces_with_low_confidence.begin(), sorted_surfaces_with_low_confidence.end(), 
            [](const std::shared_ptr<Surface>& a, const std::shared_ptr<Surface>& b) -> bool
            {
                return a->get_total_point_size() > b->get_total_point_size();
            });

        // add to the smallest uncertainty surface
        std::shared_ptr<Surface> smallest_surface = sorted_surfaces_with_low_confidence[0];
        std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(generic_point);
        new_vertex->reduce_reverse_radius_search_radius(radius);
        new_vertex->reduce_previous_radius(radius);
        smallest_surface->connect_by_edges_and_faces(new_vertex, neighboring_vertices);
        smallest_surface->connect(new_vertex);
    }
    else
    {
        return false;
    }

    return true;

    // // recompute sibling vertices
    // std::vector<std::shared_ptr<Vertex>> sibling_vertices_copy = sibling_vertices;
    // sibling_vertices.clear();
    // for (std::shared_ptr<Vertex> vertex : sibling_vertices_copy)
    // {
    //     // skip if the vertex is expired
    //     if (vertex->is_expired()) continue;

    //     sibling_vertices.push_back(vertex);
    // }

    // // if new vertex not in matched surface
    // bool new_vertex_not_in_matched_surface = sibling_vertices.size() == 0;
    // if (new_vertex_not_in_matched_surface)
    // {
    //     // start a new seed
    //     std::shared_ptr<Surface> new_surface = storage_->add_surface();
    //     std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(generic_point);
    //     new_surface->connect(new_vertex);
    //     sibling_vertices.push_back(new_vertex);
    //     sibling_vertices[0]->connect(new_vertex);
    // }

    // // if new_vertex is in multiple surfaces, try merge them
    // if (new_vertex->get_surfaces().size() > 1)
    // {
    //     new_vertex->try_merge_surfaces();
    // }

    // // get all edges the new vertex is connected to
    // std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> all_connected_vertices;
    // for (std::shared_ptr<Edge> edge : new_vertex->get_edges())
    // {
    //     all_connected_vertices.insert(edge->get_vertex(0));
    //     all_connected_vertices.insert(edge->get_vertex(1));
    // }
    // // get smallest radius
    // double smallest_radius = std::numeric_limits<double>::max();
    // for (std::shared_ptr<Vertex> vertex : all_connected_vertices)
    // {
    //     double radius = vertex->get_radius();
    //     if (radius < smallest_radius) smallest_radius = radius;
    // }
    // // set the smallest radius to new vertex
    // new_vertex->reduce_reverse_radius_search_radius(smallest_radius);

    // // if low confidence surface is still low confidence and have more than n points, remove it
    // for (const std::shared_ptr<Surface>& surface : surfaces_with_low_confidence)
    // {
    //     bool still_low_confidence = surface->get_average_projective_distance() > settings_.projective_std_threshold;
    //     bool more_than_n_points = surface->get_total_point_size() > settings_.remove_low_confidence_threshold;
    //     if (still_low_confidence && more_than_n_points)
    //     {
    //         if (settings_.log.add_point_by_radius_search) std::cout << ">> removing low confidence surface " << surface->get_id() << std::endl;
    //         storage_->delete_surface(surface);
    //     }
    // }
}

template <typename PointT>
void Application<PointT>::add_point_by_new_surface(const std::shared_ptr<GenericPoint>& generic_point, double& radius)
{
    std::shared_ptr<Surface> new_surface = storage_->add_surface();
    std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(generic_point);
    new_vertex->reduce_reverse_radius_search_radius(radius);
    new_vertex->reduce_previous_radius(radius);
    new_surface->connect(new_vertex);
}

template <typename PointT>
void Application<PointT>::get_lidar_data(Eigen::Vector3d& origin, Eigen::Vector3d& position)
{
    origin = this->origin;
    position = pointcloud->points[ith_point].getVector3fMap().cast<double>();
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
    generic_point->initialize_(storage_, thisPointVEC, thisPointOriginVEC);
    process_point(generic_point);
}

template <typename PointT>
void Application<PointT>::loop()
{
    // add all points from the cloud to the storage's processing queue
    for (std::size_t i = 0; i < ith_size; i++)
    {
        Eigen::Vector3d thisPointVEC = pointcloud->points[i].getVector3fMap().cast<double>();
        Eigen::Vector3d thisPointOriginVEC = this->origin;

        storage_->add_to_queue(thisPointVEC, thisPointOriginVEC);
    }

    // add points from repeated queue to queue
    storage_->add_points_in_repeated_queue_to_queue();

    // process all points in the queue
    while (storage_->get_queue_size() > 0)
    {
        unsigned int total_points_in_queue = storage_->get_queue_size();
        if (settings_.log.step) std::cout << "==================================================================== Processing 1 / " << total_points_in_queue << " in queue of " << ith_cloud << " cloud" << std::endl;

        std::shared_ptr<GenericPoint> generic_point = storage_->pop_from_queue();
        process_point(generic_point);
    }

    // print repeated queue size
    if (settings_.log.step) std::cout << "==================================================================== repeated queue size: " << storage_->get_repeated_queue_size() << std::endl;

    // check tree rebuild
    storage_->check_tree_rebuild();

    // load next cloud
    ith_cloud += 1;
    load_point_cloud();
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
        if (settings.color_mode == 0)
        {
            const std::tuple<int, int, int>& color = interior_point->get_surface()->get_color();
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == 1)
        {
            double distance = std::fabs(interior_point->buffer_compute_projected_distance() / 0.05);
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == 2)
        {
            double distance = interior_point->get_sibling_interior_points().size() / settings.siblings_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == 3)
        {
            double distance = interior_point->get_radius() / settings.radius_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (settings.color_mode == 4)
        {
            double distance = interior_point->get_surface()->get_surface_position_std_in_normal_direction() / settings.positional_uncertainty_denominator;
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
        if (setting.color_mode == 0)
        {
            const std::tuple<int, int, int>& color = vertex->get_surface()->get_color();
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == 1)
        {
            double distance = std::abs(vertex->buffer_compute_projected_distance() / 0.05);
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == 2)
        {
            double distance = vertex->get_sibling_vertices().size() / setting.siblings_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == 3)
        {
            double distance = vertex->get_radius() / setting.radius_denominator;
            std::tuple<int, int, int> color = valueToJet(distance);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
        }
        else if (setting.color_mode == 4)
        {
            double distance = vertex->get_surface()->get_surface_position_std_in_normal_direction() / setting.positional_uncertainty_denominator;
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