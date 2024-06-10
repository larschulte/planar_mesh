#include <pcl/common/transforms.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include "point_type/VilensPointT.hpp"
#include "eye_patch/DataLoader.hpp"
#include "eye_patch/TriangleBVH.hpp"
#include "eye_patch/RRSTree.hpp"
#include "eye_patch/utilities.hpp"

// queue
#include <queue>

#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/Storage.hpp"

// application class
template <typename PointT>
class Application
{
public:

    std::set<int> extract_existing_edge_between_points(const std::set<int>& candidate_point_set, const std::set<int>& candidate_edge_set)
    {
        // initialize
        std::set<int> out_existing_edge_set;

        // process
        // for each pair of candidate points
        for (int point_id1 : candidate_point_set)
        {
            for (int point_id2 : candidate_point_set)
            {
                // skip if repeated
                if (point_id1 >= point_id2) continue;

                // skip if edge does not exist
                std::array<int, 2> edge = {point_id1, point_id2};
                if (edge_to_point_map_reverse.find(edge) == edge_to_point_map_reverse.end()) continue;
                
                // skip if edge not in candidate_edge_set
                if (candidate_edge_set.find(edge_to_point_map_reverse[edge]) == candidate_edge_set.end()) continue;

                // else add to set
                out_existing_edge_set.insert(edge_to_point_map_reverse[edge]);
            }
        }

        // return
        return out_existing_edge_set;
    }

    // check if triangle contains any boundary point of set
    bool triangle_set_intersection(const std::array<int, 3>& triangle, int setID, const std::map<int, Eigen::Vector2d>& precompute_2d_map)
    {
        int vertexID1 = triangle[0];
        int vertexID2 = triangle[1];
        int vertexID3 = triangle[2];
        const Eigen::Vector2d& vertex1_2D = precompute_2d_map.at(vertexID1);
        const Eigen::Vector2d& vertex2_2D = precompute_2d_map.at(vertexID2);
        const Eigen::Vector2d& vertex3_2D = precompute_2d_map.at(vertexID3);

        for (int boundary_pointID : boundary_point_of_set.at(setID))
        {
            const Eigen::Vector2d& boundary_point_2D = precompute_2d_map.at(boundary_pointID);

            // contain at ends
            if (boundary_pointID == vertexID1 || boundary_pointID == vertexID2 || boundary_pointID == vertexID3) continue;

            // contain at middle
            if (point_in_triangle(boundary_point_2D, vertex1_2D, vertex2_2D, vertex3_2D)) return true;
        }
    
        return false;
    }

    // extract set from set
    std::set<int> intersection_of_sets(const std::set<int>& setA, const std::set<int>& setB)
    {
        // initialize
        std::set<int> setOut;

        // process
        for (int elemA : setA)
        {
            if (setB.find(elemA) != setB.end()) setOut.insert(elemA);
        }
        
        // return
        return setOut;
    }

    // // creates edges and triangles that connects the new point to the set
    // // to add a new point to mesh
    // // - form edge to boundary point of the mesh, skip if the edge intersects any existing boundary edge
    // // - form triangle if two used boundary points have a boundary edge between them, skip if the triangle contains other boundary points
    // void connect_point_to_set(int newPointID, int setID, const std::set<int>& searched_boundary_points)
    // {
    //     const std::set<int>& boundary_point_of_current_set = boundary_point_of_set.at(setID);
    //     const std::set<int>& boundary_edge_of_current_set = boundary_edge_of_set.at(setID);

    //     std::set<int> searched_boundary_points_in_current_set = intersection_of_sets(searched_boundary_points, boundary_point_of_current_set);
    //     std::set<int> searched_boundary_edge_in_current_set = extract_existing_edge_between_points(searched_boundary_points_in_current_set, boundary_edge_of_current_set);

    //     // // adjust radius
    //     // // choose the smallest radius among the boundary points
    //     // double smallest_radius = std::numeric_limits<double>::max();
    //     // for (int point_id : searched_boundary_points_in_current_set)
    //     // {
    //     //     double distance = (point_to_vector3d_map.at(newPointID) - point_to_vector3d_map.at(point_id)).norm();
    //     //     if (distance < smallest_radius) smallest_radius = distance;
    //     // }
    //     // rrstree.reduceRadius(newPointID, point_to_vector3d_map.at(newPointID), smallest_radius);

    //     // compute the smallest distance to searched boundary points that are not in the same set, and is not small set
    //     double smallest_distance = std::numeric_limits<double>::max();
    //     for (int point_id : searched_boundary_points)
    //     {
    //         if (point_to_set_map.at(point_id) == setID) continue;
    //         // skip if the set is small in size
    //         if (set_to_points_map.at(point_to_set_map.at(point_id)).size() < fit_plane_threshold) continue;
            
    //         double distance = (point_to_vector3d_map.at(newPointID) - point_to_vector3d_map.at(point_id)).norm();
    //         if (distance < smallest_distance) smallest_distance = distance;
    //     }
    //     if (smallest_distance < std::numeric_limits<double>::max()) rrstree.reduceRadius(newPointID, point_to_vector3d_map.at(newPointID), smallest_distance);

    //     // precompute 2d points
    //     std::map<int, Eigen::Vector2d> precompute_2d_map;
    //     for (int point_id : boundary_point_of_current_set) 
    //     {
    //         precompute_2d_map[point_id] = project_point_to_set_plane(point_id, setID);
    //     }

    //     // // get a queue of candidate edges
    //     // std::priority_queue<
    //     //     std::pair<double, std::array<int, 2>>, 
    //     //     std::vector<std::pair<double, std::array<int, 2>>>, 
    //     //     std::greater<std::pair<double, std::array<int, 2>>>
    //     // > edge_queue;
    //     // for (int point_id : searched_boundary_points_in_current_set)
    //     // {
    //     //     // new edge, smaller id first
    //     //     std::array<int, 2> newEdge = {std::min(newPointID, point_id), std::max(newPointID, point_id)};

    //     //     // skip if intersected with any boundary edge of the current set
    //     //     bool intersected = edge_set_intersection(newEdge, setID, precompute_2d_map);
    //     //     if (intersected) continue;

    //     //     // compute distance in 2d map
    //     //     double distance = (precompute_2d_map.at(newEdge[0]) - precompute_2d_map.at(newEdge[1])).norm();

    //     //     // add to queue
    //     //     edge_queue.push(std::make_pair(distance, newEdge));
    //     // }

    //     // // add from shortest edge, until two triangles are formed
    //     // std::set<int> searched_boundary_points_used;
    //     // int triangles_added = 0;
    //     // while (!edge_queue.empty() && triangles_added < 3)
    //     // {
    //     //     // get edge
    //     //     std::array<int, 2> newEdge = edge_queue.top().second;
    //     //     edge_queue.pop();

    //     //     // add edge
    //     //     int newEdgeID = getNewEdgeID();
    //     //     add_edge(newEdgeID, setID, newEdge);

    //     //     // add to used
    //     //     searched_boundary_points_used.insert(newEdge[0]);
    //     //     searched_boundary_points_used.insert(newEdge[1]);

    //     //     // add triangle
    //     //     for (const auto& edgeID : searched_boundary_edge_in_current_set)
    //     //     {   
    //     //         // skip if not both points are used
    //     //         int i1 = edge_to_point_map.at(edgeID)[0];
    //     //         int i2 = edge_to_point_map.at(edgeID)[1];
    //     //         bool i1_used = searched_boundary_points_used.find(i1) != searched_boundary_points_used.end();
    //     //         bool i2_used = searched_boundary_points_used.find(i2) != searched_boundary_points_used.end();
    //     //         if (!i1_used || !i2_used) continue;

    //     //         // new triangle, smaller id first
    //     //         std::array<int, 3> newTriangle = sortThreeInts(newPointID, i1, i2);

    //     //         // skip if triangle already exists
    //     //         if (triangle_to_vertices_map_reverse.find(newTriangle) != triangle_to_vertices_map_reverse.end()) continue;

    //     //         // skip if triangle contains other boundary points
    //     //         if (triangle_set_intersection(newTriangle, setID, precompute_2d_map)) continue;

    //     //         // add triangle
    //     //         int newTriangleID = getNewTriangleID();
    //     //         add_triangle(newTriangleID, setID, newTriangle);
    //     //         triangles_added ++;
    //     //     }
    //     // }

    //     // add edge
    //     std::set<int> searched_boundary_points_used;
    //     std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> boundary_edges;
    //     // Precompute boundary edges and their coordinates
    //     for (int boundary_edgeID : boundary_edge_of_set.at(setID)) 
    //     {
    //         const std::array<int, 2>& boundaryEdge = edge_to_point_map.at(boundary_edgeID);
    //         const Eigen::Vector2d& boundaryPoint1_2D = precompute_2d_map.at(boundaryEdge[0]);
    //         const Eigen::Vector2d& boundaryPoint2_2D = precompute_2d_map.at(boundaryEdge[1]);
    //         boundary_edges.emplace_back(boundaryPoint1_2D, boundaryPoint2_2D);
    //     }
    //     for (int point_id : searched_boundary_points_in_current_set) 
    //     {
    //         // new edge, smaller id first
    //         std::array<int, 2> newEdge = {std::min(newPointID, point_id), std::max(newPointID, point_id)};
    //         const Eigen::Vector2d& newPoint1_2D = precompute_2d_map.at(newEdge[0]);
    //         const Eigen::Vector2d& newPoint2_2D = precompute_2d_map.at(newEdge[1]);

    //         // skip if intersected with any boundary edge of the current set
    //         bool intersected = false;
    //         for (const auto& boundaryEdge : boundary_edges) 
    //         {
    //             const Eigen::Vector2d& boundaryPoint1_2D = boundaryEdge.first;
    //             const Eigen::Vector2d& boundaryPoint2_2D = boundaryEdge.second;

    //             // intersect at ends
    //             if (newPoint1_2D == boundaryPoint1_2D || newPoint1_2D == boundaryPoint2_2D || newPoint2_2D == boundaryPoint1_2D || newPoint2_2D == boundaryPoint2_2D) continue;

    //             // intersect at middle
    //             if (doIntersect(newPoint1_2D, newPoint2_2D, boundaryPoint1_2D, boundaryPoint2_2D)) { intersected = true; break; }
    //         }
            
    //         if (intersected) continue;

    //         // add edge
    //         int newEdgeID = getNewEdgeID();
    //         add_edge(newEdgeID, setID, newEdge);

    //         // add to used
    //         searched_boundary_points_used.insert(point_id);
    //     }

    //     // add triangle
    //     for (const auto& edgeID : searched_boundary_edge_in_current_set)
    //     {   
    //         // skip if not both points are used
    //         int i1 = edge_to_point_map.at(edgeID)[0];
    //         int i2 = edge_to_point_map.at(edgeID)[1];
    //         bool i1_used = searched_boundary_points_used.find(i1) != searched_boundary_points_used.end();
    //         bool i2_used = searched_boundary_points_used.find(i2) != searched_boundary_points_used.end();
    //         if (!i1_used || !i2_used) continue;

    //         // new triangle, smaller id first
    //         std::array<int, 3> newTriangle = sortThreeInts(newPointID, i1, i2);

    //         // skip if triangle already exists
    //         if (triangle_to_vertices_map_reverse.find(newTriangle) != triangle_to_vertices_map_reverse.end()) continue;

    //         // skip if triangle contains other boundary points
    //         if (triangle_set_intersection(newTriangle, setID, precompute_2d_map)) continue;

    //         // add triangle
    //         int newTriangleID = getNewTriangleID();
    //         add_triangle(newTriangleID, setID, newTriangle);
    //     }
    // }


    // merge setID2 into setID1
    std::weak_ptr<Surface> merge_surfaces(std::weak_ptr<Surface> surface1, std::weak_ptr<Surface> surface2) 
    {
        std::weak_ptr<Surface> new_surface = storage_->add_surface(surface1, surface2);
        storage_->delete_surface(surface1);
        storage_->delete_surface(surface2);

        return new_surface;
    }

    Eigen::Matrix3d merge_covariances_of_surfaces(std::weak_ptr<Surface> surface1, std::weak_ptr<Surface> surface2) 
    {
        Eigen::Matrix3d cov1 = surface1.lock()->get_covariance();
        Eigen::Matrix3d cov2 = surface2.lock()->get_covariance();
        Eigen::Vector3d mean1 = surface1.lock()->get_mean();
        Eigen::Vector3d mean2 = surface2.lock()->get_mean();
        int size1 = surface1.lock()->get_total_point_size();
        int size2 = surface2.lock()->get_total_point_size();
        return merge_covariances(cov1, cov2, mean1, mean2, size1, size2);
    }

    // try merge sets
    void try_merge_surfaces(std::set<std::weak_ptr<Surface>>& surfaces_to_merge)
    {
        while (true) 
        {
            // get all possible pairs to merge
            std::set<std::pair<std::weak_ptr<Surface>, std::weak_ptr<Surface>>> surface_pairs;
            for (std::weak_ptr<Surface> surface1 : surfaces_to_merge) 
            {
                for (std::weak_ptr<Surface> surface2 : surfaces_to_merge) 
                {
                    if (surface1 >= surface2) continue;
                    surface_pairs.insert(std::make_pair(surface1, surface2));
                }
            }
            
            // try to merge pairs
            bool again = false;
            for (const auto& pairs : surface_pairs) 
            {
                // skip if can't merge
                Eigen::Matrix3d covariance_matrix = merge_covariances_of_surfaces(pairs.first, pairs.second);
                double eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance_matrix).eigenvalues()[0];
                if (eigenvalue > merged_eigenvalue_threshold) continue;

                // merge sets
                surfaces_to_merge.erase(pairs.first);
                surfaces_to_merge.erase(pairs.second);
                std::weak_ptr<Surface> newSurface = merge_surfaces(pairs.first, pairs.second);
                surfaces_to_merge.insert(newSurface);

                // once merged, restart
                again = true;
                break;
            }
            if (!again) break;
        }
    }

    void add_point_by_radius_search(const Eigen::Vector3d& thisPointVEC, const Eigen::Vector3d& thisPointOriginVEC)
    {
        // if empty, can not set up radius search, add point to new set
        if (!storage_->can_reverse_radius_search())
        {
            std::weak_ptr<Surface> new_surface = storage_->add_surface();
            std::weak_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            new_surface.lock()->connect(new_vertex);
            return;
        }

        // perform rrstree radius search
        std::map<int, double> point_to_radius_map;
        std::set<std::weak_ptr<Vertex>> searched_boundary_vertices_set = storage_->reverse_radius_search(thisPointVEC);

        // if no searched results, add point to new set
        if (searched_boundary_vertices_set.size() == 0)
        {
            std::weak_ptr<Surface> new_surface = storage_->add_surface();
            std::weak_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            new_surface.lock()->connect(new_vertex);
            return;
        }

        // from searched points identify neighboring sets
        std::set<std::weak_ptr<Surface>> neighboring_surfaces; 
        for (std::weak_ptr<Vertex> vertex : searched_boundary_vertices_set)
        {
            neighboring_surfaces.insert(vertex.lock()->get_surface());
        }

        // try merge neighboring sets
        try_merge_surfaces(neighboring_surfaces);

        // after merging, we are left with sets that should have different normals
        // each point in the neighboring set should reduce their search radius to the closest set
        // [todo]
        // for each point, compute the shortest distance to another point that is in a different set, and that different set have enough points
        // if the distance is less than its original radius, reduce its original radius
        for (std::weak_ptr<Vertex> vertex : searched_boundary_vertices_set)
        {
            // surface
            std::weak_ptr<Surface> surface = vertex.lock()->get_surface();

            // skip if small set
            if (surface.lock()->get_total_point_size() < fit_plane_threshold) continue;

            // smallest distance
            double smallest_distance = std::numeric_limits<double>::max();

            for (std::weak_ptr<Vertex> other_vertex : searched_boundary_vertices_set)
            {
                // surface 
                std::weak_ptr<Surface> other_surface = other_vertex.lock()->get_surface();

                // skip if same set
                if (other_surface == surface) continue;

                // skip if small set
                if (other_surface.lock()->get_total_point_size() < fit_plane_threshold) continue;

                // compute distance
                double distance = (vertex.lock()->get_position() - other_vertex.lock()->get_position()).norm();

                // update radius
                if (distance < smallest_distance) smallest_distance = distance;
            }

            // adjust radius
            if (smallest_distance < vertex.lock()->get_radius()) 
            {
                vertex.lock()->set_reverse_radius_search_radius(smallest_distance);
            }
        }



        // split neighboring sets into sets with plane and sets without plane (by size)
        std::set<std::weak_ptr<Surface>> surfaces_with_plane;
        std::set<std::weak_ptr<Surface>> surfaces_without_plane;
        for (std::weak_ptr<Surface> surface : neighboring_surfaces)
        {
            if (surface.lock()->get_total_point_size() > fit_plane_threshold) surfaces_with_plane.insert(surface);
            else surfaces_without_plane.insert(surface);
        }
        
        // for sets with plane, compute the point to set intersection distance
        std::map<std::weak_ptr<Surface>, double> surface_distance_map;
        for (std::weak_ptr<Surface> surface : surfaces_with_plane)
        {
            // compute (not using normal from combined points, could implement in future)
            double distance = surface.lock()->compute_point_to_surface_distance(thisPointOriginVEC, thisPointVEC);

            // store
            surface_distance_map[surface] = std::fabs(distance);
        }

        // extract the set within distance threshold
        std::set<std::weak_ptr<Surface>> surfaces_within_threshold;
        for (const auto& pair : surface_distance_map)
        {
            if (pair.second < distance_threshold) surfaces_within_threshold.insert(pair.first);
        }

        // from the sets within threshold, find the set that is closest to the point
        std::weak_ptr<Surface> closest_surface;
        double closest_distance = std::numeric_limits<double>::max();
        for (std::weak_ptr<Surface> surface : surfaces_within_threshold)
        {
            double distance = surface_distance_map.at(surface);

            // update if closer
            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest_surface = surface;
            }
        }
        if (!closest_surface.expired() && closest_distance < distance_threshold)
        {
            // for the sets not selected as closest, update their searched points' radius
            // for all searched points
            for (std::weak_ptr<Vertex> vertex : searched_boundary_vertices_set)
            {
                // if it is in a set with plane
                if (surfaces_with_plane.find(vertex.lock()->get_surface()) != surfaces_with_plane.end())
                {
                    // and the set with plane is not the closest set
                    if (vertex.lock()->get_surface() != closest_surface)
                    {
                        // reduce their searched points' radius
                        double reduced_radius = (thisPointVEC - vertex.lock()->get_position()).norm();

                        if (reduced_radius < vertex.lock()->get_radius())
                        {
                            vertex.lock()->set_reverse_radius_search_radius(reduced_radius);
                        }
                    }
                }
            }

            // add the point to the closest set
            std::weak_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            closest_surface.lock()->connect(new_vertex);
            return;
        }

        // else, find the set that is nearest
        std::weak_ptr<Surface> nearest_surface;
        double nearest_distance = std::numeric_limits<double>::max();
        for (std::weak_ptr<Surface> surface : surfaces_without_plane)
        {
            Eigen::Vector3d mean = surface.lock()->get_mean();
            double distance = (thisPointVEC - mean).norm();
            if (distance < nearest_distance)
            {
                nearest_distance = distance;
                nearest_surface = surface;
            }
        }
        if (!nearest_surface.expired())
        {
            std::weak_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            nearest_surface.lock()->connect(new_vertex);
            return;
        }

        // else, add the point to a new set
        std::weak_ptr<Surface> new_surface = storage_->add_surface();
        std::weak_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
        new_surface.lock()->connect(new_vertex);
        return;
    }

    void load_point_cloud()
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
            // shuffle the pointcloud
            std::random_shuffle(pointcloud->points.begin(), pointcloud->points.end());
        }
    }

    void step()
    {        
        // get point
        Eigen::Vector3d thisPointVEC = pointcloud->points[ith_point].getVector3fMap().cast<double>();
        Eigen::Vector3d thisPointOriginVEC = origin;
        ith_point ++;
    
        // ------------- add point by triangle intersection

        // get list of intersected triangle by the point
        std::set<std::weak_ptr<Face>> searched_faces = bvhRoot.intersectionSearch(thisPointOriginVEC, thisPointVEC); // may include deleted triangles

        // group the faces by surface
        std::map<std::weak_ptr<Surface>, std::set<std::weak_ptr<Face>>> surface_to_searched_faces_map;
        for (std::weak_ptr<Face> face : searched_faces)
        {
            std::weak_ptr<Surface> surface = face.lock()->get_surface();
            surface_to_searched_faces_map[surface].insert(face);
        }

        // compute the intersection distance to the sets (distance measured in plane normal direction)
        std::map<std::weak_ptr<Surface>, double> surface_to_point_distance_map;
        for (const auto& pair : surface_to_searched_faces_map)
        {
            std::weak_ptr<Surface> surface = pair.first;
            double distance = surface.lock()->compute_point_to_surface_distance(thisPointOriginVEC, thisPointVEC);
            surface_to_point_distance_map[surface] = distance;
        }

        // split the sets into three categories
        std::set<std::weak_ptr<Surface>> surface_with_point_before_it;
        std::set<std::weak_ptr<Surface>> surface_with_point_within_it;
        std::set<std::weak_ptr<Surface>> surface_with_point_behind_it;
        double split_distance_threshold = distance_threshold;
        for (const auto& pair : surface_to_point_distance_map)
        {
            std::weak_ptr<Surface> surface = pair.first;
            double distance = pair.second;
            if (distance > split_distance_threshold) 
            {
                surface_with_point_before_it.insert(surface);
            }
            else if (distance < -split_distance_threshold) 
            {
                surface_with_point_behind_it.insert(surface);
            }
            else 
            {
                surface_with_point_within_it.insert(surface);
            }
        }
        
        // process point behind set set

        // get the set of triangles that are penetrated by the point
        std::set<std::weak_ptr<Face>> penetrated_faces;
        for (std::weak_ptr<Surface> surface : surface_with_point_behind_it)
        {
            penetrated_faces.insert(surface_to_searched_faces_map.at(surface).begin(), surface_to_searched_faces_map.at(surface).end());
        }

        // collect the list of points that are within the penetrated triangles, and isolated by the triangle
        for (std::weak_ptr<Face> face : penetrated_faces)
        {
            // for now, just delete the face
            storage_->delete_face(face);
        }

        // // re add the list points by radius search
        // // todo - while avoid covering the new point
        // // process point in the queue free_points_vector3d_and_origin_vector3d
        // while (!free_points_queue.empty())
        // {
        //     std::pair<Eigen::Vector3d, Eigen::Vector3d> pair = free_points_queue.front(); free_points_queue.pop();

        //     std::cout << "---------------------------------------------------------------------------------------------------- re-adding point by radius search" << std::endl;
        //     Eigen::Vector3d pointVEC = pair.first;
        //     Eigen::Vector3d pointOriginVEC = pair.second;
        //     int pointID = getNewPointID();
        //     add_point_by_radius_search(pointID, pointVEC, pointOriginVEC);
        // }
        
        // process point within set set
        bool point_added_to_surface = false;

        // try merge them
        std::set<std::weak_ptr<Surface>> merged_surface = surface_with_point_within_it;
        try_merge_surfaces(merged_surface);        

        // find the set with the smallest distance
        std::weak_ptr<Surface> smallest_surface;
        double smallest_distance = std::numeric_limits<double>::max();
        for (std::weak_ptr<Surface> surface : merged_surface)
        {
            // update if smaller
            double distance = surface.lock()->compute_point_to_surface_distance(thisPointOriginVEC, thisPointVEC);
            if (std::abs(distance) < smallest_distance)
            {
                smallest_distance = distance;
                smallest_surface = surface;
            }
        }

        // if the smallest set is within threshold, add the point to the set
        if (!smallest_surface.expired() && std::abs(smallest_distance) < distance_threshold)
        {
            // from searched_faces find the first face that belongs to the smallest surface
            for (std::weak_ptr<Face> face : searched_faces)
            {
                if (face.lock()->get_surface() != smallest_surface) continue;
                
                // add point as interior point
                storage_->add_interior_point(face, thisPointVEC, thisPointOriginVEC);

                std::cout << ith_point << " / " << ith_size << " of pointcloud " << ith_cloud << " added to set " << smallest_surface.lock()->get_id() << std::endl;
                point_added_to_surface = true;
                break;
            }
        }
        if (!point_added_to_surface)
        {
            add_point_by_radius_search(thisPointVEC, thisPointOriginVEC);
            std::cout << ith_point << " / " << ith_size << " of pointcloud " << ith_cloud << " added by radius search" << std::endl;
        }


        // todo - process point before set set


        // check if end of point cloud
        if (ith_point == ith_size) 
        {
            // next cloud
            ith_cloud += 1;
            ith_point = 0;
            load_point_cloud();
        }
    }

    void loop()
    {
        step();

        // finish all points in step
        while (ith_point != 0)
        {
            step();
        }
    }

    Application() 
    {
        // // ----------------------- DATA
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

        // // ----------------------- PARAMETERS
        distance_threshold = 0.05;
        fit_plane_threshold = 10; // may cause error if below 3
        merged_eigenvalue_threshold = 15e-5;
        

        dataset = "room";
        // dataset = "osney";
        // dataset = "blenheim";
        ith_cloud = 50;
        ith_point = 0;
        shuffle_pointcloud = false;
        pointcloud_fraction = 1;
        distance_to_radius_ratio = tan(4 * M_PI / 180);

        storage_ = std::make_shared<Storage>();        

        // // ----------------------- INITIALIZATION
        data_loader.load_dataset(dataset_map.at(dataset).first, dataset_map.at(dataset).second);
        load_point_cloud();
    }

    std::map<int, Eigen::Vector3d> get_point_to_vector3d_map() {return point_to_vector3d_map;};
    std::map<int, int> pointID_to_cloud_index_map()
    {
        // initialize
        std::map<int, int> pointID_to_cloud_index;

        // process
        int cloud_index = 0;
        for (int point_id : point_list)
        {
            pointID_to_cloud_index[point_id] = cloud_index;
            cloud_index ++;
        }

        // return
        return pointID_to_cloud_index;
    }
    std::map<int, std::array<int, 2>> get_edge_to_cloud_indices_map() 
    {
        // initialize
        std::map<int, std::array<int, 2>> edge_to_cloud_indices_map;

        // process
        std::map<int, int> pointID_to_cloud_index = pointID_to_cloud_index_map();
        for (int edge_id : edge_list)
        {
            std::array<int, 2> edge = edge_to_point_map.at(edge_id);
            std::array<int, 2> cloud_index = {pointID_to_cloud_index.at(edge[0]), pointID_to_cloud_index.at(edge[1])};
            edge_to_cloud_indices_map[edge_id] = cloud_index;
        }

        // return
        return edge_to_cloud_indices_map;
    };

    std::map<std::weak_ptr<Vertex>, int> get_vertex_to_cloud_indices_map()
    {
        return storage_->get_vertex_to_cloud_indices_map();
    } 

    std::set<std::weak_ptr<Face>> get_faces() {return storage_->get_faces();};
    std::set<std::weak_ptr<Edge>> get_edges() {return storage_->get_edges();};
    std::vector<std::weak_ptr<Vertex>> get_rrs_vertices() {return storage_->get_rrs_vertices();};
    // std::set<std::weak_ptr<Vertex>> get_boundary_edges() 
    // {
    //     // all edges
    //     std::set<std::weak_ptr<Edge>> edges = storage_->get_edges();

    //     // boundary edges
    //     std::set<std::weak_ptr<Edge>> boundary_edges;
    //     for (std::weak_ptr<Edge> edge : edges)
    //     {
    //         if (edge.lock()->is_boundary()) 
    //         {
    //             boundary_edges.insert(edge);
    //         }
    //     }
    // }

    std::map<int, std::array<int, 3>> get_triangle_to_cloud_indices_map() 
    {
        // initialize
        std::map<int, std::array<int, 3>> triangle_to_cloud_index_map;

        // process
        std::map<int, int> pointID_to_cloud_index = pointID_to_cloud_index_map();
        for (int triangle_id : triangle_list)
        {
            std::array<int, 3> triangle = triangle_to_vertices_map.at(triangle_id);
            std::array<int, 3> cloud_index = {pointID_to_cloud_index.at(triangle[0]), pointID_to_cloud_index.at(triangle[1]), pointID_to_cloud_index.at(triangle[2])};
            triangle_to_cloud_index_map[triangle_id] = cloud_index;
        }

        // return
        return triangle_to_cloud_index_map;
    };
    std::set<int> get_boundary_edge_set() {return boundary_edge_set;};

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_colored_cloud()
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (std::weak_ptr<Vertex> vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            Eigen::Vector3d position = vertex.lock()->get_position();
            std::tuple<int, int, int> color = vertex.lock()->get_surface().lock()->get_color();
            point.x = position[0];
            point.y = position[1];
            point.z = position[2];
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_distance_cloud()
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (std::weak_ptr<Vertex> vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            Eigen::Vector3d origin = vertex.lock()->get_origin();
            Eigen::Vector3d position = vertex.lock()->get_position();
            point.x = position[0];
            point.y = position[1];
            point.z = position[2];
            double value = vertex.lock()->get_surface().lock()->compute_point_to_surface_distance(origin, position) / 0.05;
            std::tuple<int, int, int> color = valueToJet(value);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_colored_cloud()
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (std::weak_ptr<Vertex> vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            Eigen::Vector3d projected_position = vertex.lock()->get_projected_position();
            std::tuple<int, int, int> color = vertex.lock()->get_surface().lock()->get_color();
            point.x = projected_position[0];
            point.y = projected_position[1];
            point.z = projected_position[2];
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_distance_cloud()
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (std::weak_ptr<Vertex> vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            Eigen::Vector3d origin = vertex.lock()->get_origin();
            Eigen::Vector3d position = vertex.lock()->get_position();
            Eigen::Vector3d projected_position = vertex.lock()->get_projected_position();
            point.x = projected_position[0];
            point.y = projected_position[1];
            point.z = projected_position[2];
            double value = vertex.lock()->get_surface().lock()->compute_point_to_surface_distance(origin, position) / 0.05;
            std::tuple<int, int, int> color = valueToJet(value);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
        }
        return cloud;
    }

    int get_number_of_triangles()
    {
        return triangle_to_vertices_map.size();
    }

    void change_color()
    {
        std::set<std::weak_ptr<Surface>> surfaces = storage_->get_surfaces();
        for (std::weak_ptr<Surface> surface : surfaces)
        {
            surface.lock()->set_random_color();
        }
    }

    void rrstree_rebuild()
    {
        rrstree.rebuild();
    }

    void rrstree_print_tree()
    {
        rrstree.print();
    }

    void rrstree_print_size()
    {
        rrstree.print_size();
    }

    int ith_cloud;
    std::size_t ith_point = 0;
    std::size_t ith_size = 0;

private:
    // storage
    std::shared_ptr<Storage> storage_;


    // data
    DataLoader<VilensPointT> data_loader;
    typename pcl::PointCloud<VilensPointT>::Ptr pointcloud;
    Eigen::Vector3d origin;
    

    // settings
    double distance_threshold;
    std::size_t fit_plane_threshold; // may cause error if below 3
    double merged_eigenvalue_threshold;
    bool shuffle_pointcloud;
    double pointcloud_fraction;
    std::string dataset;
    double distance_to_radius_ratio;

        // free points
    std::queue<std::pair<Eigen::Vector3d, Eigen::Vector3d>> free_points_queue;

        // point
    int next_point_id = 0;
    std::vector<int> point_list;
    std::map<int, Eigen::Vector3d> point_to_origin_vector3d_map;
    std::map<int, Eigen::Vector3d> point_to_intersection_vector3d_map;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, int> point_to_set_map;
    std::map<int, std::set<int>> point_to_edges_map;
    std::map<int, int> point_to_storing_triangle_map;

        // triangle
    int next_triangle_id = 0;
    std::vector<int> triangle_list;
    std::map<int, std::array<int, 3>> triangle_to_vertices_map;
    std::map<std::array<int, 3>, int> triangle_to_vertices_map_reverse;
    std::map<int, int> triangle_to_set_map;
    std::map<int, std::set<int>> triangle_to_points_map;
    
        // set
    int next_set_id = 0;
    std::vector<int> set_list;
    std::map<int, std::array<int, 3>> set_to_color_map; // map to intensity
    std::map<int, std::set<int>> set_to_points_map; // each set contains id to points
    std::map<int, std::set<int>> set_to_edges_map; // each set contains id to edges
    std::map<int, std::set<int>> set_to_triangles_map; // each set contains id to triangles
    

        // plane fitting
    std::map<int, Eigen::Vector3d> set_to_mean_map;
    std::map<int, Eigen::Matrix3d> set_to_covariance_matrix_map;
    std::map<int, Eigen::Matrix3d> set_to_eigenvectors_map;
    std::map<int, Eigen::Vector3d> set_to_eigenvalues_map;
    std::map<int, Eigen::Vector3d> set_to_normal_map;
    
        // edge
    int next_edge_id = 0;
    std::vector<int> edge_list;
    std::map<int, std::array<int, 2>> edge_to_point_map;
    std::map<std::array<int, 2>, int> edge_to_point_map_reverse;
    std::map<int, int> edge_to_set_map;

        // boundary
    RRSTree rrstree;
    std::set<int> boundary_point_set;
    std::set<int> boundary_edge_set;
    std::map<int, std::set<int>> boundary_point_of_set;
    std::map<int, std::set<int>> boundary_edge_of_set;
    std::map<int, int> edge_to_edge_count_map; // each map contains edge count

        // triangle intersection
    TriangleBVH bvhRoot;
    std::set<int> global_triangle_set;

        // projected points
    std::map<int, Eigen::Vector3d> projected_points_to_vector3d_map;
    std::map<int, double> projected_points_distance_map;
};


// interactive viewer class
template <typename PointT>
class InteractiveViewer 
{
public:
    InteractiveViewer(Application<PointT>& app) 
        : 
        app_(app),
        viewer_(new pcl::visualization::PCLVisualizer ("3D Viewer"))
    {   
        // turn off warning
        viewer_->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
        
        // set up viewports
        viewer_->setBackgroundColor (0, 0, 0);
   
        // set up coordinate system
        viewer_->initCameraParameters();
        viewer_->addCoordinateSystem(1);

        // register keyboard callback
        viewer_->registerKeyboardCallback(&InteractiveViewer::keyboard_callback, *this, nullptr);

        // // ------------------------------ parameters
        number_of_spheres_to_display = 60;

        // spin
        viewer_->spin();


        
    }

private:
    Application<PointT>& app_;

    pcl::visualization::PCLVisualizer::Ptr viewer_;

    bool show_pointcloud = true;
    bool show_triangle = true;
    bool show_edge = true;

    bool show_projected_point = false;
    bool show_error_color = false;

    bool show_wireframe = true;
    
    std::vector<std::string> sphere_name_list;
    bool show_sphere = false;

    int number_of_spheres_to_display;

    void update_display()
    {
        // data
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_cloud;
        if (show_projected_point)
        {
            if (show_error_color) 
            {
                point_cloud = app_.projected_point_to_vector3d_set_distance_cloud();
            }
            else
            {
                point_cloud = app_.projected_point_to_vector3d_set_colored_cloud();
            }
        }
        else
        {
            if (show_error_color) 
            {
                point_cloud = app_.point_to_vector3d_set_distance_cloud();
            }
            else
            {
                point_cloud = app_.point_to_vector3d_set_colored_cloud();
            }
        }
        std::map<std::weak_ptr<Vertex>, int> vertex_to_cloud_indices_map = app_.get_vertex_to_cloud_indices_map();
        std::set<std::weak_ptr<Face>> faces = app_.get_faces();
        std::set<std::weak_ptr<Edge>> edges = app_.get_edges();
        // std::set<std::weak_ptr<Edge>> boundary_edges = app_.get_boundary_edges();

        // point cloud
        viewer_->removeShape("point_cloud");
        if (show_pointcloud)
        {
            pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> color_handler(point_cloud);
            viewer_->addPointCloud<pcl::PointXYZRGB>(point_cloud, color_handler, "point_cloud");
            viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "point_cloud");
        }
        
        // triangle mesh
        viewer_->removeShape("triangle_mesh");
        if (show_triangle)
        {
            pcl::PolygonMesh triangle_mesh;
            pcl::toPCLPointCloud2(*point_cloud, triangle_mesh.cloud);
            for (const std::weak_ptr<Face>& face : faces)
            {
                pcl::Vertices triangle;
                triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face.lock()->get_vertex(0)));
                triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face.lock()->get_vertex(1)));
                triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face.lock()->get_vertex(2)));
                triangle_mesh.polygons.push_back(triangle);
            }
            viewer_->addPolygonMesh(triangle_mesh, "triangle_mesh");
        }

        // // boundary edges
        // viewer_->removeShape("boundary_edges");        
        // if (show_edge)
        // {
        //     pcl::PolygonMesh boundary_mesh;
        //     pcl::toPCLPointCloud2(*point_cloud, boundary_mesh.cloud);
        //     for (const std::weak_ptr<Edge>& edge : boundary_edges)
        //     {
        //         pcl::Vertices boundary_edge;
        //         boundary_edge.vertices.push_back(vertex_to_cloud_indices_map.at(edge.lock()->get_vertex(0)));
        //         boundary_edge.vertices.push_back(vertex_to_cloud_indices_map.at(edge.lock()->get_vertex(1)));
        //         boundary_mesh.polygons.push_back(boundary_edge);
        //     }
        //     viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");
        //     viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "boundary_edges");
        //     viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "boundary_edges");
        // }

        // boundary points spheres
        for (const std::string& sphere_name : sphere_name_list) viewer_->removeShape(sphere_name);
        sphere_name_list.clear();
        if (show_sphere)
        {
            std::vector<std::weak_ptr<Vertex>> boundary_vertices = app_.get_rrs_vertices();
            // sort
            std::sort(boundary_vertices.begin(), boundary_vertices.end());

            // for the last 20
            for (int i = 0; i < std::min(number_of_spheres_to_display, (int)boundary_vertices.size()); i++)
            {
                std::weak_ptr<Vertex> boundary_vertex = boundary_vertices[boundary_vertices.size() - 1 - i];
                std::string sphere_name = "boundary_point_" + std::to_string(boundary_vertex.lock()->get_id());
                sphere_name_list.push_back(sphere_name);
                Eigen::Vector3d position = boundary_vertex.lock()->get_position();
                viewer_->addSphere(pcl::PointXYZ(position[0], position[1], position[2]), boundary_vertex.lock()->get_radius(), 1, 1, 1, sphere_name);
            }
        }


        // display mode
        if (show_wireframe)
        {
            viewer_->setRepresentationToWireframeForAllActors();
        }
        else
        {
            viewer_->setRepresentationToSurfaceForAllActors();
        }
    }

    void keyboard_callback(const pcl::visualization::KeyboardEvent &event, void*) 
    {
        // std::cout << "key pressed: [" << event.getKeySym() << "]" << std::endl;
        

        if (event.getKeySym() == "space" && event.keyDown())
        {
            app_.loop();
            update_display();
        }
        if (event.getKeySym() == "KP_Insert" && event.keyDown())
        {
            for (int i = 0; i < 100; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "KP_Delete" && event.keyDown())
        {
            for (int i = 0; i < 1000; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "KP_Enter" && event.keyDown())
        {
            app_.loop();
            update_display();
        }
        if (event.getKeySym() == "1" && event.keyDown())
        {
            app_.step();
            update_display();
        }
        if (event.getKeySym() == "2" && event.keyDown())
        {
            for (int i = 0; i < 10; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "3" && event.keyDown())
        {
            for (int i = 0; i < 100; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "4" && event.keyDown())
        {
            for (int i = 0; i < 1000; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "0" && event.keyDown())
        {
            app_.loop();
            update_display();
        }
        if (event.getKeySym() == "Tab" && event.keyDown())
        {
            app_.change_color();
            update_display();
        }
        if (event.getKeySym() == "comma" && event.keyDown())
        {
            // toggle show point cloud
            show_pointcloud = !show_pointcloud;
            update_display();
        }
        if (event.getKeySym() == "period" && event.keyDown())
        {
            // toggle show edge
            show_edge = !show_edge;
            update_display();
        }
        if (event.getKeySym() == "slash" && event.keyDown())
        {
            // toggle show triangle
            show_triangle = !show_triangle;
            update_display();
        }
        if (event.getKeySym() == "a" && event.keyDown())
        {
            // toggle projected point 
            show_projected_point = !show_projected_point;
            update_display();
        }
        if (event.getKeySym() == "z" && event.keyDown())
        {
            // toggle set color and error color
            show_error_color = !show_error_color;
            update_display();
        }
        if (event.getKeySym() == "v" && event.keyDown())
        {
            // toggle wireframe
            show_wireframe = !show_wireframe;
            update_display();
        }
        if (event.getKeySym() == "b" && event.keyDown())
        {
            // rebuild rrstree
            app_.rrstree_rebuild();
            std::cout << "rrstree rebuilt" << std::endl;
        }
        if (event.getKeySym() == "n" && event.keyDown())
        {
            // print tree
            app_.rrstree_print_tree();
        }
        if (event.getKeySym() == "KP_Next" && event.keyDown())
        {
            app_.ith_cloud += 1;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_End" && event.keyDown())
        {
            app_.ith_cloud -= 1;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Right" && event.keyDown())
        {
            app_.ith_cloud += 10;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Left" && event.keyDown())
        {
            app_.ith_cloud -= 10;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Prior" && event.keyDown())
        {
            app_.ith_cloud += 100;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Home" && event.keyDown())
        {
            app_.ith_cloud -= 100;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "m" && event.keyDown())
        {
            show_sphere = !show_sphere;
            update_display();
        }
    }  
};

using InputPointT = VilensPointT;
int main()
{
    std::srand(42); // Fixed seed
    // test by print
    std::cout << std::rand() << std::endl;


    // application
    Application<InputPointT> app;

    // interactive viewer
    InteractiveViewer<InputPointT> iviewer(app);

   return 0;
}