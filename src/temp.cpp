#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"

#include "eye_patch/BVH.hpp"

template <typename PointT>
class flann3d
{
public:
    flann3d()
        :
        flann_tree(flann::KDTreeIndexParams(1))
        {};
    
    void set_input(std::map<int, Eigen::Vector3d>& point_to_vector3d_map, std::set<int>& boundary_points_set)
    {
        // compute data storage
        for (int point_id : boundary_points_set)
        {
            flann_data_storage.push_back(point_to_vector3d_map[point_id][0]);
            flann_data_storage.push_back(point_to_vector3d_map[point_id][1]);
            flann_data_storage.push_back(point_to_vector3d_map[point_id][2]);
        }

        // record
        flann_last_id = boundary_points_set.size() - 1;

        // add to flann
        flann_tree.buildIndex(flann::Matrix<float>(flann_data_storage.data(), boundary_points_set.size(), 3));
    }
    
    void set_input(typename pcl::PointCloud<PointT>::Ptr point_cloud)
    {
        // compute data storage
        for (const auto& point : point_cloud->points)
        {
            flann_data_storage.push_back(point.x);
            flann_data_storage.push_back(point.y);
            flann_data_storage.push_back(point.z);
        }

        // record
        flann_last_id = point_cloud->size() - 1;

        // add to flann
        flann_tree.buildIndex(flann::Matrix<float>(flann_data_storage.data(), point_cloud->size(), 3));
    }
    
    void radiusSearch(Eigen::Vector3d searchPoint, std::vector<int>& search_indices, std::vector<float>& search_dists, float radius)
    {
        // convert to vector
        std::vector<float> query_point = {searchPoint[0], searchPoint[1], searchPoint[2]};

        // intialize
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>());
        std::vector<std::vector<float>> list_of_search_dists(1, std::vector<float>());

        // search
        flann_tree.radiusSearch(flann::Matrix<float>(query_point.data(), 1, 3), list_of_search_indices, list_of_search_dists, radius * radius, flann::SearchParams(-1, 0));

        // extract
        search_indices = list_of_search_indices[0];
        search_dists = list_of_search_dists[0];
    }

    void radiusSearch(PointT searchPoint, std::vector<int>& search_indices, std::vector<float>& search_dists, float radius)
    {
        // convert to vector
        std::vector<float> query_point = {searchPoint.x, searchPoint.y, searchPoint.z};

        // intialize
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>());
        std::vector<std::vector<float>> list_of_search_dists(1, std::vector<float>());

        // search
        flann_tree.radiusSearch(flann::Matrix<float>(query_point.data(), 1, 3), list_of_search_indices, list_of_search_dists, radius * radius, flann::SearchParams(-1, 0));

        // extract
        search_indices = list_of_search_indices[0];
        search_dists = list_of_search_dists[0];
    }

    void addPoints(Eigen::Vector3d new_point)
    {
        // convert to vector
        flann_data_storage.push_back(new_point[0]);
        flann_data_storage.push_back(new_point[1]);
        flann_data_storage.push_back(new_point[2]);

        // add
        flann_tree.addPoints(flann::Matrix<float>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3));

        // update id
        flann_last_id++;
    }

    void addPoints(PointT new_point)
    {
        // convert to vector
        flann_data_storage.push_back(new_point.x);
        flann_data_storage.push_back(new_point.y);
        flann_data_storage.push_back(new_point.z);

        // add
        flann_tree.addPoints(flann::Matrix<float>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3));

        // update id
        flann_last_id++;
    }

    void removePoint(int id)
    {
        flann_tree.removePoint(id);
    }

    int flann_last_id;

private:
    std::vector<float> flann_data_storage;

    flann::Index<flann::L2_Simple<float>> flann_tree;
};


// ray plane intersection
Eigen::Vector3d 
ray_plane_intersection(Eigen::Vector3d ray_origin, Eigen::Vector3d ray_direction, Eigen::Vector3d mean, Eigen::Vector3d normal)
{   
    // if parallel, return NaN
    if (normal.dot(ray_direction) == 0)
    {
        return Eigen::Vector3d(NAN, NAN, NAN);
    }

    // compute intersection
    double distance = (mean - ray_origin).dot(normal) / ray_direction.dot(normal);
    Eigen::Vector3d intersection = ray_origin + ray_direction * distance;

    // return
    return intersection;
}


// ray set intersection 
// todo - need to handle the rayOrigin and thisPoint in different coordinate
Eigen::Vector3d ray_set_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayDirection, Eigen::Vector3d thisPoint, std::set<int> point_ids, std::map<int, Eigen::Vector3d> point_to_vector3d_map)
{
    // compute points mean
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (int point_id : point_ids)
    {
        mean += point_to_vector3d_map[point_id];
    }
    mean /= point_ids.size();
    // compute covariance matrix of the points
    Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
    for (int point_id : point_ids)
    {
        Eigen::Vector3d point = point_to_vector3d_map[point_id];
        covariance_matrix += (point - mean) * (point - mean).transpose();
    }
    // decompose to get normal
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance_matrix);
    Eigen::Vector3d normal = eigensolver.eigenvectors().col(0);
    // make normal point towards ray origin
    if (normal.dot(rayDirection) > 0)
    {
        normal = -normal;
    }
    // compute point to plane distance
    Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayDirection, mean, normal);

    return rayPlaneIntersectionPoint;
}


void add_new_point_procedure(
    const Eigen::Vector3d& thisPoint,
    const Eigen::Vector3d& rayOrigin,
    const Eigen::Vector3d& rayDirection,
    std::map<int, Eigen::Vector3d>& point_to_vector3d_map,
    std::map<int, int>& point_to_set_map,
    std::map<int, std::set<int>>& set_to_points_map,
    std::set<int>& boundary_points_set,
    std::map<int, std::array<int, 2>>& edge_to_vertices_map,
    std::map<std::array<int, 2>, int>& edge_to_vertices_map_reverse,
    std::map<int, int>& edge_occurrences_count,
    std::map<int, std::set<int>>& point_to_edge_map,
    std::set<int>& boundary_edge_set,
    std::map<int, std::array<int, 3>>& triangle_to_vertices_map,
    std::map<int, std::set<int>>& set_to_triangles_map,
    int& next_point_id,
    int& next_edge_id,
    int& next_triangle_id,
    int& next_set_id,
    double search_size,
    double distance_threshold
) 
{
    // perform radius search on boundary points

    // if no boundary points, add the point as new set
    if (boundary_points_set.size() == 0)
    {
        // ----------- add point as set ----------- 
        
        // add to point
        point_to_vector3d_map[next_point_id] = thisPoint;

        // add to set
        set_to_points_map[next_set_id].insert(next_point_id);

        // add to boundary points
        boundary_points_set.insert(next_point_id);

        // add to point to set map
        point_to_set_map[next_point_id] = next_set_id;

        // increment set id
        next_set_id ++;

        // increment point id
        next_point_id ++;

        return;
    }

    // setup flann3d
    flann3d<VilensPointT> flann_tree;
    flann_tree.set_input(point_to_vector3d_map, boundary_points_set);

    // perform radius search
    std::vector<int> search_indices;
    std::vector<float> search_dists;
    flann_tree.radiusSearch(thisPoint, search_indices, search_dists, search_size);

    // if no searched results, add point as set
    if (search_indices.size() == 0)
    {
        // ----------- add point as set ----------- 
        
        // add to point
        point_to_vector3d_map[next_point_id] = thisPoint;

        // add to set
        set_to_points_map[next_set_id].insert(next_point_id);

        // add to boundary points
        boundary_points_set.insert(next_point_id);

        // add to point to set map
        point_to_set_map[next_point_id] = next_set_id;

        // increment set id
        next_set_id ++;

        // increment point id
        next_point_id ++;

        return;
    }
    
    // group the searched points by their sets
    std::map<int, std::set<int>> set_to_searched_points_map;
    for (std::size_t j = 0; j < search_indices.size(); j++)
    {
        int point_id = search_indices[j];
        int set_id = point_to_set_map[point_id];
        set_to_searched_points_map[set_id].insert(point_id);
    }

    // for each set of searched points
    for (const auto& pair : set_to_searched_points_map)
    {
        // set id
        int set_id = pair.first;

        // searched points
        std::set<int> searched_point_ids = pair.second;

        // ray set intersection
        Eigen::Vector3d rayPlaneIntersectionPoint = ray_set_intersection(rayOrigin, rayDirection, thisPoint, searched_point_ids, point_to_vector3d_map);
        double distance = (thisPoint - rayPlaneIntersectionPoint).norm();

        // skip if not within plane
        if (distance < -distance_threshold || distance_threshold < distance) // within plane
        {
            continue;
        }

        // update set to points map
        set_to_points_map[set_id].insert(next_point_id);

        // update boundary point
        boundary_points_set.insert(next_point_id);

        // update edge to vertices map
        for (int boundary_point_id : searched_point_ids)
        {
            // form edge
            std::array<int, 2> edge = {boundary_point_id, next_point_id}; // the smaller id is first

            // update map
            edge_to_vertices_map[next_edge_id] = edge;
            edge_to_vertices_map_reverse[edge] = next_edge_id;

            // update edge count
            edge_occurrences_count[next_edge_id] = 0;

            // update point to edge map
            point_to_edge_map[boundary_point_id].insert(next_edge_id);
            point_to_edge_map[next_point_id].insert(next_edge_id);

            // update boundary edge and point (whenever a new edge is added)
            boundary_edge_set.insert(next_edge_id);
            boundary_points_set.insert(boundary_point_id);

            // increment edge id
            next_edge_id ++;
        }

        // find boundary edges between the searched points
        std::set<std::array<int, 2>> existing_boundary_edges;
        // for each searched point
        for (int searched_point_id : searched_point_ids)
        {
            // for each edge of the searched point
            for (int edge_id : point_to_edge_map[searched_point_id])
            {
                // skip if not boundary edge
                if (boundary_edge_set.find(edge_id) == boundary_edge_set.end()) continue;
                
                // add to list if the boundary edge is between searched points
                std::array<int, 2> vertices = edge_to_vertices_map[edge_id]; // [potential bug] need to make sure the vertices order are correct
                if (searched_point_ids.find(vertices[0]) != searched_point_ids.end() && searched_point_ids.find(vertices[1]) != searched_point_ids.end())
                {
                    existing_boundary_edges.insert(vertices);
                }
            }
        }

        // add triangle for each exsiting boundary edge
        for (const std::array<int, 2>& vertices : existing_boundary_edges)
        {
            // form triangle
            std::array<int, 3> triangle = {vertices[0], vertices[1], next_point_id};

            // add to triangle list
            triangle_to_vertices_map[next_triangle_id] = triangle;
            set_to_triangles_map[set_id].insert(next_triangle_id);
            next_triangle_id ++;
            
            // obtain edge id
            int edge1 = edge_to_vertices_map_reverse[{vertices[0], vertices[1]}];
            int edge2 = edge_to_vertices_map_reverse[{vertices[1], next_point_id}];
            int edge3 = edge_to_vertices_map_reverse[{vertices[0], next_point_id}];

            // update edge count (whenever a new triangle is added)
            edge_occurrences_count[edge1] ++;
            edge_occurrences_count[edge2] ++;
            edge_occurrences_count[edge3] ++;

            // update boundary edge and point (whenever a new triangle is added)
            std::array<int, 3> edges = {edge1, edge2, edge3};
            for (int edge : edges)
            {
                if (edge_occurrences_count[edge] <= 1) // boundary edge if count is 0 or 1, since ++ this can't be 0
                {
                    boundary_edge_set.insert(edge);
                    boundary_points_set.insert(edge_to_vertices_map[edge][0]);
                    boundary_points_set.insert(edge_to_vertices_map[edge][1]);
                } 
                else 
                {
                    boundary_edge_set.erase(edge);

                    // if a point is not in any boundary edge, remove from boundary points
                    int vertex1 = edge_to_vertices_map[edge][0];
                    int vertex2 = edge_to_vertices_map[edge][1];
                    std::set<int> vertex1_edges = point_to_edge_map[vertex1];
                    std::set<int> vertex2_edges = point_to_edge_map[vertex2];
                    bool vertex1_boundary = false;
                    bool vertex2_boundary = false;
                    for (int edge_id : vertex1_edges)
                    {
                        if (boundary_edge_set.find(edge_id) != boundary_edge_set.end())
                        {
                            vertex1_boundary = true;
                            break;
                        }
                    }
                    for (int edge_id : vertex2_edges)
                    {
                        if (boundary_edge_set.find(edge_id) != boundary_edge_set.end())
                        {
                            vertex2_boundary = true;
                            break;
                        }
                    }
                    if (!vertex1_boundary) boundary_points_set.erase(vertex1);
                    if (!vertex2_boundary) boundary_points_set.erase(vertex2);
                }
            }

        }

    }
}

using InputPointT = VilensPointT;
int main()
{
    /* 
    - when a new point comes in, check intersected triangle
    - if have intersection
        - from the intersected triangle, retrieve the set
        - from the set, check relation of the point with the set
        - cases
            - within - add the point to the set and to the triangle
            - in front of - move to "no intersection" process
            - behind
                - remove the triangle from the set
                - recompute boundary edge and edge points
                - for each point within the triangle, re-add them to the map
    - if no intersection with any triangles
        - perform radius search on edge points
            - from the edge points identify the set / planes
            - if new point does not match found planes
                - add the new point as new set
            - if new point match found planes
                - form edge to the found edge points
                - form triangle between the new point and any two edge points that have a boundary edge between them
                - the triangle is added to the set and the boundary edge and points are recomputed
    */

    // input data
    std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
    std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";
    DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);
    int i1 = 0;
    typename pcl::PointCloud<InputPointT>::Ptr new_cloud = data_loader.get_cloud(i1);
    Eigen::Affine3d new_pose = data_loader.get_pose(i1);

    // storage data

        // point
    int next_point_id = 0;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, int> point_to_set_map;
    std::map<int, std::set<int>> point_to_edge_map;

        // triangle
    int next_triangle_id = 0;
    std::map<int, std::array<int, 3>> triangle_to_vertices_map;
    std::map<int, int> triangle_to_set_map;
    std::map<int, std::set<int>> triangle_to_points_map;
    
        // set
    int next_set_id = 0;
    std::map<int, std::set<int>> set_to_points_map; // each set contains id to points
    std::map<int, std::set<int>> set_to_triangles_map; // each set contains id to triangles
    
        // edge
    int next_edge_id = 0;
    std::map<int, std::array<int, 2>> edge_to_vertices_map;
    std::map<std::array<int, 2>, int> edge_to_vertices_map_reverse;
    std::map<int, int> edge_occurrences_count; // number of triangles that share the edge

        // boundary
    std::set<int> boundary_points_set;
    std::set<int> boundary_edge_set;
    



    // settings
    double distance_threshold = 0.05;
    float search_size = 0.1;

    // for each point in new cloud
    for (std::size_t i = 0; i < new_cloud->size(); i++)
    {
        // precompute
        std::vector<int> triangle_map_keys;
        for (const auto& pair : triangle_to_vertices_map)
        {
            triangle_map_keys.push_back(pair.first);
        }
        
        // Build the BVH
        auto bvhRoot = buildBVH(point_to_vector3d_map, triangle_to_vertices_map, triangle_map_keys, 0, triangle_map_keys.size());

        // Define this point
        Eigen::Vector3d thisPoint = new_cloud->points[i].getVector3fMap().cast<double>();

        // Define the ray
        Eigen::Vector3d rayOrigin(0, 0, 0);
        Eigen::Vector3d rayDirection = thisPoint.normalized();


        // Perform intersection test
        std::vector<int> intersectedTriangleIdList;
        std::vector<Eigen::Vector3d> intersectionPointList;
        bool intersect = intersectBVH(bvhRoot, rayOrigin, rayDirection, point_to_vector3d_map, triangle_to_vertices_map, intersectedTriangleIdList, intersectionPointList);

        // if intersect
        if (intersect)
        {
            // for each intersected triangle
            for (std::size_t j = 0; j < intersectedTriangleIdList.size(); j++)
            {
                // get the triangle id
                int triangle_id = intersectedTriangleIdList[j];

                // get the corresponding set
                int set_id = triangle_to_set_map[triangle_id];

                // get the points in the set
                std::set<int> point_ids = set_to_points_map[set_id];

                // recompute the plane normal etc for now 
                // [todo] store precomputed results

                
                // ray set intersection
                Eigen::Vector3d rayPlaneIntersectionPoint = ray_set_intersection(rayOrigin, rayDirection, thisPoint, point_ids, point_to_vector3d_map);
                double distance = (thisPoint - rayPlaneIntersectionPoint).norm();


                // check position of point relative to plane
                if (-distance_threshold < distance && distance < distance_threshold) // within plane
                {
                    // add the point to the set
                    set_to_points_map[set_id].insert(i);
                    // add the point to the triangle
                    triangle_to_points_map[triangle_id].insert(i);
                }
                else if (distance > distance_threshold) // in front of plane
                {
                    add_new_point_procedure(
                        thisPoint,
                        rayOrigin,
                        rayDirection,
                        point_to_vector3d_map,
                        point_to_set_map,
                        set_to_points_map,
                        boundary_points_set,
                        edge_to_vertices_map,
                        edge_to_vertices_map_reverse,
                        edge_occurrences_count,
                        point_to_edge_map,
                        boundary_edge_set,
                        triangle_to_vertices_map,
                        set_to_triangles_map,
                        next_point_id,
                        next_edge_id,
                        next_triangle_id,
                        next_set_id,
                        search_size,
                        distance_threshold
                    );
                }
                else // behind plane
                {
                    // --------------------- remove triangle ---------------------------
                    
                        // remove from set to triangle map
                    set_to_triangles_map[set_id].erase(triangle_id);
                    
                        // update edge count
                    std::array<int, 3> vertices = triangle_to_vertices_map[triangle_id];
                    int edge1 = edge_to_vertices_map_reverse[{vertices[0], vertices[1]}];
                    int edge2 = edge_to_vertices_map_reverse[{vertices[1], vertices[2]}];
                    int edge3 = edge_to_vertices_map_reverse[{vertices[0], vertices[2]}];
                    edge_occurrences_count[edge1] --;
                    edge_occurrences_count[edge2] --;
                    edge_occurrences_count[edge3] --;

                        // update boundary edge and point
                    std::array<int, 3> edges = {edge1, edge2, edge3};
                    for (int edge : edges)
                    {
                        if (edge_occurrences_count[edge] <= 1) // boundary edge if count is 0 or 1
                        {
                            boundary_edge_set.insert(edge);
                            boundary_points_set.insert(edge_to_vertices_map[edge][0]);
                            boundary_points_set.insert(edge_to_vertices_map[edge][1]);
                        } 
                        else 
                        {
                            // no edge in the first place, no need to remove from boundary set
                        }
                    }

                        // remove from the triangle to vertices map
                    triangle_to_vertices_map.erase(triangle_id);

                        // remove from the triangle to set map
                    triangle_to_set_map.erase(triangle_id);

                        // remove from the triangle to points map
                    std::set<int> points = triangle_to_points_map[triangle_id];
                    triangle_to_points_map.erase(triangle_id);
                    
                        // for each point previously within the triangle, re-add them to the map
                    for (int point_id : points)
                    {
                        add_new_point_procedure(
                            thisPoint,
                            rayOrigin,
                            rayDirection,
                            point_to_vector3d_map,
                            point_to_set_map,
                            set_to_points_map,
                            boundary_points_set,
                            edge_to_vertices_map,
                            edge_to_vertices_map_reverse,
                            edge_occurrences_count,
                            point_to_edge_map,
                            boundary_edge_set,
                            triangle_to_vertices_map,
                            set_to_triangles_map,
                            next_point_id,
                            next_edge_id,
                            next_triangle_id,
                            next_set_id,
                            search_size,
                            distance_threshold
                        );
                    }
                    // --------------------- remove triangle ---------------------------
                }
            }
        }
        else // no intersection to known triangles, perform add point to map procedure
        {
            add_new_point_procedure(
                thisPoint,
                rayOrigin,
                rayDirection,
                point_to_vector3d_map,
                point_to_set_map,
                set_to_points_map,
                boundary_points_set,
                edge_to_vertices_map,
                edge_to_vertices_map_reverse,
                edge_occurrences_count,
                point_to_edge_map,
                boundary_edge_set,
                triangle_to_vertices_map,
                set_to_triangles_map,
                next_point_id,
                next_edge_id,
                next_triangle_id,
                next_set_id,
                search_size,
                distance_threshold
            );
        }
    }


    // display edge mesh
    edge_to_vertices_map;

    // make mesh
    pcl::PolygonMesh mesh;

    // add point to cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr mesh_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (const auto& pair : point_to_vector3d_map)
    {
        pcl::PointXYZ point;
        point.x = pair.second[0];
        point.y = pair.second[1];
        point.z = pair.second[2];
        mesh_cloud->push_back(point);
    }
    // convert to PCLPointCloud2
    pcl::toPCLPointCloud2(*mesh_cloud, mesh.cloud); 

    // add edges to mesh
    for (const auto& pair : edge_to_vertices_map)
    {
        pcl::Vertices edge;
        edge.vertices.push_back(pair.second[0]);
        edge.vertices.push_back(pair.second[1]);
        mesh.polygons.push_back(edge);
    }
    

    // add viewer
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1);

    // add mesh
    viewer->addPolylineFromPolygonMesh(mesh, "polyline");

    // spin
    viewer->spin();


   return 0;
}

// int main() {
//     // Define your mesh vertices and triangles
//     // vertices id to vector3d map
//     std::map<int, Eigen::Vector3d> vertex_to_vector3d_map;
//     vertex_to_vector3d_map[0] = Eigen::Vector3d(1, 1, 1);
//     vertex_to_vector3d_map[1] = Eigen::Vector3d(-1, 1, 1);
//     vertex_to_vector3d_map[2] = Eigen::Vector3d(-1, -1, 1);
//     vertex_to_vector3d_map[3] = Eigen::Vector3d(1, -1, 1);
//     vertex_to_vector3d_map[4] = Eigen::Vector3d(1, 1, -1);
//     vertex_to_vector3d_map[5] = Eigen::Vector3d(-1, 1, -1);
//     vertex_to_vector3d_map[6] = Eigen::Vector3d(-1, -1, -1);
//     vertex_to_vector3d_map[7] = Eigen::Vector3d(1, -1, -1);

//     // triangle id to indices map
//     std::map<int, std::array<int, 3>> triangle_to_indices_map;
//     triangle_to_indices_map[0] = {0, 1, 2};
//     triangle_to_indices_map[1] = {0, 3, 2};
//     triangle_to_indices_map[2] = {0, 3, 7};
//     triangle_to_indices_map[3] = {0, 4, 7};
//     triangle_to_indices_map[4] = {0, 1, 5};
//     triangle_to_indices_map[5] = {0, 4, 5};
//     triangle_to_indices_map[6] = {6, 5, 1};
//     triangle_to_indices_map[7] = {6, 2, 1};
//     triangle_to_indices_map[8] = {6, 5, 4};
//     triangle_to_indices_map[9] = {6, 7, 4};
//     triangle_to_indices_map[10] = {6, 2, 3};
//     triangle_to_indices_map[11] = {6, 7, 3};

//     // triangle id list
//     std::vector<int> triangle_id_list;
//     for (const auto& pair : triangle_to_indices_map)
//     {
//         triangle_id_list.push_back(pair.first);
//     }

//     // Build the BVH
//     auto bvhRoot = buildBVH(vertex_to_vector3d_map, triangle_to_indices_map, triangle_id_list, 0, triangle_id_list.size());

//     // Define the ray
//     Eigen::Vector3d rayOrigin(0, 0, 0);
//     Eigen::Vector3d rayDirection(0, 0, 1);
//     std::vector<int> intersectedTriangleIdList;
//     std::vector<Eigen::Vector3d> intersectionPointList;

//     // Perform intersection test
//     bool intersect = intersectBVH(bvhRoot, rayOrigin, rayDirection, vertex_to_vector3d_map, triangle_to_indices_map, intersectedTriangleIdList, intersectionPointList);

//     if (intersect) {
//         for (std::size_t i = 0; i < intersectedTriangleIdList.size(); ++i) {
//             std::cout << "Intersection with triangle id: " << intersectedTriangleIdList[i] << " at point: " << intersectionPointList[i].transpose() << std::endl;
//         }
//     }


//     // viewer
//     pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
//     viewer->initCameraParameters();
//     viewer->addCoordinateSystem(1);
    
//     // make mesh
//     pcl::PolygonMesh mesh;
//     // add points to mesh.cloud
//     pcl::PointCloud<pcl::PointXYZ>::Ptr new_cloud(new pcl::PointCloud<pcl::PointXYZ>);
//     for (const auto& pair : vertex_to_vector3d_map)
//     {
//         pcl::PointXYZ point;
//         point.x = pair.second[0];
//         point.y = pair.second[1];
//         point.z = pair.second[2];
//         new_cloud->push_back(point);
//     }


//     // add points
//     pcl::toPCLPointCloud2(*new_cloud, mesh.cloud); 

//     // add triangle to mesh
//     for (const auto& pair : triangle_to_indices_map)
//     {
//         pcl::Vertices triangle;
//         triangle.vertices.push_back(pair.second[0]);
//         triangle.vertices.push_back(pair.second[1]);
//         triangle.vertices.push_back(pair.second[2]);
//         mesh.polygons.push_back(triangle);
//     }
//     // add mesh
//     viewer->addPolylineFromPolygonMesh(mesh, "polyline");

//     // add ray
//     pcl::PointXYZ ray_origin;
//     ray_origin.x = rayOrigin[0];
//     ray_origin.y = rayOrigin[1];
//     ray_origin.z = rayOrigin[2];
//     pcl::PointXYZ ray_end;
//     double length = 2;
//     ray_end.x = rayOrigin[0] + rayDirection[0] * length;
//     ray_end.y = rayOrigin[1] + rayDirection[1] * length;
//     ray_end.z = rayOrigin[2] + rayDirection[2] * length;
//     viewer->addLine(ray_origin, ray_end, 1, 1, 0, "ray");

//     // add ray origin
//     viewer->addSphere(ray_origin, 0.01, 1, 0, 0, "ray_origin");

//     // add intersection
//     if (intersect)
//     {
//         pcl::PointXYZ intersection;
//         intersection.x = intersectionPointList[0][0];
//         intersection.y = intersectionPointList[0][1];
//         intersection.z = intersectionPointList[0][2];
//         viewer->addSphere(intersection, 0.01, 0, 1, 0, "intersection");
//     }

//     // spin
//     viewer->spin();

//     return 0;
// }







// int main()
// {
//     // test ray triangle intersection
//     Eigen::Vector3d orig(0.5, 0.25, 0);
//     Eigen::Vector3d dir(0, 0, 1);
//     Eigen::Vector3d v0(0, 0, 1);
//     Eigen::Vector3d v1(1, 0, 1);
//     Eigen::Vector3d v2(0, 1, 1);
//     Eigen::Vector3d outIntersection;
//     bool intersect = rayTriangleIntersect(orig, dir, v0, v1, v2, outIntersection);
//     std::cout << "intersect: " << intersect << std::endl;
//     std::cout << "outIntersection: " << outIntersection << std::endl;
// }

// using InputPointT = VilensPointT;
// int main()
// {
//     std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
//     std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";

//     DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);

//     // ------------------------------ data
//     // new cloud
//     int i1 = 0;
//     typename pcl::PointCloud<InputPointT>::Ptr new_cloud = data_loader.get_cloud(i1);
//     Eigen::Affine3d new_pose = data_loader.get_pose(i1);


//     // old cloud
//     typename pcl::PointCloud<InputPointT>::Ptr old_cloud(new pcl::PointCloud<InputPointT>);
//     old_cloud->push_back(new_cloud->points[0]);

//     // flann3d
//     flann3d<InputPointT> flann_tree;
//     flann_tree.set_input(old_cloud);


//     // for each point in new cloud, find the triangles in old cloud that contains it
//     // each triangle is in a set
    
    
    
//     // // 1. for each point in new cloud, find the points within radius r in old cloud, forms potential edges to those points, add the edge to the mesh object
//     // // 2. check if the found point is in a set
    
//     // // initialize
//     // std::map<int, std::vector<int>> edge_map; // for each ith point, it is connected to the points in the vector<int>

//     // for (std::size_t i = 1; i < new_cloud->size(); i++)
//     // {
//     //     std::cout << "i: " << i << std::endl;

//     //     // current point
//     //     InputPointT current_point = new_cloud->points[i];
        
//     //     // search point
//     //     std::vector<int> search_indices;
//     //     std::vector<float> search_dists;
//     //     float search_radius = 0.1;
//     //     flann_tree.radiusSearch(current_point, search_indices, search_dists, search_radius);
        
//     //     // add searched edges to map
//     //     for (std::size_t j = 0; j < search_indices.size(); j++)
//     //     {
//     //         edge_map[i].push_back(search_indices[j]);
//     //     }

//     //     // add current point to flann
//     //     flann_tree.addPoints(current_point);
//     // }

//     // // make mesh
//     // pcl::PolygonMesh mesh;
//     // // add points
//     // pcl::toPCLPointCloud2(*new_cloud, mesh.cloud); 
//     // // add edges
//     // for (const auto& pair : edge_map) // source and targets pair
//     // { 
//     //     // source
//     //     int source_i = pair.first;

//     //     // targets
//     //     for (const auto& target_j : pair.second)
//     //     {
//     //         pcl::Vertices edge;
//     //         edge.vertices.push_back(source_i);
//     //         edge.vertices.push_back(target_j);
//     //         mesh.polygons.push_back(edge);
//     //     }
//     // }




//     // // transform old cloud to global coordinate
//     // typename pcl::PointCloud<InputPointT>::Ptr new_cloud_local = transform_cloud_to_global<InputPointT>(new_cloud, new_pose);


//     // viewer
//     pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
//     viewer->initCameraParameters();
//     viewer->addCoordinateSystem(1);

//     // // add point cloud to viewer
//     // pcl::visualization::PointCloudColorHandlerCustom<InputPointT> new_cloud_color_handler(new_cloud_local, 0, 255, 0);
//     // viewer->addPointCloud<InputPointT> (new_cloud_local, new_cloud_color_handler, "new_cloud");
//     // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "new_cloud");

//     // // add mesh
//     // viewer->addPolylineFromPolygonMesh(mesh, "polyline");

//     // spin
//     viewer->spin();

//     return 0;
// }



// using InputPointT = VilensPointT;
// int main()
// {
//     // dataloader
//     std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
//     std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";
//     DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);



//     // load pointcloud and pose
//     int i1 = 100;
//     typename pcl::PointCloud<InputPointT>::Ptr old_cloud = data_loader.get_cloud(i1);
//     Eigen::Affine3d old_pose = data_loader.get_pose(i1);

//     // compute direction and origin
//     std::vector<Eigen::Vector3f> old_cloud_direction = compute_point_directions<InputPointT>(old_cloud);
//     Eigen::Vector3f origin = Eigen::Vector3f::Zero();

//     // transform cloud, direction and origin
//     typename pcl::PointCloud<InputPointT>::Ptr old_cloud_transformed = transform_cloud_to_global<InputPointT> (old_cloud, old_pose);
//     std::vector<Eigen::Vector3f> old_cloud_direction_transformed = transform_direction_to_global(old_cloud_direction, old_pose);
//     Eigen::Vector3f origin_transformed = old_pose.cast<float>() * origin;

//     // initialize viewer
//     pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
//     viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
//     viewer->initCameraParameters();

//     // convert to viewer cloud
//     typename pcl::PointCloud<InputPointT>::Ptr cloud_to_use = old_cloud_transformed;
//     std::vector<Eigen::Vector3f> direction_to_use = old_cloud_direction_transformed;
//     Eigen::Vector3f origin_to_use = origin_transformed;
//     // typename pcl::PointCloud<InputPointT>::Ptr cloud_to_use = old_cloud;
//     // std::vector<Eigen::Vector3f> direction_to_use = old_cloud_direction;
//     // Eigen::Vector3f origin_to_use = origin;

//     pcl::PointCloud<pcl::PointXYZINormal>::Ptr viewer_cloud(new pcl::PointCloud<pcl::PointXYZINormal>);
//     viewer_cloud->resize(cloud_to_use->size());
//     for (std::size_t i = 0; i < cloud_to_use->size(); i++)
//     {
//         pcl::PointXYZINormal viewer_point;
//         viewer_point.x = cloud_to_use->points[i].x;
//         viewer_point.y = cloud_to_use->points[i].y;
//         viewer_point.z = cloud_to_use->points[i].z;
//         viewer_point.intensity = 0;
//         viewer_point.normal_x = -direction_to_use[i][0];
//         viewer_point.normal_y = -direction_to_use[i][1];
//         viewer_point.normal_z = -direction_to_use[i][2];
//         viewer_cloud->points[i] = viewer_point;
//     }

//     // add to viewer
//     pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZINormal> viewer_cloud_color(viewer_cloud, 255, 0, 0);
//     viewer->addPointCloud<pcl::PointXYZINormal> (viewer_cloud, viewer_cloud_color, "pointcloud");
//     viewer->addPointCloudNormals<pcl::PointXYZINormal> (viewer_cloud, 1, 0.05, "normals");
//     viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "pointcloud");
//     viewer->addCoordinateSystem(0.5);
//     viewer->addCoordinateSystem(1, origin_to_use[0], origin_to_use[1], origin_to_use[2], "pose");

//     // spin
//     viewer->spin();

//     return 0;
// }