#include "eye_patch/DataLoader.hpp"

// #include "eye_patch/Algorithm.hpp"
#include <pcl/common/transforms.h>
#include <pcl/kdtree/kdtree_flann.h>
#include "eye_patch/BVH.hpp"

// #include "eye_patch/Visualization.hpp"
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

// #include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"

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
        flann_tree.buildIndex(flann::Matrix<double>(flann_data_storage.data(), boundary_points_set.size(), 3));
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
        flann_tree.buildIndex(flann::Matrix<double>(flann_data_storage.data(), point_cloud->size(), 3));
    }
    
    void radiusSearch(Eigen::Vector3d searchPoint, std::vector<int>& search_indices, std::vector<double>& search_dists, double radius)
    {
        // convert to vector
        std::vector<double> query_point = {searchPoint[0], searchPoint[1], searchPoint[2]}; 

        // intialize
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>());
        std::vector<std::vector<double>> list_of_search_dists(1, std::vector<double>());

        // search
        flann_tree.radiusSearch(flann::Matrix<double>(query_point.data(), 1, 3), list_of_search_indices, list_of_search_dists, radius * radius, flann::SearchParams(-1, 0));

        // extract
        search_indices = list_of_search_indices[0];
        search_dists = list_of_search_dists[0];
    }

    void radiusSearch(PointT searchPoint, std::vector<int>& search_indices, std::vector<double>& search_dists, double radius)
    {
        // convert to vector
        std::vector<double> query_point = {searchPoint.x, searchPoint.y, searchPoint.z};

        // intialize
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>());
        std::vector<std::vector<double>> list_of_search_dists(1, std::vector<double>());

        // search
        flann_tree.radiusSearch(flann::Matrix<double>(query_point.data(), 1, 3), list_of_search_indices, list_of_search_dists, radius * radius, flann::SearchParams(-1, 0));

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
        flann_tree.addPoints(flann::Matrix<double>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3));

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
        flann_tree.addPoints(flann::Matrix<double>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3));

        // update id
        flann_last_id++;
    }

    void removePoint(int id)
    {
        flann_tree.removePoint(id);
    }

    int flann_last_id;

private:
    std::vector<double> flann_data_storage;

    flann::Index<flann::L2_Simple<double>> flann_tree;
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
Eigen::Vector3d ray_set_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayDirection, std::set<int> point_ids, std::map<int, Eigen::Vector3d> point_to_vector3d_map)
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






// application class
template <typename PointT>
class Application
{
public:
    Application() 
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
        DataLoader<VilensPointT> data_loader(pcd_file_folder, pose_file_path);
        int i1 = 0;
        new_cloud = data_loader.get_cloud(i1);


        // loop
        i = 0;
        
    }

    // want boundary edge and points to prevent overlapping triangles
    // to add a new point to mesh
    // 1. extract mesh boundary points
    // 2. perform delaunay triangulation using the boundary points and the new point
    // 3. add the edge and triangles connected to the new point to the original mesh

    void add_point_to_set(int pointID, int setID)
    {
        // add to set
        set_to_points_map[setID].insert(pointID);
        point_to_set_map[pointID] = setID;

        // // add to boundary points
        // boundary_points_set.insert(newPointID);
    }

    void add_point_coordinate(int newPointID, Eigen::Vector3d thisPoint)
    {
        // add point list
        point_list.push_back(newPointID);

        // add eigen
        point_to_vector3d_map[newPointID] = thisPoint;
    }

    void add_edge(int newEdgeID, int smaller_id, int larger_id)
    {
        // add edge list
        edge_list.push_back(newEdgeID);

        // edge
        std::array<int, 2> edge = {smaller_id, larger_id}; // the smaller id is first

        // get smaller id set
        int smaller_set_id = point_to_set_map[smaller_id];
        // add larger id to smaller id set
        set_to_points_map[smaller_set_id].insert(larger_id);
        point_to_set_map[larger_id] = smaller_set_id;

        // add point to edge to vertices map
        edge_to_point_map[newEdgeID] = edge;
        edge_to_point_map_reverse[edge] = newEdgeID;

        // // initialize edge count (only increment when a new triangle is added)
        // edge_occurrences_count[newEdgeID] = 0;
        // // new edge are always boundary edge, if only search between boundary points
        // boundary_edge_set.insert(next_edge_id);
        // boundary_points_set.insert(edge[0]);
        // boundary_points_set.insert(edge[1]);

        // update point to edge map
        point_to_edge_map[smaller_id].insert(newEdgeID);
        point_to_edge_map[larger_id].insert(newEdgeID);
    }

    void add_triangle(int newTriangleID, int newSetID, std::array<int, 3> vertices)
    {
        // add triangle to vertices map
        triangle_to_vertices_map[newTriangleID] = vertices;

        // add triangle to set map
        set_to_triangles_map[newSetID].insert(newTriangleID);
    }

    int getNewPointID()
    {
        int newPointID = next_point_id;
        next_point_id ++;

        return newPointID;
    }

    int getNewSetID()
    {
        int newSetID = next_set_id;
        next_set_id ++;

        return newSetID;
    }

    int getNewEdgeID()
    {
        int newEdgeID = next_edge_id;
        next_edge_id ++;

        return newEdgeID;
    }
    int getNewTriangleID()
    {
        int newTriangleID = next_triangle_id;
        next_triangle_id ++;

        return newTriangleID;
    }

    std::map<int, std::set<int>> group_points_by_set(std::vector<int> point_indices)
    {
        std::map<int, std::set<int>> points_grouped_by_set;
        for (int point_id : point_indices)
        {
            points_grouped_by_set[point_to_set_map[point_id]].insert(point_id);
        }

        return points_grouped_by_set;
    }

    std::vector<int> compute_boundary_edge_list()
    {
        // intialize edge occurrences count
        std::map<int, int> edge_occurrences_count;
        for (int edge_id : edge_list)
        {
            edge_occurrences_count[edge_id] = 0;
        }

        // iterate through all triangles to count edge
        for (const auto& pair : triangle_to_vertices_map)
        {
            std::array<int, 3> vertices = pair.second;
            for (int i = 0; i < 3; i++)
            {
                int smaller_id = std::min(vertices[i], vertices[(i + 1) % 3]);
                int larger_id = std::max(vertices[i], vertices[(i + 1) % 3]);
                std::array<int, 2> edge = {smaller_id, larger_id};
                edge_occurrences_count[edge_to_point_map_reverse[edge]] ++;
            }
        }
        
        // identify boundary edges (edges that are shared by only 0 or 1 triangle)
        std::set<int> boundary_edge_set;
        for (const auto& pair : edge_occurrences_count)
        {
            if (pair.second <= 1)
            {
                boundary_edge_set.insert(pair.first);
            }
        }

        // convert to list
        std::vector<int> boundary_edge_list;
        for (int edge_id : boundary_edge_set)
        {
            boundary_edge_list.push_back(edge_id);
        }

        // return
        return boundary_edge_list;
    }

    std::set<int> compute_boundary_edge_set()
    {
        // get boundary edge list
        std::vector<int> boundary_edge_list = compute_boundary_edge_list();

        // return
        return std::set<int>(boundary_edge_list.begin(), boundary_edge_list.end());
    }

    std::vector<int> compute_boundary_point_list()
    {
        // get boundary edge list
        std::vector<int> boundary_edge_list = compute_boundary_edge_list();

        // identify boundary points
        // add points from boundary edges
        std::set<int> boundary_points_set;
        for (int edge_id : boundary_edge_list)
        {
            std::array<int, 2> vertices = edge_to_point_map[edge_id];
            boundary_points_set.insert(vertices[0]);
            boundary_points_set.insert(vertices[1]);
        }
        // add points that do not have edge 
        for (int point_id : point_list)
        {
            if (point_to_edge_map.find(point_id) == point_to_edge_map.end())
            {
                boundary_points_set.insert(point_id);
            }
        }

        // convert to list
        std::vector<int> boundary_points_list;
        for (int point_id : boundary_points_set)
        {
            boundary_points_list.push_back(point_id);
        }

        // return 
        return boundary_points_list;
    }

    // compute boundary point set
    std::set<int> compute_boundary_point_set()
    {
        // get boundary point list
        std::vector<int> boundary_point_list = compute_boundary_point_list();

        // return
        return std::set<int>(boundary_point_list.begin(), boundary_point_list.end());
    }

    void fit_plane_to_points(std::set<int> point_ids, Eigen::Vector3d& mean, Eigen::Vector3d& normal)
    {
        // compute mean
        mean = Eigen::Vector3d::Zero();
        for (int point_id : point_ids)
        {
            mean += point_to_vector3d_map[point_id];
        }
        mean /= point_ids.size();

        // compute covariance matrix
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (int point_id : point_ids)
        {
            Eigen::Vector3d point = point_to_vector3d_map[point_id];
            covariance_matrix += (point - mean) * (point - mean).transpose();
        }

        // decompose to get normal
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance_matrix);
        normal = eigensolver.eigenvectors().col(0);
    }



    // Function to compute the mean of the points
    Eigen::Vector3d computeMean(const std::vector<Eigen::Vector3d>& points) 
    {
        // compute
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (const auto& point : points) 
        {
            mean += point;
        }
        mean /= points.size();

        // return
        return mean;
    }

    // Function to center the points
    std::vector<Eigen::Vector3d> centerPoints(const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& mean) 
    {
        // compute
        std::vector<Eigen::Vector3d> centered_points;
        for (const auto& point : points) 
        {
            centered_points.push_back(point - mean);
        }

        // return
        return centered_points;
    }

    // Function to compute the covariance matrix
    Eigen::Matrix3d computeCovarianceMatrix(const std::vector<Eigen::Vector3d>& centered_points) 
    {
        // comptute
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (const auto& point : centered_points) 
        {
            covariance_matrix += point * point.transpose();
        }
        covariance_matrix /= centered_points.size();

        // return
        return covariance_matrix;
    }

    // Function to project points onto the plane
    std::vector<Eigen::Vector2d> projectPoints(const std::vector<Eigen::Vector3d>& centered_points, const Eigen::Matrix3d& eigenvectors) 
    {
        // compute
        Eigen::Matrix<double, 3, 2> plane_basis = eigenvectors.rightCols<2>();  // Use the last two eigenvectors
        std::vector<Eigen::Vector2d> projections;
        for (const auto& point : centered_points) 
        {
            projections.push_back((plane_basis.transpose() * point).head<2>());
        }

        // return
        return projections;
    }
    
    Eigen::Matrix3d  get_set_eigenvectors(int setID)
    {
        // mean
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (int point_id : set_to_points_map[setID])
        {
            mean += point_to_vector3d_map[point_id];
        }
        mean /= set_to_points_map[setID].size();

        // covariance matrix
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (int point_id : set_to_points_map[setID])
        {
            Eigen::Vector3d point = point_to_vector3d_map[point_id];
            covariance_matrix += (point - mean) * (point - mean).transpose();
        }
        covariance_matrix /= set_to_points_map[setID].size();

        // eigen decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(covariance_matrix);
        Eigen::Matrix3d eigenvectors = eigen_solver.eigenvectors();
    
        // return
        return eigenvectors;
    }

    // Function to find the orientation of the ordered triplet (p, q, r).
    int orientation(const Eigen::Vector2d &p, const Eigen::Vector2d &q, const Eigen::Vector2d &r) {
        double val = (q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y());
        if (val == 0) return 0;           // collinear
        return (val > 0) ? 1 : 2;         // clock or counterclockwise
    }

    // Function to check if point q lies on segment pr excluding endpoints.
    bool onSegment(const Eigen::Vector2d &p, const Eigen::Vector2d &q, const Eigen::Vector2d &r) {
        if (q.x() < std::max(p.x(), r.x()) && q.x() > std::min(p.x(), r.x()) &&
            q.y() < std::max(p.y(), r.y()) && q.y() > std::min(p.y(), r.y()))
            return true;
        return false;
    }

    // Function to check if two segments (p1, q1) and (p2, q2) intersect.
    bool doIntersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &q1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q2) {
        // Find the four orientations needed for the general and special cases
        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);

        // General case
        if (o1 != o2 && o3 != o4)
            return true;

        // Special cases
        // p1, q1 and p2 are collinear and p2 lies on segment p1q1
        if (o1 == 0 && onSegment(p1, p2, q1)) return false;

        // p1, q1 and q2 are collinear and q2 lies on segment p1q1
        if (o2 == 0 && onSegment(p1, q2, q1)) return false;

        // p2, q2 and p1 are collinear and p1 lies on segment p2q2
        if (o3 == 0 && onSegment(p2, p1, q2)) return false;

        // p2, q2 and q1 are collinear and q1 lies on segment p2q2
        if (o4 == 0 && onSegment(p2, q1, q2)) return false;

        return false; // Doesn't fall in any of the above cases
    }


    // point_in_triangle - generated by copilot
    bool point_in_triangle(Eigen::Vector2d point, Eigen::Vector2d v0, Eigen::Vector2d v1, Eigen::Vector2d v2) {
        Eigen::Vector2d v0v1 = v1 - v0;
        Eigen::Vector2d v0v2 = v2 - v0;
        Eigen::Vector2d v0p = point - v0;

        double dot00 = v0v2.dot(v0v2);
        double dot01 = v0v2.dot(v0v1);
        double dot02 = v0v2.dot(v0p);
        double dot11 = v0v1.dot(v0v1);
        double dot12 = v0v1.dot(v0p);

        double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        // Return true if point is in triangle
        return (u >= 0) && (v >= 0) && (u + v <= 1);
    }

    std::set<int> extract_existing_edge_between_points(std::set<int> candidate_point_set, std::set<int> candidate_edge_set)
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

    bool edge_edges_intersection(std::array<int, 2> edgeA, std::set<std::array<int, 2>> edgeB_set, std::map<int, Eigen::Vector2d> points_to_vector2d_map)
    {
        for (const auto& edgeB : edgeB_set)
        {
            // if intersected at end points, don't count
            if (edgeA[0] == edgeB[0] || edgeA[0] == edgeB[1] || edgeA[1] == edgeB[0] || edgeA[1] == edgeB[1]) continue;

            // intersection check
            Eigen::Vector2d pointA0 = points_to_vector2d_map[edgeA[0]];
            Eigen::Vector2d pointA1 = points_to_vector2d_map[edgeA[1]];
            Eigen::Vector2d pointB0 = points_to_vector2d_map[edgeB[0]];
            Eigen::Vector2d pointB1 = points_to_vector2d_map[edgeB[1]];
            if (doIntersect(pointA0, pointA1, pointB0, pointB1)) return true;
        }

        return false;
    }

    bool triangle_contains_point(std::array<int, 3> triangle, std::set<int> point_set, std::map<int, Eigen::Vector2d> points_to_vector2d_map)
    {
        for (int point_id : point_set)
        {
            // if contained at end points, don't count
            if (point_id == triangle[0] || point_id == triangle[1] || point_id == triangle[2]) continue;

            // containment check
            Eigen::Vector2d point = points_to_vector2d_map[point_id];
            Eigen::Vector2d a = points_to_vector2d_map[triangle[0]];
            Eigen::Vector2d b = points_to_vector2d_map[triangle[1]];
            Eigen::Vector2d c = points_to_vector2d_map[triangle[2]];
            if (point_in_triangle(point, a, b, c)) return true;
        }
    
        return false;
    }
    


    // convert to edge vertices set
    std::set<std::array<int, 2>> convert_to_edge_vertices_set(std::set<int> edge_set)
    {
        // initialize
        std::set<std::array<int, 2>> edge_vertices_set;

        // process
        for (int edge_id : edge_set)
        {
            edge_vertices_set.insert(edge_to_point_map[edge_id]);
        }

        // return
        return edge_vertices_set;
    }

    // project points to plane
    std::map<int, Eigen::Vector2d> project_points_to_plane(std::set<int> point_set, Eigen::Matrix<double, 3, 2> projection_matrix)
    {
        // initialize
        std::map<int, Eigen::Vector2d> points_to_vector2d_map;

        // process
        for (int point_id : point_set)
        {
            Eigen::Vector3d point = point_to_vector3d_map[point_id];
            Eigen::Vector2d projection = (projection_matrix.transpose() * point).head<2>();
            points_to_vector2d_map[point_id] = projection;
        }

        // return
        return points_to_vector2d_map;
    }

    // extract set from set
    std::set<int> intersection_of_sets(std::set<int> setA, std::set<int> setB)
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
    

    void step()
    {
        // for each point in new cloud
        if (i >= new_cloud->size()) return;
        std::cout << "Processing point " << i << " / " << new_cloud->size() << std::endl;
        Eigen::Vector3d thisPointVEC = new_cloud->points[i].getVector3fMap().cast<double>();
        i ++;
        
        // get new point id
        int newPointID = getNewPointID();
        
        // setup kdtreeflann
        pcl::PointCloud<pcl::PointXYZ>::Ptr kdcloud = point_to_vector3d_cloud();

        // if empty cloud, can not set up radius search, add point to new set
        if (kdcloud->size() == 0)
        {
            add_point_coordinate(newPointID, thisPointVEC);

            int newSetID = getNewSetID();
            add_point_to_set(newPointID, newSetID);
            return;
        }

        // radius search size
        double search_size = 0.2;

        // perform radius search
        pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
        kdtree.setInputCloud(kdcloud);
        std::vector<int> search_indices;
        std::vector<float> search_dists;
        pcl::PointXYZ thisPointPCL;
        thisPointPCL.x = thisPointVEC[0];
        thisPointPCL.y = thisPointVEC[1];
        thisPointPCL.z = thisPointVEC[2];
        kdtree.radiusSearch(thisPointPCL, search_size, search_indices, search_dists, 0);


        // if no searched results, add point to new set
        if (search_indices.size() == 0)
        {
            add_point_coordinate(newPointID, thisPointVEC);

            int newSetID = getNewSetID();
            add_point_to_set(newPointID, newSetID);
            return;
        }

        // assuming all searched point are from the same set
        add_point_coordinate(newPointID, thisPointVEC);

        // find boundary points in the searched points
        std::set<int> searched_point_set = std::set<int>(search_indices.begin(), search_indices.end());
        std::set<int> boundary_points_set = compute_boundary_point_set();
        std::set<int> searched_boundary_points_set = intersection_of_sets(searched_point_set, boundary_points_set);    

        // get the current set
        int closet_searched_point_id = *searched_boundary_points_set.begin();
        int closet_set_id = point_to_set_map[closet_searched_point_id];
        std::cout << "searched_set_id: " << closet_set_id << std::endl;

        // compute set eigenvectors
        // if set size too small, can't compute eigenvectors?
        // [todo]
        Eigen::Matrix3d eigenvectors = get_set_eigenvectors(closet_set_id);

        // projection matrix : (projection_matrix.transpose() * point).head<2>()
        Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();  // Use the last two eigenvectors
        
        // points_to_vector2d_map
        std::set<int> points_to_project = searched_boundary_points_set; points_to_project.insert(newPointID);
        std::map<int, Eigen::Vector2d> points_to_vector2d_map = project_points_to_plane(points_to_project, projection_matrix);

        // existing edges between searched points (boundary)
        std::set<int> existing_boundary_edge_set = extract_existing_edge_between_points(searched_boundary_points_set, compute_boundary_edge_set());
        std::set<std::array<int, 2>> existing_boundary_edge_vertices_set = convert_to_edge_vertices_set(existing_boundary_edge_set);


        // to add a new point to mesh
        // - form edge to boundary point of the mesh, skip if the edge intersects any existing boundary edge
        // - form triangle if two used boundary points have a boundary edge between them, skip if the triangle contains other boundary points

        // add edge
        std::set<int> searched_boundary_points_used;
        for (int point_id : searched_boundary_points_set)
        {
            // new edge
            std::array<int, 2> newEdge = {point_id, newPointID};

            // skip if intersected
            if (edge_edges_intersection(newEdge, existing_boundary_edge_vertices_set, points_to_vector2d_map)) continue;

            // add edge
            int newEdgeID = getNewEdgeID();
            add_edge(newEdgeID, newEdge[0], newEdge[1]);

            // add to used
            searched_boundary_points_used.insert(point_id);
        }

        // add triangle
        for (const auto& existing_edge_vertices : existing_boundary_edge_vertices_set)
        {   
            // skip if not both points are used
            int i1 = existing_edge_vertices[0];
            int i2 = existing_edge_vertices[1];
            bool i1_used = searched_boundary_points_used.find(i1) != searched_boundary_points_used.end();
            bool i2_used = searched_boundary_points_used.find(i2) != searched_boundary_points_used.end();
            if (!i1_used || !i2_used) continue;

            // new triangle
            std::array<int, 3> newTriangle = {i1, i2, newPointID};

            // skip if triangle contains other boundary points
            if (triangle_contains_point(newTriangle, searched_boundary_points_set, points_to_vector2d_map)) continue;

            // add triangle
            int newTriangleID = getNewTriangleID();
            add_triangle(newTriangleID, closet_set_id, newTriangle);
        }


        // std::map<int, std::set<int>> set_to_searched_points_map = group_points_by_set(search_indices);

        // // add edge for each set of searched points
        // for (const auto& pair : set_to_searched_points_map)
        // {
        //     int set_id = pair.first;
        //     std::set<int> point_ids = pair.second;

        //     for (int point_id : point_ids)
        //     {
        //         // add edge
        //         int newEdgeID = getNewEdgeID();
        //         add_edge(newEdgeID, point_id, newPointID);
        //     }
        // }

        // // for each set of searched points
        // bool added_to_a_set = false;
        // for (const auto& pair : set_to_searched_points_map)
        // {
        //     // set id
        //     int set_id = pair.first;

        //     // searched points
        //     std::set<int> searched_point_ids = pair.second;

        //     // check plane if at least 3 points in the set, else consider the point is added to the set
        //     std::set<int> points_in_set = set_to_points_map[set_id];
        //     if (points_in_set.size() >= 3)
        //     {
        //         // ray set intersection
        //         Eigen::Vector3d rayPlaneIntersectionPoint = ray_set_intersection(rayOrigin, rayDirection, points_in_set, point_to_vector3d_map);
        //         double distance = (thisPointVEC - rayPlaneIntersectionPoint).norm();

        //         // skip if not within plane
        //         if (distance < -distance_threshold || distance_threshold < distance) // not within plane
        //         {
        //             continue;
        //         }
        //     }

        //     // update flag
        //     added_to_a_set = true;

        //     // update set to points map
        //     set_to_points_map[set_id].insert(next_point_id);

        //     // update boundary point
        //     boundary_points_set.insert(next_point_id);

        //     // update edge to vertices map
        //     for (int boundary_point_id : searched_point_ids)
        //     {
        //         // form edge
        //         std::array<int, 2> edge = {boundary_point_id, next_point_id}; // the smaller id is first

        //         // update map
        //         edge_to_point_map[next_edge_id] = edge;
        //         edge_to_point_map_reverse[edge] = next_edge_id;

        //         // update edge count
        //         edge_occurrences_count[next_edge_id] = 0;

        //         // update point to edge map
        //         point_to_edge_map[boundary_point_id].insert(next_edge_id);
        //         point_to_edge_map[next_point_id].insert(next_edge_id);

        //         // update boundary edge and point (whenever a new edge is added)
        //         boundary_edge_set.insert(next_edge_id);
        //         boundary_points_set.insert(boundary_point_id);

        //         // increment edge id
        //         next_edge_id ++;
        //     }

        //     // find boundary edges between the searched points
        //     std::set<std::array<int, 2>> existing_boundary_edges;
        //     // for each searched point
        //     for (int searched_point_id : searched_point_ids)
        //     {
        //         // for each edge of the searched point
        //         for (int edge_id : point_to_edge_map[searched_point_id])
        //         {
        //             // skip if not boundary edge
        //             if (boundary_edge_set.find(edge_id) == boundary_edge_set.end()) continue;
                    
        //             // add to list if the boundary edge is between searched points
        //             std::array<int, 2> vertices = edge_to_point_map[edge_id]; // [potential bug] need to make sure the vertices order are correct
        //             if (searched_point_ids.find(vertices[0]) != searched_point_ids.end() && searched_point_ids.find(vertices[1]) != searched_point_ids.end())
        //             {
        //                 existing_boundary_edges.insert(vertices);
        //             }
        //         }
        //     }

        //     // add triangle for each exsiting boundary edge
        //     for (const std::array<int, 2>& vertices : existing_boundary_edges)
        //     {
        //         // form triangle
        //         std::array<int, 3> triangle = {vertices[0], vertices[1], next_point_id};

        //         // add to triangle list
        //         triangle_to_vertices_map[next_triangle_id] = triangle;
        //         triangle_to_set_map[next_triangle_id] = set_id;
        //         set_to_triangles_map[set_id].insert(next_triangle_id);
        //         next_triangle_id ++;
                
        //         // obtain edge id
        //         int edge1 = edge_to_point_map_reverse[{vertices[0], vertices[1]}];
        //         int edge2 = edge_to_point_map_reverse[{vertices[1], next_point_id}];
        //         int edge3 = edge_to_point_map_reverse[{vertices[0], next_point_id}];

        //         // update edge count (whenever a new triangle is added)
        //         edge_occurrences_count[edge1] ++;
        //         edge_occurrences_count[edge2] ++;
        //         edge_occurrences_count[edge3] ++;

        //         // update boundary edge and point (whenever a new triangle is added)
        //         std::array<int, 3> edges = {edge1, edge2, edge3};
        //         for (int edge : edges)
        //         {
        //             if (edge_occurrences_count[edge] <= 1) // boundary edge if count is 0 or 1, since ++ this can't be 0
        //             {
        //                 boundary_edge_set.insert(edge);
        //                 boundary_points_set.insert(edge_to_point_map[edge][0]);
        //                 boundary_points_set.insert(edge_to_point_map[edge][1]);
        //             } 
        //             else 
        //             {
        //                 boundary_edge_set.erase(edge);

        //                 // if a point is not in any boundary edge, remove from boundary points
        //                 int vertex1 = edge_to_point_map[edge][0];
        //                 int vertex2 = edge_to_point_map[edge][1];
        //                 std::set<int> vertex1_edges = point_to_edge_map[vertex1];
        //                 std::set<int> vertex2_edges = point_to_edge_map[vertex2];
        //                 bool vertex1_boundary = false;
        //                 bool vertex2_boundary = false;
        //                 for (int edge_id : vertex1_edges)
        //                 {
        //                     if (boundary_edge_set.find(edge_id) != boundary_edge_set.end())
        //                     {
        //                         vertex1_boundary = true;
        //                         break;
        //                     }
        //                 }
        //                 for (int edge_id : vertex2_edges)
        //                 {
        //                     if (boundary_edge_set.find(edge_id) != boundary_edge_set.end())
        //                     {
        //                         vertex2_boundary = true;
        //                         break;
        //                     }
        //                 }
        //                 if (!vertex1_boundary) boundary_points_set.erase(vertex1);
        //                 if (!vertex2_boundary) boundary_points_set.erase(vertex2);
        //             }
        //         }

        //     }

        //     // increment next point
        //     next_point_id ++;

        // }

        // // if not added to any set, add the point as new set
        // if (!added_to_a_set)
        // {
        //     // ----------- add point as set ----------- 

        //     // add to set
        //     set_to_points_map[next_set_id].insert(next_point_id);

        //     // add to boundary points
        //     boundary_points_set.insert(next_point_id);

        //     // add to point to set map
        //     point_to_set_map[next_point_id] = next_set_id;

        //     // increment set id
        //     next_set_id ++;

        //     // increment point id
        //     next_point_id ++;
        // }
    }






    std::map<int, Eigen::Vector3d> get_point_to_vector3d_map() {return point_to_vector3d_map;};
    std::map<int, std::array<int, 2>> get_edge_to_vertices_map() {return edge_to_point_map;};
    
    pcl::PointCloud<pcl::PointXYZ>::Ptr point_to_vector3d_cloud()
    {
        return point_to_vector3d_cloud(point_list);
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr point_to_vector3d_cloud(std::vector<int> point_indices)
    {
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        for (int point_id : point_indices)
        {
            pcl::PointXYZ point;
            point.x = point_to_vector3d_map[point_id][0];
            point.y = point_to_vector3d_map[point_id][1];
            point.z = point_to_vector3d_map[point_id][2];
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr boundary_point_cloud()
    {
        return point_to_vector3d_cloud(compute_boundary_point_list());
    }

    int get_number_of_triangles()
    {
        return triangle_to_vertices_map.size();
    }

    
private:
    std::size_t i;
    typename pcl::PointCloud<VilensPointT>::Ptr new_cloud;

    // settings
    double distance_threshold = 0.05;

        // point
    int next_point_id = 0;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, int> point_to_set_map;
    std::map<int, std::set<int>> point_to_edge_map;
    std::vector<int> point_list;

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
    std::map<int, std::array<int, 2>> edge_to_point_map;
    std::map<std::array<int, 2>, int> edge_to_point_map_reverse;
    // std::map<int, int> edge_occurrences_count; // number of triangles that share the edge
    std::vector<int> edge_list;

    //     // boundary
    // std::set<int> boundary_points_set;
    // std::set<int> boundary_edge_set;
    
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

        // set up initial mesh
        pcl::PolygonMesh mesh;
        viewer_->addPolylineFromPolygonMesh(mesh, "polyline");
        pcl::PolygonMesh boundary_mesh;
        viewer_->addPolylineFromPolygonMesh(mesh, "boundary_edges");

        // set up initial point cloud
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        viewer_->addPointCloud(cloud, "cloud");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
        // set up boundary point cloud
        pcl::PointCloud<pcl::PointXYZ>::Ptr boundary_point_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        viewer_->addPointCloud(boundary_point_cloud, "boundary point cloud");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "boundary point cloud");


        // register keyboard callback
        viewer_->registerKeyboardCallback(&InteractiveViewer::keyboard_callback, *this, nullptr);

        // spin
        viewer_->spin();
    }

private:
    Application<PointT>& app_;

    pcl::visualization::PCLVisualizer::Ptr viewer_;
    
    void keyboard_callback(const pcl::visualization::KeyboardEvent &event, void*) 
    {
        bool space_down = event.getKeySym() == "space" && event.keyDown();
        if (space_down)
        {
            // step application
            app_.step();
            std::cout << "Number of triangles: " << app_.get_number_of_triangles() << std::endl;

            // get data from app
            std::map<int, Eigen::Vector3d> point_to_vector3d_map = app_.get_point_to_vector3d_map();
            std::map<int, std::array<int, 2>> edge_to_vertices_map = app_.get_edge_to_vertices_map();
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = app_.point_to_vector3d_cloud();
            pcl::PointCloud<pcl::PointXYZ>::Ptr boundary_point_cloud = app_.boundary_point_cloud();

            // mesh
            pcl::PolygonMesh mesh;
            // add points
            pcl::toPCLPointCloud2(*cloud, mesh.cloud); 
            // add edges
            for (const auto& pair : edge_to_vertices_map)
            {
                pcl::Vertices edge;
                edge.vertices.push_back(pair.second[0]);
                edge.vertices.push_back(pair.second[1]);
                mesh.polygons.push_back(edge);
            }

            // boundary mesh
            pcl::PolygonMesh boundary_mesh;
            // add points
            pcl::toPCLPointCloud2(*cloud, boundary_mesh.cloud);
            // get boundary edges
            std::vector<int> boundary_edge_list = app_.compute_boundary_edge_list();
            // add edges
            for (int edge_id : boundary_edge_list)
            {
                pcl::Vertices edge;
                edge.vertices.push_back(edge_to_vertices_map[edge_id][0]);
                edge.vertices.push_back(edge_to_vertices_map[edge_id][1]);
                boundary_mesh.polygons.push_back(edge);
            }

            // remove old mesh
            viewer_->removeShape("polyline");
            viewer_->addPolylineFromPolygonMesh(mesh, "polyline");

            viewer_->removeShape("boundary_edges");
            viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");
            // modify edge color
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 0.0, "boundary_edges");

            // update cloud
            pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> color_cloud(cloud, 0, 255, 0);
            viewer_->updatePointCloud<pcl::PointXYZ>(cloud, color_cloud, "cloud");
            pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> color_boundary_point_cloud(cloud, 255, 255, 0);
            viewer_->updatePointCloud<pcl::PointXYZ>(boundary_point_cloud, color_boundary_point_cloud, "boundary point cloud");
        }
    }  
};



using InputPointT = VilensPointT;
int main()
{
    // application
    Application<InputPointT> app;

    // interactive viewer
    InteractiveViewer<InputPointT> iviewer(app);

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