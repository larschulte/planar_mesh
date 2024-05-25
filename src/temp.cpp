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


class flann3d
{
public:
    flann3d() : flann_tree(flann::KDTreeIndexParams(1))
    {
        // placeholder
        std::set<int> point_set = {-1};
        std::map<int, Eigen::Vector3d> point_to_vector3d_map = {{-1, Eigen::Vector3d(0, 0, 0)}};
        set_input(point_set, point_to_vector3d_map);
    };
    
    void set_input(std::set<int> point_set, std::map<int, Eigen::Vector3d> point_to_vector3d_map)
    {
        // reserve
        flann_data_storage.reserve(point_set.size() * 3);
        index_to_pointID.reserve(point_set.size());

        // process
        for (int point_id : point_set)
        {
            // add to data storage
            flann_data_storage.push_back(point_to_vector3d_map[point_id][0]);
            flann_data_storage.push_back(point_to_vector3d_map[point_id][1]);
            flann_data_storage.push_back(point_to_vector3d_map[point_id][2]);

            // add to index
            index_to_pointID.push_back(point_id);
        }

        // add to flann
        flann_tree.buildIndex(flann::Matrix<double>(flann_data_storage.data(), point_set.size(), 3));

        // update id
        flann_last_id = point_set.size() - 1;
    }

    void addPoint(Eigen::Vector3d new_point, int point_id)
    {
        // add to data storage
        flann_data_storage.push_back(new_point[0]);
        flann_data_storage.push_back(new_point[1]);
        flann_data_storage.push_back(new_point[2]);

        // add to index
        index_to_pointID.push_back(point_id);

        // add to flann
        flann_tree.addPoints(flann::Matrix<double>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3));

        // update id
        flann_last_id++;
    }

    // void deletePoint(int point_id)
    // {
    //     // find index
    //     int index = -1;
    //     for (auto i = 0; i < index_to_pointID.size(); i++)
    //     {
    //         if (index_to_pointID[i] == point_id)
    //         {
    //             index = i;
    //             break;
    //         }
    //     }

    //     // if not found, return
    //     if (index == -1) return;

    //     // delete from data storage
    //     flann_data_storage.erase(flann_data_storage.begin() + index * 3, flann_data_storage.begin() + index * 3 + 3);

    //     // delete from index
    //     index_to_pointID.erase(index_to_pointID.begin() + index);

    //     // delete from flann
    //     flann_tree.removePoint(index);

    //     // update id
    //     flann_last_id--;
    // }
    
    std::set<int> radiusSearch(Eigen::Vector3d searchPoint, double radius)
    {
        // initialize
        std::vector<int> searched_indices;
        std::vector<double> searched_dists;

        // convert to vector
        std::vector<double> query_point = {searchPoint[0], searchPoint[1], searchPoint[2]}; 

        // intialize
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>());
        std::vector<std::vector<double>> list_of_search_dists(1, std::vector<double>());

        // search
        flann_tree.radiusSearch(flann::Matrix<double>(query_point.data(), 1, 3), list_of_search_indices, list_of_search_dists, radius * radius, flann::SearchParams(-1, 0));

        // extract
        searched_indices = list_of_search_indices[0];
        searched_dists = list_of_search_dists[0];

        // convert to id
        std::set<int> searched_ids;
        for (const int& index : searched_indices)
        {
            searched_ids.insert(index_to_pointID[index]);
        }

        // remove placeholder
        searched_ids.erase(-1);
        
        // return
        return searched_ids;
    }

    int flann_last_id;

private:
    std::vector<double> flann_data_storage;
    flann::Index<flann::L2_Simple<double>> flann_tree;

    std::vector<int> index_to_pointID; // index to point id correspondence
};


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
    }

    // want boundary edge and points to prevent overlapping triangles
    // to add a new point to mesh
    // 1. extract mesh boundary points
    // 2. perform delaunay triangulation using the boundary points and the new point
    // 3. add the edge and triangles connected to the new point to the original mesh

    void add_point(int newPointID, int setID, Eigen::Vector3d thisPoint, Eigen::Vector3d origin)
    {
        // add eigen
        point_to_vector3d_map[newPointID] = thisPoint;
        point_to_origin_vector3d_map[newPointID] = origin;

        // update plane estimation https://stats.stackexchange.com/questions/26123/efficient-method-technique-to-update-covariance-matrix
            // compute
        double old_size = static_cast<double>(set_to_points_map[setID].size()); // use double since will involve division later
        double new_size = old_size + 1;
        Eigen::Vector3d old_mean = set_to_mean_map[setID];
        Eigen::Vector3d new_mean = (old_mean * old_size + thisPoint) / new_size;
        Eigen::Matrix3d old_covariance_matrix = set_to_covariance_matrix_map[setID];
        Eigen::Matrix3d new_covariance_matrix = (old_size / new_size) * old_covariance_matrix + (old_size / (new_size * new_size)) * (thisPoint - new_mean) * (thisPoint - new_mean).transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(new_covariance_matrix);
        Eigen::Matrix3d eigenvectors = solver.eigenvectors();
        Eigen::Vector3d normal = eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal
            // store
        point_to_mean_used_map[newPointID] = old_mean;
        set_to_mean_map[setID] = new_mean;
        set_to_covariance_matrix_map[setID] = new_covariance_matrix;
        set_to_eigenvectors_map[setID] = eigenvectors;
        set_to_normal_map[setID] = normal;

        // add to set
        set_to_points_map[setID].insert(newPointID);
        point_to_sets_map[newPointID].insert(setID);

        // add boundary point
        add_boundary_point(newPointID, setID);
    }

    void add_edge(int newEdgeID, int setID, std::array<int, 2> newEdge)
    {
        // update edge to point 
        edge_to_point_map[newEdgeID] = newEdge;
        edge_to_point_map_reverse[newEdge] = newEdgeID;

        // update point to edge
        point_to_edges_map[newEdge[0]].insert(newEdgeID);
        point_to_edges_map[newEdge[1]].insert(newEdgeID);

        // update set to edge
        set_to_edges_map[setID].insert(newEdgeID);

        // update set boundary edge
        set_to_edge_count_map[setID][newEdgeID] = 0;
        set_to_boundary_edge_set_map[setID].insert(newEdgeID);
    
        // update boundary points
        add_boundary_point(newEdge[0], setID);
        add_boundary_point(newEdge[1], setID);
    }

    void add_triangle(int newTriangleID, int newSetID, std::array<int, 3> vertices)
    {
        // add triangle to vertices map
        triangle_to_vertices_map[newTriangleID] = vertices;

        // add triangle to set map
        set_to_triangles_map[newSetID].insert(newTriangleID);

        // update boundary edge and points
        for (int i = 0; i < 3; i++)
        {
            // get correct edge order
            int smaller_id = std::min(vertices[i], vertices[(i + 1) % 3]);
            int larger_id = std::max(vertices[i], vertices[(i + 1) % 3]);
            std::array<int, 2> edge = {smaller_id, larger_id};

            // increment count
            int& count = set_to_edge_count_map[newSetID][edge_to_point_map_reverse[edge]];
            count ++;

            // cases
            if (count <= 1)
            {
                // udpate set boundary edge
                set_to_boundary_edge_set_map[newSetID].insert(edge_to_point_map_reverse[edge]);
                
                // update boundary points
                add_boundary_point(edge[0], newSetID);
                add_boundary_point(edge[1], newSetID);
            } 
            else
            {
                // remove from set boundary edge
                set_to_boundary_edge_set_map[newSetID].erase(edge_to_point_map_reverse[edge]);

                // remove from boundary points (if the edge point no longer connects to a boundary edge, it is removed from boundary points)
                remove_boundary_point(edge[0], newSetID);
                remove_boundary_point(edge[1], newSetID);
            }
        }
    }

    void add_boundary_point(int pointID, int setID)
    {
        // update local boundary points
        set_to_boundary_point_set_map[setID].insert(pointID);

        // update global boundary points
        bool inserted = global_boundary_point_set.insert(pointID).second;
        if (inserted) flann.addPoint(point_to_vector3d_map[pointID], pointID);

        // if enough points, update intersection point
        if (set_to_points_map[setID].size() >= fit_plane_threshold)
        {
            point_to_intersection_vector3d_map[pointID] = point_set_intersection(pointID, setID);
        }
    }

    void remove_boundary_point(int pointID, int setID)
    {
        // skip if still boundary
        for (int edge_id : point_to_edges_map[pointID])
        {
            if (set_to_boundary_edge_set_map[setID].find(edge_id) != set_to_boundary_edge_set_map[setID].end()) return;
        }

        // remove from both local and global boundary points (a point can only be in one set)
        set_to_boundary_point_set_map[setID].erase(pointID);
        global_boundary_point_set.erase(pointID);

        // // skip for now
        // bool erased = global_boundary_point_set.erase(pointID);
        // if (erased) flann.deletePoint(pointID);
        
    }

    int getNewPointID()
    {
        // get id
        int newPointID = next_point_id;
        next_point_id ++;

        // update point list
        point_list.push_back(newPointID);

        // return
        return newPointID;
    }

    int getNewSetID()
    {
        // get id
        int newSetID = next_set_id;
        next_set_id ++;

        // update set list
        set_list.push_back(newSetID);

        // update set color map
        set_to_color_map[newSetID] = {rand() % 256, rand() % 256, rand() % 256};

        // return
        return newSetID;
    }

    int getNewEdgeID()
    {
        // get id
        int newEdgeID = next_edge_id;
        next_edge_id ++;

        // update edge list
        edge_list.push_back(newEdgeID);

        // return
        return newEdgeID;
    }
    int getNewTriangleID()
    {
        // get id
        int newTriangleID = next_triangle_id;
        next_triangle_id ++;

        // update triangle list
        triangle_list.push_back(newTriangleID);

        // return
        return newTriangleID;
    }

    std::map<int, std::set<int>> group_points_by_set(std::set<int> point_set)
    {
        // initialize
        std::map<int, std::set<int>> points_grouped_by_set;

        // process
        for (int point_id : point_set)
        {
            for (int set_id : point_to_sets_map[point_id])
            {
                points_grouped_by_set[set_id].insert(point_id);
            }
        }

        // return
        return points_grouped_by_set;
    }

    std::set<int> compute_boundary_edge_set_of_set(int setID)
    {
        // intialize edge occurrences count
        std::map<int, int> edge_count;
        for (int edge_id : set_to_edges_map[setID])
        {
            edge_count[edge_id] = 0;
        }

        // for each triangle, increment edge count
        for (int triangle_id : set_to_triangles_map[setID])
        {
            std::array<int, 3> vertices = triangle_to_vertices_map[triangle_id];
            for (int i = 0; i < 3; i++)
            {
                // get correct edge order
                int smaller_id = std::min(vertices[i], vertices[(i + 1) % 3]);
                int larger_id = std::max(vertices[i], vertices[(i + 1) % 3]);
                std::array<int, 2> edge = {smaller_id, larger_id};

                // increment count
                edge_count[edge_to_point_map_reverse[edge]] ++;
            }
        }
        
        // identify boundary edges (edges that are shared by only 0 or 1 triangle)
        std::set<int> boundary_edge_set;
        for (const auto& pair : edge_count)
        {
            if (pair.second <= 1) 
            {
                boundary_edge_set.insert(pair.first);
            }
        }

        // return
        return boundary_edge_set;
    }

    std::set<int> compute_boundary_point_set_of_set(int setID)
    {
        // get boundary edge set
        std::set<int> boundary_edge_set = compute_boundary_edge_set_of_set(setID);

        // identify boundary points
        // add points from boundary edges
        std::set<int> boundary_points_set;
        for (int edge_id : boundary_edge_set)
        {
            std::array<int, 2> vertices = edge_to_point_map[edge_id];
            boundary_points_set.insert(vertices[0]);
            boundary_points_set.insert(vertices[1]);
        }
        // add points that do not have edge in the same set
        for (int point_id : set_to_points_map[setID])
        {
            // if point have no edge, add
            if (point_to_edges_map.find(point_id) == point_to_edges_map.end())
            {
                boundary_points_set.insert(point_id);
            }
            // if point have edge, non is in the same set, add
            else
            {
                bool one_edge_in_same_set = false;
                for (int edge_id : point_to_edges_map[point_id])
                {
                    if (set_to_edges_map[setID].find(edge_id) != set_to_edges_map[setID].end())
                    {
                        one_edge_in_same_set = true;
                        break;
                    }
                }
                if (!one_edge_in_same_set) boundary_points_set.insert(point_id);
            }
        }

        return boundary_points_set;
    }

    std::set<int> compute_boundary_edge_set()
    {
        // initialize
        std::set<int> boundary_edge_set;

        // process
        for (int set_id : set_list)
        {
            // get boundary edge set of set
            std::set<int> boundary_edge_set_of_set = compute_boundary_edge_set_of_set(set_id);

            // add to boundary edge set
            boundary_edge_set.insert(boundary_edge_set_of_set.begin(), boundary_edge_set_of_set.end());
        }

        // return
        return boundary_edge_set;
    }

    std::vector<int> compute_boundary_edge_list()
    {
        // get boundary edge set
        std::set<int> boundary_edge_set = compute_boundary_edge_set();

        // convert to list
        std::vector<int> boundary_edge_list;
        for (int edge_id : boundary_edge_set)
        {
            boundary_edge_list.push_back(edge_id);
        }

        // return
        return boundary_edge_list;
    }

    std::set<int> compute_boundary_point_set()
    {
        // get boundary edge set
        std::set<int> boundary_edge_set = compute_boundary_edge_set();

        // identify boundary points
        // add points from boundary edges
        std::set<int> boundary_points_set;
        for (int edge_id : boundary_edge_set)
        {
            std::array<int, 2> vertices = edge_to_point_map[edge_id];
            boundary_points_set.insert(vertices[0]);
            boundary_points_set.insert(vertices[1]);
        }
        // add points that do not have edge 
        for (int point_id : point_list)
        {
            if (point_to_edges_map.find(point_id) == point_to_edges_map.end())
            {
                boundary_points_set.insert(point_id);
            }
        }

        return boundary_points_set;
    }

    // compute boundary point list
    std::vector<int> compute_boundary_point_list()
    {
        // get boundary point set
        std::set<int> boundary_points_set = compute_boundary_point_set();

        // convert to list
        std::vector<int> boundary_points_list;
        for (int point_id : boundary_points_set)
        {
            boundary_points_list.push_back(point_id);
        }

        // return
        return boundary_points_list;
    }

    std::set<int> get_boundary_edge_set_of_set(int setID)
    {
        return set_to_boundary_edge_set_map[setID];
    }
        
    std::set<int> get_boundary_edge_set()
    {
        // initialize
        std::set<int> boundary_edge_set;

        // process
        for (int set_id : set_list)
        {
            // get boundary edge set of set
            std::set<int> boundary_edge_set_of_set = get_boundary_edge_set_of_set(set_id);

            // add to boundary edge set
            boundary_edge_set.insert(boundary_edge_set_of_set.begin(), boundary_edge_set_of_set.end());
        }

        // return
        return boundary_edge_set;
    }

    std::vector<int> get_boundary_edge_list()
    {
        // get boundary edge set
        std::set<int> boundary_edge_set = get_boundary_edge_set();

        // convert to list
        std::vector<int> boundary_edge_list;
        for (int edge_id : boundary_edge_set)
        {
            boundary_edge_list.push_back(edge_id);
        }

        // return
        return boundary_edge_list;
    }

    std::set<int> get_boundary_point_set_of_set(int setID)
    {
        return set_to_boundary_point_set_map[setID];
    }
    
    std::set<int> get_boundary_point_set()
    {
        // initialize
        std::set<int> boundary_points_set;

        // process
        for (int set_id : set_list)
        {
            // get boundary point set of set
            std::set<int> boundary_point_set_of_set = get_boundary_point_set_of_set(set_id);

            // add to boundary point set
            boundary_points_set.insert(boundary_point_set_of_set.begin(), boundary_point_set_of_set.end());
        }

        // return
        return boundary_points_set;
    }

    std::vector<int> get_boundary_point_list()
    {
        // get boundary point set
        std::set<int> boundary_points_set = get_boundary_point_set();

        // convert to list
        std::vector<int> boundary_points_list;
        for (int point_id : boundary_points_set)
        {
            boundary_points_list.push_back(point_id);
        }

        // return
        return boundary_points_list;
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
    
    Eigen::Vector3d get_set_mean(int setID)
    {
        // initialize
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();

        // process
        for (int point_id : set_to_points_map[setID])
        {
            mean += point_to_vector3d_map[point_id];
        }
        mean /= set_to_points_map[setID].size();

        // return
        return mean;
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

    // project points to set
    std::map<int, Eigen::Vector2d> project_points_to_set(std::set<int> point_set, int setID)
    {
        // get set eigenvectors
        Eigen::Matrix3d eigenvectors = get_set_eigenvectors(setID);

        // project points to plane
        Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();
        std::map<int, Eigen::Vector2d> points_to_vector2d_map = project_points_to_plane(point_set, projection_matrix);

        // return
        return points_to_vector2d_map;
    }

    // project points to set with range correction
    std::map<int, Eigen::Vector2d> project_points_to_set_with_range_correction(std::set<int> point_set, int setID)
    {
        // initialize
        std::map<int, Eigen::Vector2d> points_to_vector2d_map;

        // get set eigenvectors
        Eigen::Matrix3d eigenvectors = get_set_eigenvectors(setID);
        Eigen::Vector3d mean = get_set_mean(setID);
        Eigen::Vector3d normal = eigenvectors.col(0);

        // for each point
        for (int point_id : point_set)
        {
            // intersection point
            Eigen::Vector3d rayOrigin = point_to_origin_vector3d_map[point_id];
            Eigen::Vector3d rayEndPoint = point_to_vector3d_map[point_id];
            Eigen::Vector3d rayDirection = rayEndPoint.normalized();
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayDirection, mean, normal);

            // project intersected points to plane
            Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();
            Eigen::Vector2d projected_point = (projection_matrix.transpose() * rayPlaneIntersectionPoint).head<2>();
            points_to_vector2d_map[point_id] = projected_point;
        }

        // return
        return points_to_vector2d_map;
    }
    
    std::map<int, Eigen::Vector2d> project_boundary_points_of_set_to_set_plane(int setID)
    {
        // initialize
        std::map<int, Eigen::Vector2d> points_to_vector2d_map;

        // process
        for (int point_id : set_to_boundary_point_set_map[setID])
        {
            // update if not available (todo - update if not accurate enough)
            if (point_to_intersection_vector3d_map.find(point_id) == point_to_intersection_vector3d_map.end())
            {
                point_to_intersection_vector3d_map[point_id] = point_set_intersection(point_id, setID);
            }

            // get intersection point
            Eigen::Vector3d rayPlaneIntersectionPoint = point_to_intersection_vector3d_map[point_id];

            // project intersection points to plane
            Eigen::Matrix3d eigenvectors = set_to_eigenvectors_map[setID];
            Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();
            Eigen::Vector2d projected_point = (projection_matrix.transpose() * rayPlaneIntersectionPoint).head<2>();

            // store
            points_to_vector2d_map[point_id] = projected_point;
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
    

    // ray plane intersection
    Eigen::Vector3d ray_plane_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayEndPoint, Eigen::Vector3d planeMean, Eigen::Vector3d planeNormal)
    {   
        // if perpendicular, return NaN
        Eigen::Vector3d rayDirection = rayEndPoint.normalized();
        if (planeNormal.dot(rayDirection) == 0)
        {
            return Eigen::Vector3d(NAN, NAN, NAN);
        }

        // compute intersection
        double distance = (planeMean - rayOrigin).dot(planeNormal) / rayDirection.dot(planeNormal);
        Eigen::Vector3d intersection = rayOrigin + rayDirection * distance;

        // return
        return intersection;
    }

    Eigen::Vector3d point_set_intersection(int pointID, int setID)
    {
        // compute
        Eigen::Vector3d rayOrigin = point_to_origin_vector3d_map[pointID];
        Eigen::Vector3d rayEndPoint = point_to_vector3d_map[pointID];
        Eigen::Vector3d mean = set_to_mean_map[setID];
        Eigen::Vector3d normal = set_to_normal_map[setID];
        Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);

        // return
        return rayPlaneIntersectionPoint;
    }

    // ray set intersection 
    Eigen::Vector3d ray_set_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayEndPoint, int setID)
    {
        // mean
        std::set<int> point_ids = set_to_points_map[setID];
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (int point_id : point_ids)
        {
            mean += point_to_vector3d_map[point_id];
        }
        mean /= point_ids.size();

        // covariance matrix
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (int point_id : point_ids)
        {
            Eigen::Vector3d point = point_to_vector3d_map[point_id];
            covariance_matrix += (point - mean) * (point - mean).transpose();
        }

        // normal
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance_matrix);
        Eigen::Vector3d normal = eigensolver.eigenvectors().col(0);

        // correct normal direction
        Eigen::Vector3d rayDirection = rayEndPoint.normalized();
        if (normal.dot(rayDirection) > 0)
        {
            normal = -normal;
        }

        // intersection point
        Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayDirection, mean, normal);

        // return
        return rayPlaneIntersectionPoint;
    }


    void step()
    {
        // for each point in new cloud
        if (i >= new_cloud->size()) return;
        std::cout << "Processing point " << i << " / " << new_cloud->size() << std::endl;
        Eigen::Vector3d thisPointVEC = new_cloud->points[i].getVector3fMap().cast<double>();
        Eigen::Vector3d thisPointOriginVEC = Eigen::Vector3d(0, 0, 0);
        i ++;
        
        // get new point id
        int newPointID = getNewPointID();
        
        // if empty, can not set up radius search, add point to new set
        if (global_boundary_point_set.size() == 0)
        {
            int newSetID = getNewSetID();
            add_point(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
            return;
        }

        // radius search size
        double search_size = 0.2;

        // perform flann3d radius search
        std::set<int> searched_boundary_points_set = flann.radiusSearch(thisPointVEC, search_size);
        searched_boundary_points_set = intersection_of_sets(searched_boundary_points_set, global_boundary_point_set);

        // if no searched results, add point to new set
        if (searched_boundary_points_set.size() == 0)
        {
            int newSetID = getNewSetID();
            add_point(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
            return;
        }

        // for each sets of searched points
        std::map<int, std::set<int>> set_to_searched_points_map = group_points_by_set(searched_boundary_points_set);

        // record size and distance
        std::map<int, std::size_t> small_set_to_size_map;
        std::map<int, double> set_to_distance_map;
        for (const auto& pair : set_to_searched_points_map)
        {
            // set id and points
            int setID = pair.first;
            std::set<int> searched_boundary_points_in_current_set = pair.second;
            std::size_t set_size = set_to_points_map[setID].size(); 

            // if large set, record distance
            if (set_size > fit_plane_threshold)
            {
                Eigen::Vector3d mean = set_to_mean_map[setID];
                Eigen::Vector3d normal = set_to_normal_map[setID];
                Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(thisPointOriginVEC, thisPointVEC, mean, normal);

                double distance = (thisPointVEC - rayPlaneIntersectionPoint).norm();
                set_to_distance_map[setID] = abs(distance);
            }
            // if small set, record size
            else 
            {
                small_set_to_size_map[setID] = set_size;
            }
        }

        // find the smallest distance and largest size
        int smallest_distance_key;
        double smallest_distance_value = std::numeric_limits<double>::max();
        for (const auto& pair : set_to_distance_map)
        {
            if (pair.second < smallest_distance_value) 
            {
                smallest_distance_key = pair.first;
                smallest_distance_value = pair.second;
            }
        }
        int largest_size_key;
        std::size_t largest_size_value = 0;
        for (const auto& pair : small_set_to_size_map)
        {
            if (pair.second > largest_size_value) 
            {
                largest_size_key = pair.first;
                largest_size_value = pair.second;
            }
        }

        // compute the set to add to
        int setID;
        if (smallest_distance_value < distance_threshold)
        {
            setID = smallest_distance_key;
            
        }
        else if (largest_size_value > 0)
        {
            setID = largest_size_key;
        }
        else
        {
            setID = getNewSetID(); 
            add_point(newPointID, setID, thisPointVEC, thisPointOriginVEC);
            return;
        }

        // get the searched points in the set
        std::set<int> searched_boundary_points_in_current_set = set_to_searched_points_map[setID];

        // add point
        add_point(newPointID, setID, thisPointVEC, thisPointOriginVEC);

        // points_to_vector2d_map
        std::map<int, Eigen::Vector2d> points_to_vector2d_map = project_boundary_points_of_set_to_set_plane(setID);

        // existing edges between searched points (boundary)
        std::set<int> existing_boundary_edge_set = extract_existing_edge_between_points(searched_boundary_points_in_current_set, get_boundary_edge_set_of_set(setID));
        std::set<std::array<int, 2>> existing_boundary_edge_vertices_set = convert_to_edge_vertices_set(existing_boundary_edge_set);


        // to add a new point to mesh
        // - form edge to boundary point of the mesh, skip if the edge intersects any existing boundary edge
        // - form triangle if two used boundary points have a boundary edge between them, skip if the triangle contains other boundary points

        // add edge
        std::set<int> searched_boundary_points_used;
        for (int point_id : searched_boundary_points_in_current_set)
        {
            // new edge
            std::array<int, 2> newEdge = {point_id, newPointID};

            // skip if intersected with any boundary edge of the current set
            if (edge_edges_intersection(newEdge, convert_to_edge_vertices_set(get_boundary_edge_set_of_set(setID)), points_to_vector2d_map)) continue;

            // add edge
            int newEdgeID = getNewEdgeID();
            add_edge(newEdgeID, setID, newEdge);

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
            if (triangle_contains_point(newTriangle, searched_boundary_points_in_current_set, points_to_vector2d_map)) continue;

            // add triangle
            int newTriangleID = getNewTriangleID();
            add_triangle(newTriangleID, setID, newTriangle);
        }
    }

    void loop()
    {
        // finish all points in step
        while (i < new_cloud->size())
        {
            step();
        }
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

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_colored_cloud(std::vector<int> point_indices)
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (int point_id : point_indices)
        {
            pcl::PointXYZRGB point;
            point.x = point_to_vector3d_map[point_id][0];
            point.y = point_to_vector3d_map[point_id][1];
            point.z = point_to_vector3d_map[point_id][2];
            point.r = set_to_color_map[*point_to_sets_map[point_id].begin()].at(0);
            point.g = set_to_color_map[*point_to_sets_map[point_id].begin()].at(1);
            point.b = set_to_color_map[*point_to_sets_map[point_id].begin()].at(2);
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_colored_cloud()
    {
        return point_to_vector3d_set_colored_cloud(point_list);
    }
    

    pcl::PointCloud<pcl::PointXYZ>::Ptr boundary_point_cloud()
    {
        return point_to_vector3d_cloud(get_boundary_point_list());
    }

    int get_number_of_triangles()
    {
        return triangle_to_vertices_map.size();
    }

    
private:
    std::size_t i = 0;
    typename pcl::PointCloud<VilensPointT>::Ptr new_cloud;

    // settings
    double distance_threshold = 0.05;
    std::size_t fit_plane_threshold = 10; // may cause error if below 3

        // point
    int next_point_id = 0;
    std::vector<int> point_list;
    std::map<int, Eigen::Vector3d> point_to_origin_vector3d_map;
    std::map<int, Eigen::Vector3d> point_to_intersection_vector3d_map;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, std::set<int>> point_to_sets_map;
    std::map<int, std::set<int>> point_to_edges_map;

        // triangle
    int next_triangle_id = 0;
    std::vector<int> triangle_list;
    std::map<int, std::array<int, 3>> triangle_to_vertices_map;
    std::map<int, int> triangle_to_set_map;
    std::map<int, std::set<int>> triangle_to_points_map;
    
        // set
    int next_set_id = 0;
    std::vector<int> set_list;
    std::map<int, std::array<int, 3>> set_to_color_map; // map to intensity
    std::map<int, std::set<int>> set_to_points_map; // each set contains id to points
    std::map<int, std::set<int>> set_to_edges_map; // each set contains id to edges
    std::map<int, std::set<int>> set_to_triangles_map; // each set contains id to triangles
    std::map<int, std::set<int>> set_to_boundary_edge_set_map; // each set contains id to boundary edges
    std::map<int, std::set<int>> set_to_boundary_point_set_map; // each set contains id to boundary points
    std::map<int, std::map<int, int>> set_to_edge_count_map; // each map contains edge count

        // plane fitting
    std::map<int, Eigen::Vector3d> point_to_mean_used_map;
    std::map<int, Eigen::Vector3d> set_to_mean_map;
    std::map<int, Eigen::Matrix3d> set_to_covariance_matrix_map;
    std::map<int, Eigen::Matrix3d> set_to_eigenvectors_map;
    std::map<int, Eigen::Vector3d> set_to_normal_map;
    
        // edge
    int next_edge_id = 0;
    std::vector<int> edge_list;
    std::map<int, std::array<int, 2>> edge_to_point_map;
    std::map<std::array<int, 2>, int> edge_to_point_map_reverse;

        // boundary
    flann3d flann;
    std::set<int> global_boundary_point_set;
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

        // mesh
        pcl::PolygonMesh mesh;
        viewer_->addPolylineFromPolygonMesh(mesh, "polyline");

        // boundary edges
        pcl::PolygonMesh boundary_mesh;
        viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");

        // colored cloud
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        viewer_->addPointCloud(cloud, "cloud");

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
            app_.loop();
            // std::cout << "Number of triangles: " << app_.get_number_of_triangles() << std::endl;

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
            std::vector<int> boundary_edge_list = app_.get_boundary_edge_list();
            // add edges
            for (int edge_id : boundary_edge_list)
            {
                pcl::Vertices edge;
                edge.vertices.push_back(edge_to_vertices_map[edge_id][0]);
                edge.vertices.push_back(edge_to_vertices_map[edge_id][1]);
                boundary_mesh.polygons.push_back(edge);
            }

            // mesh
            viewer_->removeShape("polyline");
            viewer_->addPolylineFromPolygonMesh(mesh, "polyline");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.1, 0.1, 0.1, "polyline");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 0.3, "polyline");

            // boundary edges
            viewer_->removeShape("boundary_edges");
            viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.2, 0.2, 0.2, "boundary_edges");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "boundary_edges");

            // set colored cloud
            pcl::PointCloud<pcl::PointXYZRGB>::Ptr intensity_cloud = app_.point_to_vector3d_set_colored_cloud();
            pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> intensity_handler(intensity_cloud);
            viewer_->updatePointCloud<pcl::PointXYZRGB>(intensity_cloud, intensity_handler, "cloud");
            viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "cloud");
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