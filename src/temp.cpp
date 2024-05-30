#include <pcl/common/transforms.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include "point_type/VilensPointT.hpp"
#include "eye_patch/DataLoader.hpp"
#include "eye_patch/TriangleBVH.hpp"
#include "eye_patch/KDTree.hpp"
#include "eye_patch/utilities.hpp"

// queue
#include <queue>

// application class
template <typename PointT>
class Application
{
public:

    int getNewSetID()
    {
        // get id
        int newSetID = next_set_id;
        next_set_id ++;

        // initialize
        set_to_color_map[newSetID] = {};
        set_to_points_map[newSetID] = {};
        set_to_edges_map[newSetID] = {};
        set_to_triangles_map[newSetID] = {};
        set_to_mean_map[newSetID] = Eigen::Vector3d::Zero();
        set_to_covariance_matrix_map[newSetID] = Eigen::Matrix3d::Zero();
        set_to_eigenvectors_map[newSetID] = Eigen::Matrix3d::Identity();
        set_to_eigenvalues_map[newSetID] = Eigen::Vector3d::Zero();
        set_to_normal_map[newSetID] = Eigen::Vector3d(0, 0, 1);

        // add set
        set_list.push_back(newSetID);
        set_to_color_map.at(newSetID) = {rand() % 256, rand() % 256, rand() % 256};

        // return
        return newSetID;
    }

    int getNewPointID()
    {
        // get id
        int newPointID = next_point_id;
        next_point_id ++;

        // initialize
        point_to_origin_vector3d_map[newPointID] = Eigen::Vector3d(999, 999, 999);
        point_to_intersection_vector3d_map[newPointID] = Eigen::Vector3d(999, 999, 999);
        point_to_vector3d_map[newPointID] = Eigen::Vector3d(999, 999, 999);
        point_to_set_map[newPointID] = -999;
        point_to_edges_map[newPointID] = {};

        // return
        return newPointID;
    }

    int getNewEdgeID()
    {
        // get id
        int newEdgeID = next_edge_id;
        next_edge_id ++;

        // initialize (skip edge_to_point_map_reverse as can't be initialized)
        edge_to_point_map[newEdgeID] = {};
        edge_to_set_map[newEdgeID] = -999;

        // return
        return newEdgeID;
    }

    int getNewTriangleID()
    {
        // get id
        int newTriangleID = next_triangle_id;
        next_triangle_id ++;

        // initialize (skip triangle_to_vertices_map_reverse as can't be initialized)
        triangle_to_vertices_map[newTriangleID] = {};
        triangle_to_set_map[newTriangleID] = -999;
        triangle_to_points_map[newTriangleID] = {};

        // return
        return newTriangleID;
    }

    void add_point(int newPointID, int setID, Eigen::Vector3d thisPoint, Eigen::Vector3d origin)
    {
        point_list.push_back(newPointID);
        point_to_vector3d_map.at(newPointID) = thisPoint;
        point_to_origin_vector3d_map.at(newPointID) = origin;
        point_to_set_map.at(newPointID) = setID;
        set_to_points_map.at(setID).insert(newPointID);
        add_to_plane_estimate(newPointID, setID);
        update_point_type(newPointID);
    }

    void add_edge(int newEdgeID, int setID, std::array<int, 2> vertices)
    {
        vertices = sortTwoInts(vertices[0], vertices[1]);
        int pointID1 = vertices[0];
        int pointID2 = vertices[1];

        edge_list.push_back(newEdgeID);
        edge_to_point_map.at(newEdgeID) = vertices;
        edge_to_point_map_reverse[vertices] = newEdgeID;
        point_to_edges_map.at(pointID1).insert(newEdgeID);
        point_to_edges_map.at(pointID2).insert(newEdgeID);
        set_to_edges_map.at(setID).insert(newEdgeID);
        edge_to_set_map.at(newEdgeID) = setID;
        edge_to_edge_count_map[newEdgeID] = 0;
        update_edge_type(newEdgeID);
        update_point_type(pointID1);
        update_point_type(pointID2);
    }

    void add_triangle(int newTriangleID, int setID, std::array<int, 3> vertices)
    {
        vertices = sortThreeInts(vertices[0], vertices[1], vertices[2]);
        int pointID1 = vertices[0];
        int pointID2 = vertices[1];
        int pointID3 = vertices[2];
        std::array<int, 2> edge1 = {pointID1, pointID2};
        std::array<int, 2> edge2 = {pointID2, pointID3};
        std::array<int, 2> edge3 = {pointID1, pointID3};
        int edgeID1 = edge_to_point_map_reverse.at(edge1);
        int edgeID2 = edge_to_point_map_reverse.at(edge2);
        int edgeID3 = edge_to_point_map_reverse.at(edge3);

        triangle_list.push_back(newTriangleID);
        triangle_to_vertices_map.at(newTriangleID) = vertices;
        triangle_to_vertices_map_reverse[vertices] = newTriangleID;
        set_to_triangles_map.at(setID).insert(newTriangleID);
        triangle_to_set_map.at(newTriangleID) = setID;

        edge_to_edge_count_map.at(edgeID1) ++;
        edge_to_edge_count_map.at(edgeID2) ++;
        edge_to_edge_count_map.at(edgeID3) ++;
        update_edge_type(edgeID1);
        update_edge_type(edgeID2);
        update_edge_type(edgeID3);
        update_point_type(pointID1);
        update_point_type(pointID2);
        update_point_type(pointID3);

        bool inserted = global_triangle_set.insert(newTriangleID).second;
        if (inserted) bvhRoot.addTriangle(newTriangleID, vertices, point_to_vector3d_map.at(pointID1), point_to_vector3d_map.at(pointID2), point_to_vector3d_map.at(pointID3));
    }

    void remove_point(int pointID)
    {
        // original info
        Eigen::Vector3d thisPoint = point_to_vector3d_map.at(pointID);
        Eigen::Vector3d origin = point_to_origin_vector3d_map.at(pointID);
        int setID = point_to_set_map.at(pointID);

        // remove
        point_list.erase(std::remove(point_list.begin(), point_list.end(), pointID), point_list.end());
        point_to_vector3d_map.erase(pointID);
        point_to_origin_vector3d_map.erase(pointID);
        point_to_set_map.erase(pointID);
        set_to_points_map.at(setID).erase(pointID);
        remove_from_plane_estimate(pointID, setID);

        // update boundary
        remove_from_point_type(pointID);
        
        // add to free points
        free_points_queue.push(std::make_pair(thisPoint, origin));
    }

    void remove_edge(int edgeID)
    {
        // original info
        std::array<int, 2> vertices = edge_to_point_map.at(edgeID);
        int pointID1 = vertices[0];
        int pointID2 = vertices[1];
        int setID = edge_to_set_map.at(edgeID);

        // remove
        edge_list.erase(std::remove(edge_list.begin(), edge_list.end(), edgeID), edge_list.end());
        edge_to_point_map.erase(edgeID);
        edge_to_point_map_reverse.erase(vertices);
        point_to_edges_map.at(pointID1).erase(edgeID);
        point_to_edges_map.at(pointID2).erase(edgeID);
        set_to_edges_map.at(setID).erase(edgeID);
        edge_to_set_map.erase(edgeID);
        edge_to_edge_count_map.erase(edgeID);

        // update boundary
        remove_from_edge_type(edgeID);
        update_point_type(pointID1);
        update_point_type(pointID2);
        
        // if the removal of edge causes a point to be isolated, remove the point
        if (point_to_edges_map.at(pointID1).empty()) remove_point(pointID1);
        if (point_to_edges_map.at(pointID2).empty()) remove_point(pointID2);
    }

    void remove_triangle(int triangleID)
    {
        // original info
        std::array<int, 3> vertices = triangle_to_vertices_map.at(triangleID);
        int pointID1 = vertices[0];
        int pointID2 = vertices[1];
        int pointID3 = vertices[2];
        int edgeID1 = edge_to_point_map_reverse.at({pointID1, pointID2});
        int edgeID2 = edge_to_point_map_reverse.at({pointID2, pointID3});
        int edgeID3 = edge_to_point_map_reverse.at({pointID1, pointID3});
        int setID = triangle_to_set_map.at(triangleID);
        std::set<int> points_within_triangles = triangle_to_points_map.at(triangleID);

        // remove
        triangle_list.erase(std::remove(triangle_list.begin(), triangle_list.end(), triangleID), triangle_list.end());
        triangle_to_vertices_map.erase(triangleID);
        triangle_to_vertices_map_reverse.erase(vertices);
        triangle_to_set_map.erase(triangleID);
        set_to_triangles_map.at(setID).erase(triangleID);
        triangle_to_points_map.erase(triangleID);
        bool erased = global_triangle_set.erase(triangleID);
        if (erased) bvhRoot.deleteTriangle(triangleID);

        // update boundary
        edge_to_edge_count_map.at(edgeID1) --;
        edge_to_edge_count_map.at(edgeID2) --;
        edge_to_edge_count_map.at(edgeID3) --;
        update_edge_type(edgeID1);
        update_edge_type(edgeID2);
        update_edge_type(edgeID3);
        update_point_type(pointID1);
        update_point_type(pointID2);
        update_point_type(pointID3);

        // if the removal of triangle causes an edge to be isolated, remove the edge
        if (edge_to_edge_count_map.at(edgeID1) == 0) remove_edge(edgeID1);
        if (edge_to_edge_count_map.at(edgeID2) == 0) remove_edge(edgeID2);
        if (edge_to_edge_count_map.at(edgeID3) == 0) remove_edge(edgeID3);

        // if the removal of triangle causes a point to be isolated, remove the point
        for (int pointID : points_within_triangles) remove_point(pointID);
    }

    void update_point_type(int pointID)
    {
        // is isolated
        bool isolated = point_to_edges_map.at(pointID).empty();
        // is connected to boundary edge
        bool connected_to_boundary_edge = false;
        for (int edge_id : point_to_edges_map.at(pointID))
        {
            if (boundary_edge_set.find(edge_id) != boundary_edge_set.end())
            {
                connected_to_boundary_edge = true;
                break;
            }
        }

        // update
        if (isolated || connected_to_boundary_edge)
        {
            bool inserted = boundary_point_set.insert(pointID).second;
            if (inserted) kdtree.addPoint(point_to_vector3d_map.at(pointID), pointID);
        }
        else
        {
            bool erased = boundary_point_set.erase(pointID);
            if (erased) kdtree.deletePoint(pointID);
        }
    }

    void update_edge_type(int edgeID)
    {
        int count = edge_to_edge_count_map.at(edgeID);

        // is a boundary edge
        bool is_boundary_edge = count == 0 || count == 1;
        if (is_boundary_edge)
        {
            boundary_edge_set.insert(edgeID);            
        }
        else
        {
            boundary_edge_set.erase(edgeID);
        }
    }

    void remove_from_point_type(int pointID)
    {
        bool erased = boundary_point_set.erase(pointID);
        if (erased) kdtree.deletePoint(pointID);
    }

    void remove_from_edge_type(int edgeID)
    {
        boundary_edge_set.erase(edgeID);
    }

    void add_to_plane_estimate(int pointID, int setID)
    {
        if (set_to_points_map.at(setID).size() == 1)
        {
            // initialize
            set_to_mean_map.at(setID) = point_to_vector3d_map.at(pointID);
            set_to_covariance_matrix_map.at(setID) = Eigen::Matrix3d::Zero();
            set_to_eigenvectors_map.at(setID) = Eigen::Matrix3d::Identity();
            set_to_eigenvalues_map.at(setID) = Eigen::Vector3d::Zero();
            set_to_normal_map.at(setID) = Eigen::Vector3d(0, 0, 1);
        }
        else
        {
            // set
            int size1 = set_to_points_map.at(setID).size();
            Eigen::Vector3d mean1 = set_to_mean_map.at(setID);
            Eigen::Matrix3d cov1 = set_to_covariance_matrix_map.at(setID);

            // point
            int size2 = 1;
            Eigen::Vector3d mean2 = point_to_vector3d_map.at(pointID);
            Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();

            // set + point
            Eigen::Vector3d new_mean = merge_means(mean1, mean2, size1, size2);
            Eigen::Matrix3d new_cov = merge_covariances(cov1, cov2, mean1, mean2, size1, size2);

            // plane estimate
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(new_cov);
            Eigen::Matrix3d new_eigenvectors = solver.eigenvectors();
            Eigen::Vector3d new_eigenvalues = solver.eigenvalues();
            Eigen::Vector3d new_normal = new_eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal
            Eigen::Vector3d vector_towards_origin = point_to_origin_vector3d_map.at(pointID) - point_to_vector3d_map.at(pointID);
            if (new_normal.dot(vector_towards_origin) < 0) new_normal *= -1; // normal should points towards the origin

            // store
            set_to_mean_map.at(setID) = new_mean;
            set_to_covariance_matrix_map.at(setID) = new_cov;
            set_to_eigenvectors_map.at(setID) = new_eigenvectors;
            set_to_eigenvalues_map.at(setID) = new_eigenvalues;
            set_to_normal_map.at(setID) = new_normal;
        }
    }

    void remove_from_plane_estimate(int pointID, int setID)
    {
        // pass
    }

    std::set<int> get_boundary_point_of_set(int setID)
    {
        // initialize
        std::set<int> boundary_point_of_set;

        // process
        for (int pointID : set_to_points_map.at(setID))
        {
            if (boundary_point_set.find(pointID) != boundary_point_set.end()) 
            {
                boundary_point_of_set.insert(pointID);
            }
        }

        // return
        return boundary_point_of_set;
    }

    std::set<int> get_boundary_edge_of_set(int setID)
    {
        // initialize
        std::set<int> boundary_edge_of_set;

        // process
        for (int edgeID : set_to_edges_map.at(setID))
        {
            if (boundary_edge_set.find(edgeID) != boundary_edge_set.end()) 
            {
                boundary_edge_of_set.insert(edgeID);
            }
        }

        // return
        return boundary_edge_of_set;
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
            Eigen::Vector2d pointA0 = points_to_vector2d_map.at(edgeA[0]);
            Eigen::Vector2d pointA1 = points_to_vector2d_map.at(edgeA[1]);
            Eigen::Vector2d pointB0 = points_to_vector2d_map.at(edgeB[0]);
            Eigen::Vector2d pointB1 = points_to_vector2d_map.at(edgeB[1]);
            if (doIntersect(pointA0, pointA1, pointB0, pointB1)) return true;
        }

        return false;
    }

    bool edge_edges_intersection(std::array<int, 2> edgeA, std::set<int> edgeB_set, std::map<int, Eigen::Vector2d> points_to_vector2d_map)
    {
        for (const auto& edgeB_ID : edgeB_set)
        {
            // if intersected at end points, don't count
            std::array<int, 2> edgeB = edge_to_point_map.at(edgeB_ID);
            if (edgeA[0] == edgeB[0] || edgeA[0] == edgeB[1] || edgeA[1] == edgeB[0] || edgeA[1] == edgeB[1]) continue;

            if (points_to_vector2d_map.find(edgeB[0]) == points_to_vector2d_map.end() || points_to_vector2d_map.find(edgeB[1]) == points_to_vector2d_map.end() || points_to_vector2d_map.find(edgeA[0]) == points_to_vector2d_map.end() || points_to_vector2d_map.find(edgeA[1]) == points_to_vector2d_map.end())
            {
                std::cout << point_to_set_map.at(edgeA[0]) << " " << point_to_set_map.at(edgeA[1]) << std::endl;
                std::cout << point_to_set_map.at(edgeB[0]) << " " << point_to_set_map.at(edgeB[1]) << std::endl;
            }

            // intersection check
            Eigen::Vector2d pointA0 = points_to_vector2d_map.at(edgeA[0]);
            Eigen::Vector2d pointA1 = points_to_vector2d_map.at(edgeA[1]);
            Eigen::Vector2d pointB0 = points_to_vector2d_map.at(edgeB[0]);
            Eigen::Vector2d pointB1 = points_to_vector2d_map.at(edgeB[1]);
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
            Eigen::Vector2d point = points_to_vector2d_map.at(point_id);
            Eigen::Vector2d a = points_to_vector2d_map.at(triangle[0]);
            Eigen::Vector2d b = points_to_vector2d_map.at(triangle[1]);
            Eigen::Vector2d c = points_to_vector2d_map.at(triangle[2]);
            if (point_in_triangle(point, a, b, c)) return true;
        }
    
        return false;
    }
    
    std::map<int, Eigen::Vector2d> project_boundary_points_of_set_to_set_plane(int setID)
    {
        // initialize
        std::map<int, Eigen::Vector2d> points_to_vector2d_map;

        // process
        for (int point_id : get_boundary_point_of_set(setID))
        {
            // // update if not available (todo - update if not accurate enough)
            // if (point_to_intersection_vector3d_map.find(point_id) == point_to_intersection_vector3d_map.end())
            // {
            //     point_to_intersection_vector3d_map.at(point_id) = point_set_intersection(point_id, setID);
            // }

            // get intersection point
            Eigen::Vector3d rayOrigin = point_to_origin_vector3d_map.at(point_id);
            Eigen::Vector3d rayEndPoint = point_to_vector3d_map.at(point_id);
            Eigen::Vector3d mean = set_to_mean_map.at(setID);
            Eigen::Vector3d normal = set_to_normal_map.at(setID);
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);

            // project intersection points to plane
            Eigen::Matrix3d eigenvectors = set_to_eigenvectors_map.at(setID);
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

    // creates edges and triangles that connects the new point to the set
    void connect_point_to_set(int newPointID, int setID, std::set<int> searched_boundary_points_in_current_set)
    {
        // add point as boundary point
        update_point_type(newPointID);

        // check if newPointID is boundary point
        if (boundary_point_set.find(newPointID) != boundary_point_set.end()) 
        {
            std::cout << "IS boundary point" << std::endl;
        }
        else
        {
            std::cout << "NOT boundary point" << std::endl;
        }

        // points_to_vector2d_map
        std::map<int, Eigen::Vector2d> points_to_vector2d_map = project_boundary_points_of_set_to_set_plane(setID);

        // existing edges between searched points (boundary)
        std::set<int> existing_boundary_edge_set = extract_existing_edge_between_points(searched_boundary_points_in_current_set, get_boundary_edge_of_set(setID));


        // to add a new point to mesh
        // - form edge to boundary point of the mesh, skip if the edge intersects any existing boundary edge
        // - form triangle if two used boundary points have a boundary edge between them, skip if the triangle contains other boundary points

        // add edge
        std::set<int> searched_boundary_points_used;
        for (int point_id : searched_boundary_points_in_current_set)
        {
            // new edge, smaller id first
            std::array<int, 2> newEdge = {std::min(newPointID, point_id), std::max(newPointID, point_id)};

            // skip if edge already exists
            if (edge_to_point_map_reverse.find(newEdge) != edge_to_point_map_reverse.end()) continue;

            // skip if intersected with any boundary edge of the current set
            if (edge_edges_intersection(newEdge, get_boundary_edge_of_set(setID), points_to_vector2d_map)) continue;

            // add edge
            int newEdgeID = getNewEdgeID();
            add_edge(newEdgeID, setID, newEdge);

            // add to used
            searched_boundary_points_used.insert(point_id);
        }

        // add triangle
        for (const auto& edgeID : existing_boundary_edge_set)
        {   
            // skip if not both points are used
            int i1 = edge_to_point_map.at(edgeID)[0];
            int i2 = edge_to_point_map.at(edgeID)[1];
            bool i1_used = searched_boundary_points_used.find(i1) != searched_boundary_points_used.end();
            bool i2_used = searched_boundary_points_used.find(i2) != searched_boundary_points_used.end();
            if (!i1_used || !i2_used) continue;

            // new triangle, smaller id first
            std::array<int, 3> newTriangle = sortThreeInts(newPointID, i1, i2);

            // skip if triangle already exists
            if (triangle_to_vertices_map_reverse.find(newTriangle) != triangle_to_vertices_map_reverse.end()) continue;

            // skip if triangle contains other boundary points
            if (triangle_contains_point(newTriangle, searched_boundary_points_in_current_set, points_to_vector2d_map)) continue;

            // add triangle
            int newTriangleID = getNewTriangleID();
            add_triangle(newTriangleID, setID, newTriangle);
        }
    }

    Eigen::Vector3d merge_means_of_sets(int setID1, int setID2) 
    {
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID1);
        Eigen::Vector3d mean2 = set_to_mean_map.at(setID2);
        int size1 = set_to_points_map.at(setID1).size();
        int size2 = set_to_points_map.at(setID2).size();
        return merge_means(mean1, mean2, size1, size2);
    }

    Eigen::Vector3d merge_means_between_set_and_point(int setID, Eigen::Vector3d point)
    {
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID);
        Eigen::Vector3d mean2 = point;
        int size1 = set_to_points_map.at(setID).size();
        int size2 = 1;
        return merge_means(mean1, mean2, size1, size2);
    }

    Eigen::Matrix3d merge_covariances_of_sets(int setID1, int setID2) 
    {
        Eigen::Matrix3d cov1 = set_to_covariance_matrix_map.at(setID1);
        Eigen::Matrix3d cov2 = set_to_covariance_matrix_map.at(setID2);
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID1);
        Eigen::Vector3d mean2 = set_to_mean_map.at(setID2);
        int size1 = set_to_points_map.at(setID1).size();
        int size2 = set_to_points_map.at(setID2).size();
        return merge_covariances(cov1, cov2, mean1, mean2, size1, size2);
    }

    Eigen::Matrix3d merge_covariances_between_set_and_point(int setID, Eigen::Vector3d point)
    {
        Eigen::Matrix3d cov1 = set_to_covariance_matrix_map.at(setID);
        Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero(); 
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID);
        Eigen::Vector3d mean2 = point;
        int size1 = set_to_points_map.at(setID).size();
        int size2 = 1;
        return merge_covariances(cov1, cov2, mean1, mean2, size1, size2);
    }

    // merge setID2 into setID1
    void merge_sets(int setID1, int setID2, int newSetID) 
    {
        // initialize merged statistics
        std::set<int> combined_points;
        std::set<int> combined_edges;
        std::set<int> combined_triangles;
        std::set<int> combined_boundary_points;
        std::set<int> combined_boundary_edges;
        std::map<int, int> combined_edge_count_map;
        Eigen::Vector3d combined_mean;
        Eigen::Matrix3d combined_covariance_matrix;
        Eigen::Matrix3d combined_eigenvectors;
        Eigen::Vector3d combined_eigenvalues;
        Eigen::Vector3d combined_normal;
        std::array<int, 3> combined_color;
        
        // compute merged statistics
        // insert if setID exist
        if (set_to_points_map.find(setID1) != set_to_points_map.end()) combined_points.insert(set_to_points_map.at(setID1).begin(), set_to_points_map.at(setID1).end());
        if (set_to_points_map.find(setID2) != set_to_points_map.end()) combined_points.insert(set_to_points_map.at(setID2).begin(), set_to_points_map.at(setID2).end());
        if (set_to_edges_map.find(setID1) != set_to_edges_map.end()) combined_edges.insert(set_to_edges_map.at(setID1).begin(), set_to_edges_map.at(setID1).end());
        if (set_to_edges_map.find(setID2) != set_to_edges_map.end()) combined_edges.insert(set_to_edges_map.at(setID2).begin(), set_to_edges_map.at(setID2).end());
        if (set_to_triangles_map.find(setID1) != set_to_triangles_map.end()) combined_triangles.insert(set_to_triangles_map.at(setID1).begin(), set_to_triangles_map.at(setID1).end());
        if (set_to_triangles_map.find(setID2) != set_to_triangles_map.end()) combined_triangles.insert(set_to_triangles_map.at(setID2).begin(), set_to_triangles_map.at(setID2).end());
        combined_mean = merge_means_of_sets(setID1, setID2);
        combined_covariance_matrix = merge_covariances_of_sets(setID1, setID2);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(combined_covariance_matrix);
        combined_eigenvectors = solver.eigenvectors();
        combined_eigenvalues = solver.eigenvalues();
        combined_normal = combined_eigenvectors.col(0);
        combined_color = {rand() % 256, rand() % 256, rand() % 256};

        // remove old statistics
        set_to_points_map.erase(setID1);
        set_to_edges_map.erase(setID1);
        set_to_triangles_map.erase(setID1);
        set_to_mean_map.erase(setID1);
        set_to_covariance_matrix_map.erase(setID1);
        set_to_eigenvectors_map.erase(setID1);
        set_to_eigenvalues_map.erase(setID1);
        set_to_normal_map.erase(setID1);
        set_to_color_map.erase(setID1);
        set_list.erase(std::remove(set_list.begin(), set_list.end(), setID1), set_list.end());

        set_to_points_map.erase(setID2);
        set_to_edges_map.erase(setID2);
        set_to_triangles_map.erase(setID2);
        set_to_mean_map.erase(setID2);
        set_to_covariance_matrix_map.erase(setID2);
        set_to_eigenvectors_map.erase(setID2);
        set_to_eigenvalues_map.erase(setID2);
        set_to_normal_map.erase(setID2);
        set_to_color_map.erase(setID2);
        set_list.erase(std::remove(set_list.begin(), set_list.end(), setID2), set_list.end());

        // store merged statistics
        set_to_points_map.at(newSetID) = combined_points;
        set_to_edges_map.at(newSetID) = combined_edges;
        set_to_triangles_map.at(newSetID) = combined_triangles;
        set_to_mean_map.at(newSetID) = combined_mean;
        set_to_covariance_matrix_map.at(newSetID) = combined_covariance_matrix;
        set_to_eigenvectors_map.at(newSetID) = combined_eigenvectors;
        set_to_eigenvalues_map.at(newSetID) = combined_eigenvalues;
        set_to_normal_map.at(newSetID) = combined_normal;
        set_to_color_map.at(newSetID) = combined_color;
        for (int point_id : combined_points) point_to_set_map.at(point_id) = newSetID;
        for (int edge_id : combined_edges) edge_to_set_map.at(edge_id) = newSetID;
        for (int triangle_id : combined_triangles) triangle_to_set_map.at(triangle_id) = newSetID;
    }

    // try merge sets
    std::map<int, int> try_merge_sets(std::set<int>& sets_to_merge)
    {
        std::map<int, int> merge_map;

        while (true) 
        {
            // get all possible pairs to merge
            std::set<std::pair<int, int>> pairs_set;
            for (int setID1 : sets_to_merge) 
            {
                for (int setID2 : sets_to_merge) 
                {
                    if (setID1 >= setID2) continue;
                    pairs_set.insert(std::make_pair(setID1, setID2));
                }
            }
            
            // try to merge pairs
            bool again = false;
            for (const auto& pairs : pairs_set) 
            {
                // skip if can't merge
                Eigen::Matrix3d covariance_matrix = merge_covariances_of_sets(pairs.first, pairs.second);
                double eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance_matrix).eigenvalues()[0];
                if (eigenvalue > merged_eigenvalue_threshold) continue;

                // merge sets
                int newSetID = getNewSetID();
                merge_sets(pairs.first, pairs.second, newSetID);
                sets_to_merge.erase(pairs.first);
                sets_to_merge.erase(pairs.second);
                sets_to_merge.insert(newSetID);
                merge_map[pairs.first] = newSetID;
                merge_map[pairs.second] = newSetID;

                // once merged, restart
                again = true;
                break;
            }
            if (!again) break;
        }

        // return
        return merge_map;
    }

    void add_point_by_radius_search(int newPointID, Eigen::Vector3d thisPointVEC, Eigen::Vector3d thisPointOriginVEC)
    {
        // if empty, can not set up radius search, add point to new set
        if (boundary_point_set.size() == 0)
        {
            int newSetID = getNewSetID();
            add_point(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
            return;
        }

        // radius search size
        double search_size = 0.2;

        // perform kdtree radius search
        std::set<int> searched_boundary_points_set = kdtree.radiusSearch(thisPointVEC, search_size);

        // if no searched results, add point to new set
        if (searched_boundary_points_set.size() == 0)
        {
            int newSetID = getNewSetID();
            add_point(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
            return;
        }

        // from searched points identify neighboring sets
        std::set<int> neighboring_sets; 
        for (int point_id : searched_boundary_points_set)
        {
            neighboring_sets.insert(point_to_set_map.at(point_id));
        }

        // try merge neighboring sets
        std::map<int, int> merge_map = try_merge_sets(neighboring_sets);

        // split neighboring sets into sets with plane and sets without plane (by size)
        std::set<int> sets_with_plane;
        std::set<int> sets_without_plane;
        for (int setID : neighboring_sets)
        {
            if (set_to_points_map.at(setID).size() > fit_plane_threshold) sets_with_plane.insert(setID);
            else sets_without_plane.insert(setID);
        }
        
        // for sets with plane, compute the point to set intersection distance
        std::map<int, double> set_distance_map;
        for (int setID : sets_with_plane)
        {
            // use mean and normal of the set with the point added
            Eigen::Vector3d mean = merge_means_between_set_and_point(setID, thisPointVEC);
            Eigen::Matrix3d cov = merge_covariances_between_set_and_point(setID, thisPointVEC);

            // use eigen solver to get normal
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(cov);
            Eigen::Vector3d normal = eigensolver.eigenvectors().col(0);

            // compute distance
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(thisPointOriginVEC, thisPointVEC, mean, normal);
            double distance = (thisPointVEC - rayPlaneIntersectionPoint).norm();

            // store
            set_distance_map[setID] = distance;
        }

        // extract the set within distance threshold
        std::set<int> sets_within_threshold;
        for (const auto& pair : set_distance_map)
        {
            if (pair.second < distance_threshold) sets_within_threshold.insert(pair.first);
        }

        // from the sets within threshold, find the set that is closest to the point
        int closest_setID = -1;
        double closest_distance = std::numeric_limits<double>::max();
        for (int setID : sets_within_threshold)
        {
            // use mean and normal of the set with the point added
            Eigen::Vector3d mean = merge_means_between_set_and_point(setID, thisPointVEC);
            Eigen::Matrix3d cov = merge_covariances_between_set_and_point(setID, thisPointVEC);

            // use eigen solver to get normal
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(cov);
            Eigen::Vector3d normal = eigensolver.eigenvectors().col(0);

            // compute distance
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(thisPointOriginVEC, thisPointVEC, mean, normal);
            double distance = (thisPointVEC - rayPlaneIntersectionPoint).norm();

            // update if closer
            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest_setID = setID;
            }
        }
        // somewhere above, the global and local boundary points sets are not in sync
        if (closest_setID != -1 && closest_distance < distance_threshold)
        {
            std::set<int> searched_boundary_points_in_current_set = intersection_of_sets(searched_boundary_points_set, get_boundary_point_of_set(closest_setID));
            add_point(newPointID, closest_setID, thisPointVEC, thisPointOriginVEC);
            connect_point_to_set(newPointID, closest_setID, searched_boundary_points_in_current_set);
            return;
        }

        // else, find the set with the largest size
        int largest_setID = -1;
        std::size_t largest_size = 0;
        for (int setID : sets_without_plane)
        {
            if (set_to_points_map.at(setID).size() > largest_size)
            {
                largest_size = set_to_points_map.at(setID).size();
                largest_setID = setID;
            }
        }
        if (largest_setID != -1)
        {
            std::set<int> searched_boundary_points_in_current_set = intersection_of_sets(searched_boundary_points_set, get_boundary_point_of_set(largest_setID));
            add_point(newPointID, largest_setID, thisPointVEC, thisPointOriginVEC);
            connect_point_to_set(newPointID, largest_setID, searched_boundary_points_in_current_set);
            return;
        }

        // else, add the point to a new set
        int newSetID = getNewSetID();
        add_point(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
    }



    void step()
    {
        // if all points are processed
        if (ith_point >= ith_size) 
        {
            // next cloud
            ith_cloud += 1;
            ith_point = 0;

            // load data
            typename pcl::PointCloud<PointT>::Ptr pointcloud_local = data_loader.get_cloud(ith_cloud);
            Eigen::Affine3d pose = data_loader.get_pose(ith_cloud);
            pointcloud = transform_cloud_to_global<PointT>(pointcloud_local, pose);
            origin = pose.translation();

            ith_size = pointcloud->size();

            // shuffle the pointcloud
            std::random_shuffle(pointcloud->points.begin(), pointcloud->points.end());
        }

        // get point
        Eigen::Vector3d thisPointVEC = pointcloud->points[ith_point].getVector3fMap().cast<double>();
        Eigen::Vector3d thisPointOriginVEC = origin;
        ith_point ++;
        
        // get new point id
        int newPointID = getNewPointID();


        // ------------- add point by triangle intersection

        // get list of intersected triangle by the point
        std::set<int> searched_triangles = bvhRoot.intersectionSearch(thisPointOriginVEC, thisPointVEC); // may include deleted triangles
        // std::set<int> searched_triangles = intersection_of_sets(candidate_searched_triangles, global_triangle_set);

        // group the triangles by set
        std::map<int, std::set<int>> set_to_searched_triangle_map;
        for (int triangleID : searched_triangles)
        {
            int setID = triangle_to_set_map.at(triangleID);
            set_to_searched_triangle_map[setID].insert(triangleID);
        }

        // compute the intersection distance to the sets (distance measured in plane normal direction)
        std::map<int, double> set_distance_map;
        for (const auto& pair : set_to_searched_triangle_map)
        {
            int setID = pair.first;
            Eigen::Vector3d mean = set_to_mean_map.at(setID);
            Eigen::Vector3d normal = set_to_normal_map.at(setID);
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(thisPointOriginVEC, thisPointVEC, mean, normal);
            double distance = (thisPointVEC - rayPlaneIntersectionPoint).dot(normal);
            set_distance_map[setID] = distance;
        }

        // split the sets into three categories
        std::set<int> set_with_point_before_it;
        std::set<int> set_with_point_within_it;
        std::set<int> set_with_point_behind_it;
        double split_distance_threshold = distance_threshold;
        for (const auto& pair : set_distance_map)
        {
            int setID = pair.first;
            double distance = pair.second;
            if (distance > split_distance_threshold) 
            {
                set_with_point_before_it.insert(setID);
            }
            else if (distance < -split_distance_threshold) 
            {
                set_with_point_behind_it.insert(setID);
            }
            else 
            {
                set_with_point_within_it.insert(setID);
            }
        }
        
        // process point behind set set

        // get the set of triangles that are penetrated by the point
        std::set<int> penetrated_triangles;
        for (int setID : set_with_point_behind_it)
        {
            penetrated_triangles.insert(set_to_searched_triangle_map.at(setID).begin(), set_to_searched_triangle_map.at(setID).end());
        }

        // collect the list of points that are within the penetrated triangles, and isolated by the triangle
        for (int triangleID : penetrated_triangles)
        {
            // penetrate triangles
            remove_triangle(triangleID);
        }

        // re add the list points by radius search
        // todo - while avoid covering the new point
        // process point in the queue free_points_vector3d_and_origin_vector3d
        while (!free_points_queue.empty())
        {
            std::pair<Eigen::Vector3d, Eigen::Vector3d> pair = free_points_queue.front(); free_points_queue.pop();

            std::cout << "---------------------------------------------------------------------------------------------------- re-adding point by radius search" << std::endl;
            Eigen::Vector3d pointVEC = pair.first;
            Eigen::Vector3d pointOriginVEC = pair.second;
            int pointID = getNewPointID();
            add_point_by_radius_search(pointID, pointVEC, pointOriginVEC);
        }
        
        // process point within set set
        bool point_added_to_set = false;

        // try merge them
        std::set<int> merged_set = set_with_point_within_it;
        std::map<int, int> merge_map = try_merge_sets(merged_set);

        // update the set_to_searched_triangle_map using the merge map
        while (true)
        {
            bool restart = false;

            for (const auto& pair : merge_map)
            {
                // skip if the key is already handled
                if (set_to_searched_triangle_map.find(pair.first) == set_to_searched_triangle_map.end()) continue;

                // change the key
                std::set<int> content = set_to_searched_triangle_map.at(pair.first);
                set_to_searched_triangle_map.erase(pair.first);
                set_to_searched_triangle_map[pair.second] = content;

                // restart
                restart = true;
                break;
            }

            if (!restart) break;
        }
        

        // find the set with the smallest distance
        int smallest_setID = -1;
        double smallest_distance = std::numeric_limits<double>::max();
        for (int setID : merged_set)
        {
            // if the set does not have distance, it is new, compute its distance
            if (set_distance_map.find(setID) == set_distance_map.end())
            {
                Eigen::Vector3d mean = set_to_mean_map.at(setID);
                Eigen::Vector3d normal = set_to_normal_map.at(setID);
                Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(thisPointOriginVEC, thisPointVEC, mean, normal);
                double distance = (thisPointVEC - rayPlaneIntersectionPoint).dot(normal);
                set_distance_map[setID] = distance;
            }

            // update if smaller
            if (std::abs(set_distance_map.at(setID)) < smallest_distance)
            {
                smallest_distance = set_distance_map.at(setID);
                smallest_setID = setID;
            }
        }

        // if the smallest set is within threshold, add the point to the set
        if (smallest_setID != -1 && std::abs(smallest_distance) < distance_threshold)
        {
            add_point(newPointID, smallest_setID, thisPointVEC, thisPointOriginVEC);
            std::cout << ith_point << " / " << pointcloud->size() << " of pointcloud " << ith_cloud << " added to set " << smallest_setID << std::endl;

            // add the point to the triangle (if multiple, add to first)
            int triangleID = *set_to_searched_triangle_map.at(smallest_setID).begin();
            triangle_to_points_map.at(triangleID).insert(newPointID);

            point_added_to_set = true;
        }

        if (!point_added_to_set)
        {
            add_point_by_radius_search(newPointID, thisPointVEC, thisPointOriginVEC);
            std::cout << ith_point << " / " << pointcloud->size() << " of pointcloud " << ith_cloud << " added by radius search" << std::endl;
        }


        // todo - process point before set set

    }

    void loop()
    {
        step();

        // finish all points in step
        while (ith_point < ith_size)
        {
            step();
        }
    }

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
        data_loader.load_dataset(pcd_file_folder, pose_file_path);
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
        for (int point_id : point_list)
        {
            pcl::PointXYZRGB point;
            point.x = point_to_vector3d_map.at(point_id)[0];
            point.y = point_to_vector3d_map.at(point_id)[1];
            point.z = point_to_vector3d_map.at(point_id)[2];
            point.r = set_to_color_map.at(point_to_set_map.at(point_id)).at(0);
            point.g = set_to_color_map.at(point_to_set_map.at(point_id)).at(1);
            point.b = set_to_color_map.at(point_to_set_map.at(point_id)).at(2);
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_distance_cloud()
    {
        compute_projected_point_to_vector3d_map();

        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (int point_id : point_list)
        {
            pcl::PointXYZRGB point;
            point.x = point_to_vector3d_map.at(point_id)[0];
            point.y = point_to_vector3d_map.at(point_id)[1];
            point.z = point_to_vector3d_map.at(point_id)[2];
            double value = projected_points_distance_map.at(point_id) / 0.05;
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
        compute_projected_point_to_vector3d_map();

        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (int point_id : point_list)
        {
            pcl::PointXYZRGB point;
            point.x = projected_points_to_vector3d_map.at(point_id)[0];
            point.y = projected_points_to_vector3d_map.at(point_id)[1];
            point.z = projected_points_to_vector3d_map.at(point_id)[2];
            point.r = set_to_color_map.at(point_to_set_map.at(point_id)).at(0);
            point.g = set_to_color_map.at(point_to_set_map.at(point_id)).at(1);
            point.b = set_to_color_map.at(point_to_set_map.at(point_id)).at(2);
            cloud->push_back(point);
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_distance_cloud()
    {
        compute_projected_point_to_vector3d_map();

        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (int point_id : point_list)
        {
            pcl::PointXYZRGB point;
            point.x = projected_points_to_vector3d_map.at(point_id)[0];
            point.y = projected_points_to_vector3d_map.at(point_id)[1];
            point.z = projected_points_to_vector3d_map.at(point_id)[2];
            double value = projected_points_distance_map.at(point_id) / 0.05;
            std::tuple<int, int, int> color = valueToJet(value);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
        }
        return cloud;
    }
    
    void compute_projected_point_to_vector3d_map()
    {
        // for each point
        for (int point_id : point_list)
        {
            // get set id
            int setID = point_to_set_map.at(point_id);

            // get projected point
            Eigen::Vector3d rayOrigin = point_to_origin_vector3d_map.at(point_id);
            Eigen::Vector3d rayEndPoint = point_to_vector3d_map.at(point_id);
            Eigen::Vector3d mean = set_to_mean_map.at(setID);
            Eigen::Vector3d normal = set_to_normal_map.at(setID);
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);

            // compute distance
            double distance = (point_to_vector3d_map.at(point_id) - rayPlaneIntersectionPoint).norm();

            // store
            projected_points_to_vector3d_map[point_id] = rayPlaneIntersectionPoint;   
            projected_points_distance_map[point_id] = distance;
        }
    }

    int get_number_of_triangles()
    {
        return triangle_to_vertices_map.size();
    }

    void change_color()
    {
        for (int setID : set_list)
        {
            set_to_color_map.at(setID) = {rand() % 256, rand() % 256, rand() % 256};
        }
    }

    void kdtree_rebuild()
    {
        kdtree.rebuild();
    }

    void kdtree_print_tree()
    {
        kdtree.print_tree();
    }
    
private:
    // data
    DataLoader<VilensPointT> data_loader;
    typename pcl::PointCloud<VilensPointT>::Ptr pointcloud;
    Eigen::Vector3d origin;

    int ith_cloud = -1;
    std::size_t ith_point = -1;
    std::size_t ith_size = -2;
    

    // settings
    double distance_threshold = 0.05;
    std::size_t fit_plane_threshold = 10; // may cause error if below 3
    double merged_eigenvalue_threshold = 15e-5;

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
    KDTree kdtree;
    std::set<int> boundary_point_set;
    std::set<int> boundary_edge_set;
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
        std::map<int, std::array<int, 3>> triangle_to_cloud_indices_map = app_.get_triangle_to_cloud_indices_map();
        std::map<int, std::array<int, 2>> edge_to_cloud_indices_map = app_.get_edge_to_cloud_indices_map();
        std::set<int> boundary_edge_set = app_.get_boundary_edge_set();

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
            for (const auto& pair : triangle_to_cloud_indices_map)
            {
                pcl::Vertices triangle;
                triangle.vertices.push_back(pair.second[0]);
                triangle.vertices.push_back(pair.second[1]);
                triangle.vertices.push_back(pair.second[2]);
                triangle_mesh.polygons.push_back(triangle);
            }
            viewer_->addPolygonMesh(triangle_mesh, "triangle_mesh");
        }

        // boundary edges
        viewer_->removeShape("boundary_edges");        
        if (show_edge)
        {
            pcl::PolygonMesh boundary_mesh;
            pcl::toPCLPointCloud2(*point_cloud, boundary_mesh.cloud);
            for (int edge_id : boundary_edge_set)
            {
                pcl::Vertices edge;
                edge.vertices.push_back(edge_to_cloud_indices_map.at(edge_id)[0]);
                edge.vertices.push_back(edge_to_cloud_indices_map.at(edge_id)[1]);
                boundary_mesh.polygons.push_back(edge);
            }
            viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "boundary_edges");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "boundary_edges");
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
        // std::cout << "key pressed: " << event.getKeySym() << std::endl;

        if (event.getKeySym() == "space" && event.keyDown())
        {
            app_.step();
            app_.step();
            app_.step();
            app_.step();
            app_.step();
            app_.step();
            app_.step();
            app_.step();
            app_.step();
            app_.step();
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
            // rebuild kdtree
            app_.kdtree_rebuild();
            std::cout << "kdtree rebuilt" << std::endl;
        }
        if (event.getKeySym() == "n" && event.keyDown())
        {
            // print tree
            app_.kdtree_print_tree();
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