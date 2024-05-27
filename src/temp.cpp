#include "eye_patch/DataLoader.hpp"

// #include "eye_patch/Algorithm.hpp"
#include <pcl/common/transforms.h>
#include <pcl/kdtree/kdtree_flann.h>

// #include "eye_patch/Visualization.hpp"
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

// #include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"

class TriangleBVH
{
private:

    // https://chatgpt.com/share/96c43118-6cb4-4549-9478-2725dee3b44d
    struct BoundingBox 
    {
        Eigen::Vector3d min = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
        Eigen::Vector3d max = Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());

        BoundingBox(){};

        void expand(const Eigen::Vector3d& point) 
        {
            min = min.cwiseMin(point);
            max = max.cwiseMax(point);
        }

        bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir, double& tMin, double& tMax) const 
        {
            for (int i = 0; i < 3; ++i) 
            {
                double invD = 1.0 / dir[i];
                double t0 = (min[i] - orig[i]) * invD;
                double t1 = (max[i] - orig[i]) * invD;
                if (invD < 0.0) std::swap(t0, t1);
                tMin = std::max(tMin, t0);
                tMax = std::min(tMax, t1);
                if (tMax <= tMin) return false;
            }
            return true;
        }

        bool intersect(const Eigen::Vector3d& orig, const Eigen::Vector3d& dir) const 
        {
            double tMin = -std::numeric_limits<double>::infinity();
            double tMax = std::numeric_limits<double>::infinity();
            return intersect(orig, dir, tMin, tMax);
        }

        int get_longest_axis()
        {
            // compute
            Eigen::Vector3d diagonal_line = max - min;
            int axis = 0;
            if (diagonal_line[1] > diagonal_line[axis]) axis = 1;
            if (diagonal_line[2] > diagonal_line[axis]) axis = 2;

            // return
            return axis;
        }
    };

    struct Node 
    {
        // bounding box, for intersection test
        BoundingBox box;

        // as branch node
        double split_value;
        int split_axis;
        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;
        
        // as leaf node
        bool isLeaf() const {return !left && !right;}
        std::set<int> triangleIDs;
    };

    Eigen::Vector3d compute_triangle_center(int triangleID)
    {
        int point0_id = triangle_to_indices_map.at(triangleID)[0];
        int point1_id = triangle_to_indices_map.at(triangleID)[1];
        int point2_id = triangle_to_indices_map.at(triangleID)[2];
        Eigen::Vector3d v0 = point_to_vector3d_map.at(point0_id);
        Eigen::Vector3d v1 = point_to_vector3d_map.at(point1_id);
        Eigen::Vector3d v2 = point_to_vector3d_map.at(point2_id);
        return (v0 + v1 + v2) / 3.0;
    }

    // nth_element sort
    double sort_triangle_list_in_axis(std::vector<int>& triangle_list, int axis, int start, int mid, int end)
    {
        // partial sort
        std::nth_element(triangle_list.begin() + start, triangle_list.begin() + mid, triangle_list.begin() + end, 
            [&](const int& triangle_a, const int& triangle_b) 
            {
                // compare center in axis
                Eigen::Vector3d centerA = triangle_to_center_vector3d_map.at(triangle_a);
                Eigen::Vector3d centerB = triangle_to_center_vector3d_map.at(triangle_b);
                return centerA[axis] < centerB[axis];
            });

        // return split value
        return triangle_to_center_vector3d_map.at(triangle_list[mid])[axis];
    }

    void expand_node_box(std::shared_ptr<Node> node, int triangle_id)
    {
        // expand
        Eigen::Vector3d v0 = point_to_vector3d_map.at(triangle_to_indices_map.at(triangle_id)[0]);
        Eigen::Vector3d v1 = point_to_vector3d_map.at(triangle_to_indices_map.at(triangle_id)[1]);
        Eigen::Vector3d v2 = point_to_vector3d_map.at(triangle_to_indices_map.at(triangle_id)[2]);
        node->box.expand(v0);
        node->box.expand(v1);
        node->box.expand(v2);
    }
    
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
    bool rayTriangleIntersect(
        const Eigen::Vector3d& orig, const Eigen::Vector3d& dir,
        const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
        Eigen::Vector3d& outIntersection)
    {
        const double EPSILON = 1e-8;
        Eigen::Vector3d edge1 = v1 - v0;
        Eigen::Vector3d edge2 = v2 - v0;
        
        // compute determinant
        Eigen::Vector3d pvec = dir.cross(edge2);
        double det = edge1.dot(pvec);
        if (std::fabs(det) < EPSILON) {
            return false; // This ray is parallel to this triangle.
        }

        // compute inverse determinant
        double invDet = 1.0 / det;

        // compute u
        Eigen::Vector3d tvec = orig - v0;
        double u = tvec.dot(pvec) * invDet;
        if (u < 0.0 || u > 1.0) return false;
        
        // compute v
        Eigen::Vector3d qvec = tvec.cross(edge1);
        double v = dir.dot(qvec) * invDet;
        if (v < 0.0 || u + v > 1.0) return false;

        // compute t
        double t = edge2.dot(qvec) * invDet;
        if (t < EPSILON) return false; // This means that there is a line intersection but not a ray intersection.

        // comptue intersection
        outIntersection = orig + dir * t;
        return true;
    }

    // return a list of triangle ids that could intersect with the ray
    std::set<int> intersectHierarchy(const std::shared_ptr<Node>& node, Eigen::Vector3d orig, Eigen::Vector3d dir) 
    {
        // check intersection with current box
        bool intersected = node->box.intersect(orig, dir);
        if (!intersected) return std::set<int>();
        
        // if leaf node
        if (node->isLeaf())
        {
            return node->triangleIDs;
        }
        // if branch node
        else
        {
            std::set<int> triangles_left = intersectHierarchy(node->left, orig, dir);
            std::set<int> triangles_right = intersectHierarchy(node->right, orig, dir);
            triangles_left.insert(triangles_right.begin(), triangles_right.end());
            return triangles_left;
        }
    }

    void convert_leaf_to_branch(std::shared_ptr<Node> node)
    {
        // sort current triangle list along longest axis
        std::vector<int> triangle_list(node->triangleIDs.begin(), node->triangleIDs.end());
        int start = 0;
        int end = triangle_list.size();
        int mid = (start + end) / 2;
        int axis = node->box.get_longest_axis();
        double split_value = sort_triangle_list_in_axis(triangle_list, axis, start, mid, end);
        
        // build children nodes
        node->triangleIDs.clear();
        node->split_axis = axis;
        node->split_value = split_value;
        node->left = build_node(std::vector<int>(triangle_list.begin(), triangle_list.begin() + mid));
        node->right = build_node(std::vector<int>(triangle_list.begin() + mid, triangle_list.end()));
    }

    std::shared_ptr<Node> build_node(std::vector<int> triangle_list)
    {
        // build current node
        auto node = std::make_shared<Node>();
        for (int triangle_id : triangle_list) 
        {
            expand_node_box(node, triangle_id);
            node->triangleIDs.insert(triangle_id);
        }

        // convert to branch node if too many triangles
        if (node->triangleIDs.size() > 4) convert_leaf_to_branch(node);

        // return
        return node;
    }
    
    void addTriangleToNode(std::shared_ptr<Node> node, int triangleID)
    {
        // leaf node
        if (node->isLeaf())
        {
            // expand box
            expand_node_box(node, triangleID);

            // add to leaf
            node->triangleIDs.insert(triangleID);

            // convert to branch if too many triangles
            if (node->triangleIDs.size() > 4) convert_leaf_to_branch(node);
        }
        // branch node
        else
        {
            // expand box
            expand_node_box(node, triangleID);

            // process to left or right 
            // (since using stored split axis and value instead of rebalancing, tree efficency is lowered)
            if (triangle_to_center_vector3d_map.at(triangleID)[node->split_axis] < node->split_value)
            {
                addTriangleToNode(node->left, triangleID);
            }
            else 
            {
                addTriangleToNode(node->right, triangleID);
            }
        }
    }

    // settings
    double rebuild_threshold = 2;
    
    // data
    std::vector<int> triangle_list;
    std::map<int, std::array<int, 3>> triangle_to_indices_map;
    std::map<int, Eigen::Vector3d> point_to_vector3d_map;
    std::map<int, Eigen::Vector3d> triangle_to_center_vector3d_map;
    std::shared_ptr<Node> root;
    int size_at_last_rebuild = 0;


public:
    TriangleBVH()
    {
        rebuild();
    }

    void addData(std::vector<int> _triangle_list, std::map<int, std::array<int, 3>> _triangle_to_indices_map, std::map<int, Eigen::Vector3d> _point_to_vector3d_map)
    {
        // store
        triangle_list = _triangle_list;
        triangle_to_indices_map = _triangle_to_indices_map;
        point_to_vector3d_map = _point_to_vector3d_map;
        for (int triangleID : triangle_list)
        {
            triangle_to_center_vector3d_map[triangleID] = compute_triangle_center(triangleID);
        }
    }

    void rebuild()
    {
        root = build_node(triangle_list);
    }
    
    void addTriangle(int triangleID, std::array<int, 3> indices, Eigen::Vector3d v0, Eigen::Vector3d v1, Eigen::Vector3d v2)
    {
        // add triangle data
        triangle_list.push_back(triangleID);
        triangle_to_indices_map[triangleID] = indices;
        point_to_vector3d_map[indices[0]] = v0;
        point_to_vector3d_map[indices[1]] = v1;
        point_to_vector3d_map[indices[2]] = v2;
        triangle_to_center_vector3d_map[triangleID] = (v0 + v1 + v2) / 3.0;

        // if rebuild threshold reached
        if (triangle_list.size() > size_at_last_rebuild * rebuild_threshold)
        {
            // rebuild
            rebuild();
            size_at_last_rebuild = triangle_list.size();
        }
        else
        {
            // add to existing
            addTriangleToNode(root, triangleID);
        }
    }

    std::set<int> intersectionSearch(Eigen::Vector3d origin, Eigen::Vector3d endPoint)
    {
        // initialize
        std::set<int> intersected_triangleIDs;

        // process
        Eigen::Vector3d dir = (endPoint - origin).normalized();
        std::set<int> triangleIDs = intersectHierarchy(root, origin, dir);
        for (int triangleID : triangleIDs)
        {
            int point0_id = triangle_to_indices_map.at(triangleID)[0];
            int point1_id = triangle_to_indices_map.at(triangleID)[1];
            int point2_id = triangle_to_indices_map.at(triangleID)[2];
            Eigen::Vector3d v0 = point_to_vector3d_map.at(point0_id);
            Eigen::Vector3d v1 = point_to_vector3d_map.at(point1_id);
            Eigen::Vector3d v2 = point_to_vector3d_map.at(point2_id);
            Eigen::Vector3d intersection;
            if (rayTriangleIntersect(origin, dir, v0, v1, v2, intersection)) intersected_triangleIDs.insert(triangleID);
        }

        // return
        return intersected_triangleIDs;
    }
};

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
            flann_data_storage.push_back(point_to_vector3d_map.at(point_id)[0]);
            flann_data_storage.push_back(point_to_vector3d_map.at(point_id)[1]);
            flann_data_storage.push_back(point_to_vector3d_map.at(point_id)[2]);

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
        flann_tree.addPoints(flann::Matrix<double>(flann_data_storage.data() + flann_data_storage.size() - 3, 1, 3)); // addPoints introduces randomness to radiusSearch, since flann tree randomize something
        // flann_tree.buildIndex(flann::Matrix<double>(flann_data_storage.data(), flann_data_storage.size() / 3, 3));

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

    // want boundary edge and points to prevent overlapping triangles
    // to add a new point to mesh
    // 1. extract mesh boundary points
    // 2. perform delaunay triangulation using the boundary points and the new point
    // 3. add the edge and triangles connected to the new point to the original mesh


    Eigen::Vector3d merge_means(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2) 
    {
        return (size1 * mean1 + size2 * mean2) / (size1 + size2);
    }

    Eigen::Matrix3d merge_covariances(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                    const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                    int size1, int size2) 
    {
        Eigen::Vector3d combined_mean = merge_means(mean1, mean2, size1, size2);
        Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
        Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
        Eigen::Matrix3d combined_covariance = (size1 * cov1 + size2 * cov2 + size1 * mean_diff1 + size2 * mean_diff2) / (size1 + size2);

        return combined_covariance;
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

    std::array<int, 3> sortThreeInts(int a, int b, int c) {
        std::array<int, 3> sorted = {a, b, c};
        if (sorted[0] > sorted[1]) std::swap(sorted[0], sorted[1]);
        if (sorted[0] > sorted[2]) std::swap(sorted[0], sorted[2]);
        if (sorted[1] > sorted[2]) std::swap(sorted[1], sorted[2]);
        return sorted;
    }

    std::array<int, 2> sortTwoInts(int a, int b) {
        std::array<int, 2> sorted = {a, b};
        if (sorted[0] > sorted[1]) std::swap(sorted[0], sorted[1]);
        return sorted;
    }

    void add_point_to_set(int newPointID, int setID, Eigen::Vector3d thisPoint, Eigen::Vector3d origin)
    {
        // add eigen
        point_to_vector3d_map.at(newPointID) = thisPoint;
        point_to_origin_vector3d_map.at(newPointID) = origin;

        // add point to set
        point_to_set_map.at(newPointID) = setID;

        // if first point to a set
        if (set_to_points_map.find(setID) == set_to_points_map.end())
        {
            // plane statistics
            set_to_mean_map.at(setID) = thisPoint;
            set_to_covariance_matrix_map.at(setID) = Eigen::Matrix3d::Zero();
            set_to_eigenvectors_map.at(setID) = Eigen::Matrix3d::Identity();
            set_to_eigenvalues_map.at(setID) = Eigen::Vector3d::Zero();
            set_to_normal_map.at(setID) = Eigen::Vector3d(0, 0, 1);

            // set to point
            set_to_points_map.at(setID) = {newPointID};
        }
        else
        {
            int size1 = set_to_points_map.at(setID).size();
            int size2 = 1;
            Eigen::Vector3d mean1 = set_to_mean_map.at(setID);
            Eigen::Vector3d mean2 = thisPoint;
            Eigen::Vector3d new_mean = merge_means(mean1, mean2, size1, size2);
            Eigen::Matrix3d cov1 = set_to_covariance_matrix_map.at(setID);
            Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d new_cov = merge_covariances(cov1, cov2, mean1, mean2, size1, size2);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(new_cov);
            Eigen::Matrix3d eigenvectors = solver.eigenvectors();
            Eigen::Vector3d eigenvalues = solver.eigenvalues();
            Eigen::Vector3d normal = eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal
            // normal should points towards the origin
            if (normal.dot(origin - thisPoint) < 0) normal *= -1;
            set_to_mean_map.at(setID) = new_mean;
            set_to_covariance_matrix_map.at(setID) = new_cov;
            set_to_eigenvectors_map.at(setID) = eigenvectors;
            set_to_eigenvalues_map.at(setID) = eigenvalues;
            set_to_normal_map.at(setID) = normal;

            // set to point
            set_to_points_map.at(setID).insert(newPointID);
        }
    }

    void add_edge(int newEdgeID, int setID, std::array<int, 2> newEdge)
    {
        // sort edge
        newEdge = sortTwoInts(newEdge[0], newEdge[1]);

        // update edge to point 
        edge_to_point_map.at(newEdgeID) = newEdge;
        edge_to_point_map_reverse[newEdge] = newEdgeID;

        // update point to edge
        point_to_edges_map.at(newEdge[0]).insert(newEdgeID);
        point_to_edges_map.at(newEdge[1]).insert(newEdgeID);

        // update set to edge
        set_to_edges_map.at(setID).insert(newEdgeID);
        edge_to_set_map.at(newEdgeID) = setID;

        // update set boundary edge
        set_to_edge_count_map.at(setID)[newEdgeID] = 0;
        set_to_boundary_edge_set_map.at(setID).insert(newEdgeID);
    
        // update boundary points
        add_boundary_point(newEdge[0], setID);
        add_boundary_point(newEdge[1], setID);
    }

    void add_triangle(int newTriangleID, int setID, std::array<int, 3> vertices)
    {
        // sort vertices
        vertices = sortThreeInts(vertices[0], vertices[1], vertices[2]);

        // add triangle to vertices map
        triangle_to_vertices_map.at(newTriangleID) = vertices;
        triangle_to_vertices_map_reverse[vertices] = newTriangleID;

        // add triangle to set map
        set_to_triangles_map.at(setID).insert(newTriangleID);
        triangle_to_set_map.at(newTriangleID) = setID;

        // update boundary edge and points
        for (int i = 0; i < 3; i++)
        {
            // get correct edge order
            int smaller_id = std::min(vertices[i], vertices[(i + 1) % 3]);
            int larger_id = std::max(vertices[i], vertices[(i + 1) % 3]);
            std::array<int, 2> edge = {smaller_id, larger_id};

            // increment count
            int& count = set_to_edge_count_map.at(setID)[edge_to_point_map_reverse[edge]];
            count ++;

            // cases
            if (count <= 1)
            {
                // udpate set boundary edge
                set_to_boundary_edge_set_map.at(setID).insert(edge_to_point_map_reverse[edge]);
                
                // update boundary points
                add_boundary_point(edge[0], setID);
                add_boundary_point(edge[1], setID);
            } 
            else
            {
                // remove from set boundary edge
                set_to_boundary_edge_set_map.at(setID).erase(edge_to_point_map_reverse[edge]);

                // remove from boundary points (if the edge point no longer connects to a boundary edge, it is removed from boundary points)
                remove_boundary_point(edge[0], setID);
                remove_boundary_point(edge[1], setID);
            }
        }

        // add to bvh
        bvhRoot.addTriangle(newTriangleID, vertices, point_to_vector3d_map.at(vertices[0]), point_to_vector3d_map.at(vertices[1]), point_to_vector3d_map.at(vertices[2]));
    }

    void add_boundary_point(int pointID, int setID)
    {
        // update local boundary points
        set_to_boundary_point_set_map.at(setID).insert(pointID);

        // update global boundary points
        bool inserted = global_boundary_point_set.insert(pointID).second;
        if (inserted) flann.addPoint(point_to_vector3d_map.at(pointID), pointID);
    }

    void remove_boundary_point(int pointID, int setID)
    {
        // skip if still boundary
        for (int edge_id : point_to_edges_map.at(pointID))
        {
            if (set_to_boundary_edge_set_map.at(setID).find(edge_id) != set_to_boundary_edge_set_map.at(setID).end()) return;
        }

        // remove from both local and global boundary points (a point can only be in one set)
        set_to_boundary_point_set_map.at(setID).erase(pointID);
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

        // initialize
        point_list.push_back(newPointID);
        point_to_origin_vector3d_map[newPointID] = Eigen::Vector3d(999, 999, 999);
        point_to_intersection_vector3d_map[newPointID] = Eigen::Vector3d(999, 999, 999);
        point_to_vector3d_map[newPointID] = Eigen::Vector3d(999, 999, 999);
        point_to_set_map[newPointID] = -999;
        point_to_edges_map[newPointID] = {};

        // return
        return newPointID;
    }

    int getNewSetID()
    {
        // get id
        int newSetID = next_set_id;
        next_set_id ++;

        // initialize
        set_list.push_back(newSetID);
        set_to_color_map[newSetID] = {rand() % 256, rand() % 256, rand() % 256};
        set_to_points_map[newSetID] = {};
        set_to_edges_map[newSetID] = {};
        set_to_triangles_map[newSetID] = {};
        set_to_boundary_edge_set_map[newSetID] = {};
        set_to_boundary_point_set_map[newSetID] = {};
        set_to_edge_count_map[newSetID] = {};
        set_to_mean_map[newSetID] = Eigen::Vector3d::Zero();
        set_to_covariance_matrix_map[newSetID] = Eigen::Matrix3d::Zero();
        set_to_eigenvectors_map[newSetID] = Eigen::Matrix3d::Identity();
        set_to_eigenvalues_map[newSetID] = Eigen::Vector3d::Zero();
        set_to_normal_map[newSetID] = Eigen::Vector3d(0, 0, 1);

        // return
        return newSetID;
    }

    int getNewEdgeID()
    {
        // get id
        int newEdgeID = next_edge_id;
        next_edge_id ++;

        // initialize (skip edge_to_point_map_reverse as can't be initialized)
        edge_list.push_back(newEdgeID);
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
        triangle_list.push_back(newTriangleID);
        triangle_to_vertices_map[newTriangleID] = {};
        triangle_to_set_map[newTriangleID] = -999;
        triangle_to_points_map[newTriangleID] = {};

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
            points_grouped_by_set[point_to_set_map.at(point_id)].insert(point_id);
        }

        // return
        return points_grouped_by_set;
    }

    std::set<int> compute_boundary_edge_set_of_set(int setID)
    {
        // intialize edge occurrences count
        std::map<int, int> edge_count;
        for (int edge_id : set_to_edges_map.at(setID))
        {
            edge_count[edge_id] = 0;
        }

        // for each triangle, increment edge count
        for (int triangle_id : set_to_triangles_map.at(setID))
        {
            std::array<int, 3> vertices = triangle_to_vertices_map.at(triangle_id);
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
            std::array<int, 2> vertices = edge_to_point_map.at(edge_id);
            boundary_points_set.insert(vertices[0]);
            boundary_points_set.insert(vertices[1]);
        }
        // add points that do not have edge in the same set
        for (int point_id : set_to_points_map.at(setID))
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
                for (int edge_id : point_to_edges_map.at(point_id))
                {
                    if (set_to_edges_map.at(setID).find(edge_id) != set_to_edges_map.at(setID).end())
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
            std::array<int, 2> vertices = edge_to_point_map.at(edge_id);
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
        if (set_to_boundary_edge_set_map.find(setID) == set_to_boundary_edge_set_map.end())
        {
            return std::set<int>();
        }
        else
        {
            return set_to_boundary_edge_set_map.at(setID);
        }
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
        return set_to_boundary_point_set_map.at(setID);
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
            mean += point_to_vector3d_map.at(point_id);
        }
        mean /= point_ids.size();

        // compute covariance matrix
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (int point_id : point_ids)
        {
            Eigen::Vector3d point = point_to_vector3d_map.at(point_id);
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
        for (int point_id : set_to_points_map.at(setID))
        {
            mean += point_to_vector3d_map.at(point_id);
        }
        mean /= set_to_points_map.at(setID).size();

        // return
        return mean;
    }

    Eigen::Matrix3d  get_set_eigenvectors(int setID)
    {
        // mean
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (int point_id : set_to_points_map.at(setID))
        {
            mean += point_to_vector3d_map.at(point_id);
        }
        mean /= set_to_points_map.at(setID).size();

        // covariance matrix
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (int point_id : set_to_points_map.at(setID))
        {
            Eigen::Vector3d point = point_to_vector3d_map.at(point_id);
            covariance_matrix += (point - mean) * (point - mean).transpose();
        }
        covariance_matrix /= set_to_points_map.at(setID).size();

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
    


    // convert to edge vertices set
    std::set<std::array<int, 2>> convert_to_edge_vertices_set(std::set<int> edge_set)
    {
        // initialize
        std::set<std::array<int, 2>> edge_vertices_set;

        // process
        for (int edge_id : edge_set)
        {
            edge_vertices_set.insert(edge_to_point_map.at(edge_id));
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
            Eigen::Vector3d point = point_to_vector3d_map.at(point_id);
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
            Eigen::Vector3d rayOrigin = point_to_origin_vector3d_map.at(point_id);
            Eigen::Vector3d rayEndPoint = point_to_vector3d_map.at(point_id);
            Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);

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
        for (int point_id : set_to_boundary_point_set_map.at(setID))
        {
            // // update if not available (todo - update if not accurate enough)
            // if (point_to_intersection_vector3d_map.find(point_id) == point_to_intersection_vector3d_map.end())
            // {
            //     point_to_intersection_vector3d_map.at(point_id) = point_set_intersection(point_id, setID);
            // }

            // get intersection point
            Eigen::Vector3d rayPlaneIntersectionPoint = point_set_intersection(point_id, setID);

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
    

    // ray plane intersection
    Eigen::Vector3d ray_plane_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayEndPoint, Eigen::Vector3d planeMean, Eigen::Vector3d planeNormal)
    {   
        // if perpendicular, return NaN
        Eigen::Vector3d rayDirection = (rayEndPoint - rayOrigin).normalized();
        if (planeNormal.dot(rayDirection) == 0)
        {
            // terminate programe by throwing an error
            throw std::invalid_argument("Ray and plane are perpendicular");
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
        Eigen::Vector3d rayOrigin = point_to_origin_vector3d_map.at(pointID);
        Eigen::Vector3d rayEndPoint = point_to_vector3d_map.at(pointID);
        Eigen::Vector3d mean = set_to_mean_map.at(setID);
        Eigen::Vector3d normal = set_to_normal_map.at(setID);
        Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);

        // return
        return rayPlaneIntersectionPoint;
    }

    // ray set intersection 
    Eigen::Vector3d ray_set_intersection(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayEndPoint, int setID)
    {
        // mean
        std::set<int> point_ids = set_to_points_map.at(setID);
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (int point_id : point_ids)
        {
            mean += point_to_vector3d_map.at(point_id);
        }
        mean /= point_ids.size();

        // covariance matrix
        Eigen::Matrix3d covariance_matrix = Eigen::Matrix3d::Zero();
        for (int point_id : point_ids)
        {
            Eigen::Vector3d point = point_to_vector3d_map.at(point_id);
            covariance_matrix += (point - mean) * (point - mean).transpose();
        }

        // normal
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance_matrix);
        Eigen::Vector3d normal = eigensolver.eigenvectors().col(0);

        // correct normal direction
        Eigen::Vector3d rayDirection = (rayEndPoint - rayOrigin).normalized();
        if (normal.dot(rayDirection) > 0)
        {
            normal = -normal;
        }

        // intersection point
        Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);

        // return
        return rayPlaneIntersectionPoint;
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
        if (set_to_boundary_point_set_map.find(setID1) != set_to_boundary_point_set_map.end()) combined_boundary_points.insert(set_to_boundary_point_set_map.at(setID1).begin(), set_to_boundary_point_set_map.at(setID1).end());
        if (set_to_boundary_point_set_map.find(setID2) != set_to_boundary_point_set_map.end()) combined_boundary_points.insert(set_to_boundary_point_set_map.at(setID2).begin(), set_to_boundary_point_set_map.at(setID2).end());
        if (set_to_boundary_edge_set_map.find(setID1) != set_to_boundary_edge_set_map.end()) combined_boundary_edges.insert(set_to_boundary_edge_set_map.at(setID1).begin(), set_to_boundary_edge_set_map.at(setID1).end());
        if (set_to_boundary_edge_set_map.find(setID2) != set_to_boundary_edge_set_map.end()) combined_boundary_edges.insert(set_to_boundary_edge_set_map.at(setID2).begin(), set_to_boundary_edge_set_map.at(setID2).end());
        if (set_to_edge_count_map.find(setID1) != set_to_edge_count_map.end()) for (const auto& pair : set_to_edge_count_map.at(setID1)) combined_edge_count_map[pair.first] = pair.second;
        if (set_to_edge_count_map.find(setID2) != set_to_edge_count_map.end()) for (const auto& pair : set_to_edge_count_map.at(setID2)) combined_edge_count_map[pair.first] = pair.second;
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
        set_to_edge_count_map.erase(setID1);
        set_to_boundary_edge_set_map.erase(setID1);
        set_to_boundary_point_set_map.erase(setID1);
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
        set_to_edge_count_map.erase(setID2);
        set_to_boundary_edge_set_map.erase(setID2);
        set_to_boundary_point_set_map.erase(setID2);
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
        set_to_edge_count_map.at(newSetID) = combined_edge_count_map;
        set_to_boundary_edge_set_map.at(newSetID) = combined_boundary_edges;
        set_to_boundary_point_set_map.at(newSetID) = combined_boundary_points;
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

    std::set<int> extract_points_by_setID(std::set<int> point_set, int setID)
    {
        // initialize
        std::set<int> out_point_set;

        // process
        for (int point_id : point_set)
        {
            if (point_to_set_map.at(point_id) == setID) out_point_set.insert(point_id);
        }

        // return
        return out_point_set;
    }


    // creates edges and triangles that connects the new point to the set
    void connect_point_to_set(int newPointID, int setID, std::set<int> searched_boundary_points_in_current_set)
    {
        // add point as boundary point
        add_boundary_point(newPointID, setID);

        // points_to_vector2d_map
        std::map<int, Eigen::Vector2d> points_to_vector2d_map = project_boundary_points_of_set_to_set_plane(setID);

        // existing edges between searched points (boundary)
        std::set<int> existing_boundary_edge_set = extract_existing_edge_between_points(searched_boundary_points_in_current_set, get_boundary_edge_set_of_set(setID));


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
            if (edge_edges_intersection(newEdge, get_boundary_edge_set_of_set(setID), points_to_vector2d_map)) continue;

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

    void add_point_by_radius_search(int newPointID, Eigen::Vector3d thisPointVEC, Eigen::Vector3d thisPointOriginVEC)
    {
        // if empty, can not set up radius search, add point to new set
        if (global_boundary_point_set.size() == 0)
        {
            int newSetID = getNewSetID();
            add_point_to_set(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
            add_boundary_point(newPointID, newSetID);
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
            add_point_to_set(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
            add_boundary_point(newPointID, newSetID);
            return;
        }

        // from searched points identify neighboring sets
        std::set<int> neighboring_sets; 
        for (int point_id : searched_boundary_points_set)
        {
            neighboring_sets.insert(point_to_set_map.at(point_id));
        }

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

        // try merge the sets_within_threshold (only merge if the new point can connect to both sets)
        std::set<int> sets_to_merge = sets_within_threshold;
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

                // once merged, restart
                again = true;
                break;
            }
            if (!again) break;
        }

        // from the merged sets, find the set that is closest to the point
        int closest_setID = -1;
        double closest_distance = std::numeric_limits<double>::max();
        for (int setID : sets_to_merge)
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
        if (closest_setID != -1 && closest_distance < distance_threshold)
        {
            add_point_to_set(newPointID, closest_setID, thisPointVEC, thisPointOriginVEC);            
            connect_point_to_set(newPointID, closest_setID, extract_points_by_setID(searched_boundary_points_set, closest_setID));
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
            add_point_to_set(newPointID, largest_setID, thisPointVEC, thisPointOriginVEC);
            connect_point_to_set(newPointID, largest_setID, extract_points_by_setID(searched_boundary_points_set, largest_setID));
            return;
        }

        // else, add the point to a new set
        int newSetID = getNewSetID();
        add_point_to_set(newPointID, newSetID, thisPointVEC, thisPointOriginVEC);
        add_boundary_point(newPointID, newSetID);
    }


    typename pcl::PointCloud<PointT>::Ptr 
    transform_cloud_to_global(typename pcl::PointCloud<PointT>::Ptr cloud, Eigen::Affine3d pose_eigen)
    {
        // transform point cloud
        typename pcl::PointCloud<PointT>::Ptr transformed_cloud (new pcl::PointCloud<PointT> ());
        pcl::transformPointCloud (*cloud, *transformed_cloud, pose_eigen);
        return transformed_cloud;
    }

    void step()
    {
        // if all points are processed
        if (ith_point >= ith_size) 
        {
            // next cloud
            ith_cloud++;
            ith_point = 0;

            // load data
            typename pcl::PointCloud<PointT>::Ptr pointcloud_local = data_loader.get_cloud(ith_cloud);
            Eigen::Affine3d pose = data_loader.get_pose(ith_cloud);
            pointcloud = transform_cloud_to_global(pointcloud_local, pose);
            origin = pose.translation();

            ith_size = pointcloud->size();
        }

        // get point
        Eigen::Vector3d thisPointVEC = pointcloud->points[ith_point].getVector3fMap().cast<double>();
        Eigen::Vector3d thisPointOriginVEC = origin;
        ith_point ++;
        
        // get new point id
        int newPointID = getNewPointID();


        // ------------- add point by triangle intersection

        // get list of intersected triangle by the point
        std::set<int> searched_triangles = bvhRoot.intersectionSearch(thisPointOriginVEC, thisPointVEC);

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
        for (const auto& pair : set_distance_map)
        {
            int setID = pair.first;
            double distance = pair.second;
            if (distance > distance_threshold) set_with_point_before_it.insert(setID);
            else if (distance < -distance_threshold) set_with_point_behind_it.insert(setID);
            else set_with_point_within_it.insert(setID);
        }
        
        bool point_added_to_set = false;

        // process point within set set
        // todo - if there are multiple sets with point within it, merge them
        for (int setID : set_with_point_within_it)
        {
            // add the point to the set
            add_point_to_set(newPointID, setID, thisPointVEC, thisPointOriginVEC);
            std::cout << ith_point << " / " << pointcloud->size() << " of pointcloud " << ith_cloud << " added to set " << setID << std::endl;
            point_added_to_set = true;

            // add the point to the triangle
            for (int triangleID : set_to_searched_triangle_map.at(setID))
            {
                triangle_to_points_map.at(triangleID).insert(newPointID);
            }
        }

        // process point behind set set

        // for each set that the point is behind

        // find the triangle penetrated by the point

        // remove the triangle

        // get the list of points inside the triangle

        // get the list of points that are cut off by the triangle

        // re add the list points by radius search, while avoid covering the new point



        // ------------- add point by radius search
        if (!point_added_to_set)
        {
            add_point_by_radius_search(newPointID, thisPointVEC, thisPointOriginVEC);
            std::cout << ith_point << " / " << pointcloud->size() << " of pointcloud " << ith_cloud << " added by radius search" << std::endl;
        }
    }

    void loop()
    {
        step();

        // finish all points in step
        while (ith_point < pointcloud->size())
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
    std::map<int, std::array<int, 2>> get_edge_to_vertices_map() {return edge_to_point_map;};
    std::map<int, std::array<int, 3>> get_triangle_to_vertices_map() {return triangle_to_vertices_map;};

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
    
    void compute_projected_point_to_vector3d_map()
    {
        // for each point
        for (int point_id : point_list)
        {
            // get set id
            int setID = point_to_set_map.at(point_id);

            // get projected point
            Eigen::Vector3d rayPlaneIntersectionPoint = point_set_intersection(point_id, setID);

            // store
            projected_points_to_vector3d_map[point_id] = rayPlaneIntersectionPoint;   
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
    std::map<int, std::set<int>> set_to_boundary_edge_set_map; // each set contains id to boundary edges
    std::map<int, std::set<int>> set_to_boundary_point_set_map; // each set contains id to boundary points
    std::map<int, std::map<int, int>> set_to_edge_count_map; // each map contains edge count

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
    flann3d flann;
    std::set<int> global_boundary_point_set;

        // triangle intersection
    TriangleBVH bvhRoot;

        // projected points
    std::map<int, Eigen::Vector3d> projected_points_to_vector3d_map;
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

        // triangle mesh
        pcl::PolygonMesh triangle_mesh;
        viewer_->addPolygonMesh(triangle_mesh, "triangle_mesh");

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

    bool show_triangle = false;
    bool show_projected_point = false;
    
    void update_display()
    {
        // get data from app
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr intensity_cloud;
        if (show_projected_point)
        {
            intensity_cloud = app_.projected_point_to_vector3d_set_colored_cloud();
        }
        else
        {
            intensity_cloud = app_.point_to_vector3d_set_colored_cloud();
        }
        std::map<int, std::array<int, 2>> edge_to_vertices_map = app_.get_edge_to_vertices_map();
        std::map<int, std::array<int, 3>> triangle_to_vertices_map = app_.get_triangle_to_vertices_map();
        

        // triangle faces
        pcl::PolygonMesh triangle_mesh;
        // add points
        pcl::toPCLPointCloud2(*intensity_cloud, triangle_mesh.cloud);
        // add triangles
        for (const auto& pair : triangle_to_vertices_map)
        {
            pcl::Vertices triangle;
            triangle.vertices.push_back(pair.second[0]);
            triangle.vertices.push_back(pair.second[1]);
            triangle.vertices.push_back(pair.second[2]);
            triangle_mesh.polygons.push_back(triangle);
        }

        // mesh
        pcl::PolygonMesh mesh;
        // add points
        pcl::toPCLPointCloud2(*intensity_cloud, mesh.cloud); 
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
        pcl::toPCLPointCloud2(*intensity_cloud, boundary_mesh.cloud);
        // get boundary edges
        std::vector<int> boundary_edge_list = app_.get_boundary_edge_list();
        // add edges
        for (int edge_id : boundary_edge_list)
        {
            pcl::Vertices edge;
            edge.vertices.push_back(edge_to_vertices_map.at(edge_id)[0]);
            edge.vertices.push_back(edge_to_vertices_map.at(edge_id)[1]);
            boundary_mesh.polygons.push_back(edge);
        }

        // triangle mesh
        if (show_triangle)
        {
            // add shape
            if (!viewer_->updatePolygonMesh(triangle_mesh, "triangle_mesh"))
            {
                viewer_->addPolygonMesh(triangle_mesh, "triangle_mesh");
            }
        }
        else
        {
            // remove shape
            if (viewer_->updatePolygonMesh(triangle_mesh, "triangle_mesh"))
            {
                viewer_->removeShape("triangle_mesh");
            }
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
        pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> intensity_handler(intensity_cloud);
        viewer_->updatePointCloud<pcl::PointXYZRGB>(intensity_cloud, intensity_handler, "cloud");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "cloud");
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
        if (event.getKeySym() == "t" && event.keyDown())
        {
            // toggle triangle visibility
            show_triangle = !show_triangle;
            update_display();
        }
        if (event.getKeySym() == "a" && event.keyDown())
        {
            // toggle projected point 
            show_projected_point = !show_projected_point;
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