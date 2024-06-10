void add_set(int newSetID)
{
    set_list.push_back(newSetID);
    set_to_color_map[newSetID] = {rand() % 256, rand() % 256, rand() % 256};
    set_to_points_map[newSetID] = {};
    set_to_edges_map[newSetID] = {};
    set_to_triangles_map[newSetID] = {};
    set_to_mean_map[newSetID] = Eigen::Vector3d::Zero();
    set_to_covariance_matrix_map[newSetID] = Eigen::Matrix3d::Zero();
    set_to_eigenvectors_map[newSetID] = Eigen::Matrix3d::Identity();
    set_to_eigenvalues_map[newSetID] = Eigen::Vector3d::Zero();
    set_to_normal_map[newSetID] = Eigen::Vector3d(0, 0, 1);
    boundary_edge_of_set[newSetID] = {};
    boundary_point_of_set[newSetID] = {};
}

void add_point(int newPointID, int setID, const Eigen::Vector3d& thisPoint, const Eigen::Vector3d& origin)
{
    point_list.push_back(newPointID);
    point_to_vector3d_map[newPointID] = thisPoint;
    point_to_origin_vector3d_map[newPointID] = origin;
    point_to_edges_map[newPointID] = {};
    point_to_set_map[newPointID] = setID;
    set_to_points_map.at(setID).insert(newPointID);
    add_to_plane_estimate(newPointID, setID);

    update_boundary_point_record(newPointID, setID);
}

void store_point_in_triangle(int newPointID, int triangleID, const Eigen::Vector3d& thisPoint, const Eigen::Vector3d& origin)
{
    int setID = triangle_to_set_map.at(triangleID);

    point_list.push_back(newPointID);
    point_to_vector3d_map[newPointID] = thisPoint;
    point_to_origin_vector3d_map[newPointID] = origin;
    point_to_edges_map[newPointID] = {};
    point_to_storing_triangle_map[newPointID] = triangleID;
    point_to_set_map[newPointID] = setID;
    set_to_points_map.at(setID).insert(newPointID);
    add_to_plane_estimate(newPointID, setID);
    // update_boundary_point_record(newPointID, setID);
    triangle_to_points_map.at(triangleID).insert(newPointID);
}

void add_edge(int newEdgeID, int setID, std::array<int, 2> vertices)
{
    vertices = sortTwoInts(vertices[0], vertices[1]);
    int pointID1 = vertices[0];
    int pointID2 = vertices[1];

    edge_list.push_back(newEdgeID);
    edge_to_point_map[newEdgeID] = vertices;
    edge_to_point_map_reverse[vertices] = newEdgeID;
    point_to_edges_map.at(pointID1).insert(newEdgeID);
    point_to_edges_map.at(pointID2).insert(newEdgeID);
    set_to_edges_map.at(setID).insert(newEdgeID);
    edge_to_set_map[newEdgeID] = setID;

    edge_to_edge_count_map[newEdgeID] = 0;
    update_boundary_edge_record(newEdgeID, setID);
    update_boundary_point_record(pointID1, setID);
    update_boundary_point_record(pointID2, setID);
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
    triangle_to_vertices_map[newTriangleID] = vertices;
    triangle_to_vertices_map_reverse[vertices] = newTriangleID;
    triangle_to_points_map[newTriangleID] = {};
    triangle_to_set_map[newTriangleID] = setID;
    set_to_triangles_map.at(setID).insert(newTriangleID);

    edge_to_edge_count_map.at(edgeID1) ++;
    edge_to_edge_count_map.at(edgeID2) ++;
    edge_to_edge_count_map.at(edgeID3) ++;
    update_boundary_edge_record(edgeID1, setID);
    update_boundary_edge_record(edgeID2, setID);
    update_boundary_edge_record(edgeID3, setID);
    update_boundary_point_record(pointID1, setID);
    update_boundary_point_record(pointID2, setID);
    update_boundary_point_record(pointID3, setID);

    bool inserted = global_triangle_set.insert(newTriangleID).second;
    if (inserted) bvhRoot.addTriangle(newTriangleID, vertices, point_to_vector3d_map.at(pointID1), point_to_vector3d_map.at(pointID2), point_to_vector3d_map.at(pointID3));
}

void remove_set(int setID)
{
    // remove
    set_list.erase(std::remove(set_list.begin(), set_list.end(), setID), set_list.end());
    set_to_color_map.erase(setID);
    set_to_points_map.erase(setID);
    set_to_edges_map.erase(setID);
    set_to_triangles_map.erase(setID);
    set_to_mean_map.erase(setID);
    set_to_covariance_matrix_map.erase(setID);
    set_to_eigenvectors_map.erase(setID);
    set_to_eigenvalues_map.erase(setID);
    set_to_normal_map.erase(setID);
    boundary_edge_of_set.erase(setID);
    boundary_point_of_set.erase(setID);
}

void remove_point(int pointID)
{
    // original info
    Eigen::Vector3d thisPoint = point_to_vector3d_map.at(pointID);
    Eigen::Vector3d origin = point_to_origin_vector3d_map.at(pointID);
    int setID = point_to_set_map.at(pointID);

    // remove
    point_list.erase(std::remove(point_list.begin(), point_list.end(), pointID), point_list.end());
    point_to_set_map.erase(pointID);
    update_boundary_point_record(pointID, setID); // update boundary before point to vector3d map is erased
    point_to_vector3d_map.erase(pointID);
    point_to_origin_vector3d_map.erase(pointID);
    set_to_points_map.at(setID).erase(pointID);
    point_to_storing_triangle_map.erase(pointID);
    remove_from_plane_estimate(pointID, setID);        
    
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
    update_boundary_edge_record(edgeID, setID);
    update_boundary_point_record(pointID1, setID);
    update_boundary_point_record(pointID2, setID);
    
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
    if (erased) bvhRoot.delete_face(triangleID, vertices, point_to_vector3d_map.at(pointID1), point_to_vector3d_map.at(pointID2), point_to_vector3d_map.at(pointID3));

    // update boundary
    edge_to_edge_count_map.at(edgeID1) --;
    edge_to_edge_count_map.at(edgeID2) --;
    edge_to_edge_count_map.at(edgeID3) --;
    update_boundary_edge_record(edgeID1, setID);
    update_boundary_edge_record(edgeID2, setID);
    update_boundary_edge_record(edgeID3, setID);
    update_boundary_point_record(pointID1, setID);
    update_boundary_point_record(pointID2, setID);
    update_boundary_point_record(pointID3, setID);

    // if the removal of triangle causes an edge to be isolated, remove the edge
    if (edge_to_edge_count_map.at(edgeID1) == 0) remove_edge(edgeID1);
    if (edge_to_edge_count_map.at(edgeID2) == 0) remove_edge(edgeID2);
    if (edge_to_edge_count_map.at(edgeID3) == 0) remove_edge(edgeID3);

    // // if the removal of triangle causes a point to be isolated, remove the point
    for (int pointID : points_within_triangles) remove_point(pointID);
}

void update_boundary_point_record(int pointID, int setID)
{
    // if point is removed
    if (point_to_set_map.find(pointID) == point_to_set_map.end())
    {
        bool erased = boundary_point_set.erase(pointID);
        boundary_point_of_set.at(setID).erase(pointID);
        if (erased) rrstree.deleteBoundaryPoint(pointID, point_to_vector3d_map.at(pointID));
        return;
    }

    // if point is stored in triangle
    if (point_to_storing_triangle_map.find(pointID) != point_to_storing_triangle_map.end())
    {
        bool erased = boundary_point_set.erase(pointID);
        boundary_point_of_set.at(setID).erase(pointID);
        if (erased) rrstree.deleteBoundaryPoint(pointID, point_to_vector3d_map.at(pointID));
        return;
    }

    // if point to edge map is updated
    bool isolated = point_to_edges_map.at(pointID).empty();
    bool connected_to_boundary_edge = false;
    for (int edge_id : point_to_edges_map.at(pointID))
    {
        if (boundary_edge_set.find(edge_id) != boundary_edge_set.end()) {connected_to_boundary_edge = true; break;}
    }
    if (isolated || connected_to_boundary_edge)
    {
        bool inserted = boundary_point_set.insert(pointID).second;
        boundary_point_of_set.at(setID).insert(pointID);
        if (inserted) 
        {
            // origin to point distance 
            double distance = (point_to_vector3d_map.at(pointID) - point_to_origin_vector3d_map.at(pointID)).norm();
            double radius = distance * distance_to_radius_ratio;
            rrstree.addBoundaryPoint(pointID, point_to_vector3d_map.at(pointID), radius);
        }
    }
    else
    {
        bool erased = boundary_point_set.erase(pointID);
        boundary_point_of_set.at(setID).erase(pointID);
        if (erased) rrstree.deleteBoundaryPoint(pointID, point_to_vector3d_map.at(pointID));
    }
}

void update_boundary_edge_record(int edgeID, int setID)
{
    // if edge is removed
    if (edge_to_set_map.find(edgeID) == edge_to_set_map.end())
    {
        boundary_edge_set.erase(edgeID);
        boundary_edge_of_set.at(setID).erase(edgeID);
        return;
    }

    // if edge count is updated
    int count = edge_to_edge_count_map.at(edgeID);
    bool is_boundary_edge = count == 0 || count == 1;
    if (is_boundary_edge)
    {
        boundary_edge_set.insert(edgeID);
        boundary_edge_of_set.at(setID).insert(edgeID);
    }
    else
    {
        boundary_edge_set.erase(edgeID);
        boundary_edge_of_set.at(setID).erase(edgeID);
    }
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


    Eigen::Vector2d project_point_to_set_plane(int pointID, int setID)
    {
        const Eigen::Vector3d& rayOrigin = point_to_origin_vector3d_map.at(pointID);
        const Eigen::Vector3d& rayEndPoint = point_to_vector3d_map.at(pointID);
        const Eigen::Vector3d& mean = set_to_mean_map.at(setID);
        const Eigen::Vector3d& normal = set_to_normal_map.at(setID);
        Eigen::Matrix3d eigenvectors = set_to_eigenvectors_map.at(setID);
        Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();
        Eigen::Vector3d rayPlaneIntersectionPoint = ray_plane_intersection(rayOrigin, rayEndPoint, mean, normal);
        Eigen::Vector2d projected_point = (projection_matrix.transpose() * rayPlaneIntersectionPoint).head<2>();
        return projected_point;
    }

        Eigen::Vector3d merge_means_of_sets(int setID1, int setID2) 
    {
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID1);
        Eigen::Vector3d mean2 = set_to_mean_map.at(setID2);
        int size1 = set_to_points_map.at(setID1).size();
        int size2 = set_to_points_map.at(setID2).size();
        return merge_means(mean1, mean2, size1, size2);
    }

    Eigen::Vector3d merge_means_between_set_and_point(int setID, const Eigen::Vector3d& point)
    {
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID);
        Eigen::Vector3d mean2 = point;
        int size1 = set_to_points_map.at(setID).size();
        int size2 = 1;
        return merge_means(mean1, mean2, size1, size2);
    }

    Eigen::Matrix3d merge_covariances_between_set_and_point(int setID, const Eigen::Vector3d& point)
    {
        Eigen::Matrix3d cov1 = set_to_covariance_matrix_map.at(setID);
        Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero(); 
        Eigen::Vector3d mean1 = set_to_mean_map.at(setID);
        Eigen::Vector3d mean2 = point;
        int size1 = set_to_points_map.at(setID).size();
        int size2 = 1;
        return merge_covariances(cov1, cov2, mean1, mean2, size1, size2);
    }

        std::vector<BoundaryPoint> rrstree_get_boundary_points()
    {
        return rrstree.getBoundaryPoints();
    }
