#include "./main.cpp"

// test for eye_patch intersection 

int main()
{
    // input
    Eigen::Vector3f v1(0, 0, 0);
    Eigen::Vector3f v2(1, 0, 0);
    Eigen::Vector3f v3(0, 2, 0);
    Eigen::Vector3f current_point(0.5, 0.5, 1);
    Eigen::Vector3f current_point_direction(0, 0, -1);
    double vertex_variance = std::pow(0.01, 2);
    float sphere_radius = 3;

    // output
    Eigen::Vector3f likelihood_point;
    double likelihood_variance;

    // compute likelihood point
    bool sphere_exists = eye_patch_intersection(v1, v2, v3, current_point, current_point_direction, vertex_variance, sphere_radius, likelihood_point, likelihood_variance);

    // print
    if (!sphere_exists)
    {
        std::cout << "no sphere exists" << std::endl;
        return 0;
    }
    else
    {
    std::cout << "likelihood point: " << likelihood_point.transpose() << std::endl;
    std::cout << "likelihood variance: " << likelihood_variance << std::endl;
    }

    return 0;
}

// test for is_inside_triangle
int main()
{
    Eigen::Vector3f p0(0, 0, 0);
    Eigen::Vector3f p1(1, 0, 0);
    Eigen::Vector3f p2(0.5, 1, 0);

    // plot triangles given three points
    std::vector<float> x_list = {p0(0), p1(0), p2(0), p0(0)};
    std::vector<float> y_list = {p0(1), p1(1), p2(1), p0(1)};
    plt::plot(x_list, y_list, std::map<std::string, std::string>{{"linewidth", "0.1"}, {"color", "black"}});

    // scatter plot random points and color by if they are inside the triangle
    std::vector<float> x_list_random;
    std::vector<float> y_list_random;
    std::vector<int> inside_list;
    for (int i = 0; i < 1000; i++)
    {
        float x = (rand() % 100) / 100.0;
        float y = (rand() % 100) / 100.0;
        Eigen::Vector3f point(x, y, 0);
        bool inside = is_inside_triangle(p0, p1, p2, point);
        x_list_random.push_back(x);
        y_list_random.push_back(y);
        inside_list.push_back(inside);
    }
    // red if not inside, green if inside
    plt::scatter_colored(x_list_random, y_list_random, inside_list, 1, std::map<std::string, std::string>{{"cmap", "viridis"}});
    
    plt::show();

}


int main()
{
    // given index number, add pointcloud to display
    std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";
    std::string pcd_file2 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd";
    pcl::PointCloud<InputPointT>::Ptr cloud2 = load_pointcloud<InputPointT>(pcd_file2);

    // triangulate cloud
    delaunator::Delaunator d = obtain_triangulation<InputPointT>(cloud2);

    // cloud point to triangle map
    // point: cloud index of the point
    // triangle: list index of the first triangle point in d.triangles, which are stored as cloud index
    std::map<int, std::vector<int>> pt_map = point_to_triangle_map(d);

    // update old points using triangles
    // todo


    // ------------------------------- search
    // 2d cloud from coords
    pcl::PointCloud<pcl::PointXY>::Ptr cloud2d (new pcl::PointCloud<pcl::PointXY>);
    cloud2d->resize(d.coords.size() / 2);
    for (std::size_t i = 0; i < d.coords.size(); i+=2)
    {
        pcl::PointXY point;
        point.x = d.coords[i];
        point.y = d.coords[i + 1];
        cloud2d->points[i / 2] = point;
    }

    // kd tree search
    pcl::KdTreeFLANN<pcl::PointXY> kdtree;
    kdtree.setInputCloud (cloud2d);
    pcl::PointXY searchPoint;
    searchPoint.x = 0;
    searchPoint.y = 0;

    // K nearest neighbor search
    int K = 1;
    std::vector<int> pointIdxKNNSearch(K);
    std::vector<float> pointKNNSquaredDistance(K);
    if ( kdtree.nearestKSearch (searchPoint, K, pointIdxKNNSearch, pointKNNSquaredDistance) > 0 )
    {
        for (std::size_t i = 0; i < pointIdxKNNSearch.size (); ++i)
        std::cout << "    "  <<   cloud2d->points[pointIdxKNNSearch[i]].x 
                    << " " << cloud2d->points[pointIdxKNNSearch[i]].y 
                    << " (squared distance: " << pointKNNSquaredDistance[i] << ")" << std::endl;
    }

    // triangles that contains the search point
    std::vector<int> triangles_containing_search_point = pt_map[pointIdxKNNSearch[0]];

    // plt - plot searched triangles
    for (std::size_t i = 0; i < triangles_containing_search_point.size(); i++)
    {
        int i1 = d.triangles[triangles_containing_search_point[i]];
        int i2 = d.triangles[triangles_containing_search_point[i] + 1];
        int i3 = d.triangles[triangles_containing_search_point[i] + 2];

        // plt - plot the triangle
        float tx0 = d.coords[2 * i1];
        float ty0 = d.coords[2 * i1 + 1];
        float tx1 = d.coords[2 * i2];
        float ty1 = d.coords[2 * i2 + 1];
        float tx2 = d.coords[2 * i3];
        float ty2 = d.coords[2 * i3 + 1];
        plt::plot({tx0, tx1, tx2, tx0}, {ty0, ty1, ty2, ty0}, std::map<std::string, std::string>{{"linewidth", "1"}, {"color", "black"}});
    }

    // plt - plot all points
    std::vector<double> x_list, y_list;
    for (std::size_t i = 0; i < cloud2d->size(); i++)
    {
        x_list.push_back((*cloud2d)[i].x);
        y_list.push_back((*cloud2d)[i].y);
    }
    plt::scatter(x_list, y_list, 10);

    // plt - plot search point
    std::vector<double> search_x_list = {searchPoint.x};
    std::vector<double> search_y_list = {searchPoint.y};
    plt::scatter(search_x_list, search_y_list, 10, std::map<std::string, std::string>{{"color", "red"}});

    // plt - plot K nearest neighbors
    std::vector<double> knn_x_list, knn_y_list;
    for (std::size_t i = 0; i < pointIdxKNNSearch.size(); i++)
    {
        knn_x_list.push_back((*cloud2d)[pointIdxKNNSearch[i]].x);
        knn_y_list.push_back((*cloud2d)[pointIdxKNNSearch[i]].y);
    }
    plt::scatter(knn_x_list, knn_y_list, 10, std::map<std::string, std::string>{{"color", "green"}});

    // plt show
    plt::show();

    return (0);
}


int main() {
    /* x0, y0, x1, y1, ... */
    std::vector<double> coords = {-1, 1, 1, 1, 1, -1, -1, -1};

    //triangulation happens here
    delaunator::Delaunator d(coords);

    // print content of hull_tri
    std::cout << "hull_tri: ";
    for (std::size_t i = 0; i < d.hull_tri.size(); i++)
    {
        std::cout << d.hull_tri[i] << " ";
    }

}




int main()
{
    pcl::PointCloud<pcl::PointXY>::Ptr cloud (new pcl::PointCloud<pcl::PointXY>);

    // Generate pointcloud data
    cloud->width = 1000;
    cloud->height = 1;
    cloud->points.resize (cloud->width * cloud->height);

    for (std::size_t i = 0; i < cloud->size (); ++i)
    {
        (*cloud)[i].x = 1024.0f * rand () / (RAND_MAX + 1.0f);
        (*cloud)[i].y = 1024.0f * rand () / (RAND_MAX + 1.0f);
    }

    


    pcl::KdTreeFLANN<pcl::PointXY> kdtree;
    kdtree.setInputCloud (cloud);
    pcl::PointXY searchPoint;
    searchPoint.x = 1024.0f * rand () / (RAND_MAX + 1.0f);
    searchPoint.y = 1024.0f * rand () / (RAND_MAX + 1.0f);

    // K nearest neighbor search

    int K = 10;

    std::vector<int> pointIdxKNNSearch(K);
    std::vector<float> pointKNNSquaredDistance(K);

    std::cout << "K nearest neighbor search at (" << searchPoint.x 
                << " " << searchPoint.y 
                << ") with K=" << K << std::endl;

    if ( kdtree.nearestKSearch (searchPoint, K, pointIdxKNNSearch, pointKNNSquaredDistance) > 0 )
    {
        for (std::size_t i = 0; i < pointIdxKNNSearch.size (); ++i)
        std::cout << "    "  <<   (*cloud)[ pointIdxKNNSearch[i] ].x 
                    << " " << (*cloud)[ pointIdxKNNSearch[i] ].y 
                    << " (squared distance: " << pointKNNSquaredDistance[i] << ")" << std::endl;
    }

    // plt - plot points
    std::vector<double> x_list, y_list;
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        x_list.push_back((*cloud)[i].x);
        y_list.push_back((*cloud)[i].y);
    }
    plt::scatter(x_list, y_list, 10);

    // plt - plot search point
    std::vector<double> search_x_list = {searchPoint.x};
    std::vector<double> search_y_list = {searchPoint.y};
    plt::scatter(search_x_list, search_y_list, 10, std::map<std::string, std::string>{{"color", "red"}});

    // plt - plot K nearest neighbors
    std::vector<double> knn_x_list, knn_y_list;
    for (std::size_t i = 0; i < pointIdxKNNSearch.size(); i++)
    {
        knn_x_list.push_back((*cloud)[pointIdxKNNSearch[i]].x);
        knn_y_list.push_back((*cloud)[pointIdxKNNSearch[i]].y);
    }
    plt::scatter(knn_x_list, knn_y_list, 10, std::map<std::string, std::string>{{"color", "green"}});

    plt::show();

}


int main()
{
    // given index number, add pointcloud to display
    std::string pose_file = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam";

    // std::string pcd_file1 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460854_847538000.pcd";
    // pcl::PointCloud<InputPointT>::Ptr cloud1 = load_pointcloud<InputPointT>(pcd_file1);
    // Eigen::Affine3d pose_eigen1 = find_pose(pcd_file1, pose_file);
    // pcl::PointCloud<InputPointT>::Ptr transformed_cloud1 = transform_to_global<InputPointT>(cloud1, pose_eigen1);

    std::string pcd_file2 = "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/cloud_1711460869_333305000.pcd";
    pcl::PointCloud<InputPointT>::Ptr cloud2 = load_pointcloud<InputPointT>(pcd_file2);
    Eigen::Affine3d pose_eigen2 = find_pose(pcd_file2, pose_file);
    pcl::PointCloud<InputPointT>::Ptr transformed_cloud2 = transform_to_global<InputPointT>(cloud2, pose_eigen2);


    // update old points using triangles

    // // transform cloud1 to cloud2 frame 
    // Eigen::Affine3d pose1_to_pose2 = pose_eigen2.inverse() * pose_eigen1;
    // pcl::PointCloud<InputPointT>::Ptr cloud1_in_cloud2_frame (new pcl::PointCloud<InputPointT> ());
    // pcl::transformPointCloud (*cloud1, *cloud1_in_cloud2_frame, pose1_to_pose2);

    // triangulate cloud2
    delaunator::Delaunator d = obtain_triangulation<InputPointT>(cloud2);


    // // ------------------------------------------------------ plt
    // // plot black background
    // plt_plot_black_background();
    // // plot triangles
    // plt_plot_triangles(d);
    // // display
    // plt::show();



    // // ------------------------------------------------------ pclvisuliazer    
    // // set up viewer
    // pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    // viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    
    // // set up viewports
    // int port1(0);
    // viewer->createViewPort (0.0, 0.0, 0.5, 1.0, port1);
    // viewer->setBackgroundColor (0, 0, 0, port1);
    // int port2(0);
    // viewer->createViewPort (0.5, 0.0, 1.0, 1.0, port2);
    // viewer->setBackgroundColor (0, 0, 0, port2);

    // // set up coordinate system
    // viewer->initCameraParameters();
    // viewer->addCoordinateSystem(1);

    // // // add triangles to viewer
    // // viewer_add_triangles<InputPointT>(viewer, port1, cloud2, transformed_cloud2, d2);
    // // display pointclouds
    // add_to_viewer<InputPointT>(viewer, port1, transformed_cloud1, "transformed cloud", color_tuple(0, 255, 0), 1);
    // add_to_viewer<InputPointT>(viewer, port1, transformed_cloud2, "transformed cloud2", color_tuple(255, 0, 0), 1);
    // add_to_viewer<InputPointT>(viewer, port2, cloud1_in_cloud2_frame, "cloud1 in cloud2 frame", color_tuple(0, 255, 0), 1);
    // add_to_viewer<InputPointT>(viewer, port2, cloud2, "cloud2", color_tuple(255, 0, 0), 1);
    
    // // display
    // viewer->spin();


    // ------------------------------------------------------ pclvisuliazer    
    // set up viewer
    pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    
    // set up coordinate system
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1);

    // add triangles to viewer
    viewer_add_triangles<InputPointT>(viewer, 0, cloud2, transformed_cloud2, d);
    // display pointclouds
    add_to_viewer<InputPointT>(viewer, 0, transformed_cloud2, "cloud", color_tuple(0, 255, 0), 1);

    // display
    viewer->spin();

    return (0);
}