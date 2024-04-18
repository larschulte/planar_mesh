#include "./main.cpp"

// test for eye_patch intersection 

// viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "control cloud");
// viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_OPACITY, 0.7, "oldcloud mean");


// test curvature computation
// compute vertex to link map
std::map<int, std::vector<int>> compute_vertex_to_link_map(delaunator::Delaunator d)
{
    // initialize
    std::map<int, std::vector<int>> vertex_to_link_map;

    // compute link
    for (std::size_t i = 0; i < d.triangles.size(); i+=3)
    {
        // vertcies index
        int v1_index = d.triangles[i];
        int v2_index = d.triangles[i + 1];
        int v3_index = d.triangles[i + 2];
        
        // push link
        vertex_to_link_map[v1_index].push_back(v2_index);
        vertex_to_link_map[v1_index].push_back(v3_index);
        vertex_to_link_map[v2_index].push_back(v1_index);
        vertex_to_link_map[v2_index].push_back(v3_index);
        vertex_to_link_map[v3_index].push_back(v1_index);
        vertex_to_link_map[v3_index].push_back(v2_index);
    }
    
    for (auto& [vertex_index, links] : vertex_to_link_map)
    {
        unique_vector(links);
    }

    // return
    return vertex_to_link_map;
}
// compute vertex curvature 
template <typename PointT>
std::map<int, float> compute_vertex_to_curvature_map(typename pcl::PointCloud<PointT>::Ptr cloud, std::map<int, std::vector<int>> vertex_to_link_map, std::map<int, Eigen::Vector3f> vertex_to_normal_map)
{
    // initialize
    std::map<int, float> vertex_to_curvature_map;

    // compute curvature
    for (auto& [vertex_index, links] : vertex_to_link_map)
    {
        Eigen::Vector3f p_c = cloud->at(vertex_index).getVector3fMap();
        Eigen::Vector3f n_c = vertex_to_normal_map[vertex_index];

        // compute mean curvature
        float total_curvature = 0;
        int count = 0;
        for (auto& link : links)
        {
            // individual link curvature
            Eigen::Vector3f p_l = cloud->at(link).getVector3fMap();
            Eigen::Vector3f n_l = vertex_to_normal_map[link];
            float link_curvature = compute_edge_curvature(p_c, p_l, n_c, n_l);

            // mean
            total_curvature += link_curvature;
            count++;
        }
        float mean_curvature = total_curvature / count;

        // store
        vertex_to_curvature_map[vertex_index] = mean_curvature;
    }

    // return
    return vertex_to_curvature_map;
}   
using InputPointT = VilensPointT;
int main()
{
    // parameters
    std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
    std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";

    // data loader
    DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);
    pcl::PointCloud<InputPointT>::Ptr cloud (new pcl::PointCloud<InputPointT>);
    cloud = data_loader.get_cloud(0);

    // triangulate the cloud
    delaunator::Delaunator d = obtain_triangulation<InputPointT>(cloud);

    // plt_plot_black_background();
    // plt_plot_triangles(d);
    // // scatter color
    // std::vector<double> scatter_color;
    // for (std::size_t i = 0; i < cloud->size(); i++)
    // {
    //     float distance = cloud->points[i].getVector3fMap().norm();
    //     float distance_threshold = 6;
    //     if (distance > distance_threshold)
    //     {
    //         distance = distance_threshold;
    //     }
        
    //     scatter_color.push_back(distance);
    // }
    // plt_scatter_plot_coords(d, scatter_color); 
    // plt::show();

    

    // compute vertex normal
    std::map<int, Eigen::Vector3f> vertex_to_normal_map = compute_vertex_to_normal_map<InputPointT>(cloud, d);
    
    // compute vertex link
    std::map<int, std::vector<int>> vertex_to_link_map = compute_vertex_to_link_map(d);

    // compute vertex curvature
    std::map<int, float> vertex_to_mean_curvature_map = compute_vertex_to_curvature_map<InputPointT>(cloud, vertex_to_link_map, vertex_to_normal_map);


    // the current update assume planar surface within each triangle, and does not filter the planar surface even if the triangle is very large
    // this will be solved when introducing eye patch

    // ------------------------------------------------------ pclvisuliazer    
    // set up viewer
    pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    
    // set up viewports
    int port1(0);
    viewer->createViewPort (0.0, 0.0, 1, 1.0, port1);
    viewer->setBackgroundColor (0, 0, 0, port1);
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1);


    // convert to viewer cloud
    pcl::PointCloud<pcl::PointXYZINormal>::Ptr viewer_cloud(new pcl::PointCloud<pcl::PointXYZINormal>);
    viewer_cloud->resize(cloud->size());
    for (std::size_t i = 0; i < cloud->size(); i++)
    {
        pcl::PointXYZINormal viewer_point;
        viewer_point.x = cloud->points[i].x;
        viewer_point.y = cloud->points[i].y;
        viewer_point.z = cloud->points[i].z;
        float radius = abs(1.0 / vertex_to_mean_curvature_map[i]);
        std::cout << "radius is " << radius << std::endl;
        // cap the radius
        float radius_cap = 0.8;
        if (radius > radius_cap)
        {
            radius = radius_cap;
        }
        viewer_point.intensity = radius;
        viewer_point.normal_x = vertex_to_normal_map[i][0];
        viewer_point.normal_y = vertex_to_normal_map[i][1];
        viewer_point.normal_z = vertex_to_normal_map[i][2];
        viewer_cloud->points[i] = viewer_point;
    }

    // color
    // std::tuple<int, int, int> color(0, 255, 0);
    // pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZINormal> color_handler(xyz_cloud, std::get<0>(color), std::get<1>(color), std::get<2>(color));
    

    // add to viewer
    pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZINormal> color_handler(viewer_cloud, "intensity");
    viewer->addPointCloud<pcl::PointXYZINormal> (viewer_cloud, color_handler, "pointcloud", port1);
    viewer->addPointCloudNormals<pcl::PointXYZINormal> (viewer_cloud, 1, 0.05, "normals", port1);
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "pointcloud");
    
    

    // add spheres
    for (std::size_t i = 0; i < cloud->size(); i = i + 10)
    {
        pcl::PointXYZINormal point = viewer_cloud->points[i];
        float radius = point.intensity;
        if (radius > 0.1)
        {
            continue;
        }
        Eigen::Vector3f point_vector = point.getVector3fMap();
        Eigen::Vector3f point_normal = point.getNormalVector3fMap();
        Eigen::Vector3f sphere_center = point_vector - radius * point_normal;

        viewer->addSphere(pcl::PointXYZ(sphere_center[0], sphere_center[1], sphere_center[2]), radius, 1, 0, 0, "sphere" + std::to_string(i), port1);
    }


    // display
    viewer->spin();

    return (0);
}



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