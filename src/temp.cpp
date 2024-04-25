#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"

using InputPointT = VilensPointT;
int main()
{
    std::string pcd_file_folder = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/";
    std::string pose_file_path = "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam";

    DataLoader<InputPointT> data_loader(pcd_file_folder, pose_file_path);

    // // ------------------------------------------------------- viewer
    // // viewer
    // pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
    // viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    // viewer->initCameraParameters();
    // viewer->addCoordinateSystem(1);

    // // load cloud and pose
    // int i1 = 0;
    // typename pcl::PointCloud<InputPointT>::Ptr cloud1 = data_loader.get_cloud(i1);
    // Eigen::Affine3d pose1 = data_loader.get_pose(i1);
    // typename pcl::PointCloud<InputPointT>::Ptr cloud1_transformed = transform_to_global<InputPointT> (cloud1, pose1);
    // pcl::visualization::PointCloudColorHandlerCustom<InputPointT> color1(cloud1_transformed, 0, 255, 0);
    // viewer->addPointCloud<InputPointT> (cloud1_transformed, color1, "cloud1");
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud1");

    // int i2 = 60;
    // typename pcl::PointCloud<InputPointT>::Ptr cloud2 = data_loader.get_cloud(i2);
    // Eigen::Affine3d pose2 = data_loader.get_pose(i2);
    // typename pcl::PointCloud<InputPointT>::Ptr cloud2_transformed = transform_to_global<InputPointT> (cloud2, pose2);
    // pcl::visualization::PointCloudColorHandlerCustom<InputPointT> color2(cloud2_transformed, 255, 0, 0);
    // viewer->addPointCloud<InputPointT> (cloud2_transformed, color2, "cloud2");
    // viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud2");

    // // spin
    // viewer->spin();


    // // // ------------------------------------------------------- PLT
    // // plt triangulation
    // plt_plot_black_background();

    // // plt plot triangles
    // for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
    //     // indices
    //     int i1 = d.triangles[i];
    //     int i2 = d.triangles[i + 1];
    //     int i3 = d.triangles[i + 2];

    //     // plt - plot the triangle
    //     float tx0 = d.coords[2 * i1];
    //     float ty0 = d.coords[2 * i1 + 1];
    //     float tx1 = d.coords[2 * i2];
    //     float ty1 = d.coords[2 * i2 + 1];
    //     float tx2 = d.coords[2 * i3];
    //     float ty2 = d.coords[2 * i3 + 1];
    //     plt::plot({tx0, tx1, tx2, tx0}, {ty0, ty1, ty2, ty0}, std::map<std::string, std::string>{{"linewidth", "0.1"}, {"color", "white"}});
        
    //     // message
    //     std::cout << "adding triange " << i <<  " out of " << d.triangles.size() << std::endl;
    // }

    // // reverse x axis order
    // plt::xlim(180, -180);

    // // show
    // plt::show();

    
    // old cloud
    int i1 = 0;
    typename pcl::PointCloud<InputPointT>::Ptr old_cloud = data_loader.get_cloud(i1);
    Eigen::Affine3d old_pose = data_loader.get_pose(i1);

    // new cloud
    int i2 = 50;
    typename pcl::PointCloud<InputPointT>::Ptr new_cloud = data_loader.get_cloud(i2);
    Eigen::Affine3d new_pose = data_loader.get_pose(i2);

    // obtain old cloud triangulation in old cloud frame
    delaunator::Delaunator old_d = obtain_triangulation<InputPointT> (old_cloud);
    std::map<int, std::vector<int>> triangle_map = d_to_triangle_map(old_d);

    // transform old cloud to new cloud frame
    typename pcl::PointCloud<InputPointT>::Ptr old_cloud_local = transform_cloud_to_frame<InputPointT>(old_cloud, old_pose, new_pose);


    
    
    // // ------------- flann ---------------
    pcl::PointCloud<pcl::PointXYZ>::Ptr triangle_center_cloud = compute_triangle_center_cloud<InputPointT>(old_cloud_local, triangle_map);
    pcl::PointCloud<pcl::PointXY>::Ptr triangle_center_cloud_polar = compute_2d_polar_cloud<pcl::PointXYZ>(triangle_center_cloud);
    std::vector<float> flann_data_storage;
    for (const auto& point : triangle_center_cloud_polar->points)
    {
        flann_data_storage.push_back(point.x);
        flann_data_storage.push_back(point.y);
    }
    flann::Matrix<float> flann_data_matrix(flann_data_storage.data(), triangle_center_cloud_polar->size(), 2);
    flann::Index<flann::L2<float>> flann_tree(flann_data_matrix, flann::KDTreeIndexParams(1));
    flann_tree.buildIndex();
    int flann_last_id = triangle_center_cloud_polar->size() - 1;


    pcl::PointCloud<pcl::PointXY>::Ptr new_point_polar_cloud = compute_2d_polar_cloud<InputPointT>(new_cloud);
    for (std::size_t i = 0; i < new_point_polar_cloud->size(); i++)
    {
        // get new point (vector)
        InputPointT p_new_point = new_cloud->points[i];
        Eigen::Vector3f v_new_point = p_new_point.getVector3fMap();

        // get new point polar
        pcl::PointXY searchPoint = new_point_polar_cloud->points[i];

        // output stat
        std::cout << "processing point " << i << " out of " << new_point_polar_cloud->size() << std::endl;

        // search
        int K = 4;
        std::vector<float> query_point = {searchPoint.x, searchPoint.y};
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>(K));
        std::vector<std::vector<float>> list_of_search_dists(1, std::vector<float>(K));
        flann_tree.knnSearch(flann::Matrix<float>(query_point.data(), 1, 2), list_of_search_indices, list_of_search_dists, K, flann::SearchParams(-1, 0));
        std::vector<int> search_indices = list_of_search_indices[0];
        std::vector<float> search_dists = list_of_search_dists[0];
 
        // for each searched triangle center
        std::vector<int> intersected_triangle_indices;
        std::vector<float> intersected_triangle_distances;
        for (std::size_t j = 0; j < search_indices.size(); j++)
        {
            int index = search_indices[j];

            // get vertices index
            std::vector<int> vertices_index = triangle_map[index];

            // get vertices point
            InputPointT p0 = old_cloud_local->points[vertices_index[0]];
            InputPointT p1 = old_cloud_local->points[vertices_index[1]];
            InputPointT p2 = old_cloud_local->points[vertices_index[2]];

            // get vertices vector
            Eigen::Vector3f v0 = p0.getVector3fMap();
            Eigen::Vector3f v1 = p1.getVector3fMap();
            Eigen::Vector3f v2 = p2.getVector3fMap();

            // compute ray triangle intersection
            Eigen::Vector3f ray_origin = v_new_point;
            Eigen::Vector3f ray_direction = v_new_point.normalized();
            Eigen::Vector3f intersection = ray_triangle_intersection(ray_origin, ray_direction, v0, v1, v2);
            float distance = (intersection - v_new_point).norm();

            // check if intersection is inside the triangle
            bool inside = is_inside_triangle(v0, v1, v2, intersection);
            if (!inside) continue;
            
            // store the intersection
            intersected_triangle_indices.push_back(index);
            intersected_triangle_distances.push_back(distance);
        }

        // skip if no intersection
        if (intersected_triangle_indices.size() == 0) continue;

        // --- point --- 
        // add 
        old_cloud_local->push_back(p_new_point);
        int new_point_index = old_cloud_local->size() - 1;


        // find the triangle to remove (has the smallest distance to the new point)
        int min_index = std::distance(intersected_triangle_distances.begin(), std::min_element(intersected_triangle_distances.begin(), intersected_triangle_distances.end()));
        // std::cout << "min_index: " << min_index << std::endl;
        int min_triangle_index = intersected_triangle_indices[min_index];
        
        
        // --- triangle --- 
        // get vertices index
        std::vector<int> vertices_index = triangle_map[min_triangle_index];
        // compute
        std::vector<int> new_triangle1 = {vertices_index[0], vertices_index[1], new_point_index};
        std::vector<int> new_triangle2 = {vertices_index[1], vertices_index[2], new_point_index};
        std::vector<int> new_triangle3 = {vertices_index[2], vertices_index[0], new_point_index};
        // remove
        triangle_map.erase(min_triangle_index);
        // add
        triangle_map[flann_last_id+1] = new_triangle1;
        triangle_map[flann_last_id+2] = new_triangle2;
        triangle_map[flann_last_id+3] = new_triangle3;
        flann_last_id += 3;

        // --- flann ---
        // remove
        flann_tree.removePoint(min_triangle_index);
        // compute azimuth and altitude for new triangle center
        // get vertices point
        InputPointT p0 = old_cloud_local->points[vertices_index[0]];
        InputPointT p1 = old_cloud_local->points[vertices_index[1]];
        InputPointT p2 = old_cloud_local->points[vertices_index[2]];
        // get vertices vector
        Eigen::Vector3f v0 = p0.getVector3fMap();
        Eigen::Vector3f v1 = p1.getVector3fMap();
        Eigen::Vector3f v2 = p2.getVector3fMap();
        // compute azimuth and altitude
        float azimuth1, altitude1;
        float azimuth2, altitude2;
        float azimuth3, altitude3;
        Eigen::Vector3f new_triangle_center1 = (v0 + v1 + v_new_point) / 3;
        Eigen::Vector3f new_triangle_center2 = (v1 + v2 + v_new_point) / 3;
        Eigen::Vector3f new_triangle_center3 = (v2 + v0 + v_new_point) / 3;
        compute_azimuth_and_altitude(new_triangle_center1, azimuth1, altitude1);
        compute_azimuth_and_altitude(new_triangle_center2, azimuth2, altitude2);
        compute_azimuth_and_altitude(new_triangle_center3, azimuth3, altitude3);
        // generate polar dataset
        std::vector<float> new_triangle_center1_polar = {azimuth1, altitude1};
        std::vector<float> new_triangle_center2_polar = {azimuth2, altitude2};
        std::vector<float> new_triangle_center3_polar = {azimuth3, altitude3};
        // add
        flann_tree.addPoints(flann::Matrix<float>(new_triangle_center1_polar.data(), 1, 2), 10);
        flann_tree.addPoints(flann::Matrix<float>(new_triangle_center2_polar.data(), 1, 2), 10);
        flann_tree.addPoints(flann::Matrix<float>(new_triangle_center3_polar.data(), 1, 2), 10);
    }



    // todo
    // add each new point to triangles
    // find triangles in old cloud that contains the new point

    // compute polar coordinates for old_cloud_local
    

    
    
    
    

    
    // add point to mesh        
    // treat each new point as individual points

    // when new points are inside old triangles
    // need to find new point to old triangle correspondence
    // since new point have direction aligned with the view direction, if a new point has azimuth and altitude inside an old triangle, 
    // the intersection of the new point to the old triangle will be inside the old triangle
    // don't compute the delaunay triangles for new cloud, as we treat each new point as individual points
    // 
    // project old triangle centers to 


    

    // viewer
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1);

    // create mesh
    pcl::PolygonMesh old_mesh = convert_to_mesh<InputPointT>(old_cloud_local, triangle_map);
    // // add mesh
    // viewer->addPolygonMesh(old_mesh, "mesh");
    // add polyline
    viewer->addPolylineFromPolygonMesh(old_mesh, "polyline");

    // spin
    viewer->spin();

    return 0;
}



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