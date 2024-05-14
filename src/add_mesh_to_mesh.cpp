#include "eye_patch/DataLoader.hpp"
#include "eye_patch/Algorithm.hpp"
#include "eye_patch/Visualization.hpp"

#include "point_type/BagPointT.hpp"
#include "point_type/VilensPointT.hpp"


class flann2d
{
public:
    flann2d()
        :
        flann_tree(flann::KDTreeIndexParams(1))
        {};
    
    void set_input(pcl::PointCloud<pcl::PointXY>::Ptr triangle_center_cloud_polar)
    {
        // compute data storage
        for (const auto& point : triangle_center_cloud_polar->points)
        {
            flann_data_storage.push_back(point.x);
            flann_data_storage.push_back(point.y);
        }

        // record
        flann_last_id = triangle_center_cloud_polar->size() - 1;

        // add to flann
        flann_tree.buildIndex(flann::Matrix<float>(flann_data_storage.data(), triangle_center_cloud_polar->size(), 2));
    }
    
    void knnSearch(pcl::PointXY searchPoint, std::vector<int>& search_indices, std::vector<float>& search_dists, int K)
    {
        // convert to vector
        std::vector<float> query_point = {searchPoint.x, searchPoint.y};

        // intialize
        std::vector<std::vector<int>> list_of_search_indices(1, std::vector<int>(K));
        std::vector<std::vector<float>> list_of_search_dists(1, std::vector<float>(K));

        // search
        flann_tree.knnSearch(flann::Matrix<float>(query_point.data(), 1, 2), list_of_search_indices, list_of_search_dists, K, flann::SearchParams(-1, 0));

        // extract
        search_indices = list_of_search_indices[0];
        search_dists = list_of_search_dists[0];
    }

    void addPoints(pcl::PointXY new_point_polar)
    {
        // convert to vector
        std::vector<float> new_triangle_center_polar;
        new_triangle_center_polar.push_back(new_point_polar.x);
        new_triangle_center_polar.push_back(new_point_polar.y);

        // add
        flann_tree.addPoints(flann::Matrix<float>(new_triangle_center_polar.data(), 1, 2));

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

    flann::Index<flann::L2<float>> flann_tree;
};



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

    // ------------------------------ dataloader
    // old cloud
    int i1 = 0;
    typename pcl::PointCloud<InputPointT>::Ptr old_cloud = data_loader.get_cloud(i1);
    Eigen::Affine3d old_pose = data_loader.get_pose(i1);
    // new cloud
    int i2 = 50;
    typename pcl::PointCloud<InputPointT>::Ptr new_cloud = data_loader.get_cloud(i2);
    Eigen::Affine3d new_pose = data_loader.get_pose(i2);

    // ------------------------------ algorithm
    // delaunay triangulation old cloud
    delaunator::Delaunator old_d = obtain_triangulation<InputPointT> (old_cloud);
    std::map<int, std::vector<int>> triangle_map = d_to_triangle_map(old_d);

    // compute triangle center cloud and polar
    typename pcl::PointCloud<InputPointT>::Ptr old_cloud_local = transform_cloud_to_frame<InputPointT>(old_cloud, old_pose, new_pose);
    pcl::PointCloud<pcl::PointXYZ>::Ptr triangle_center_cloud = compute_triangle_center_cloud<InputPointT>(old_cloud_local, triangle_map);
    pcl::PointCloud<pcl::PointXY>::Ptr triangle_center_cloud_polar = compute_2d_polar_cloud<pcl::PointXYZ>(triangle_center_cloud);

    // initialize flann
    flann2d flann_tree;
    flann_tree.set_input(triangle_center_cloud_polar);

    // compute polar coordinate of each new point
    pcl::PointCloud<pcl::PointXY>::Ptr new_point_polar_cloud = compute_2d_polar_cloud<InputPointT>(new_cloud);

    // search for each new point
    for (std::size_t i = 0; i < new_point_polar_cloud->size(); i++)
    {
        // get new point (vector)
        InputPointT p_new_point = new_cloud->points[i];
        Eigen::Vector3f v_new_point = p_new_point.getVector3fMap();

        // get new point polar
        pcl::PointXY searchPoint = new_point_polar_cloud->points[i];

        // output stat
        std::cout << "processing point " << i << " out of " << new_point_polar_cloud->size() << std::endl;

        // search for nearest triangles
        int K = 4;
        std::vector<int> knn_indices;
        std::vector<float> knn_dists;
        flann_tree.knnSearch(searchPoint, knn_indices, knn_dists, K);
 
        // for each searched triangle center
        std::vector<int> intersected_triangle_indices;
        std::vector<float> intersected_triangle_distances;
        for (int knn_index : knn_indices)
        {
            // get vertices index
            std::vector<int> vertex_indices = triangle_map[knn_index];

            // get vertices vector
            Eigen::Vector3f v0 = old_cloud_local->points[vertex_indices[0]].getVector3fMap();
            Eigen::Vector3f v1 = old_cloud_local->points[vertex_indices[1]].getVector3fMap();
            Eigen::Vector3f v2 = old_cloud_local->points[vertex_indices[2]].getVector3fMap();

            // compute ray triangle intersection
            Eigen::Vector3f intersection = ray_triangle_intersection(v_new_point, v_new_point.normalized(), v0, v1, v2);
            float distance = (intersection - v_new_point).norm();

            // check if intersection is inside the triangle
            bool inside = is_inside_triangle(v0, v1, v2, intersection);
            if (!inside) continue;
            
            // store the intersection
            intersected_triangle_indices.push_back(knn_index);
            intersected_triangle_distances.push_back(distance);
        }

        // skip if no intersection
        if (intersected_triangle_indices.size() == 0) continue;

        // find the closest triangle (out of the k searched ones)
        int min_index = std::distance(intersected_triangle_distances.begin(), std::min_element(intersected_triangle_distances.begin(), intersected_triangle_distances.end()));
        int min_triangle_index = intersected_triangle_indices[min_index];

        // --- update point --- 
        // add 
        old_cloud_local->push_back(p_new_point);
        int new_point_index = old_cloud_local->size() - 1;

        // --- update triangle map --- 
        // vertices index
        std::vector<int> vertices_index = triangle_map[min_triangle_index];
        // new triangles
        std::vector<int> new_triangle1 = {vertices_index[0], vertices_index[1], new_point_index};
        std::vector<int> new_triangle2 = {vertices_index[1], vertices_index[2], new_point_index};
        std::vector<int> new_triangle3 = {vertices_index[2], vertices_index[0], new_point_index};
        // remove
        triangle_map.erase(min_triangle_index);
        // add
        triangle_map[flann_tree.flann_last_id+1] = new_triangle1;
        triangle_map[flann_tree.flann_last_id+2] = new_triangle2;
        triangle_map[flann_tree.flann_last_id+3] = new_triangle3;

        // --- update flann triangle polar data ---
        // vertices point
        InputPointT p0 = old_cloud_local->points[vertices_index[0]];
        InputPointT p1 = old_cloud_local->points[vertices_index[1]];
        InputPointT p2 = old_cloud_local->points[vertices_index[2]];
        // center point
        pcl::PointXYZ new_triangle_center1 = compute_triangle_center_point(p0, p1, p_new_point);
        pcl::PointXYZ new_triangle_center2 = compute_triangle_center_point(p1, p2, p_new_point);
        pcl::PointXYZ new_triangle_center3 = compute_triangle_center_point(p2, p0, p_new_point);
        // center polar point
        pcl::PointXY new_triangle_center1_polar = compute_2d_polar_point(new_triangle_center1);
        pcl::PointXY new_triangle_center2_polar = compute_2d_polar_point(new_triangle_center2);
        pcl::PointXY new_triangle_center3_polar = compute_2d_polar_point(new_triangle_center3);
        // remove
        flann_tree.removePoint(min_triangle_index);
        // add
        flann_tree.addPoints(new_triangle_center1_polar);
        flann_tree.addPoints(new_triangle_center2_polar);
        flann_tree.addPoints(new_triangle_center3_polar);
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