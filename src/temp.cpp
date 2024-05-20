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


int main() {
    // Define your mesh vertices and triangles
    // vertices id to vector3d map
    std::map<int, Eigen::Vector3d> vertex_to_vector3d_map;
    vertex_to_vector3d_map[0] = Eigen::Vector3d(1, 1, 1);
    vertex_to_vector3d_map[1] = Eigen::Vector3d(-1, 1, 1);
    vertex_to_vector3d_map[2] = Eigen::Vector3d(-1, -1, 1);
    vertex_to_vector3d_map[3] = Eigen::Vector3d(1, -1, 1);
    vertex_to_vector3d_map[4] = Eigen::Vector3d(1, 1, -1);
    vertex_to_vector3d_map[5] = Eigen::Vector3d(-1, 1, -1);
    vertex_to_vector3d_map[6] = Eigen::Vector3d(-1, -1, -1);
    vertex_to_vector3d_map[7] = Eigen::Vector3d(1, -1, -1);

    // triangle id to indices map
    std::map<int, std::array<int, 3>> triangle_to_indices_map;
    triangle_to_indices_map[0] = {0, 1, 2};
    triangle_to_indices_map[1] = {0, 3, 2};
    triangle_to_indices_map[2] = {0, 3, 7};
    triangle_to_indices_map[3] = {0, 4, 7};
    triangle_to_indices_map[4] = {0, 1, 5};
    triangle_to_indices_map[5] = {0, 4, 5};
    triangle_to_indices_map[6] = {6, 5, 1};
    triangle_to_indices_map[7] = {6, 2, 1};
    triangle_to_indices_map[8] = {6, 5, 4};
    triangle_to_indices_map[9] = {6, 7, 4};
    triangle_to_indices_map[10] = {6, 2, 3};
    triangle_to_indices_map[11] = {6, 7, 3};

    // triangle id list
    std::vector<int> triangle_id_list;
    for (const auto& pair : triangle_to_indices_map)
    {
        triangle_id_list.push_back(pair.first);
    }

    // Build the BVH
    auto bvhRoot = buildBVH(vertex_to_vector3d_map, triangle_to_indices_map, triangle_id_list, 0, triangle_id_list.size());

    // Define the ray
    Eigen::Vector3d rayOrigin(0.5, 0.5, 0.5);
    Eigen::Vector3d rayDirection(1, 1, 1);
    std::vector<Eigen::Vector3d> intersectionPointList;

    // Perform intersection test
    bool intersect = intersectBVH(bvhRoot, rayOrigin, rayDirection, vertex_to_vector3d_map, triangle_to_indices_map, intersectionPointList);

    if (intersect) {
        for (const auto& intersectionPoint : intersectionPointList) {
            std::cout << "Intersection point: " << intersectionPoint.transpose() << std::endl;
        }
    }


    // viewer
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1);
    
    // make mesh
    pcl::PolygonMesh mesh;
    // add points to mesh.cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr new_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    for (const auto& pair : vertex_to_vector3d_map)
    {
        pcl::PointXYZ point;
        point.x = pair.second[0];
        point.y = pair.second[1];
        point.z = pair.second[2];
        new_cloud->push_back(point);
    }


    // add points
    pcl::toPCLPointCloud2(*new_cloud, mesh.cloud); 

    // add triangle to mesh
    for (const auto& pair : triangle_to_indices_map)
    {
        pcl::Vertices triangle;
        triangle.vertices.push_back(pair.second[0]);
        triangle.vertices.push_back(pair.second[1]);
        triangle.vertices.push_back(pair.second[2]);
        mesh.polygons.push_back(triangle);
    }
    // add mesh
    viewer->addPolylineFromPolygonMesh(mesh, "polyline");

    // add ray
    pcl::PointXYZ ray_origin;
    ray_origin.x = rayOrigin[0];
    ray_origin.y = rayOrigin[1];
    ray_origin.z = rayOrigin[2];
    pcl::PointXYZ ray_end;
    double length = 2;
    ray_end.x = rayOrigin[0] + rayDirection[0] * length;
    ray_end.y = rayOrigin[1] + rayDirection[1] * length;
    ray_end.z = rayOrigin[2] + rayDirection[2] * length;
    viewer->addLine(ray_origin, ray_end, 1, 1, 0, "ray");

    // add ray origin
    viewer->addSphere(ray_origin, 0.01, 1, 0, 0, "ray_origin");

    // add intersection
    if (intersect)
    {
        pcl::PointXYZ intersection;
        intersection.x = intersectionPointList[0][0];
        intersection.y = intersectionPointList[0][1];
        intersection.z = intersectionPointList[0][2];
        viewer->addSphere(intersection, 0.01, 0, 1, 0, "intersection");
    }

    // spin
    viewer->spin();

    return 0;
}







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