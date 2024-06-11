#include <pcl/common/transforms.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

#include "point_type/VilensPointT.hpp"
#include "eye_patch/DataLoader.hpp"
#include "eye_patch/utilities.hpp"

#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/Storage.hpp"
#include "MeshObject/TriangleBVH.hpp"
#include "MeshObject/RRSTree.hpp"

// application class
template <typename PointT>
class Application
{
public:

    Eigen::Matrix3d merge_covariances_of_surfaces(std::shared_ptr<Surface> surface1, std::shared_ptr<Surface> surface2) 
    {
        const Eigen::Matrix3d& cov1 = surface1->get_covariance();
        const Eigen::Matrix3d& cov2 = surface2->get_covariance();
        const Eigen::Vector3d& mean1 = surface1->get_mean();
        const Eigen::Vector3d& mean2 = surface2->get_mean();
        int size1 = surface1->get_total_point_size();
        int size2 = surface2->get_total_point_size();
        return merge_covariances(cov1, cov2, mean1, mean2, size1, size2);
    }

    // try merge sets
    void try_merge_surfaces(std::set<std::shared_ptr<Surface>>& surfaces_to_merge)
    {
        while (true) 
        {
            // get all possible pairs to merge
            std::set<std::pair<std::shared_ptr<Surface>, std::shared_ptr<Surface>>> surface_pairs;
            for (std::shared_ptr<Surface> surface1 : surfaces_to_merge) 
            {
                for (std::shared_ptr<Surface> surface2 : surfaces_to_merge) 
                {
                    if (surface1 >= surface2) continue;
                    surface_pairs.insert(std::make_pair(surface1, surface2));
                }
            }
            
            // try to merge pairs
            bool again = false;
            for (const auto& pairs : surface_pairs) 
            {
                // skip if can't merge
                Eigen::Matrix3d covariance_matrix = merge_covariances_of_surfaces(pairs.first, pairs.second);
                double eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(covariance_matrix).eigenvalues()[0];
                if (eigenvalue > merged_eigenvalue_threshold) continue;

                // merge surfaces
                surfaces_to_merge.erase(pairs.second);
                pairs.first->merge_surface(pairs.second);

                // once merged, restart
                again = true;
                break;
            }
            if (!again) break;
        }
    }

    void add_point_by_radius_search(const Eigen::Vector3d& thisPointVEC, const Eigen::Vector3d& thisPointOriginVEC)
    {
        // if empty, can not set up radius search, add point to new set
        if (!storage_->can_reverse_radius_search())
        {
            std::shared_ptr<Surface> new_surface = storage_->add_surface();
            std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            new_surface->connect(new_vertex);
            return;
        }

        // perform rrstree radius search
        std::map<int, double> point_to_radius_map;
        std::set<std::shared_ptr<Vertex>> searched_boundary_vertices_set = storage_->reverse_radius_search(thisPointVEC);

        // if no searched results, add point to new set
        if (searched_boundary_vertices_set.size() == 0)
        {
            std::shared_ptr<Surface> new_surface = storage_->add_surface();
            std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            new_surface->connect(new_vertex);
            return;
        }

        // from searched points identify neighboring sets
        std::set<std::shared_ptr<Surface>> neighboring_surfaces; 
        for (std::shared_ptr<Vertex> vertex : searched_boundary_vertices_set)
        {
            neighboring_surfaces.insert(vertex->get_surface());
        }

        // try merge neighboring sets
        try_merge_surfaces(neighboring_surfaces);

        // after merging, we are left with sets that should have different normals
        // each point in the neighboring set should reduce their search radius to the closest set
        // [todo]
        // for each point, compute the shortest distance to another point that is in a different set, and that different set have enough points
        // if the distance is less than its original radius, reduce its original radius
        // group searched vertex by surface
        std::map<std::shared_ptr<Surface>, std::set<std::shared_ptr<Vertex>>> surface_to_searched_vertices_map;
        for (std::shared_ptr<Vertex> vertex : searched_boundary_vertices_set)
        {
            surface_to_searched_vertices_map[vertex->get_surface()].insert(vertex);
        }
        // for each surfaces
        for (const auto& pair : surface_to_searched_vertices_map)
        {
            std::shared_ptr<Surface> surface = pair.first;
            const std::set<std::shared_ptr<Vertex>>& searched_boundary_vertices_set = pair.second;

            // skip if small surface
            if (surface->get_total_point_size() < fit_plane_threshold) continue;

            // for each vertex in the surface
            for (std::shared_ptr<Vertex> this_vertex : searched_boundary_vertices_set)
            {
                // smallest distance
                double smallest_distance = std::numeric_limits<double>::max();

                // for each other surface
                for (const auto& other_pair : surface_to_searched_vertices_map)
                {
                    // skip if same surface
                    if (other_pair.first == surface) continue;

                    // skip if small set
                    if (other_pair.first->get_total_point_size() < fit_plane_threshold) continue;

                    // for each point in the other surface
                    for (std::shared_ptr<Vertex> other_vertex : other_pair.second)
                    {
                        // compute distance
                        double distance = (other_vertex->get_position() - this_vertex->get_position()).norm();

                        // update radius
                        if (distance < smallest_distance) smallest_distance = distance;
                    }
                }

                // adjust radius of this vertex
                if (smallest_distance < this_vertex->get_radius()) this_vertex->set_reverse_radius_search_radius(smallest_distance);
            }
        }

        // split neighboring sets into sets with plane and sets without plane (by size)
        std::set<std::shared_ptr<Surface>> surfaces_with_plane;
        std::set<std::shared_ptr<Surface>> surfaces_without_plane;
        for (std::shared_ptr<Surface> surface : neighboring_surfaces)
        {
            if (surface->get_total_point_size() > fit_plane_threshold) surfaces_with_plane.insert(surface);
            else surfaces_without_plane.insert(surface);
        }
        
        // for sets with plane, compute the point to set intersection distance
        std::map<std::shared_ptr<Surface>, double> surface_distance_map;
        for (std::shared_ptr<Surface> surface : surfaces_with_plane)
        {
            // compute (not using normal from combined points, could implement in future)
            double distance = surface->compute_point_to_surface_distance(thisPointOriginVEC, thisPointVEC);

            // store
            surface_distance_map[surface] = std::fabs(distance);
        }

        // extract the set within distance threshold
        std::set<std::shared_ptr<Surface>> surfaces_within_threshold;
        for (const auto& pair : surface_distance_map)
        {
            if (pair.second < distance_threshold) surfaces_within_threshold.insert(pair.first);
        }

        // from the sets within threshold, find the set that is closest to the point
        std::shared_ptr<Surface> closest_surface = std::make_shared<Surface>();
        double closest_distance = std::numeric_limits<double>::max();
        for (std::shared_ptr<Surface> surface : surfaces_within_threshold)
        {
            double distance = surface_distance_map.at(surface);

            // update if closer
            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest_surface = surface;
            }
        }
        if (!closest_surface->is_expired() && closest_distance < distance_threshold)
        {
            // for the sets not selected as closest, update their searched points' radius
            // for all searched points
            for (std::shared_ptr<Vertex> vertex : searched_boundary_vertices_set)
            {
                // if it is in a set with plane
                if (surfaces_with_plane.find(vertex->get_surface()) != surfaces_with_plane.end())
                {
                    // and the set with plane is not the closest set
                    if (vertex->get_surface() != closest_surface)
                    {
                        // reduce their searched points' radius
                        double reduced_radius = (thisPointVEC - vertex->get_position()).norm();

                        if (reduced_radius < vertex->get_radius())
                        {
                            vertex->set_reverse_radius_search_radius(reduced_radius);
                        }
                    }
                }
            }

            // add the point to the closest set
            std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            closest_surface->connect(new_vertex, searched_boundary_vertices_set);
            return;
        }

        // else, find the set that is nearest
        std::shared_ptr<Surface> nearest_surface = std::make_shared<Surface>();
        double nearest_distance = std::numeric_limits<double>::max();
        for (std::shared_ptr<Surface> surface : surfaces_without_plane)
        {
            const Eigen::Vector3d& mean = surface->get_mean();
            double distance = (thisPointVEC - mean).norm();
            if (distance < nearest_distance)
            {
                nearest_distance = distance;
                nearest_surface = surface;
            }
        }
        if (!nearest_surface->is_expired())
        {
            std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
            nearest_surface->connect(new_vertex, searched_boundary_vertices_set);
            return;
        }

        // else, add the point to a new set
        std::shared_ptr<Surface> new_surface = storage_->add_surface();
        std::shared_ptr<Vertex> new_vertex = storage_->add_vertex(thisPointOriginVEC, thisPointVEC);
        new_surface->connect(new_vertex);
        return;
    }

    void load_point_cloud()
    {
        if (ith_cloud < 0)
        {
            std::cout << "reached the first pointcloud" << std::endl;
            ith_cloud = 0;
        }
        if (ith_cloud >= data_loader.size())
        {
            std::cout << "reached the last pointcloud" << std::endl;
            ith_cloud = data_loader.size() - 1;
        }
        
        typename pcl::PointCloud<PointT>::Ptr pointcloud_local = data_loader.get_cloud(ith_cloud);
        Eigen::Affine3d pose = data_loader.get_pose(ith_cloud);
        pointcloud = transform_cloud_to_global<PointT>(pointcloud_local, pose);
        origin = pose.translation();
        ith_size = pointcloud->size() * pointcloud_fraction;

        std::cout << "loaded pointcloud " << ith_cloud << " with " << pointcloud->size() << " points" << std::endl;

        if (shuffle_pointcloud) 
        {
            // shuffle the pointcloud
            std::random_shuffle(pointcloud->points.begin(), pointcloud->points.end());
        }
    }

    void step()
    {        
        // get point
        Eigen::Vector3d thisPointVEC = pointcloud->points[ith_point].getVector3fMap().cast<double>();
        Eigen::Vector3d thisPointOriginVEC = origin;
        ith_point ++;
    
        // ------------- add point by triangle intersection

        // get list of intersected triangle by the point
        std::set<std::shared_ptr<Face>> searched_faces = storage_->face_intersection_search(thisPointOriginVEC, thisPointVEC); // may include deleted triangles

        // group the faces by surface
        std::map<std::shared_ptr<Surface>, std::set<std::shared_ptr<Face>>> surface_to_searched_faces_map;
        for (std::shared_ptr<Face> face : searched_faces)
        {
            surface_to_searched_faces_map[face->get_surface()].insert(face);
        }

        // compute the intersection distance to the sets (distance measured in plane normal direction)
        std::map<std::shared_ptr<Surface>, double> surface_to_point_distance_map;
        for (const auto& pair : surface_to_searched_faces_map)
        {
            std::shared_ptr<Surface> surface = pair.first;
            double distance = surface->compute_point_to_surface_distance(thisPointOriginVEC, thisPointVEC);
            surface_to_point_distance_map[surface] = distance;
        }

        // split the sets into three categories
        std::set<std::shared_ptr<Surface>> surface_with_point_before_it;
        std::set<std::shared_ptr<Surface>> surface_with_point_within_it;
        std::set<std::shared_ptr<Surface>> surface_with_point_behind_it;
        double split_distance_threshold = distance_threshold;
        for (const auto& pair : surface_to_point_distance_map)
        {
            std::shared_ptr<Surface> surface = pair.first;
            double distance = pair.second;
            if (distance > split_distance_threshold) 
            {
                surface_with_point_before_it.insert(surface);
            }
            else if (distance < -split_distance_threshold) 
            {
                surface_with_point_behind_it.insert(surface);
            }
            else 
            {
                surface_with_point_within_it.insert(surface);
            }
        }
        
        // process point behind set set

        // get the set of triangles that are penetrated by the point
        std::set<std::shared_ptr<Face>> penetrated_faces;
        for (std::shared_ptr<Surface> surface : surface_with_point_behind_it)
        {
            penetrated_faces.insert(surface_to_searched_faces_map.at(surface).begin(), surface_to_searched_faces_map.at(surface).end());
        }

        // collect the list of points that are within the penetrated triangles, and isolated by the triangle
        for (std::shared_ptr<Face> face : penetrated_faces)
        {
            // for now, just delete the face
            storage_->delete_face(face);
        }

        // // re add the list points by radius search
        // // todo - while avoid covering the new point
        // // process point in the queue free_points_vector3d_and_origin_vector3d
        // while (!free_points_queue.empty())
        // {
        //     std::pair<Eigen::Vector3d, Eigen::Vector3d> pair = free_points_queue.front(); free_points_queue.pop();

        //     std::cout << "---------------------------------------------------------------------------------------------------- re-adding point by radius search" << std::endl;
        //     Eigen::Vector3d pointVEC = pair.first;
        //     Eigen::Vector3d pointOriginVEC = pair.second;
        //     int pointID = getNewPointID();
        //     add_point_by_radius_search(pointID, pointVEC, pointOriginVEC);
        // }
        
        // process point within set set
        bool point_added_to_surface = false;

        // try merge them
        std::set<std::shared_ptr<Surface>> merged_surface = surface_with_point_within_it;
        try_merge_surfaces(merged_surface);        

        // find the set with the smallest distance
        std::shared_ptr<Surface> smallest_surface = std::make_shared<Surface>();
        double smallest_distance = std::numeric_limits<double>::max();
        for (std::shared_ptr<Surface> surface : merged_surface)
        {
            // update if smaller
            double distance = surface->compute_point_to_surface_distance(thisPointOriginVEC, thisPointVEC);
            if (std::abs(distance) < smallest_distance)
            {
                smallest_distance = distance;
                smallest_surface = surface;
            }
        }

        // if the smallest set is within threshold, add the point to the set
        if (!smallest_surface->is_expired() && std::abs(smallest_distance) < distance_threshold)
        {
            // from searched_faces find the first face that belongs to the smallest surface
            for (const std::shared_ptr<Face>& face : searched_faces)
            {
                // if face is expired
                if (face->is_expired()) continue;

                // if not the same surface
                if (face->get_surface() != smallest_surface) continue;
                
                // add point as interior point
                storage_->add_interior_point(face, thisPointVEC, thisPointOriginVEC);

                std::cout << ith_point << " / " << ith_size << " of pointcloud " << ith_cloud << " added to set " << smallest_surface->get_id() << std::endl;
                point_added_to_surface = true;
                break;
            }
        }
        if (!point_added_to_surface)
        {
            add_point_by_radius_search(thisPointVEC, thisPointOriginVEC);
            std::cout << ith_point << " / " << ith_size << " of pointcloud " << ith_cloud << " added by radius search" << std::endl;
        }


        // todo - process point before set set


        // check if end of point cloud
        if (ith_point == ith_size) 
        {
            // next cloud
            ith_cloud += 1;
            ith_point = 0;
            load_point_cloud();
        }
    }

    void loop()
    {
        step();

        // finish all points in step
        while (ith_point != 0)
        {
            step();
        }
    }

    Application() 
    {
        // // ----------------------- DATA
        std::map<std::string, std::pair<std::string, std::string>> dataset_map;
        dataset_map["room"] = std::make_pair(
            "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_clouds/",
            "/home/jiahao/datasets/bag2pcd_output/mission2_reverse/slam_poses/slam_poss_graph.slam"
        );
        dataset_map["osney"] = std::make_pair(
            "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_clouds/",
            "/home/jiahao/datasets/osney power station/2024-03-26_13-47-27_rec004_osney_power_station/slam_pose_graph.slam"
        );
        dataset_map["blenheim"] = std::make_pair(
            "/home/jiahao/datasets/2024-03-14-09-09-02-lenord-walk-for-lintong/individual_clouds/",
            "/home/jiahao/datasets/2024-03-14-09-09-02-lenord-walk-for-lintong/slam_pose_graph.g2o"
        );

        // // ----------------------- PARAMETERS
        distance_threshold = 0.05;
        fit_plane_threshold = 10; // may cause error if below 3
        merged_eigenvalue_threshold = 15e-5;
        

        dataset = "room";
        // dataset = "osney";
        // dataset = "blenheim";
        ith_cloud = 50;
        ith_point = 0;
        shuffle_pointcloud = false;
        pointcloud_fraction = 1;
        distance_to_radius_ratio = tan(4 * M_PI / 180);

        storage_ = std::make_shared<Storage>();        

        // // ----------------------- INITIALIZATION
        data_loader.load_dataset(dataset_map.at(dataset).first, dataset_map.at(dataset).second);
        load_point_cloud();
    }

    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map()
    {
        return vertex_to_cloud_indices_map;
    } 

    const std::set<std::shared_ptr<Face>>& get_faces() {return storage_->get_faces();};
    const std::set<std::shared_ptr<Edge>>& get_edges() {return storage_->get_edges();};
    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices() {return storage_->get_rrs_vertices();};
    
    std::set<std::shared_ptr<Edge>> get_boundary_edges() 
    {
        // initialize
        std::set<std::shared_ptr<Edge>> boundary_edges;

        // process
        for (const std::shared_ptr<Edge>& edge : storage_->get_edges())
        {
            if (edge->is_boundary()) 
            {
                boundary_edges.insert(edge);
            }
        }

        // return
        return boundary_edges;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_colored_cloud()
    {
        vertex_to_cloud_indices_map.clear();

        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (const std::shared_ptr<Vertex>& vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            const Eigen::Vector3d& position = vertex->get_position();
            const std::tuple<int, int, int>& color = vertex->get_surface()->get_color();
            point.x = position[0];
            point.y = position[1];
            point.z = position[2];
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
            vertex_to_cloud_indices_map[vertex] = cloud->size() - 1;
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr point_to_vector3d_set_distance_cloud()
    {
        vertex_to_cloud_indices_map.clear();

        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (const std::shared_ptr<Vertex>& vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            const Eigen::Vector3d& origin = vertex->get_origin();
            const Eigen::Vector3d& position = vertex->get_position();
            point.x = position[0];
            point.y = position[1];
            point.z = position[2];
            double value = vertex->get_surface()->compute_point_to_surface_distance(origin, position) / 0.05;
            std::tuple<int, int, int> color = valueToJet(value);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
            vertex_to_cloud_indices_map[vertex] = cloud->size() - 1;
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_colored_cloud()
    {
        vertex_to_cloud_indices_map.clear();

        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (std::shared_ptr<Vertex> vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            Eigen::Vector3d projected_position = vertex->get_projected_position();
            const std::tuple<int, int, int>& color = vertex->get_surface()->get_color();
            point.x = projected_position[0];
            point.y = projected_position[1];
            point.z = projected_position[2];
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
            vertex_to_cloud_indices_map[vertex] = cloud->size() - 1;
        }
        return cloud;
    }

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr projected_point_to_vector3d_set_distance_cloud()
    {
        vertex_to_cloud_indices_map.clear();
        
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
        for (const std::shared_ptr<Vertex>& vertex : storage_->get_vertices())
        {
            pcl::PointXYZRGB point;
            const Eigen::Vector3d& origin = vertex->get_origin();
            const Eigen::Vector3d& position = vertex->get_position();
            Eigen::Vector3d projected_position = vertex->get_projected_position();
            point.x = projected_position[0];
            point.y = projected_position[1];
            point.z = projected_position[2];
            double value = vertex->get_surface()->compute_point_to_surface_distance(origin, position) / 0.05;
            std::tuple<int, int, int> color = valueToJet(value);
            point.r = std::get<0>(color);
            point.g = std::get<1>(color);
            point.b = std::get<2>(color);
            cloud->push_back(point);
            vertex_to_cloud_indices_map[vertex] = cloud->size() - 1;
        }
        return cloud;
    }

    void change_color()
    {
        for (const std::shared_ptr<Surface>& surface : storage_->get_surfaces())
        {
            surface->set_random_color();
        }
    }

    int ith_cloud;
    std::size_t ith_point = 0;
    std::size_t ith_size = 0;

private:
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map;

    // storage
    std::shared_ptr<Storage> storage_;

    // data
    DataLoader<VilensPointT> data_loader;
    typename pcl::PointCloud<VilensPointT>::Ptr pointcloud;
    Eigen::Vector3d origin;
    

    // settings
    double distance_threshold;
    std::size_t fit_plane_threshold; // may cause error if below 3
    double merged_eigenvalue_threshold;
    bool shuffle_pointcloud;
    double pointcloud_fraction;
    std::string dataset;
    double distance_to_radius_ratio;
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

        // // ------------------------------ parameters
        number_of_spheres_to_display = 60;

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
    
    std::vector<std::string> sphere_name_list;
    bool show_sphere = false;

    int number_of_spheres_to_display;

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
        std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map = app_.get_vertex_to_cloud_indices_map();
        std::set<std::shared_ptr<Face>> faces = app_.get_faces();
        std::set<std::shared_ptr<Edge>> boundary_edges = app_.get_boundary_edges();

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
            for (const std::shared_ptr<Face>& face : faces)
            {
                pcl::Vertices triangle;
                triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face->get_vertex(0)));
                triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face->get_vertex(1)));
                triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face->get_vertex(2)));
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
            for (const std::shared_ptr<Edge>& edge : boundary_edges)
            {
                pcl::Vertices boundary_edge;
                boundary_edge.vertices.push_back(vertex_to_cloud_indices_map.at(edge->get_vertex(0)));
                boundary_edge.vertices.push_back(vertex_to_cloud_indices_map.at(edge->get_vertex(1)));
                boundary_mesh.polygons.push_back(boundary_edge);
            }
            viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "boundary_edges");
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "boundary_edges");
        }

        // boundary points spheres
        for (const std::string& sphere_name : sphere_name_list) viewer_->removeShape(sphere_name);
        sphere_name_list.clear();
        if (show_sphere)
        {
            std::vector<std::shared_ptr<Vertex>> boundary_vertices = app_.get_rrs_vertices();
            // sort
            std::sort(boundary_vertices.begin(), boundary_vertices.end());

            // for the last 20
            for (int i = 0; i < std::min(number_of_spheres_to_display, (int)boundary_vertices.size()); i++)
            {
                std::shared_ptr<Vertex> boundary_vertex = boundary_vertices[boundary_vertices.size() - 1 - i];
                std::string sphere_name = "boundary_point_" + std::to_string(boundary_vertex->get_id());
                sphere_name_list.push_back(sphere_name);
                const Eigen::Vector3d& position = boundary_vertex->get_position();
                viewer_->addSphere(pcl::PointXYZ(position[0], position[1], position[2]), boundary_vertex->get_radius(), 1, 1, 1, sphere_name);
            }
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
        // std::cout << "key pressed: [" << event.getKeySym() << "]" << std::endl;
        

        if (event.getKeySym() == "space" && event.keyDown())
        {
            app_.loop();
            update_display();
        }
        if (event.getKeySym() == "KP_Insert" && event.keyDown())
        {
            for (int i = 0; i < 100; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "KP_Delete" && event.keyDown())
        {
            for (int i = 0; i < 1000; i++) app_.step();
            update_display();
        }
        if (event.getKeySym() == "KP_Enter" && event.keyDown())
        {
            app_.loop();
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
        if (event.getKeySym() == "KP_Next" && event.keyDown())
        {
            app_.ith_cloud += 1;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_End" && event.keyDown())
        {
            app_.ith_cloud -= 1;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Right" && event.keyDown())
        {
            app_.ith_cloud += 10;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Left" && event.keyDown())
        {
            app_.ith_cloud -= 10;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Prior" && event.keyDown())
        {
            app_.ith_cloud += 100;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "KP_Home" && event.keyDown())
        {
            app_.ith_cloud -= 100;
            app_.load_point_cloud();
        }
        if (event.getKeySym() == "m" && event.keyDown())
        {
            show_sphere = !show_sphere;
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