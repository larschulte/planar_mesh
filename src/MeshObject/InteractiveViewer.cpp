#include "MeshObject/InteractiveViewer.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include <unordered_set>

template class InteractiveViewer<VilensPointT>;

template <typename PointT>
InteractiveViewer<PointT>::InteractiveViewer(Application<PointT>& app) 
    : 
    app_(app),
    viewer_(new pcl::visualization::PCLVisualizer ("3D Viewer"))
{   
    viewer_->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    viewer_->setBackgroundColor (0, 0, 0);
    viewer_->initCameraParameters();
    viewer_->addCoordinateSystem(1);
    viewer_->registerKeyboardCallback(&InteractiveViewer::keyboard_callback, *this, nullptr);
    viewer_->spin();
}

template <typename PointT>
void InteractiveViewer<PointT>::update_display()
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr vertex_pointcloud = app_.compute_vertex_point_pointcloud(settings_.show_projected_point, settings_.show_error_color);
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map = app_.get_vertex_to_cloud_indices_map();
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = app_.get_faces();
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> boundary_edges = app_.get_boundary_edges();

    // pointcloud
    viewer_->removeShape("point_cloud");
    if (settings_.show_pointcloud)
    {
        viewer_->addPointCloud<pcl::PointXYZRGB>(vertex_pointcloud, "point_cloud");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "point_cloud");
    }

    // interior points
    viewer_->removeShape("interior_points");
    if (settings_.show_interior_points)
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr interior_point_cloud = app_.compute_interior_point_pointcloud(settings_.show_projected_point, settings_.show_error_color);
        viewer_->addPointCloud<pcl::PointXYZRGB>(interior_point_cloud, "interior_points");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "interior_points");
    }

    // generic points
    viewer_->removeShape("generic_points");
    if (settings_.show_generic_points)
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr generic_point_cloud = app_.compute_generic_point_pointcloud();
        viewer_->addPointCloud<pcl::PointXYZRGB>(generic_point_cloud, "generic_points");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "generic_points");
    }

    // triangle
    viewer_->removeShape("triangle_mesh");
    if (settings_.show_triangle)
    {
        pcl::PolygonMesh triangle_mesh;
        pcl::toPCLPointCloud2(*vertex_pointcloud, triangle_mesh.cloud);
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

    // boundary edge
    viewer_->removeShape("boundary_edges");        
    if (settings_.show_edge)
    {
        pcl::PolygonMesh boundary_mesh;
        pcl::toPCLPointCloud2(*vertex_pointcloud, boundary_mesh.cloud);
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

    // spheres
    for (const std::string& sphere_name : sphere_name_list) viewer_->removeShape(sphere_name);
    sphere_name_list.clear();
    if (settings_.show_sphere)
    {
        std::vector<std::shared_ptr<Vertex>> boundary_vertices = app_.get_rrs_vertices();
        std::sort(boundary_vertices.begin(), boundary_vertices.end());
        for (int i = 0; i < std::min(settings_.number_of_spheres_to_display, (int)boundary_vertices.size()); i++)
        {
            std::shared_ptr<Vertex> boundary_vertex = boundary_vertices[boundary_vertices.size() - 1 - i];
            std::string sphere_name = "boundary_point_" + std::to_string(boundary_vertex->get_id());
            sphere_name_list.push_back(sphere_name);
            const Eigen::Vector3d& position = boundary_vertex->get_position();
            viewer_->addSphere(pcl::PointXYZ(position[0], position[1], position[2]), boundary_vertex->get_radius(), 1, 1, 1, sphere_name);
        }
    }

    // wireframe
    if (settings_.show_wireframe)
    {
        viewer_->setRepresentationToWireframeForAllActors();
    }
    else
    {
        viewer_->setRepresentationToSurfaceForAllActors();
    }
}

template <typename PointT>
void InteractiveViewer<PointT>::keyboard_callback(const pcl::visualization::KeyboardEvent &event, void*) 
{
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
        settings_.show_pointcloud = !settings_.show_pointcloud;
        update_display();
    }
    if (event.getKeySym() == "period" && event.keyDown())
    {
        settings_.show_edge = !settings_.show_edge;
        update_display();
    }
    if (event.getKeySym() == "slash" && event.keyDown())
    {
        settings_.show_triangle = !settings_.show_triangle;
        update_display();
    }
    if (event.getKeySym() == "a" && event.keyDown())
    {
        settings_.show_projected_point = !settings_.show_projected_point;
        update_display();
    }
    if (event.getKeySym() == "z" && event.keyDown())
    {
        settings_.show_error_color = !settings_.show_error_color;
        update_display();
    }
    if (event.getKeySym() == "v" && event.keyDown())
    {
        settings_.show_wireframe = !settings_.show_wireframe;
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
        settings_.show_sphere = !settings_.show_sphere;
        update_display();
    }
    if (event.getKeySym() == "b" && event.keyDown())
    {
        app_.refine_surfaces();
        update_display();
    }
    if (event.getKeySym() == "n" && event.keyDown())
    {
        app_.add_back_generic_points();
        update_display();
    }
    if (event.getKeySym() == "k" && event.keyDown())
    {
        // toggle generic points
        settings_.show_generic_points = !settings_.show_generic_points;
        update_display();
    }
    if (event.getKeySym() == "l" && event.keyDown())
    {
        // toggle generic points
        settings_.show_interior_points = !settings_.show_interior_points;
        update_display();
    }
    if (event.getKeySym() == "r" && event.keyDown())
    {
        // restart
        app_.restart();
        update_display();
    }
}