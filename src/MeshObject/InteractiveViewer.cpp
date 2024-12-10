#include "MeshObject/InteractiveViewer.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include <unordered_set>
#include <pcl/io/ply_io.h>

#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"

template class InteractiveViewer<VilensPointT>;
template class InteractiveViewer<BagPointT>;

template <typename PointT>
InteractiveViewer<PointT>::InteractiveViewer(Application<PointT>& app) 
    : 
    app_(app),
    viewer_(new pcl::visualization::PCLVisualizer ("3D Viewer"))
{   
    viewer_->getRenderWindow()->GlobalWarningDisplayOff(); // Add This Line
    if (settings_.flip_color) 
    {
        viewer_->setBackgroundColor (1, 1, 1);
    }
    else
    {
        viewer_->setBackgroundColor (0, 0, 0);
    }
    viewer_->initCameraParameters();
    viewer_->addCoordinateSystem(1);
    viewer_->registerKeyboardCallback(&InteractiveViewer::keyboard_callback, *this, nullptr);
    viewer_->spin();
}

template <typename PointT>
void InteractiveViewer<PointT>::update_display(bool export_ply)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr vertex_pointcloud = app_.compute_vertex_point_pointcloud(settings_);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr interior_point_cloud = app_.compute_interior_point_pointcloud(settings_);
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map = app_.get_vertex_to_cloud_indices_map();
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = app_.get_faces();
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> boundary_edges = app_.get_boundary_edges();

    // vertex points
    viewer_->removeShape("point_cloud");
    if (settings_.show_pointcloud)
    {
        viewer_->addPointCloud<pcl::PointXYZRGB>(vertex_pointcloud, "point_cloud");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "point_cloud");

        if (export_ply)
        {
            pcl::io::savePLYFile("vertex.ply", *vertex_pointcloud);
            std::cout << "exported vertex.ply" << std::endl;
        }
    }

    // interior points
    viewer_->removeShape("interior_points");
    if (settings_.show_interior_points)
    {
        viewer_->addPointCloud<pcl::PointXYZRGB>(interior_point_cloud, "interior_points");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "interior_points");

        if (export_ply) 
        {
            pcl::io::savePLYFile("interior.ply", *interior_point_cloud);
            std::cout << "exported interior.ply" << std::endl;
        }
    }

    // generic points
    viewer_->removeShape("generic_points");
    if (settings_.show_generic_points)
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr generic_point_cloud = app_.compute_generic_point_pointcloud();
        viewer_->addPointCloud<pcl::PointXYZRGB>(generic_point_cloud, "generic_points");
        viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "generic_points");

        if (export_ply) 
        {
            pcl::io::savePLYFile("generic.ply", *generic_point_cloud);
            std::cout << "exported generic.ply" << std::endl;
        }
    }

    // triangle
    viewer_->removeShape("triangle_mesh");
    if (settings_.show_triangle)
    {
        pcl::PolygonMesh triangle_mesh;
        pcl::toPCLPointCloud2(*vertex_pointcloud, triangle_mesh.cloud);
        for (const std::shared_ptr<Face>& face : faces)
        {
            // skip if not confirmed
            if (settings_.show_confirmed_only && !face->is_confirmed()) continue;

            pcl::Vertices triangle;
            triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face->get_vertex(0)));
            triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face->get_vertex(1)));
            triangle.vertices.push_back(vertex_to_cloud_indices_map.at(face->get_vertex(2)));
            triangle_mesh.polygons.push_back(triangle);
        }
        viewer_->addPolygonMesh(triangle_mesh, "triangle_mesh");

        if (export_ply) 
        {
            pcl::io::savePLYFile("triangle.ply", triangle_mesh);
            std::cout << "exported triangle.ply" << std::endl;
        }
    }

    // boundary edge
    viewer_->removeShape("boundary_edges");        
    if (settings_.show_edge)
    {
        pcl::PolygonMesh boundary_mesh;
        pcl::toPCLPointCloud2(*vertex_pointcloud, boundary_mesh.cloud);
        for (const std::shared_ptr<Edge>& edge : boundary_edges)
        {
            // skip if not confirmed
            if (settings_.show_confirmed_only && !edge->is_confirmed()) continue;

            // skip if singular
            if (!settings_.show_singular_edge && edge->is_singular()) continue;
            
            pcl::Vertices boundary_edge;
            boundary_edge.vertices.push_back(vertex_to_cloud_indices_map.at(edge->get_vertex(0)));
            boundary_edge.vertices.push_back(vertex_to_cloud_indices_map.at(edge->get_vertex(1)));
            boundary_mesh.polygons.push_back(boundary_edge);
        }
        viewer_->addPolylineFromPolygonMesh(boundary_mesh, "boundary_edges");
        if (settings_.flip_color) 
        {
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 0, 0, "boundary_edges");
        }
        else
        {
            viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 1, 1, "boundary_edges");
        }
        viewer_->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "boundary_edges");

        if (export_ply) 
        {
            pcl::io::savePLYFile("boundary.ply", boundary_mesh);
            std::cout << "exported boundary.ply" << std::endl;
        }
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
    if (settings_.show_keycode)
    {
        std::cout << "keySym: " << event.getKeySym() << " | keyCode: " << event.getKeyCode() << std::endl;
    }
    if (event.getKeySym() == "Num_Lock" && event.keyDown())
    {
        settings_.show_keycode = !settings_.show_keycode;
        
        // log
        std::cout << "show_keycode: " << settings_.show_keycode << std::endl;
    }
    if (event.getKeySym() == "Return" && event.keyDown())
    {
        // update display
        update_display();
        
        // log
        std::cout << "update display" << std::endl;
    }
    if (event.getKeySym() == "space" && event.keyDown())
    {
        app_.loop();
        if (settings_.update_display) update_display();

        // log
        std::cout << "loop" << std::endl;
    }
    // kp number 0
    if (event.getKeySym() == "KP_Insert" && event.keyDown())
    {
        settings_.color_mode = ColorMode::ID;
        update_display();

        // log
        std::cout << "color_mode: original" << std::endl;
    }
    // kp number 1
    if (event.getKeySym() == "KP_End" && event.keyDown())
    {
        settings_.color_mode = ColorMode::POSITIONAL_UNCERTAINTY;
        update_display();

        // log
        std::cout << "color_mode: projected distance" << std::endl;
    }
    // kp number 2
    if (event.getKeySym() == "KP_Down" && event.keyDown())
    {
        settings_.color_mode = ColorMode::POSITIONAL_UNCERTAINTY_NORMALIZED;
        update_display();

        // log
        std::cout << "color_mode: positional uncertainty normalized" << std::endl;
    }
    // kp number 3
    if (event.getKeySym() == "KP_Next" && event.keyDown())
    {
        settings_.color_mode = ColorMode::RADIUS;
        update_display();

        // log
        std::cout << "color_mode: radius" << std::endl;
    }
    // kp number 4
    if (event.getKeySym() == "KP_Left" && event.keyDown())
    {
        settings_.color_mode = ColorMode::SURFACE_UNCERTAINTY;
        update_display();

        // log
        std::cout << "color_mode: plane positional uncertainty" << std::endl;
    }
    // kp number 5
    if (event.getKeySym() == "KP_Begin" && event.keyDown())
    {
        settings_.color_mode = ColorMode::CONTENTION;
        update_display();

        // log
        std::cout << "color_mode: contention" << std::endl;
    }
    // kp number 6
    if (event.getKeySym() == "KP_Right" && event.keyDown())
    {
        settings_.color_mode = ColorMode::DISTANCE_TRAVELLED;
        update_display();

        // log
        std::cout << "color_mode: distance travelled" << std::endl;
    }
    // kp number 7
    if (event.getKeySym() == "KP_Home" && event.keyDown())
    { 
        settings_.color_mode = ColorMode::PROJECTED_UNCERTAINTY;
        update_display();

        // log
        std::cout << "color_mode: projected uncertainty" << std::endl;
    }
    // kp numebr 8
    if (event.getKeySym() == "KP_Up" && event.keyDown())
    {
    }
    // kp numebr 9
    if (event.getKeySym() == "KP_Prior" && event.keyDown())
    {
        // singular edge
        settings_.show_singular_edge = !settings_.show_singular_edge;
        settings_.show_singular_vertex = !settings_.show_singular_vertex;
        update_display();

        // log
        std::cout << "show_singular_edge: " << settings_.show_singular_edge << std::endl;
        std::cout << "show_singular_vertex: " << settings_.show_singular_vertex << std::endl;
    }
    
    if (event.getKeySym() == "1" && event.keyDown())
    {
        app_.step();
        update_display();

        // log
        std::cout << "step" << std::endl;
    }
    if (event.getKeySym() == "2" && event.keyDown())
    {
        for (int i = 0; i < 10; i++) app_.step();
        update_display();

        // log
        std::cout << "step 10" << std::endl;
    }
    if (event.getKeySym() == "3" && event.keyDown())
    {
        for (int i = 0; i < 100; i++) app_.step();
        update_display();

        // log
        std::cout << "step 100" << std::endl;
    }
    if (event.getKeySym() == "4" && event.keyDown())
    {
        for (int i = 0; i < 1000; i++) app_.step();
        update_display();

        // log
        std::cout << "step 1000" << std::endl;
    }
    if (event.getKeySym() == "9" && event.keyDown())
    {
        update_display(true);

        // log
        std::cout << "exported complete" << std::endl;
    }
    if (event.getKeySym() == "0" && event.keyDown())
    {
        app_.process_the_rest();
        update_display();

        // log
        std::cout << "process the rest" << std::endl;
    }
    if (event.getKeySym() == "Tab" && event.keyDown())
    {
        app_.change_color();
        update_display();

        // log
        std::cout << "changed color" << std::endl;
    }
    if (event.getKeySym() == "comma" && event.keyDown())
    {
        settings_.show_pointcloud = !settings_.show_pointcloud;
        update_display();

        // log
        std::cout << "show_verticies: " << settings_.show_pointcloud << std::endl;
    }
    if (event.getKeySym() == "period" && event.keyDown())
    {
        settings_.show_edge = !settings_.show_edge;
        update_display();

        // log
        std::cout << "show_edges: " << settings_.show_edge << std::endl;
    }
    if (event.getKeySym() == "slash" && event.keyDown())
    {
        settings_.show_triangle = !settings_.show_triangle;
        update_display();

        // log
        std::cout << "show_faces: " << settings_.show_triangle << std::endl;
    }
    if (event.getKeySym() == "a" && event.keyDown())
    {
        settings_.show_projected_point = !settings_.show_projected_point;
        update_display();

        // log
        std::cout << "show_projected_point: " << settings_.show_projected_point << std::endl;
    }
    if (event.getKeySym() == "z" && event.keyDown())
    {
        // show confirmed only
        settings_.show_confirmed_only = !settings_.show_confirmed_only;
        update_display();

        // log
        std::cout << "show_confirmed_only: " << settings_.show_confirmed_only << std::endl;
    }
    if (event.getKeySym() == "v" && event.keyDown())
    {
        settings_.show_wireframe = !settings_.show_wireframe;
        update_display();

        // log
        std::cout << "show_wireframe: " << settings_.show_wireframe << std::endl;
    }
    if (event.getKeySym() == "m" && event.keyDown())
    {
        settings_.show_sphere = !settings_.show_sphere;
        update_display();

        // log
        std::cout << "show_sphere: " << settings_.show_sphere << std::endl;
    }
    if (event.getKeySym() == "b" && event.keyDown())
    {
        app_.refine_surfaces();
        update_display();

        // log
        std::cout << "refined surfaces" << std::endl;
    }
    if (event.getKeySym() == "k" && event.keyDown())
    {
        // toggle generic points
        settings_.show_generic_points = !settings_.show_generic_points;
        update_display();

        // log
        std::cout << "show_generic_points: " << settings_.show_generic_points << std::endl;
    }
    if (event.getKeySym() == "l" && event.keyDown())
    {
        // toggle generic points
        settings_.show_interior_points = !settings_.show_interior_points;
        update_display();

        // log
        std::cout << "show_interior_points: " << settings_.show_interior_points << std::endl;
    }
    if (event.getKeySym() == "r" && event.keyDown())
    {
        // restart
        app_.restart();
        update_display();

        // log
        std::cout << "restarted" << std::endl;
    }
    if (event.getKeySym() == "t" && event.keyDown())
    {
        // restart
        app_.rebuild_tree();

        // log
        std::cout << "rebuild tree" << std::endl;
    }
}