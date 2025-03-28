#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PolygonMesh.h>
#include <pcl/io/ply_io.h>
#include "MeshObject/Surface.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Settings.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>

#include "utilities/simplified_mesh.hpp"

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef K::Point_2 Point_2;
typedef K::Triangle_2 Triangle_2;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> KdTree;

// AABB Tree definitions
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_triangle_primitive_2.h>
typedef std::vector<Triangle_2>::iterator TriangleIterator;
typedef CGAL::AABB_triangle_primitive_2<K, TriangleIterator> Primitive;
typedef CGAL::AABB_traits_2<K, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;

bool is_point_inside_triangle(const Triangle_2& triangle, const KdTree& kd_tree) 
{
    // Get triangle bounding box
    double min_x = std::min({triangle[0].x(), triangle[1].x(), triangle[2].x()});
    double max_x = std::max({triangle[0].x(), triangle[1].x(), triangle[2].x()});
    double min_y = std::min({triangle[0].y(), triangle[1].y(), triangle[2].y()});
    double max_y = std::max({triangle[0].y(), triangle[1].y(), triangle[2].y()});

    // Query points within bounding box using KD-tree
    CGAL::Fuzzy_iso_box<Traits> bbox_query(Point_2(min_x, min_y), Point_2(max_x, max_y));
    std::vector<Point_2> nearby_points;
    kd_tree.search(back_inserter(nearby_points), bbox_query);

    // Check only the nearby points
    unsigned int num_points_inside = 0;
    for (const auto& point : nearby_points) 
    {
        if (triangle.has_on_bounded_side(point)) 
        {
            num_points_inside++;
        }
    }

    // compute point density inside
    double triangle_area = CGAL::to_double(triangle.area());
    double point_density_inside = num_points_inside / triangle_area;
    
    // threashold density 
    Settings settings_;
    double threshold_density = settings_.simplify_surfaces_density_threshold;

    // return true if point density inside is greater than threshold
    if (point_density_inside > threshold_density) 
    {
        return true;
    }
    else
    {
        return false;
    }
}

pcl::PolygonMesh create_simplified_mesh(const std::shared_ptr<Surface>& surface, bool with_color)
{
    if (with_color)
    {
        return create_simplified_mesh_impl<pcl::PointXYZRGB>(surface);
    }
    else
    {
        return create_simplified_mesh_impl<pcl::PointXYZ>(surface);
    }
}

template <typename PointT>
pcl::PolygonMesh create_simplified_mesh_impl(const std::shared_ptr<Surface>& surface)
{
    // skip if surface is seed
    if (surface->is_seed()) return pcl::PolygonMesh();

    // get list of vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices = surface->get_vertices();

    Settings settings_;

    // filter the vertices
    std::vector<std::shared_ptr<Vertex>> vertices_filtered;
    {
        // cap radius of boundary vertices
        std::unordered_map<std::shared_ptr<Vertex>, double, MeshObjectHash> vertex_radii;
        for (const std::shared_ptr<Vertex>& vertex : vertices)
        {
            const double radius_original = vertex->get_radius() - settings_.extra_radius;

            vertex_radii[vertex] = radius_original;
        }

        // sort the vertices by radius
        std::vector<std::shared_ptr<Vertex>> vertices_sorted(vertices.begin(), vertices.end());
        std::sort(vertices_sorted.begin(), vertices_sorted.end(), [&vertex_radii](const std::shared_ptr<Vertex>& a, const std::shared_ptr<Vertex>& b) 
        {
            return vertex_radii[a] < vertex_radii[b];
        });
        
        // filter the vertices
        for (const std::shared_ptr<Vertex>& vertex : vertices_sorted)
        {
            // skip if new vertex is too close to existing vertices
            bool too_close = false;
            for (const std::shared_ptr<Vertex>& vertex_filtered : vertices_filtered)
            {
                // radius original
                const double radius_original = vertex_radii[vertex_filtered]; // this checks if existing vertex is within current vertex's radius
                
                // edge length
                const double edge_length = (vertex->get_position() - vertex_filtered->get_position()).norm();

                // skip if too close
                if (edge_length < radius_original) 
                {
                    too_close = true;
                    break;
                }
            }

            // skip if too close
            if (too_close) continue;

            // store
            vertices_filtered.push_back(vertex);
        }
    }

    // cloud 3d 
    typename pcl::PointCloud<PointT>::Ptr cloud_3d(new pcl::PointCloud<PointT>);
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_3d_indices_map;
    {
        unsigned int index = 0;

        // for each vertex
        for (const std::shared_ptr<Vertex>& vertex : vertices_filtered)
        {
            // get 3d point
            PointT point;
            Eigen::Vector3d projected_position = vertex->compute_projected_position();
            point.x = projected_position[0];
            point.y = projected_position[1];
            point.z = projected_position[2];
            // Eigen::Vector2d surface_coordinate = vertex->get_surface_coordinate();
            // point.x = surface_coordinate[0];
            // point.y = surface_coordinate[1];
            // point.z = 0;

            // If the point type supports color, assign it
            if constexpr (std::is_same<PointT, pcl::PointXYZRGB>::value)
            {
                const std::tuple<int, int, int>& color = vertex->get_surface()->get_color();
                point.r = std::get<0>(color);
                point.g = std::get<1>(color);
                point.b = std::get<2>(color);
            }
            
            // store point
            cloud_3d->push_back(point);

            // store mapping
            vertex_to_cloud_3d_indices_map[vertex] = index++;
        }
    }    
    
    // delaunay triangulation
    Delaunay dt;

    // add point to delaunay triangulation
    std::map<Delaunay::Vertex_handle, std::shared_ptr<Vertex>> vertex_handle_to_vertex_map;
    std::map<std::shared_ptr<Vertex>, Delaunay::Vertex_handle> vertex_to_vertex_handle_map;
    for (const std::shared_ptr<Vertex>& vertex : vertices_filtered)
    {
        // get 2d point
        Eigen::Vector2d surface_coordinate = vertex->get_surface_coordinate();
        Point_2 p2d(surface_coordinate[0], surface_coordinate[1]);

        // add point to delaunay triangulation
        Delaunay::Vertex_handle vertex_handle = dt.insert(p2d);

        // store mapping
        vertex_handle_to_vertex_map[vertex_handle] = vertex;
        vertex_to_vertex_handle_map[vertex] = vertex_handle;
    }

    // extract faces from triangulation
    std::vector<pcl::Vertices> mesh_faces;
    for (auto delaunay_face_iterator = dt.finite_faces_begin(); delaunay_face_iterator != dt.finite_faces_end(); ++delaunay_face_iterator) 
    {
        Triangle_2 triangle
        (
            delaunay_face_iterator->vertex(0)->point(),
            delaunay_face_iterator->vertex(1)->point(),
            delaunay_face_iterator->vertex(2)->point()
        );
        
        // skip if center of delaunay face triangle is not inside original triangulation
        Point_2 center = CGAL::centroid(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));

        // check if center is within any existing circles from points
        bool within_circle = false;
        for (const std::shared_ptr<Vertex>& vertex : vertices_filtered)
        {
            // get 2d point
            Eigen::Vector2d surface_coordinate = vertex->get_surface_coordinate();
            Point_2 p2d(surface_coordinate[0], surface_coordinate[1]);

            // check if center is within circle
            double radius = vertex->get_radius();
            if (CGAL::squared_distance(center, p2d) < radius * radius) 
            {
                within_circle = true;
            }
        }

        if (!within_circle) continue;

        // vertex handle
        Delaunay::Vertex_handle vertex_handle_0 = delaunay_face_iterator->vertex(0);
        Delaunay::Vertex_handle vertex_handle_1 = delaunay_face_iterator->vertex(1);
        Delaunay::Vertex_handle vertex_handle_2 = delaunay_face_iterator->vertex(2);

        // vertices
        std::shared_ptr<Vertex> vertex_0 = vertex_handle_to_vertex_map[vertex_handle_0];
        std::shared_ptr<Vertex> vertex_1 = vertex_handle_to_vertex_map[vertex_handle_1];
        std::shared_ptr<Vertex> vertex_2 = vertex_handle_to_vertex_map[vertex_handle_2];

        // skip if any vertices is nullptr (created due to insertion of constraint)
        if (vertex_0 == nullptr || vertex_1 == nullptr || vertex_2 == nullptr) continue;

        // cloud 3d indices
        unsigned int index_0 = vertex_to_cloud_3d_indices_map[vertex_0];
        unsigned int index_1 = vertex_to_cloud_3d_indices_map[vertex_1];
        unsigned int index_2 = vertex_to_cloud_3d_indices_map[vertex_2];

        // store to face
        pcl::Vertices face;
        face.vertices.push_back(index_0);
        face.vertices.push_back(index_1);
        face.vertices.push_back(index_2);

        // store to mesh
        mesh_faces.push_back(face);
    }

    // Convert point cloud to PCL PolygonMesh format
    pcl::PolygonMesh triangle_mesh;
    pcl::toPCLPointCloud2(*cloud_3d, triangle_mesh.cloud);
    triangle_mesh.polygons = mesh_faces;

    // return
    return triangle_mesh;    
}

void merge_polygon_mesh(pcl::PolygonMesh& mesh1, pcl::PolygonMesh& mesh2)
{
    // modify the cloud indices
    for (pcl::Vertices& face : mesh2.polygons)
    {
        for (int& index : face.vertices)
        {
            index += mesh1.cloud.width;
        }
    }

    // add to all meshes
    mesh1.polygons.insert(mesh1.polygons.end(), mesh2.polygons.begin(), mesh2.polygons.end());
    mesh1.cloud += mesh2.cloud;
}