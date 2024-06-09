#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "MeshObject/InteriorPoint.hpp"

#include "eye_patch/RRSTree.hpp"
#include "eye_patch/TriangleBVH.hpp"

std::weak_ptr<Vertex> Storage::add_vertex(Eigen::Vector3d origin, Eigen::Vector3d position) 
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), position, origin);

    // store
    vertices_.insert(vertex);

    // return
    return vertex;
}

std::weak_ptr<Vertex> Storage::add_vertex(Eigen::Vector3d origin, Eigen::Vector3d position, double radius)
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), position, origin, radius);

    // store
    vertices_.insert(vertex);

    // return
    return vertex;
}

std::weak_ptr<Edge> Storage::add_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2)
{    
    // create
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), vertex1, vertex2);

    // store
    edges_.insert(edge);

    // return
    return edge;

}

std::weak_ptr<Face> Storage::add_face(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2, std::weak_ptr<Vertex> vertex3) 
{
    // create
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), vertex1, vertex2, vertex3);

    // store
    faces_.insert(face);

    // return
    return face;
}

std::weak_ptr<Surface> Storage::add_surface() 
{
    // create
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this());

    // store
    surfaces_.insert(surface);

    // return
    return surface;
}

std::weak_ptr<Surface> Storage::add_surface(std::weak_ptr<Surface> surface1, std::weak_ptr<Surface> surface2) 
{
    // create
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this(), surface1, surface2);

    // store
    surfaces_.insert(surface);

    // return
    return surface;
}

std::weak_ptr<GenericPoint> Storage::add_generic_point(Eigen::Vector3d position, Eigen::Vector3d origin) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), position, origin);

    // store
    genertic_points_.insert(genertic_point);

    // return
    return genertic_point;
}

std::weak_ptr<InteriorPoint> Storage::add_interior_point(std::weak_ptr<Face> face, Eigen::Vector3d position, Eigen::Vector3d origin) 
{
    // create
    std::shared_ptr<InteriorPoint> interior_point = std::make_shared<InteriorPoint>();
    interior_point->initialize_(shared_from_this(), face, position, origin);

    // store
    interior_points_.insert(interior_point);

    // return
    return interior_point;
}

// need to ensure the vertex/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vertex(std::weak_ptr<Vertex> vertex) 
{
    // check input
    if (vertex.expired()) throw std::runtime_error("Attempts to delete expired vertex.");

    // member delete
    vertex.lock()->delete_();

    // storage delete
    vertices_.erase(vertex.lock());
}

void Storage::delete_edge(std::weak_ptr<Edge> edge) 
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to delete expired edge.");

    // member delete
    edge.lock()->delete_();

    // storage delete
    edges_.erase(edge.lock());
}

void Storage::delete_face(std::weak_ptr<Face> face) 
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to delete expired face.");

    // member delete
    face.lock()->delete_();

    // storage delete
    faces_.erase(face.lock());
}

void Storage::delete_surface(std::weak_ptr<Surface> surface) 
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to delete expired surface.");

    // member delete
    surface.lock()->delete_();

    // storage delete
    surfaces_.erase(surface.lock());
}

void Storage::delete_genertic_point(std::weak_ptr<GenericPoint> genertic_point) 
{
    // check input
    if (genertic_point.expired()) throw std::runtime_error("Attempts to delete expired genertic point.");

    // member delete
    genertic_point.lock()->delete_();

    // storage delete
    genertic_points_.erase(genertic_point.lock());
}

void Storage::delete_interior_point(std::weak_ptr<InteriorPoint> interior_point) 
{
    // check input
    if (interior_point.expired()) throw std::runtime_error("Attempts to delete expired interior point.");

    // member delete
    interior_point.lock()->delete_();

    // storage delete
    interior_points_.erase(interior_point.lock());
}

void Storage::add_searchable_vertex(std::weak_ptr<Vertex> vertex)
{
    // add to rrs_tree
    rrs_tree_.add_vertex(vertex);
}

void Storage::remove_searchable_vertex(std::weak_ptr<Vertex> vertex)
{
    // remove from rrs_tree
    rrs_tree_.delete_vertex(vertex);
}

void Storage::add_searchable_face(std::weak_ptr<Face> face)
{
    // add to triangle_bvh
    triangle_bvh_.add_face(face);
}

void Storage::remove_searchable_face(std::weak_ptr<Face> face)
{
    // remove from triangle_bvh
    triangle_bvh_.delete_face(face);
}

// get id
int Storage::get_next_vertex_id() { return next_vertex_id_++; }
int Storage::get_next_edge_id() { return next_edge_id_++; }
int Storage::get_next_face_id() { return next_face_id_++; }
int Storage::get_next_surface_id() { return next_surface_id_++; }
int Storage::get_next_generic_point_id() { return next_genertic_point_id_++; }
int Storage::get_next_interior_point_id() { return next_interior_point_id_++; }

bool Storage::can_reverse_radius_search() 
{ 
    return rrs_tree_.can_reverse_radius_search(); 
}

// reverse radius search
std::set<std::weak_ptr<Vertex>> Storage::reverse_radius_search(Eigen::Vector3d point) 
{
    return rrs_tree_.reverse_radius_search(point);
}

// face intersection search
std::set<std::weak_ptr<Face>> Storage::face_intersection_search(Eigen::Vector3d origin, Eigen::Vector3d point) 
{
    return triangle_bvh_.intersectionSearch(origin, point);
}

std::set<std::weak_ptr<Vertex>> Storage::get_vertices() const
{
    std::set<std::weak_ptr<Vertex>> vertices;
    for (auto vertex : vertices_) vertices.insert(vertex);
    return vertices;
}

std::set<std::weak_ptr<Edge>> Storage::get_edges() const
{
    std::set<std::weak_ptr<Edge>> edges;
    for (auto edge : edges_) edges.insert(edge);
    return edges;
}

std::set<std::weak_ptr<Face>> Storage::get_faces() const
{
    std::set<std::weak_ptr<Face>> faces;
    for (auto face : faces_) faces.insert(face);
    return faces;
}

std::map<std::weak_ptr<Vertex>, int> Storage::get_vertex_to_cloud_indices_map() const
{
    // initialize
    std::map<std::weak_ptr<Vertex>, int> vertex_to_cloud_indices_map;

    // fill
    int id = 0;
    for (auto vertex : vertices_)
    {
        vertex_to_cloud_indices_map[vertex] = id;
        id++;
    }

    // return
    return vertex_to_cloud_indices_map;
} 


// get edge
std::weak_ptr<Edge> Storage::get_edge(std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2) const
{
    // check input
    if (vertex1.expired() || vertex2.expired()) throw std::runtime_error("Attempts to get edge with invalid vertex.");

    // search
    for (auto edge : edges_)
    {
        if (edge->has_vertex(vertex1) && edge->has_vertex(vertex2)) return edge;
    }

    // error
    throw std::runtime_error("Edge not found.");
}

void Storage::print_rrs() const
{
    rrs_tree_.print();
}

void Storage::print_bvh() const
{
    triangle_bvh_.print();
}

void Storage::rebuild_tree()
{
    rrs_tree_.rebuild();
}
