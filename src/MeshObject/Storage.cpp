#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "MeshObject/InteriorPoint.hpp"

#include "MeshObject/RRSTree.hpp"
#include "MeshObject/TriangleBVH.hpp"

Storage::Storage()
{
    is_expired_ = false;
}

Storage::~Storage()
{
    is_expired_ = true;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const Eigen::Vector3d& origin, const Eigen::Vector3d& position) 
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), position, origin);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius)
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), position, origin, radius);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Edge>& Storage::add_edge(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{    
    // create
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), surface, vertex1, vertex2);

    // store
    return *edges_.insert(edge).first;
}

const std::shared_ptr<Face>& Storage::add_face(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3) 
{
    // create
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), vertex1, vertex2, vertex3);

    // store
    return *faces_.insert(face).first;
}

const std::shared_ptr<Surface>& Storage::add_surface() 
{
    // create
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this());

    // store
    return *surfaces_.insert(surface).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), position, origin);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<InteriorPoint>& Storage::add_interior_point(const std::shared_ptr<Face>& face, const Eigen::Vector3d& position, const Eigen::Vector3d& origin) 
{
    // create
    std::shared_ptr<InteriorPoint> interior_point = std::make_shared<InteriorPoint>();
    interior_point->initialize_(shared_from_this(), face, position, origin);

    // store
    return *interior_points_.insert(interior_point).first;
}

// need to ensure the vertex/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vertex(const std::shared_ptr<Vertex>& vertex) 
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to delete expired vertex.");

    // storage delete
    vertices_.erase(vertex);

    // member delete
    vertex->delete_();    
}

void Storage::delete_edge(const std::shared_ptr<Edge>& edge) 
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to delete expired edge.");

    // storage delete
    edges_.erase(edge);
    
    // member delete
    edge->delete_();
}

void Storage::delete_face(const std::shared_ptr<Face>& face) 
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to delete expired face.");

    // storage delete
    faces_.erase(face);

    // member delete
    face->delete_();
}

void Storage::delete_surface(const std::shared_ptr<Surface>& surface) 
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to delete expired surface.");

    // storage delete
    surfaces_.erase(surface);

    // member delete
    surface->delete_();
}

void Storage::delete_genertic_point(const std::shared_ptr<GenericPoint>& genertic_point) 
{
    // check input
    if (genertic_point->is_expired()) throw std::runtime_error("Attempts to delete expired genertic point.");

    // storage delete
    genertic_points_.erase(genertic_point);

    // member delete
    genertic_point->delete_();
}

void Storage::delete_interior_point(const std::shared_ptr<InteriorPoint>& interior_point) 
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to delete expired interior point.");

    // storage delete
    interior_points_.erase(interior_point);
    
    // member delete
    interior_point->delete_();
}

void Storage::add_searchable_vertex(const std::shared_ptr<Vertex>& vertex)
{
    // add to rrs_tree
    rrs_tree_.tree_add_vertex(vertex);
}

void Storage::remove_searchable_vertex(const std::shared_ptr<Vertex>& vertex)
{
    // remove from rrs_tree
    rrs_tree_.tree_delete_vertex(vertex);
}

void Storage::add_searchable_face(const std::shared_ptr<Face>& face)
{
    // add to triangle_bvh
    triangle_bvh_.tree_add_face(face);
}

void Storage::remove_searchable_face(const std::shared_ptr<Face>& face)
{
    // remove from triangle_bvh
    triangle_bvh_.tree_delete_face(face);
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
std::set<std::shared_ptr<Vertex>> Storage::reverse_radius_search(const Eigen::Vector3d& point) 
{
    std::set<std::shared_ptr<Vertex>> result;
    rrs_tree_.tree_reverse_radius_search(point, result);
    return result;
}

// face intersection search
std::set<std::shared_ptr<Face>> Storage::face_intersection_search(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) 
{
    std::set<std::shared_ptr<Face>> result;
    triangle_bvh_.tree_intersection_search(origin, point, result);
    return result;
}

const std::set<std::shared_ptr<Vertex>>& Storage::get_vertices() const
{
    return vertices_;
}

const std::set<std::shared_ptr<Edge>>& Storage::get_edges() const
{
    return edges_;
}

const std::set<std::shared_ptr<Face>>& Storage::get_faces() const
{
    return faces_;
}

const std::set<std::shared_ptr<Surface>>& Storage::get_surfaces() const
{
    return surfaces_;
}

const std::set<std::shared_ptr<GenericPoint>>& Storage::get_generic_points() const
{
    return genertic_points_;
}

std::vector<std::shared_ptr<Vertex>> Storage::get_rrs_vertices()
{
    return rrs_tree_.compute_vertices_list();
}

std::map<std::shared_ptr<Vertex>, int> Storage::get_vertex_to_cloud_indices_map() const
{
    // initialize
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map;

    // fill
    int id = 0;
    for (const auto& vertex : vertices_)
    {
        vertex_to_cloud_indices_map[vertex] = id;
        id++;
    }

    // return
    return vertex_to_cloud_indices_map;
} 

bool Storage::is_expired() const
{
    return is_expired_;
}

// get edge
const std::shared_ptr<Edge>& Storage::get_edge(std::shared_ptr<Vertex> vertex1, std::shared_ptr<Vertex> vertex2) const
{
    // check input
    if (vertex1->is_expired() || vertex2->is_expired()) throw std::runtime_error("Attempts to get edge with invalid vertex.");

    // search
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->has_vertex(vertex1) && edge->has_vertex(vertex2)) return edge;
    }

    // not found
    throw std::runtime_error("Edge not found.");
}

void Storage::print_rrs() const
{
    rrs_tree_.tree_print();
}

void Storage::print_bvh() const
{
    triangle_bvh_.tree_print();
}

void Storage::rebuild_tree()
{
    rrs_tree_.rebuild();
}
