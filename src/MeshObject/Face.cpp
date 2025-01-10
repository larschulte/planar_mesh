#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/GenericPoint.hpp"

Settings Face::settings_;

void Face::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface> surface, const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{
    // set expired
    is_expired_ = false;

    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create face with invalid storage.");
    if (vertex0->is_expired()) throw std::runtime_error("Attempts to create face with invalid vertex0.");
    if (vertex1->is_expired()) throw std::runtime_error("Attempts to create face with invalid vertex1.");
    if (vertex2->is_expired()) throw std::runtime_error("Attempts to create face with invalid vertex2.");

    // store
    storage_ = storage;

    // get id
    id_ = storage_->get_next_face_id();

    // connect surface
    connect(surface);

    // connect
    connect(vertex0);
    connect(vertex1);
    connect(vertex2);

    // store vertices position
    v0_ = vertex0->get_position();
    v1_ = vertex1->get_position();
    v2_ = vertex2->get_position();

    // compute min and max
    BoundingBox box;
    box.expand(vertex0->get_position());
    box.expand(vertex1->get_position());
    box.expand(vertex2->get_position());
    min_ = box.min;
    max_ = box.max;

    // first vertex
    first_vertex_ = vertex0;

    // get edges
    std::shared_ptr<Edge> edge0;
    std::shared_ptr<Edge> edge1;
    std::shared_ptr<Edge> edge2;
    for (const std::shared_ptr<Edge>& edge : vertex0->get_edges())
    {
        if (edge->get_surface() != surface_) continue;
        if (edge->has_vertex(vertex1))
        {
            edge0 = edge;
            break;
        }
    }
    for (const std::shared_ptr<Edge>& edge : vertex1->get_edges())
    {
        if (edge->get_surface() != surface_) continue;
        if (edge->has_vertex(vertex2))
        {
            edge1 = edge;
            break;
        }
    }
    for (const std::shared_ptr<Edge>& edge : vertex2->get_edges())
    {
        if (edge->get_surface() != surface_) continue;
        if (edge->has_vertex(vertex0))
        {
            edge2 = edge;
            break;
        }
    }
    connect(edge0);
    connect(edge1);
    connect(edge2);

    // compute center
    const Eigen::Vector3d& pos0 = vertex0->get_position();
    const Eigen::Vector3d& pos1 = vertex1->get_position();
    const Eigen::Vector3d& pos2 = vertex2->get_position();
    center_ = (pos0 + pos1 + pos2) / 3;

    // add to search
    // storage_->add_searchable_face(shared_from_this());
    storage_->add_to_set_of_faces_to_update_bvh_tree(shared_from_this());

    // log
    if (settings_.log.initialize) std::cout << "Face " << id_ << " created between vertex " << vertex0->get_id() << ", vertex " << vertex1->get_id() << " and vertex " << vertex2->get_id() << std::endl;
}

void Face::update_radius(const std::shared_ptr<GenericPoint>& generic_point)
{
}

void Face::delete_()
{
    // lock face
    while (!omp_test_nest_lock(&face_lock)) 
    {
        std::cout << "waiting to lock face " << id_ << std::endl;
    }

    // add to affected faces set
    storage_->add_to_set_of_faces_to_update_bvh_tree(shared_from_this());

    // log
    if (settings_.log.deletion) std::cout << "Destroying face " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    // make a copy of the set to avoid iterator invalidation
    std::vector<std::shared_ptr<Vertex>> vertices = vertices_;
    std::vector<std::shared_ptr<Edge>> edges = edges_;
    std::vector<std::shared_ptr<InteriorPoint>> interior_points = interior_points_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> sibling_faces = sibling_faces_;
    for (const auto& vertex : vertices) disconnect(vertex);
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& interior_point : interior_points) disconnect(interior_point);
    if (surface_ != nullptr) disconnect(surface_);
    for (const auto& sibling_face : sibling_faces) disconnect(sibling_face);

    // log
    if (settings_.log.deletion) std::cout << "---------- face " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;

    // release face lock
    omp_unset_nest_lock(&face_lock);
}

void Face::temp_initialize(const Eigen::Vector3d& end_point)
{
    // temp initialize a vertex
    std::shared_ptr<Vertex> temp0 = std::make_shared<Vertex>();
    std::shared_ptr<Vertex> temp1 = std::make_shared<Vertex>();
    std::shared_ptr<Vertex> temp2 = std::make_shared<Vertex>();
    temp0->temp_initialize(end_point, 0);
    temp1->temp_initialize(end_point, 1);
    temp2->temp_initialize(end_point, 2);

    // set the vertex as the first vertex
    first_vertex_ = temp0;
    vertices_.push_back(temp0);
    vertices_.push_back(temp1);
    vertices_.push_back(temp2);

    // set the bounding box
    double radius = 0.001;
    min_ = end_point - Eigen::Vector3d(radius, radius, radius);
    max_ = end_point + Eigen::Vector3d(radius, radius, radius);
}

void Face::un_add_face()
{
    // get copy of vertices
    std::shared_ptr<Vertex> vertex0 = *vertices_.begin();
    std::shared_ptr<Vertex> vertex1 = *std::next(vertices_.begin(), 1);
    std::shared_ptr<Vertex> vertex2 = *std::next(vertices_.begin(), 2);

    // set can self destruct
    vertex0->set_can_self_destruct(false);
    vertex1->set_can_self_destruct(false);
    vertex2->set_can_self_destruct(false);
    vertex0->get_edge(vertex1)->set_can_self_destruct(false);
    vertex1->get_edge(vertex2)->set_can_self_destruct(false);
    vertex2->get_edge(vertex0)->set_can_self_destruct(false);

    // delete face
    storage_->delete_face(shared_from_this());

    // set can self destruct
    vertex0->set_can_self_destruct(true);
    vertex1->set_can_self_destruct(true);
    vertex2->set_can_self_destruct(true);
    vertex0->get_edge(vertex1)->set_can_self_destruct(true);
    vertex1->get_edge(vertex2)->set_can_self_destruct(true);
    vertex2->get_edge(vertex0)->set_can_self_destruct(true);
}

const int& Face::get_id() const
{
    return id_;
}

const Eigen::Vector3d& Face::get_center() const
{
    return center_;
}

const std::vector<std::shared_ptr<Vertex>>& Face::get_vertices() const
{
    return vertices_;
}

const std::vector<std::shared_ptr<InteriorPoint>>& Face::get_interior_points() const
{
    return interior_points_;
}

const std::vector<std::shared_ptr<Edge>>& Face::get_edges() const
{
    return edges_;
}

const std::shared_ptr<Vertex>& Face::get_vertex(int index) const
{
    if (index < 0 || index > 2) throw std::runtime_error("Invalid index for vertex.");
    auto it = vertices_.begin();
    for (int i = 0; i < index; i++) it++;
    return *it;
}

const std::shared_ptr<Vertex>& Face::get_first_vertex() const
{
    return first_vertex_;
}

const std::shared_ptr<Surface>& Face::get_surface() const
{
    return surface_;
}

const Eigen::Vector3d& Face::get_min() const
{
    return min_;
}

const Eigen::Vector3d& Face::get_max() const
{
    return max_;
}

const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Face::get_sibling_faces() const
{
    return sibling_faces_;
}

bool Face::is_expired() const
{
    return is_expired_;
}

bool Face::is_deleting() const
{
    return deleting_;
}

bool Face::is_searchable() const
{
    return node != nullptr;
}

bool Face::has_vertex(const std::shared_ptr<Vertex>& vertex) const
{
    for (const std::shared_ptr<Vertex>& vertex_ : vertices_) if (vertex_ == vertex) return true;
    return false;
}

bool Face::intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction) const
{    
    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1_ - v0_;
    Eigen::Vector3d edge2 = v2_ - v0_;
    
    Eigen::Vector3d pvec = direction.cross(edge2);
    double det = edge1.dot(pvec);
    if (std::fabs(det) < EPSILON) return false;

    double invDet = 1.0 / det;

    Eigen::Vector3d tvec = origin - v0_;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) return false;
    
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = direction.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return false;

    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) return false;

    return true;
}   

bool Face::intersects_point(const std::shared_ptr<GenericPoint>& generic_point)
{
    return intersects_point(generic_point->get_origin(), generic_point->get_direction());
}

Eigen::Vector3d Face::compute_intersection_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction)
{
    // get first Vertex in vertices_
    auto it = vertices_.begin();
    std::shared_ptr vertex0 = *(it++);
    std::shared_ptr vertex1 = *(it++);
    std::shared_ptr vertex2 = *(it++);
    const Eigen::Vector3d& v0 = vertex0->get_position();
    const Eigen::Vector3d& v1 = vertex1->get_position();
    const Eigen::Vector3d& v2 = vertex2->get_position();

    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    Eigen::Vector3d pvec = direction.cross(edge2);
    double det = edge1.dot(pvec);
    if (std::fabs(det) < EPSILON) throw std::runtime_error("No intersection point found.");

    double invDet = 1.0 / det;

    Eigen::Vector3d tvec = origin - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) throw std::runtime_error("No intersection point found.");
    
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = direction.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) throw std::runtime_error("No intersection point found.");

    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) throw std::runtime_error("No intersection point found.");

    return origin + direction * t;
}


void Face::connect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to connect face with invalid vertex.");

    // skip if already connected
    for (const std::shared_ptr<Vertex>& vertex_ : vertices_) if (vertex_ == vertex) return;

    // connect
    vertices_.push_back(vertex);
    vertex->connect(shared_from_this());

    // check size
    if (vertices_.size() > 3) throw std::runtime_error("Face connected to more than 3 vertices.");
}

void Face::connect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to connect face with invalid edge.");

    // skip if already connected
    for (const std::shared_ptr<Edge>& edge_ : edges_) if (edge_ == edge) return;

    // connect
    edges_.push_back(edge);
    edge->connect(shared_from_this());

    // check size
    if (edges_.size() > 3) throw std::runtime_error("Face connected to more than 3 edges.");
}

void Face::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect face with invalid surface.");

    // connect
    bool inserted = surface_ != surface;
    if (inserted) surface_ = surface;
    if (inserted) surface->connect(shared_from_this());
}

void Face::connect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to connect face with invalid interior point.");

    // skip if already connected
    for (const std::shared_ptr<InteriorPoint>& interior_point_ : interior_points_) if (interior_point_ == interior_point) return;

    // connect
    interior_points_.push_back(interior_point);
    interior_point->connect(shared_from_this());
    update_confirmed_status();
}

void Face::connect(const std::shared_ptr<Face>& sibling_face)
{
    // check input
    if (sibling_face->is_expired()) throw std::runtime_error("Attempts to connect face with invalid sibling face.");

    // skip if try to connect to itself
    if (sibling_face == shared_from_this()) return;

    // connect
    bool inserted = sibling_faces_.insert(sibling_face).second;
    if (inserted) std::cout << "Connected face " << id_ << " with face " << sibling_face->get_id() << " as sibling."<< std::endl;
    if (inserted) sibling_face->connect(shared_from_this());
    if (inserted)
    {
        for (const std::shared_ptr<Face>& sibling_face_ : sibling_faces_)
        {
            sibling_face_->connect(sibling_face);
        }
    }
}

void Face::disconnect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) return;

    // skip if not connected
    auto it = std::find(vertices_.begin(), vertices_.end(), vertex);
    if (it == vertices_.end()) return;

    // delete
    vertices_.erase(it);
    vertex->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && can_self_destruct_) storage_->delete_face(shared_from_this());
}

void Face::disconnect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) return;

    // skip if not connected
    auto it = std::find(edges_.begin(), edges_.end(), edge);
    if (it == edges_.end()) return;

    // delete
    edges_.erase(it);
    edge->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_->delete_face(shared_from_this());
}

void Face::disconnect(const std::shared_ptr<Surface>& surface)
{
    // lock node
    std::shared_ptr<Node> node_copy = node ? node : std::make_shared<Node>(); // lock if node exists
    while (!omp_test_nest_lock(&node_copy->omp_lock)) 
    {
        std::cout << "disconnect face waiting " << id_ << std::endl;
    }

    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surface_ == surface;
    if (erased) surface->disconnect(shared_from_this());
    if (erased) surface_ = nullptr;

    // self destruct
    if (!deleting_ && erased && can_self_destruct_) storage_->delete_face(shared_from_this());

    // release lock
    omp_unset_nest_lock(&node_copy->omp_lock);
}

void Face::disconnect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    // skip if not connected
    auto it = std::find(interior_points_.begin(), interior_points_.end(), interior_point);
    if (it == interior_points_.end()) return;

    // delete
    interior_points_.erase(it);
    interior_point->disconnect(shared_from_this());
    update_confirmed_status();

    // self destruct
    if (!deleting_) storage_->delete_face(shared_from_this());
}

void Face::disconnect(const std::shared_ptr<Face>& sibling_face)
{
    // check input
    if (sibling_face->is_expired()) return;

    // delete
    bool erased = sibling_faces_.erase(sibling_face);
    if (erased) sibling_face->disconnect(shared_from_this());
}

bool Face::is_connected_to_boundary_edges(std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& all_connected_faces, std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& all_connected_edges) const
{
    // check for each edge
    for (const std::shared_ptr<Edge>& edge : get_edges())
    {
        // add to visited edges
        const bool inserted = all_connected_edges.insert(edge).second;

        // skip if edge is already visited
        if (!inserted) continue;

        // return true if edge is boundary
        if (edge->is_boundary()) return true;
        
        // else, recursively check connected faces
        if (edge->is_connected_to_boundary_edges(all_connected_faces, all_connected_edges)) return true;        
    }

    return false;
}

bool Face::is_non_manifold() const
{
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> all_connected_faces;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> all_connected_edges;

    return !is_connected_to_boundary_edges(all_connected_faces, all_connected_edges);
}

void Face::swap(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{
    // if contains vertex1
    auto it = std::find(vertices_.begin(), vertices_.end(), vertex1);
    if (it != vertices_.end())
    {
        // std::cout << "Swapping face " << id_ << " vertex " << vertex1->get_id() << " with vertex " << vertex2->get_id() << std::endl;

        can_self_destruct_ = false;
        disconnect(vertex1);
        connect(vertex2);
        can_self_destruct_ = true;
    }
}

void Face::update_confirmed_status()
{
    bool previous_status = is_confirmed_;
    bool current_status = interior_points_.size() > 0;

    // if changed, update connected edges and vertices and interior points
    if (current_status != previous_status)
    {
        is_confirmed_ = current_status;
        for (const std::shared_ptr<Edge>& edge : edges_)
        {
            edge->update_confirmed_status();
        }
        for (const std::shared_ptr<Vertex>& vertex : vertices_)
        {
            vertex->update_confirmed_status();
        }
        for (const std::shared_ptr<InteriorPoint>& interior_point : interior_points_)
        {
            interior_point->update_confirmed_status();
        }
    }
}

bool Face::is_confirmed() const
{
    return is_confirmed_;
}

// swap surface1 with surface2
void Face::swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2)
{
    // if contains surfacce1    
    if (surface_ == surface1)
    {
        // std::cout << "Swapping face " << id_ << " surface " << surface1->get_id() << " with surface " << surface2->get_id() << std::endl;

        can_self_destruct_ = false;
        disconnect(surface1);
        connect(surface2);
        can_self_destruct_ = true;

        // cascade swap
        for (const std::shared_ptr<Vertex>& vertex : vertices_)
        {
            vertex->swap(surface1, surface2);
        }
        for (const std::shared_ptr<Edge>& edge : edges_)
        {
            edge->swap(surface1, surface2);
        }
        for (const std::shared_ptr<InteriorPoint>& interior_point : interior_points_)
        {
            interior_point->swap(surface1, surface2);
        }
    }
}

bool operator<(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired faces");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<Face>& lhs, const std::shared_ptr<Face>& rhs)
{
    if (!lhs && !rhs) return true; // true if both are nullptr
    if (!lhs || !rhs) return false; // false if either is nullptr
    // if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired faces");
    return lhs->get_id() == rhs->get_id();
}