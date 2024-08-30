#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

#include "MeshObject/InteriorPoint.hpp"

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
    if (!is_searchable_)
    {
        storage_->add_searchable_face(shared_from_this());
        is_searchable_ = true;
    }

    // log
    std::cout << "Face " << id_ << " created between vertex " << vertex0->get_id() << ", vertex " << vertex1->get_id() << " and vertex " << vertex2->get_id() << std::endl;
}

void Face::delete_()
{
    // log
    std::cout << "Destroying face " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // if there is penetrating point, update the radius of face vertices
    if (storage_->has_penetrating_point())
    {
        // make a copy of vertices, as reduce_reverse_radius_search_radius will modify the set
        std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> copy_vertices = vertices_;
        for (const auto& vertex : copy_vertices)
        {
            // compute distance
            const Eigen::Vector3d& vertex_position = vertex->get_position();
            const Eigen::Vector3d& penetrating_point_position = storage_->get_penetrating_point();
            double distance = (vertex_position - penetrating_point_position).norm();
            vertex->reduce_reverse_radius_search_radius(distance);
        }
    }

    // disconnect
    // make a copy of the set to avoid iterator invalidation
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices = vertices_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges = edges_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points = interior_points_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> sibling_faces = sibling_faces_;
    for (const auto& vertex : vertices) disconnect(vertex);
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& interior_point : interior_points) disconnect(interior_point);
    if (surface_ != nullptr) disconnect(surface_);
    for (const auto& sibling_face : sibling_faces) disconnect(sibling_face);

    // remove from search tree
    if (is_searchable_)
    {
        storage_->remove_searchable_face(shared_from_this());
        is_searchable_ = false;
    }

    // log
    std::cout << "---------- face " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& Face::get_id() const
{
    return id_;
}

const Eigen::Vector3d& Face::get_center() const
{
    return center_;
}

const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& Face::get_vertices() const
{
    return vertices_;
}

const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& Face::get_interior_points() const
{
    return interior_points_;
}

const std::shared_ptr<Vertex>& Face::get_vertex(int index) const
{
    if (index < 0 || index > 2) throw std::runtime_error("Invalid index for vertex.");
    auto it = vertices_.begin();
    for (int i = 0; i < index; i++) it++;
    return *it;
}

const std::shared_ptr<Surface>& Face::get_surface() const
{
    return surface_;
}

const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Face::get_sibling_faces() const
{
    return sibling_faces_;
}

bool Face::is_expired() const
{
    return is_expired_;
}

bool Face::has_vertex(const std::shared_ptr<Vertex>& vertex) const
{
    return vertices_.find(vertex) != vertices_.end();
}

bool Face::intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction)
{    
    const Eigen::Vector3d& v0 = get_vertex(0)->get_position();
    const Eigen::Vector3d& v1 = get_vertex(1)->get_position();
    const Eigen::Vector3d& v2 = get_vertex(2)->get_position();

    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    Eigen::Vector3d pvec = direction.cross(edge2);
    double det = edge1.dot(pvec);
    if (std::fabs(det) < EPSILON) return false;

    double invDet = 1.0 / det;

    Eigen::Vector3d tvec = origin - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) return false;
    
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = direction.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return false;

    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) return false;

    return true;
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

    // connect
    bool inserted = vertices_.insert(vertex).second;
    if (inserted) vertex->connect(shared_from_this());

    // check size
    if (vertices_.size() > 3) throw std::runtime_error("Face connected to more than 3 vertices.");
}

void Face::connect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to connect face with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge->connect(shared_from_this());

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

    // connect
    bool inserted = interior_points_.insert(interior_point).second;
    if (inserted) interior_point->connect(shared_from_this());
    if (inserted) update_confirmed_status();
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

    // delete
    bool erased = vertices_.erase(vertex);
    if (erased) vertex->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && can_self_destruct_) storage_->delete_face(shared_from_this());
}

void Face::disconnect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) return;

    // delete
    bool erased = edges_.erase(edge);
    if (erased) edge->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_->delete_face(shared_from_this());
}

void Face::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surface_ == surface;
    if (erased) surface->disconnect(shared_from_this());
    if (erased) surface_ = nullptr;

    // self destruct
    if (!deleting_ && erased && can_self_destruct_) storage_->delete_face(shared_from_this());
}

void Face::disconnect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    // delete
    bool erased = interior_points_.erase(interior_point);
    if (erased) interior_point->disconnect(shared_from_this());
    if (erased) update_confirmed_status();

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

void Face::swap(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{
    // if contains vertex1
    if (vertices_.find(vertex1) != vertices_.end())
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
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired faces");
    return lhs->get_id() == rhs->get_id();
}