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

    // add to vertices
    vertices_.push_back(vertex0);
    vertices_.push_back(vertex1);
    vertices_.push_back(vertex2);

    // connect vertices
    vertex0->connect(shared_from_this());
    vertex1->connect(shared_from_this());
    vertex2->connect(shared_from_this());

    // store vertices position
    v0_ = vertex0->get_position();
    v1_ = vertex1->get_position();
    v2_ = vertex2->get_position();

    // compute area
    Eigen::Vector3d v0_v1 = v1_ - v0_;
    Eigen::Vector3d v0_v2 = v2_ - v0_;
    area_ = 0.5 * v0_v1.cross(v0_v2).norm();

    // connect surface
    {
        // connect surface
        surface_ = surface;
        surface_->connect(shared_from_this());
    }

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
    std::shared_ptr<Edge> edge0 = vertex0->get_edge(vertex1);
    std::shared_ptr<Edge> edge1 = vertex1->get_edge(vertex2);
    std::shared_ptr<Edge> edge2 = vertex2->get_edge(vertex0);
    
    // add to edges
    edges_.push_back(edge0);
    edges_.push_back(edge1);
    edges_.push_back(edge2);
    
    // connect edges
    edge0->connect(shared_from_this());
    edge1->connect(shared_from_this());
    edge2->connect(shared_from_this());

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

void Face::initialize_(
    const std::shared_ptr<Storage>& storage, 
    const std::shared_ptr<Surface> surface, 
    const std::shared_ptr<Vertex>& vertex0, 
    const std::shared_ptr<Vertex>& vertex1, 
    const std::shared_ptr<Vertex>& vertex2,
    const std::shared_ptr<Edge>& edge0,
    const std::shared_ptr<Edge>& edge1,
    const std::shared_ptr<Edge>& edge2)
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

    // add to vertices
    vertices_.push_back(vertex0);
    vertices_.push_back(vertex1);
    vertices_.push_back(vertex2);

    // connect vertices
    vertex0->connect(shared_from_this());
    vertex1->connect(shared_from_this());
    vertex2->connect(shared_from_this());

    // store vertices position
    v0_ = vertex0->get_position();
    v1_ = vertex1->get_position();
    v2_ = vertex2->get_position();

    // compute area
    Eigen::Vector3d v0_v1 = v1_ - v0_;
    Eigen::Vector3d v0_v2 = v2_ - v0_;
    area_ = 0.5 * v0_v1.cross(v0_v2).norm();

    // connect surface
    {
        // connect surface
        surface_ = surface;
        surface_->connect(shared_from_this());
    }

    // compute min and max
    BoundingBox box;
    box.expand(vertex0->get_position());
    box.expand(vertex1->get_position());
    box.expand(vertex2->get_position());
    min_ = box.min;
    max_ = box.max;

    // first vertex
    first_vertex_ = vertex0;

    // add to edges
    edges_.push_back(edge0);
    edges_.push_back(edge1);
    edges_.push_back(edge2);
    
    // connect edges
    edge0->connect(shared_from_this());
    edge1->connect(shared_from_this());
    edge2->connect(shared_from_this());

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
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_lifecycle_);
    
    // skip if already deleted
    if (is_expired_) return;
    
    // set deletion flag
    deleting_ = true;

    // log
    if (settings_.log.deletion) std::cout << "Destroying face " << id_ << std::endl;

    // surface (disconnect)
    {
        // disconnect surface
        surface_->disconnect(shared_from_this());

        // clear surface
        surface_ = nullptr;
    }

    // vertices (disconnect)
    {
        // disconnect vertices
        for (const auto& vertex : vertices_) vertex->disconnect(shared_from_this());

        // clear vertices
        vertices_.clear();
    }

    // edges (disconnect)
    {
        // disconnect edges
        for (const auto& edge : edges_) edge->disconnect(shared_from_this());

        // clear edges
        edges_.clear();
    }

    // interior points (delete)
    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_points_);

        // delete interior points
        for (const auto& interior_point : interior_points_) storage_->add_interior_point_to_be_deleted(interior_point);

        // clear interior points
        interior_points_.clear();
    }

    // update tree
    {
        // add to affected faces set
        storage_->add_to_set_of_faces_to_update_bvh_tree(shared_from_this());
    }

    // log
    if (settings_.log.deletion) std::cout << "---------- face " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

Face::~Face()
{
    // std::cout << "Face " << id_ << " deleted." << std::endl;
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

std::vector<std::shared_ptr<InteriorPoint>> Face::get_interior_points() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_interior_points_);
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

double Face::get_area() const
{
    return area_;
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

void Face::connect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to connect face with invalid interior point.");

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_points_);

        // skip if already connected
        for (const std::shared_ptr<InteriorPoint>& interior_point_ : interior_points_) if (interior_point_ == interior_point) return;

        // connect
        interior_points_.push_back(interior_point);
    }
}

void Face::disconnect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_points_);

        // skip if not connected
        auto it = std::find(interior_points_.begin(), interior_points_.end(), interior_point);
        if (it == interior_points_.end()) return;

        // delete
        interior_points_.erase(it);
    }
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