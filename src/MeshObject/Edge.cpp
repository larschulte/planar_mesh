#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

Settings Edge::settings_;

void Edge::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{
    // set expired
    is_expired_ = false;

    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create edge with invalid storage.");
    if (vertex1->is_expired()) throw std::runtime_error("Attempts to create edge with invalid vertex1.");
    if (vertex2->is_expired()) throw std::runtime_error("Attempts to create edge with invalid vertex2.");
    auto storage_valid = storage;
    auto vertex1_valid = vertex1;
    auto vertex2_valid = vertex2;

    // store
    storage_ = storage;

    // get id
    id_ = storage_valid->get_next_edge_id();

    // sort
    if (vertex1_valid->get_id() > vertex2_valid->get_id()) std::swap(vertex1_valid, vertex2_valid);

    // make connections
    connect(vertex1);
    connect(vertex2);

    // center min and max are computed once, and are not updated, otherwise BVH tree will cause search error
    // compute center
    center_ = 0.5 * (vertex1_valid->get_position() + vertex2_valid->get_position());

    // compute max and min
    double margin = settings_.range_accuracy + settings_.envelope_size * settings_.range_precision;
    max_ = vertex1->get_position().cwiseMax(vertex2->get_position()) + margin * Eigen::Vector3d::Ones();
    min_ = vertex1->get_position().cwiseMin(vertex2->get_position()) - margin * Eigen::Vector3d::Ones();

    // compute length
    length_ = (vertex1->get_position() - vertex2->get_position()).norm();

    // update boundary state
    update_boundary_state();

    // log
    if (settings_.log.initialize) std::cout << "Edge " << id_ << " created between vertex " << vertex1_valid->get_id() << " and vertex " << vertex2_valid->get_id() << std::endl;
}

void Edge::delete_()
{
    // log
    if (settings_.log.deletion) std::cout << "Destroying edge " << id_ << std::endl;
    

    // set deletion flag
    deleting_ = true;

    // disconnect
    // make copy of vertices_ and faces_ to avoid iterator invalidation
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices = vertices_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> sibling_edges = sibling_edges_;
    for (const auto& vertex : vertices) disconnect(vertex);
    for (const auto& face : faces) disconnect(face);
    if (surface_ != nullptr) disconnect(surface_);
    for (const auto& sibling_edge : sibling_edges) disconnect(sibling_edge);

    // log
    if (settings_.log.deletion) std::cout << "---------- edge " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

void Edge::connect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to connect edge with invalid vertex.");

    // connect
    bool inserted = vertices_.insert(vertex).second;
    if (inserted) vertex->connect(shared_from_this());

    // check size
    if (vertices_.size() > 2) throw std::runtime_error("Edge connected to more than 2 vertices.");
}

void Edge::connect(const std::shared_ptr<Face>& face) 
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect edge with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face->connect(shared_from_this());

    // update boundary state
    if (inserted) update_boundary_state();

    // update confirmed status
    if (inserted) update_confirmed_status();
    if (inserted) update_singular_state();
}

void Edge::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect edge with invalid surface.");

    // connect
    bool inserted = surface_ != surface;
    if (inserted) surface_ = surface;
    if (inserted) surface->connect(shared_from_this());
    if (inserted) is_searchable_ = false;
    if (inserted) is_boundary_ = false;
    if (inserted) is_singular_ = true;
    if (inserted) update_singular_state();
    if (inserted) update_boundary_state();
}

void Edge::connect(const std::shared_ptr<Edge>& sibling_edge)
{
    // check input
    if (sibling_edge->is_expired()) throw std::runtime_error("Attempts to connect edge with invalid sibling edge.");

    // skip if try to connect to itself
    if (sibling_edge == shared_from_this()) return;

    // connect
    bool inserted = sibling_edges_.insert(sibling_edge).second;
    if (inserted) std::cout << "Connected edge " << id_ << " with edge " << sibling_edge->get_id() << " as sibling."<< std::endl;
    if (inserted) sibling_edge->connect(shared_from_this());
    if (inserted)
    {
        for (const std::shared_ptr<Edge>& sibling_edge_ : sibling_edges_)
        {
            sibling_edge_->connect(sibling_edge);
        }
    }
}

void Edge::disconnect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) return;

    // disconnect
    bool erased = vertices_.erase(vertex);
    if (erased) vertex->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && can_self_destruct_) storage_->delete_edge(shared_from_this());
}

void Edge::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());

    // update boundary state
    if (erased) update_boundary_state();

    // update confirmed status
    if (erased) update_confirmed_status();
    if (erased) update_singular_state();

    // do not self destruct when have no face
    // check self destruct
    if (erased && faces_.empty() && !deleting_ && can_self_destruct_) storage_->delete_edge(shared_from_this());
}

void Edge::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surface_ == surface;
    if (erased) surface->disconnect(shared_from_this());
    if (erased) remove_searchable_state();
    if (erased) is_boundary_ = false;
    if (erased) is_singular_ = false;
    if (erased) is_searchable_ = false;
    if (erased) surface_ = nullptr;

    // check self destruct
    if (!deleting_ && erased && can_self_destruct_) storage_->delete_edge(shared_from_this());
}

void Edge::disconnect(const std::shared_ptr<Edge>& sibling_edge)
{
    // check input
    if (sibling_edge->is_expired()) return;

    // disconnect
    bool erased = sibling_edges_.erase(sibling_edge);
    if (erased) sibling_edge->disconnect(shared_from_this());
}

void Edge::swap(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{
    // make sure the position of vertex 1 and 2 are identical, otherwise need to update BVH which is not yet implemented
    if ((vertex1->get_position() - vertex2->get_position()).norm() > 1e-8) throw std::runtime_error("Swapping edge with non-identical vertices.");
    
    if (vertices_.find(vertex1) != vertices_.end())
    {
        // std::cout << "Swapping edge " << id_ << " vertex " << vertex1->get_id() << " with vertex " << vertex2->get_id() << std::endl;
        
        can_self_destruct_ = false;
        disconnect(vertex1);
        connect(vertex2);
        can_self_destruct_ = true;
    }
}

void Edge::update_confirmed_status()
{
    // update number of confirmed faces
    num_confirmed_faces = 0;
    for (const std::shared_ptr<Face>& face : faces_)
    {
        if (face->is_confirmed()) num_confirmed_faces++;
    }

    // update confirmed status
    if (num_confirmed_faces >= 1) is_confirmed_ = true;
    else is_confirmed_ = false;
}

bool Edge::is_confirmed() const
{
    return is_confirmed_;
}

// swap surface1 with surface2
void Edge::swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2)
{
    if (surface_ == surface1)
    {
        // std::cout << "Swapping edge " << id_ << " surface " << surface1->get_id() << " with surface " << surface2->get_id() << std::endl;

        can_self_destruct_ = false;
        disconnect(surface1);
        connect(surface2);
        can_self_destruct_ = true;

        // cascade swap
        for (const std::shared_ptr<Vertex>& vertex : vertices_)
        {
            vertex->swap(surface1, surface2);
        }
        for (const std::shared_ptr<Face>& face : faces_)
        {
            face->swap(surface1, surface2);
        }
    }
}

const int& Edge::get_id() const 
{ 
    return id_; 
}

const std::shared_ptr<Vertex>& Edge::get_vertex(int index) const
{
    if (is_expired_) throw std::runtime_error("Accessing expired edge.");
    if (index < 0 || index > 1) throw std::runtime_error("Invalid vertex index.");
    if (vertices_.size() != 2) throw std::runtime_error("Edge does not have 2 vertices.");
    return *std::next(vertices_.begin(), index);
}

bool Edge::is_expired() const
{
    return is_expired_;
}

// has vertex
bool Edge::has_vertex(const std::shared_ptr<Vertex>& vertex) const
{
    return get_vertex(0) == vertex || get_vertex(1) == vertex;
}

bool Edge::is_boundary() const
{
    return is_boundary_;
}

bool Edge::is_singular() const
{
    return is_singular_;
}

bool Edge::is_deleting() const
{
    return deleting_;
}

void Edge::update_boundary_state()
{
    if (deleting_) return;

    // update boundary state
    bool previous_boundary_state = is_boundary_;
    // count number of faces in this surface
    int num_faces_in_this_surface = faces_.size();
    if (num_faces_in_this_surface <= 1) 
    {
        is_boundary_ = true;
    }
    else 
    {
        is_boundary_ = false;
    }
    bool changed_boundary_state = previous_boundary_state != is_boundary_;

    // if changed state
    if (changed_boundary_state) 
    {
        // update connected surfaces
        update_searchable_state();

        // update connected vertices
        for (const std::shared_ptr<Vertex>& vertex : vertices_)
        {
            vertex->check_if_update_search_tree();
        }
    }
}

void Edge::update_singular_state()
{
    // count number of faces in this surface
    int num_faces_in_this_surface = faces_.size();

    // update singular state
    if (num_faces_in_this_surface == 0) is_singular_ = true;
    else is_singular_ = false; 
}

void Edge::update_searchable_state()
{
    if (surface_ == nullptr) return;

    // upon changes of boundary state
    // add
    if (is_boundary_ && !is_searchable_)
    {
        // surface_->add_searchable_edge(shared_from_this());
        is_searchable_ = true;
    }
    // remove
    if (!is_boundary_ && is_searchable_)
    {   
        // surface_->remove_searchable_edge(shared_from_this());
        is_searchable_ = false;
    }
}

void Edge::remove_searchable_state()
{
    if (is_searchable_)
    {
        // surface_->remove_searchable_edge(shared_from_this());
        is_searchable_ = false;
    }
}

bool Edge::is_non_manifold() const
{
    // non manifold if connected to more than 2 faces
    return faces_.size() > 2;
}

const std::shared_ptr<Surface>& Edge::get_surface() const
{
    return surface_;
}

const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Edge::get_faces() const
{
    return faces_;
}

const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Edge::get_sibling_edges() const
{
    return sibling_edges_;
}

const Eigen::Vector3d& Edge::get_center() const
{
    return center_;
}

const Eigen::Vector3d& Edge::get_max() const
{
    return max_;
}

const Eigen::Vector3d& Edge::get_min() const
{
    return min_;
}

const double& Edge::get_length() const
{
    return length_;
}

bool Edge::intersects_edge(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
{
    // skip if vertices are connected
    if (has_vertex(vertex0) || has_vertex(vertex1)) return false;

    // get surface coordinates
    const Eigen::Vector2d& p1 = vertex0->get_surface_coordinate(surface);
    const Eigen::Vector2d& p2 = vertex1->get_surface_coordinate(surface);
    const Eigen::Vector2d& q1 = get_vertex(0)->get_surface_coordinate(surface);
    const Eigen::Vector2d& q2 = get_vertex(1)->get_surface_coordinate(surface);
    
    // check if edge intersects
    return segments_intersect(p1, p2, q1, q2);
}

bool operator<(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs)
{
    if (!lhs && !rhs) return true; // true if both are nullptr
    if (!lhs || !rhs) return false; // false if either is nullptr
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs->get_id() == rhs->get_id();
}

// Function to check if two 2D segments (p1, p2) and (q1, q2) intersect
bool segments_intersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2) 
{
    auto orientation = [](const Eigen::Vector2d &p, const Eigen::Vector2d &q, const Eigen::Vector2d &r) 
    {
        double val = (q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y());
        if (val == 0) return 0;  // collinear
        return (val > 0) ? 1 : 2; // clock or counterclockwise
    };

    int o1 = orientation(p1, p2, q1);
    int o2 = orientation(p1, p2, q2);
    int o3 = orientation(q1, q2, p1);
    int o4 = orientation(q1, q2, p2);
    if (o1 != o2 && o3 != o4) 
    {
        return true;
    }
    else
    {
        return false;
    }
}