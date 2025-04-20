#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

Settings Edge::settings_;

void Edge::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
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

    // put into vertices list
    {
        // put into vertices list
        vertices_.push_back(vertex1_valid);
        vertices_.push_back(vertex2_valid);
    }

    // connect to vertices
    {
        // connect to vertices
        vertex1_valid->connect(shared_from_this());
        vertex2_valid->connect(shared_from_this());
    }

    // connect to surface
    {
        // connect to surface
        surface_ = surface;
        surface->connect(shared_from_this());
    }

    // center min and max are computed once, and are not updated, otherwise BVH tree will cause search error
    // compute center
    center_ = 0.5 * (vertex1_valid->get_position() + vertex2_valid->get_position());

    // compute max and min
    double margin = 0.1f;
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
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_lifecycle_);

    // skip if already deleted
    if (is_expired_) return;
    
    // set deletion flag
    deleting_ = true;

    // log
    if (settings_.log.deletion) std::cout << "Destroying edge " << id_ << std::endl;

    // surface (disconnect)
    {
        // disconnect
        surface_.lock()->disconnect(shared_from_this());

        // add to be split
        storage_.lock()->add_surface_to_be_split(surface_.lock());
        
        // clear surface
        surface_.reset();
    }

    // vertices (disconnect)
    {
        // disconnect from vertices
        for (const auto& vertex : vertices_)
        {
            if (vertex.lock())
            {
                // check if vertex is expired
                if (vertex.lock()->is_expired()) continue;

                // disconnect
                vertex.lock()->disconnect(shared_from_this());
            }
        } 

        // clear vertices
        vertices_.clear();
    }

    // faces (delete)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock_faces(rwlock_faces_);

        // delete face
        for (const auto& face : faces_) storage_.lock()->add_face_to_be_deleted(face.lock());

        // clear faces
        faces_.clear();
    }
    
    // log
    if (settings_.log.deletion) std::cout << "---------- edge " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

Edge::~Edge()
{
    // std::cout << "Edge " << id_ << " deleted." << std::endl;
}

void Edge::connect(const std::shared_ptr<Face>& face) 
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect edge with invalid face.");

    {
        // write lock
        std::unique_lock lock(rwlock_faces_);

        // skip if already connected
        for (const std::weak_ptr<Face>& face_ : faces_) if (face_.lock() == face) return;

        // connect
        faces_.push_back(face);
        // if (faces_.size() > 2) throw std::runtime_error("Edge connected to more than 2 faces.");
    }
    
    // update boundary state
    update_boundary_state();
}

void Edge::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    {
        // write lock
        std::unique_lock lock(rwlock_faces_);

        // skip if not connected
        auto it = std::find(faces_.begin(), faces_.end(), face);
        if (it == faces_.end()) return;

        // disconnect
        faces_.erase(it);
    }

    // update boundary state
    update_boundary_state();

    // do not self destruct when have no face
    // check self destruct
    if (faces_.empty() && !deleting_ && can_self_destruct_) storage_.lock()->add_edge_to_be_deleted(shared_from_this());
}


// swap surface1 with surface2
void Edge::swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2)
{
    if (surface_.lock() == surface1)
    {
        surface_ = surface2;
        surface1->disconnect(shared_from_this());
        surface2->connect(shared_from_this());
    }
}

void Edge::set_can_self_destruct(bool can_self_destruct)
{
    can_self_destruct_ = can_self_destruct;
}

bool Edge::is_connected_to_boundary_edges(std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& all_connected_faces, std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& all_connected_edges) const
{
    // check for each face
    for (std::weak_ptr<Face> face : get_faces())
    {
        // add to visited faces
        const bool inserted = all_connected_faces.insert(face.lock()).second;

        // skip if face is already visited
        if (!inserted) continue;

        // return true if face is boundary
        if (face.lock()->is_connected_to_boundary_edges(all_connected_faces, all_connected_edges)) return true;
    }

    return false;
}

const int& Edge::get_id() const 
{ 
    return id_; 
}

std::shared_ptr<Vertex> Edge::get_vertex(int index) const
{
    if (is_expired_) throw std::runtime_error("Accessing expired edge.");
    if (index < 0 || index > 1) throw std::runtime_error("Invalid vertex index.");
    if (vertices_.size() != 2) throw std::runtime_error("Edge does not have 2 vertices.");
    return (*std::next(vertices_.begin(), index)).lock();
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

    // read lock
    std::shared_lock lock(rwlock_faces_);

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
        // update connected vertices
        for (const std::weak_ptr<Vertex>& vertex : vertices_)
        {
            vertex.lock()->check_if_update_search_tree();
        }
    }
}

bool Edge::is_non_manifold() const
{
    // non manifold if connected to more than 2 faces
    return faces_.size() > 2;
}

std::shared_ptr<Surface> Edge::get_surface() const
{
    return surface_.lock();
}

std::vector<std::weak_ptr<Face>> Edge::get_faces() const
{
    // read lock
    std::shared_lock lock(rwlock_faces_);
    return faces_;
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

bool Edge::intersects_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
{
    // edge can never be not expired without vertices

    // skip if deleting
    if (is_deleting()) return false;
    
    // skip if vertices are connected
    if (vertices_[0] == vertex0 || vertices_[0] == vertex1 || vertices_[1] == vertex0 || vertices_[1] == vertex1) return false;

    // get surface coordinates
    const Eigen::Vector2d& p1 = vertex0->get_surface_coordinate(surface_.lock());
    const Eigen::Vector2d& p2 = vertex1->get_surface_coordinate(surface_.lock());
    const Eigen::Vector2d& q1 = vertices_[0].lock()->get_surface_coordinate(surface_.lock());
    const Eigen::Vector2d& q2 = vertices_[1].lock()->get_surface_coordinate(surface_.lock());
    
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