#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

void Edge::initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create edge with invalid storage.");
    if (vertex1.expired()) throw std::runtime_error("Attempts to create edge with invalid vertex1.");
    if (vertex2.expired()) throw std::runtime_error("Attempts to create edge with invalid vertex2.");
    auto storage_valid = storage.lock();
    auto vertex1_valid = vertex1.lock();
    auto vertex2_valid = vertex2.lock();

    // store
    storage_ = storage;

    // get id
    id_ = storage_valid->get_next_edge_id();

    // sort
    if (vertex1_valid->get_id() > vertex2_valid->get_id()) std::swap(vertex1_valid, vertex2_valid);

    // make connections
    connect(vertex1);
    connect(vertex2);

    // compute center
    center_ = 0.5 * (vertex1_valid->get_position() + vertex2_valid->get_position());

    // compute max and min
    double margin = 0.05;
    max_ = vertex1.lock()->get_position().cwiseMax(vertex2.lock()->get_position()) + margin * Eigen::Vector3d::Ones();
    min_ = vertex1.lock()->get_position().cwiseMin(vertex2.lock()->get_position()) - margin * Eigen::Vector3d::Ones();

    // log
    std::cout << "Edge " << id_ << " created between vertex " << vertex1_valid->get_id() << " and vertex " << vertex2_valid->get_id() << std::endl;
}

void Edge::delete_()
{
    // log
    std::cout << "Destroying edge " << id_ << std::endl;
    

    // set deletion flag
    deleting_ = true;

    // disconnect
    while (!vertices_.empty())
    {
        disconnect(*vertices_.begin());
    }
    while (!faces_.empty())
    {
        disconnect(*faces_.begin());
    }
    while (!surfaces_.empty())
    {
        disconnect(*surfaces_.begin());
    }

    // log
    std::cout << "---------- edge " << id_ << " destroyed" << std::endl;
}

void Edge::connect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) throw std::runtime_error("Attempts to connect edge with invalid vertex.");

    // connect
    bool inserted = vertices_.insert(vertex).second;
    if (inserted) vertex.lock()->connect(shared_from_this());

    // check size
    if (vertices_.size() > 2) throw std::runtime_error("Edge connected to more than 2 vertices.");
}

void Edge::connect(std::weak_ptr<Face> face) 
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to connect edge with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face.lock()->connect(shared_from_this());

    // update boundary state
    update_boundary_state();
}

void Edge::connect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to connect edge with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface.lock()->connect(shared_from_this());
}

void Edge::disconnect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) return;

    // disconnect
    bool erased = vertices_.erase(vertex);
    if (erased) vertex.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_.lock()->delete_edge(shared_from_this());
}

void Edge::disconnect(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face.lock()->disconnect(shared_from_this());

    // update boundary state
    update_boundary_state();

    // check self destruct
    if (faces_.empty())
    {
        if (!deleting_) storage_.lock()->delete_edge(shared_from_this());
    }
}

void Edge::disconnect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface.lock()->disconnect(shared_from_this());
}

int Edge::get_id() const 
{ 
    return id_; 
}

std::weak_ptr<Vertex> Edge::get_vertex(int index) const
{
    if (index < 0 || index > 1) throw std::runtime_error("Invalid vertex index.");
    if (vertices_.size() != 2) throw std::runtime_error("Edge does not have 2 vertices.");
    return *std::next(vertices_.begin(), index);
}

// has vertex
bool Edge::has_vertex(std::weak_ptr<Vertex> vertex) const
{
    return vertices_.find(vertex) != vertices_.end();
}

bool Edge::is_boundary() const
{
    return is_boundary_;
}

void Edge::update_boundary_state()
{
    // becomes boundary when 0 or 1 face connected
    if (faces_.size() <= 1) 
    {
        is_boundary_ = true;
    }
    else 
    {
        is_boundary_ = false;
    }

    // update boundary state of connected vertices
    for (std::weak_ptr<Vertex> vertex : vertices_)
    {
        vertex.lock()->update_boundary_state();
    }
}

Eigen::Vector3d Edge::get_center() const
{
    return center_;
}

Eigen::Vector3d Edge::get_max() const
{
    return max_;
}

Eigen::Vector3d Edge::get_min() const
{
    return min_;
}

bool Edge::intersects_edge(std::weak_ptr<Vertex> vertex0, std::weak_ptr<Vertex> vertex1)
{
    // skip if vertices are connected
    if (has_vertex(vertex0) || has_vertex(vertex1)) return false;

    // get surface coordinates
    Eigen::Vector2d p1 = vertex0.lock()->get_surface_coordinate();
    Eigen::Vector2d p2 = vertex1.lock()->get_surface_coordinate();
    Eigen::Vector2d q1 = get_vertex(0).lock()->get_surface_coordinate();
    Eigen::Vector2d q2 = get_vertex(1).lock()->get_surface_coordinate();
    
    // check if edge intersects
    return segments_intersect(p1, p2, q1, q2);
}

bool operator<(const std::weak_ptr<Edge>& lhs, const std::weak_ptr<Edge>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Edge>& lhs, const std::weak_ptr<Edge>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
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