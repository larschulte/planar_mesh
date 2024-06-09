#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

#include "MeshObject/GenericPoint.hpp"

void Vertex::initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create vertex with invalid storage.");

    // get id
    id_ = storage.lock()->get_next_vertex_id();

    // store
    storage_ = storage;
    position_ = position;
    origin_ = origin;

    // compute reverse search radius based on distance
    double distance_to_radius_ratio = tan(4 * M_PI / 180);
    double distance = (position - origin).norm();
    double radius = distance * distance_to_radius_ratio;
    set_reverse_radius_search_radius(radius);

    // update boundary state
    update_boundary_state();

    // log
    std::cout << "Vertex " << id_ << " created.\n";
}

void Vertex::initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin, double radius)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create vertex with invalid storage.");

    // get id
    id_ = storage.lock()->get_next_vertex_id();

    // store
    storage_ = storage;
    position_ = position;
    origin_ = origin;

    // set reverse search radius based on input parameter
    set_reverse_radius_search_radius(radius);

    // update boundary state
    update_boundary_state();

    // log
    std::cout << "Vertex " << id_ << " created.\n";
}

void Vertex::delete_()
{
    // log
    std::cout << "Destroying vertex " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    while (!edges_.empty())
    {
        disconnect(*edges_.begin());
    }
    while (!faces_.empty())
    {
        disconnect(*faces_.begin());
    }
    while (!surfaces_.empty())
    {
        disconnect(*surfaces_.begin());
    }

    // remove from search tree
    if (is_boundary_)
    {
        storage_.lock()->remove_searchable_vertex(shared_from_this());
    }

    // add to storage as generic point
    storage_.lock()->add_generic_point(get_position(), get_origin());

    // log
    std::cout << "---------- vertex " << id_ << " destroyed" << std::endl;
}

int Vertex::get_id() const 
{ 
    return id_; 
}

Eigen::Vector3d Vertex::get_position() const 
{ 
    return position_; 
}

Eigen::Vector3d Vertex::get_projected_position() const
{
    return get_surface().lock()->compute_point_to_surface_position(get_origin(), get_position());
}

Eigen::Vector3d Vertex::get_origin() const 
{ 
    return origin_; 
}

std::weak_ptr<Surface> Vertex::get_surface() const
{    
    // return surfaces_.empty() ? std::weak_ptr<Surface>() : *surfaces_.begin();
    return *surfaces_.begin();
}

void Vertex::connect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to connect vertex with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge.lock()->connect(shared_from_this());

    // update boundary state
    update_boundary_state();
}

void Vertex::connect(std::weak_ptr<Face> face) 
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to connect vertex with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face.lock()->connect(shared_from_this());
}

void Vertex::connect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to connect vertex with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface.lock()->connect(shared_from_this());
}

void Vertex::disconnect(std::weak_ptr<Edge> edge) 
{
    // check input
    if (edge.expired()) return;

    // disconnect
    bool erased = edges_.erase(edge);
    if (erased) edge.lock()->disconnect(shared_from_this());

    // update boundary state
    update_boundary_state();

    // check self destruct
    if (!deleting_ && edges_.empty()) storage_.lock()->delete_vertex(shared_from_this());
}

void Vertex::disconnect(std::weak_ptr<Face> face)
{
    // check pointer validity
    if (face.expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face.lock()->disconnect(shared_from_this());
}

void Vertex::disconnect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface.lock()->disconnect(shared_from_this());
}

void Vertex::update_boundary_state()
{
    // becomes boundary when one of the connected edges is boundary, or when the point is alone
    is_boundary_ = false;
    for (std::weak_ptr<Edge> edge : edges_)
    {
        if (edge.lock()->is_boundary())
        {
            is_boundary_ = true;
            break;
        }
    }
    if (edges_.empty())
    {
        is_boundary_ = true;
    }

    // update search tree
    if (is_boundary_)
    {
        if (storage_.expired()) throw std::runtime_error("Storage expired in vertex update boundary state.");
        storage_.lock()->add_searchable_vertex(shared_from_this());
    }
    else
    {
        storage_.lock()->remove_searchable_vertex(shared_from_this());
    }
}

void Vertex::set_reverse_radius_search_radius(double radius)
{
    // set radius
    reverse_search_radius_ = radius;

    // update min and max
    min_ = position_ - Eigen::Vector3d(radius, radius, radius);
    max_ = position_ + Eigen::Vector3d(radius, radius, radius);

    // should update search tree if expand radius
}

Eigen::Vector3d Vertex::get_min() const
{
    return min_;
}

Eigen::Vector3d Vertex::get_max() const
{
    return max_;
}

double Vertex::get_radius() const
{
    return reverse_search_radius_;
}

bool Vertex::contains(const Eigen::Vector3d& point) const
{
    return (point - position_).norm() < reverse_search_radius_;
}

bool Vertex::approx_contains(const Eigen::Vector3d& point) const
{
    // by comparing to max and min
    return (point.x() > min_.x() && point.x() < max_.x() &&
            point.y() > min_.y() && point.y() < max_.y() &&
            point.z() > min_.z() && point.z() < max_.z());
}

bool operator<(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs)
{
    // when updating the third point's boundary state (due to first point deleted), the first point is already deleted. 
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}