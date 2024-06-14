#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

#include "MeshObject/GenericPoint.hpp"

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
{
    // set expired
    is_expired_ = false;

    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create vertex with invalid storage.");

    // get id
    id_ = storage->get_next_vertex_id();

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

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius)
{
    // set expired
    is_expired_ = false;
    
    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create vertex with invalid storage.");

    // get id
    id_ = storage->get_next_vertex_id();

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
    if (is_searchable_)
    {
        storage_->remove_searchable_vertex(shared_from_this());
        is_searchable_ = false;
    }

    // add to storage as generic point
    storage_->add_generic_point(get_position(), get_origin());

    // log
    std::cout << "---------- vertex " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& Vertex::get_id() const 
{ 
    return id_; 
}

const Eigen::Vector3d& Vertex::get_position() const 
{ 
    return position_; 
}

void Vertex::try_update_surface_projection()
{
    // update if surface changes
    if (normal_used_ != get_surface()->get_normal())
    {
        normal_used_ = get_surface()->get_normal();
        projected_position_ = get_surface()->compute_point_to_surface_position(get_origin(), get_position());
        projected_distance_ = get_surface()->compute_point_to_surface_distance(get_origin(), get_position());
    }
}

const Eigen::Vector3d& Vertex::get_projected_position()
{
    try_update_surface_projection();
    return projected_position_;
}

const double& Vertex::get_projected_distance()
{
    try_update_surface_projection();
    return projected_distance_;

    // if (std::fabs(distance) > 0.03)
    // {
    //     storage_->delete_vertex(shared_from_this());
    // }    
}

const Eigen::Vector3d& Vertex::get_origin() const 
{ 
    return origin_; 
}

const std::shared_ptr<Surface>& Vertex::get_surface() const
{    
    // return surfaces_.empty() ? std::shared_ptr<Surface>() : *surfaces_.begin();
    return *surfaces_.begin();
}

const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Vertex::get_edges() const 
{ 
    return edges_; 
}

Eigen::Vector2d Vertex::get_surface_coordinate()
{
    const Eigen::Matrix3d& eigenvectors = get_surface()->get_eigenvectors();
    if (eigenvectors_used_ == eigenvectors)
    {
        // use stored coordinate if eigenvectors are the same
        return surface_coordinate_;
    }
    else
    {
        // compute new coordinate
        Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();
        Eigen::Vector3d projected_position = get_projected_position();
        surface_coordinate_ = (projection_matrix.transpose() * projected_position).head<2>();
        eigenvectors_used_ = eigenvectors;
        return surface_coordinate_;
    }
}

bool Vertex::is_expired() const
{
    return is_expired_;
}

void Vertex::connect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge->connect(shared_from_this());

    // update boundary state
    update_boundary_state();
}

void Vertex::connect(const std::shared_ptr<Face>& face) 
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face->connect(shared_from_this());
}

void Vertex::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface->connect(shared_from_this());

    // cascade connect edges and faces to surface
    if (inserted) 
    {
        for (const std::shared_ptr<Edge>& edge : edges_)
        {
            edge->connect(surface);
        }
        for (const std::shared_ptr<Face>& face : faces_)
        {
            face->connect(surface);
        }
    }
}

void Vertex::disconnect(const std::shared_ptr<Edge>& edge) 
{
    // check input
    if (edge->is_expired()) return;

    // disconnect
    bool erased = edges_.erase(edge);
    if (erased) edge->disconnect(shared_from_this());

    // update boundary state
    update_boundary_state();

    // check self destruct
    if (!deleting_ && edges_.empty()) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Face>& face)
{
    // check pointer validity
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface->disconnect(shared_from_this());

    // cascade disconnect edges and faces from surface
    if (erased)
    {
        for (const std::shared_ptr<Edge>& edge : edges_)
        {
            edge->disconnect(surface);
        }
        for (const std::shared_ptr<Face>& face : faces_)
        {
            face->disconnect(surface);
        }
    }
}

void Vertex::update_boundary_state()
{
    if (deleting_) return;

    // becomes boundary when one of the connected edges is boundary, or when the point is alone
    is_boundary_ = false;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->is_boundary())
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
    if (is_boundary_ && !is_searchable_)
    {
        storage_->add_searchable_vertex(shared_from_this());
        is_searchable_ = true;
    }
    else if (!is_boundary_ && is_searchable_)
    {
        storage_->remove_searchable_vertex(shared_from_this());
        is_searchable_ = false;
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

bool operator<(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs)
{
    // when updating the third point's boundary state (due to first point deleted), the first point is already deleted. 
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs->get_id() < rhs->get_id();
}

bool operator<=(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs)
{
    // when updating the third point's boundary state (due to first point deleted), the first point is already deleted. 
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs->get_id() <= rhs->get_id();
}

bool operator==(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs->get_id() == rhs->get_id();
}