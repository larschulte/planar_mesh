#include "MeshObject/Storage.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Face>& face, const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
{
    // set expired
    is_expired_ = false;

    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create interior point with invalid storage.");

    // store
    storage_ = storage;

    // get id
    id_ = storage_->get_next_interior_point_id();

    // store
    position_ = position;
    origin_ = origin;

    // connect
    connect(face);
    connect(face->get_surface());

    // log
    std::cout << "InteriorPoint " << id_ << " created.\n";
}

void InteriorPoint::delete_()
{
    // log
    std::cout << "Destroying InteriorPoint " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    while (!faces_.empty())
    {
        disconnect(*faces_.begin());
    }
    while (!surfaces_.empty())
    {
        disconnect(*surfaces_.begin());
    }

    // add to storage as generic point
    storage_->add_generic_point(get_position(), get_origin());

    // log
    std::cout << "---------- InteriorPoint " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& InteriorPoint::get_id() const
{
    return id_;
}

const Eigen::Vector3d& InteriorPoint::get_position() const
{
    return position_;
}

const Eigen::Vector3d& InteriorPoint::get_origin() const
{
    return origin_;
}

const std::shared_ptr<Surface>& InteriorPoint::get_surface() const
{    
    return *surfaces_.begin();
}

void InteriorPoint::try_update_surface_projection()
{
    // update if surface changes
    if (normal_used_ != get_surface()->get_normal())
    {
        normal_used_ = get_surface()->get_normal();
        projected_position_ = get_surface()->compute_point_to_surface_position(get_origin(), get_position());
        projected_distance_ = get_surface()->compute_point_to_surface_distance(get_origin(), get_position());
    }
}

const Eigen::Vector3d& InteriorPoint::get_projected_position()
{
    try_update_surface_projection();
    return projected_position_;
}

const double& InteriorPoint::get_projected_distance()
{
    try_update_surface_projection();
    return projected_distance_;

    // if (std::fabs(distance) > 0.03)
    // {
    //     storage_->delete_vertex(shared_from_this());
    // }    
}

bool InteriorPoint::is_expired() const
{
    return is_expired_;
}

void InteriorPoint::connect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect interior point with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face->connect(shared_from_this());
}

void InteriorPoint::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect interior point with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface->connect(shared_from_this());
}

void InteriorPoint::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && faces_.empty()) storage_->delete_interior_point(shared_from_this());
}

void InteriorPoint::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && surfaces_.empty()) storage_->delete_interior_point(shared_from_this());
}

bool operator<(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired InteriorPoints");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired InteriorPoints");
    return lhs->get_id() == rhs->get_id();
}
