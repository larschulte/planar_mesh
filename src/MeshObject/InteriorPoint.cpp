#include "MeshObject/Storage.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"

void InteriorPoint::initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Face> face, Eigen::Vector3d position, Eigen::Vector3d origin)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create interior point with invalid storage.");

    // store
    storage_ = storage;

    // get id
    id_ = storage_.lock()->get_next_interior_point_id();

    // store
    position_ = position;
    origin_ = origin;

    // connect
    connect(face);
    connect(face.lock()->get_surface());

    // log
    std::cout << "InteriorPoint " << id_ << " created.\n";
}

void InteriorPoint::delete_()
{
    // log
    std::cout << "Destroying InteriorPoint " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // add to storage as generic point
    storage_.lock()->add_generic_point(get_position(), get_origin());

    // log
    std::cout << "---------- InteriorPoint " << id_ << " destroyed" << std::endl;
}

int InteriorPoint::get_id() const
{
    return id_;
}

Eigen::Vector3d InteriorPoint::get_position() const
{
    return position_;
}

Eigen::Vector3d InteriorPoint::get_origin() const
{
    return origin_;
}

void InteriorPoint::connect(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) throw std::runtime_error("Attempts to connect interior point with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face.lock()->connect(shared_from_this());
}

void InteriorPoint::connect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to connect interior point with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface.lock()->connect(shared_from_this());
}

void InteriorPoint::disconnect(std::weak_ptr<Face> face)
{
    // check input
    if (face.expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && faces_.empty()) storage_.lock()->delete_interior_point(shared_from_this());
}

void InteriorPoint::disconnect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && surfaces_.empty()) storage_.lock()->delete_interior_point(shared_from_this());
}

bool operator<(const std::weak_ptr<InteriorPoint>& lhs, const std::weak_ptr<InteriorPoint>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired InteriorPoints");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<InteriorPoint>& lhs, const std::weak_ptr<InteriorPoint>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired InteriorPoints");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}
