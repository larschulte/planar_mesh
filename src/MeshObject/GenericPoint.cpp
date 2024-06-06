#include "MeshObject/Storage.hpp"
#include "MeshObject/GenericPoint.hpp"

void GenericPoint::initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create generic point with invalid storage.");

    // get id
    id_ = storage_valid->get_next_generic_point_id();

    // store
    position_ = position;
    origin_ = origin;

    // log
    std::cout << "GenericPoint " << id_ << " created.\n";
}

void GenericPoint::delete_()
{
    // log
    std::cout << "Destroying GenericPoint " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // log
    std::cout << "---------- GenericPoint " << id_ << " destroyed" << std::endl;
}

int GenericPoint::get_id() const
{
    return id_;
}

Eigen::Vector3d GenericPoint::get_position() const
{
    return position_;
}

Eigen::Vector3d GenericPoint::get_origin() const
{
    return origin_;
}

bool operator<(const std::weak_ptr<GenericPoint>& lhs, const std::weak_ptr<GenericPoint>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired GenericPoints");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired GenericPoints");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}
