#include "MeshObject/Storage.hpp"
#include "MeshObject/GenericPoint.hpp"

void GenericPoint::initialize_(std::shared_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin)
{
    // set expired
    is_expired_ = false;

    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create generic point with invalid storage.");

    // store
    storage_ = storage;
    
    // get id
    id_ = storage_->get_next_generic_point_id();

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

    // set expired
    is_expired_ = true;
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

bool GenericPoint::is_expired() const
{
    return is_expired_;
}

bool operator<(const std::shared_ptr<GenericPoint>& lhs, const std::shared_ptr<GenericPoint>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired GenericPoints");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<GenericPoint>& lhs, const std::shared_ptr<GenericPoint>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired GenericPoints");
    return lhs->get_id() == rhs->get_id();
}
