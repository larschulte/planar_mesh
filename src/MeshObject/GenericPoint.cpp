#include "MeshObject/Storage.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/InteriorPoint.hpp"

Settings GenericPoint::settings_;

void GenericPoint::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
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
    direction_ = (position - origin).normalized();
    if (GenericPoint::settings_.use_radius_value)
    {
        radius_ = GenericPoint::settings_.radius_value;
    }
    else
    {
        radius_ = (position - origin).norm() * GenericPoint::settings_.radius_ratio;  
    }

    num_deletes_ = 0;

    // log
    if (settings_.log.initialize) std::cout << "GenericPoint " << id_ << " created.\n";
}

void GenericPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Vertex>& vertex)
{
    initialize_(storage, vertex->get_position(), vertex->get_origin());
    radius_ = vertex->get_radius();

    num_deletes_ = vertex->get_num_deletes();
}

void GenericPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<InteriorPoint>& interior_point)
{
    initialize_(storage, interior_point->get_position(), interior_point->get_origin());
    radius_ = interior_point->get_radius();

    num_deletes_ = interior_point->get_num_deletes();
}

void GenericPoint::delete_()
{
    // log
    if (settings_.log.deletion) std::cout << "Destroying GenericPoint " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // log
    if (settings_.log.deletion) std::cout << "---------- GenericPoint " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& GenericPoint::get_id() const
{
    return id_;
}

const Eigen::Vector3d& GenericPoint::get_position() const
{
    return position_;
}

const Eigen::Vector3d& GenericPoint::get_origin() const
{
    return origin_;
}

const Eigen::Vector3d& GenericPoint::get_direction() const
{
    return direction_;
}

const double& GenericPoint::get_radius() const
{
    return radius_;
}

bool GenericPoint::is_expired() const
{
    return is_expired_;
}

std::size_t GenericPoint::get_num_deletes() const
{
    return num_deletes_;
}

bool operator<(const std::shared_ptr<GenericPoint>& lhs, const std::shared_ptr<GenericPoint>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired GenericPoints");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<GenericPoint>& lhs, const std::shared_ptr<GenericPoint>& rhs)
{
    if (!lhs && !rhs) return true; // true if both are nullptr
    if (!lhs || !rhs) return false; // false if either is nullptr
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired GenericPoints");
    return lhs->get_id() == rhs->get_id();
}
