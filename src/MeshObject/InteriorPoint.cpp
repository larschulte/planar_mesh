#include "MeshObject/Storage.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "utilities/covariance_math.hpp"

#include "MeshObject/Vertex.hpp"

Settings InteriorPoint::settings_;

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius, double distance_travelled)
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
    distance_travelled_ = distance_travelled;
    direction_ = (position_ - origin_).normalized();
    radius_ = radius;

    num_deletes_ = 0;

    // log
    if (settings_.log.initialize) std::cout << "InteriorPoint " << id_ << " created.\n";
}

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled)
{
    std::shared_ptr<GenericPoint> generic_point = storage->add_generic_point(position, origin, distance_travelled);
    initialize_(storage, generic_point);
    storage->delete_generic_point(generic_point);
}

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point)
{
    initialize_(storage, generic_point->get_position(), generic_point->get_origin(), generic_point->get_radius(), generic_point->get_distance_travelled());
    num_deletes_ = generic_point->get_num_deletes();
    projected_uncertainty_ = generic_point->get_projected_uncertainty();
}

void InteriorPoint::delete_()
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_lifecycle_);

    // set deletion flag
    deleting_ = true;

    // log
    if (settings_.log.deletion) std::cout << "Destroying InteriorPoint " << id_ << std::endl;

    // subscribers and publishers
    {
        delete_subscribers();
    }

    // surface (disconnect)
    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_surface_);

        // disconnect from surface
        surface_->disconnect(shared_from_this());
    }

    // faces (disconnect)
    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

        // disconnect from faces
        for (const auto& face : faces_) face->disconnect(shared_from_this());

        // clear
        faces_.clear();
    }
    
    // create generic point
    {
        storage_->add_to_queue(shared_from_this());
    }

    // log
    if (settings_.log.deletion) std::cout << "---------- InteriorPoint " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& InteriorPoint::get_id() const
{
    return id_;
}

const Eigen::Vector3d& InteriorPoint::get_original_position() const
{
    return position_;
}

const Eigen::Vector3d& InteriorPoint::get_position() const
{
    if (projected_position_.isZero()) throw std::runtime_error("Interior projected position is not set.");

    return projected_position_;
    // return position_;
}

const Eigen::Vector3d& InteriorPoint::get_origin() const
{
    return origin_;
}

const double& InteriorPoint::get_distance_travelled() const
{
    return distance_travelled_;
}

const Eigen::Vector3d& InteriorPoint::get_direction() const
{
    return direction_;
}

std::shared_ptr<Surface> InteriorPoint::get_surface() const
{    
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_surface_);
    return surface_;
}

const double& InteriorPoint::get_radius() const
{
    return radius_;
}

Eigen::Vector3d InteriorPoint::compute_projected_position() 
{ 
    return get_surface()->compute_point_projective_position(get_origin(), get_original_position());
}

double InteriorPoint::compute_projected_distance() 
{
    return get_surface()->compute_point_projective_distance(get_origin(), get_original_position());
}

bool InteriorPoint::is_expired() const
{
    return is_expired_;
}

double& InteriorPoint::get_projected_uncertainty()
{
    return projected_uncertainty_;
}

std::size_t InteriorPoint::get_num_deletes() const
{
    return num_deletes_;
}

void InteriorPoint::connect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect interior point with invalid face.");

    {
        // wirte lock
        std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

        // skip if already exist
        if (faces_.find(face) != faces_.end()) return; // Already exists

        // connect
        faces_.insert(face);
    }

    // reverse connection
    face->connect(shared_from_this());
}

void InteriorPoint::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect interior point with invalid surface.");

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_surface_);

        // skip if already exist
        if (surface_ == surface) return; // Already exists

        // connect
        surface_ = surface;
    }

    // compute projected position
    projected_position_ = surface->compute_point_projective_position(get_origin(), get_original_position());
    
    // reverse connection
    surface->connect(shared_from_this());    
}

void InteriorPoint::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

        // skip if not exist
        if (faces_.find(face) == faces_.end()) return; // skip if not exist

        // disconnect
        faces_.erase(face);
    }

    // reverse disconnection
    face->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && faces_.empty()) storage_->add_interior_point_to_be_deleted(shared_from_this());
}

void InteriorPoint::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_surface_);

        // skip if not exist
        if (surface_ != surface) return; // skip if not exist

        // disconnect
        surface_ = nullptr;
    }

    // reset projected position
    projected_position_ = Eigen::Vector3d::Zero();

    // reverse disconnection
    surface->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && can_self_destruct_) storage_->add_interior_point_to_be_deleted(shared_from_this());
}

void InteriorPoint::delete_subscribers()
{
    // interior point subscribers
    std::vector<std::shared_ptr<Vertex>> interior_point_distance_subscribers_copy;
    {
        // lock
        std::shared_lock<std::shared_mutex> lock(rwlock_interior_point_distance_subscribers_);

        // copy
        interior_point_distance_subscribers_copy = interior_point_distance_subscribers_;
    }
    
    for (const auto& interior_point_subscriber : interior_point_distance_subscribers_copy)
    {
        // delete
        delete_interior_point_distance_subscriber(interior_point_subscriber);
    }
}

void InteriorPoint::add_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber)
{
    // check input
    if (interior_point_subscriber->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_subscribers_);

        // skip if already exist
        for (const auto& interior_point_subscriber_ : interior_point_distance_subscribers_) if (interior_point_subscriber_ == interior_point_subscriber) return; // Already exists

        // add subscriber
        interior_point_distance_subscribers_.push_back(interior_point_subscriber);
    }

    // add self to subscriber vertex as publisher
    interior_point_subscriber->add_interior_point_distance_publisher(shared_from_this());
}

void InteriorPoint::delete_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber)
{
    // check input
    if (interior_point_subscriber->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_subscribers_);

        // skip if not exist
        auto it = std::find(interior_point_distance_subscribers_.begin(), interior_point_distance_subscribers_.end(), interior_point_subscriber);
        if (it == interior_point_distance_subscribers_.end()) return; // skip if not exist

        // delete subscriber
        interior_point_distance_subscribers_.erase(it);
    }    
    
    // delete self from subscriber vertex as publisher
    interior_point_subscriber->delete_interior_point_distance_publisher(shared_from_this());
}


void InteriorPoint::set_reverse_radius_search_radius(double radius)
{
    radius_ = radius;
}

// swap surface1 with surface2
void InteriorPoint::swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2)
{
    // if contains surfacce1
    if (surface_ == surface1)
    {
        // std::cout << "Swapping interior point " << id_ << " surface " << surface1->get_id() << " with surface " << surface2->get_id() << std::endl;
        
        can_self_destruct_ = false;
        disconnect(surface1);
        connect(surface2);
        can_self_destruct_ = true;

        // cascade swap
        for (const std::shared_ptr<Face>& face : faces_)
        {
            face->swap(surface1, surface2);
        }
    }
}

bool operator<(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs)
{
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired InteriorPoints");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs)
{
    if (!lhs && !rhs) return true; // true if both are nullptr
    if (!lhs || !rhs) return false; // false if either is nullptr
    return lhs->get_id() == rhs->get_id();
}
