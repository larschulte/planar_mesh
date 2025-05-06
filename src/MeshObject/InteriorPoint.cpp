#include "MeshObject/Storage.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "utilities/covariance_math.hpp"

#include "MeshObject/Vertex.hpp"

Settings InteriorPoint::settings_;

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<Face>& face, const std::shared_ptr<GenericPoint>& generic_point)
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
    position_ = generic_point->get_position();
    origin_ = generic_point->get_origin();
    distance_travelled_ = generic_point->get_distance_travelled();
    direction_ = (position_ - origin_).normalized();
    radius_ = generic_point->get_radius();
    num_deletes_ = generic_point->get_num_deletes();
    projected_uncertainty_ = generic_point->get_projected_uncertainty();

    // connect to surface
    {
        // projected position
        projected_position_ = surface->compute_point_projective_position(origin_, position_);

        // connect to surface
        surface_ = surface;
        surface_.lock()->connect(shared_from_this());
    }

    // connect to face
    if (face)
    {
        // connect to face
        face_ = face;
        face_.lock()->connect(shared_from_this());
    }

    // log
    if (settings_.log.initialize) std::cout << "InteriorPoint " << id_ << " created.\n";
}

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex, const std::shared_ptr<GenericPoint>& generic_point)
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
    position_ = generic_point->get_position();
    origin_ = generic_point->get_origin();
    distance_travelled_ = generic_point->get_distance_travelled();
    direction_ = (position_ - origin_).normalized();
    radius_ = generic_point->get_radius();
    num_deletes_ = generic_point->get_num_deletes();
    projected_uncertainty_ = generic_point->get_projected_uncertainty();

    // connect to surface
    {
        // projected position
        projected_position_ = surface->compute_point_projective_position(origin_, position_);

        // connect to surface
        surface_ = surface;
        surface_.lock()->connect(shared_from_this());
    }

    // connect to face
    if (vertex)
    {
        // connect to face
        vertex_ = vertex;
        vertex_.lock()->connect(shared_from_this());
    }

    // log
    if (settings_.log.initialize) std::cout << "InteriorPoint " << id_ << " created.\n";
}

void InteriorPoint::delete_()
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_lifecycle_);

    // skip if already deleted
    if (is_expired_) return;

    // set deletion flag
    deleting_ = true;

    // log
    if (settings_.log.deletion) std::cout << "Destroying InteriorPoint " << id_ << std::endl;

    // subscribers (disconnect)
    {
        // lock, copy and clear
        std::unordered_set<std::weak_ptr<Vertex>, MeshObjectHash> interior_point_distance_subscribers_copy;
        {
            std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_subscribers_);
            interior_point_distance_subscribers_copy = interior_point_distance_subscribers_;
            interior_point_distance_subscribers_.clear();
        }

        // disconnect
        for (const auto& interior_point_subscriber : interior_point_distance_subscribers_copy) 
        {
            interior_point_subscriber.lock()->delete_interior_point_distance_publisher(shared_from_this());
        }

        // add vertex to update radius
        for (const auto& interior_point_subscriber : interior_point_distance_subscribers_copy)
        {
            storage_->add_vertex_that_have_deleted_publishers(interior_point_subscriber.lock());
        } 
    }

    // surface (disconnect)
    {
        // disconnect from surface
        surface_.lock()->disconnect(shared_from_this());

        // clear
        surface_.reset();
    }

    // faces (disconnect)
    if (face_.lock())
    {
        // disconnect from faces
        face_.lock()->disconnect(shared_from_this());

        // clear
        face_.reset();
    }

    // vertex (disconnect)
    if (vertex_.lock())
    {
        // disconnect from vertices
        vertex_.lock()->disconnect(shared_from_this());

        // clear
        vertex_.reset();
    }
    
    // create generic point
    {
        num_deletes_++;

        if (do_not_add_back_due_to_seed_surface_)
        {
            // do nothing
        }
        else
        {
            storage_->add_to_queue(shared_from_this());
        }
    }

    // log
    if (settings_.log.deletion) std::cout << "---------- InteriorPoint " << id_ << " destroyed" << std::endl;

    storage_.reset();
    
    // set expired
    is_expired_ = true;
}

InteriorPoint::~InteriorPoint()
{
    // std::cout << "InteriorPoint " << id_ << " deleted." << std::endl;
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
    return surface_.lock();
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

const Eigen::Vector2d& InteriorPoint::get_surface_coordinate(const std::shared_ptr<Surface>& surface)
{
    const Eigen::Matrix3d& eigenvectors = surface->get_eigenvectors();
    if (eigenvectors_used_ == eigenvectors)
    {
        // use stored coordinate if eigenvectors are the same
        return surface_coordinate_;
    }
    else
    {
        // compute new coordinate
        Eigen::Matrix<double, 3, 2> projection_matrix = eigenvectors.rightCols<2>();
        Eigen::Vector3d projected_position = surface->compute_point_projective_position(get_origin(), get_original_position());
        surface_coordinate_ = (projection_matrix.transpose() * projected_position).head<2>();
        eigenvectors_used_ = eigenvectors;
        return surface_coordinate_;
    }
}

const Eigen::Vector2d& InteriorPoint::get_surface_coordinate()
{
    return get_surface_coordinate(get_surface());
}


std::size_t InteriorPoint::get_num_deletes() const
{
    return num_deletes_;
}

void InteriorPoint::add_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber)
{
    // check input
    if (interior_point_subscriber->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_subscribers_);

        // skip if already exist
        if (interior_point_distance_subscribers_.find(interior_point_subscriber) != interior_point_distance_subscribers_.end()) return; // Already exists

        // add subscriber
        interior_point_distance_subscribers_.insert(interior_point_subscriber);
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
        auto it = interior_point_distance_subscribers_.find(interior_point_subscriber);
        if (it == interior_point_distance_subscribers_.end()) return; // skip if not exist

        // delete subscriber
        interior_point_distance_subscribers_.erase(it);
    }
}


void InteriorPoint::set_reverse_radius_search_radius(double radius)
{
    radius_ = radius;
}

void InteriorPoint::set_do_not_add_back_due_to_seed_surface(bool do_not_add_back_due_to_seed_surface)
{
    do_not_add_back_due_to_seed_surface_ = do_not_add_back_due_to_seed_surface;
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
