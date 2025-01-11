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
    // log
    if (settings_.log.deletion) std::cout << "Destroying InteriorPoint " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    delete_subscribers();
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    for (const auto& face : faces) disconnect(face);
    if (surface_) disconnect(surface_);

    // only create penetrated point / generic point if sibling is empty
    if (sibling_interior_points_.empty())
    {
        storage_->add_to_queue(shared_from_this());
    }
    
    // disconnect from sibling interior points
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> sibling_interior_points = sibling_interior_points_;
    for (const auto& sibling_interior_point : sibling_interior_points) disconnect(sibling_interior_point);

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

const std::shared_ptr<Surface>& InteriorPoint::get_surface() const
{    
    return surface_;
}

const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& InteriorPoint::get_sibling_interior_points() const
{
    return sibling_interior_points_;
}

const double& InteriorPoint::get_radius() const
{
    return radius_;
}

const Eigen::Vector3d& InteriorPoint::buffer_compute_projected_position(const std::shared_ptr<Surface> surface)
{
    // do cartersian rounding now, swtich to Locality Sensitive Hashing later

    // compute hash
    std::size_t hash = surface->get_approximate_normal_hash();

    // add to cache if not exist
    if (!buffer_projected_position_.exists(hash)) 
    {
        const Eigen::Vector3d computedResult = surface->compute_point_projective_position(get_origin(), get_original_position());
        buffer_projected_position_.put(hash, computedResult);
    }

    // return
    return buffer_projected_position_.get(hash);
}

const double& InteriorPoint::buffer_compute_projected_distance(const std::shared_ptr<Surface> surface)
{
    // do cartersian rounding now, swtich to Locality Sensitive Hashing later

    // compute hash
    std::size_t hash = surface->get_approximate_normal_hash();

    // add to cache if not exist
    if (!buffer_projected_distance_.exists(hash)) 
    {
        const double computedResult = surface->compute_point_projective_distance(get_origin(), get_original_position());
        buffer_projected_distance_.put(hash, computedResult);
    }

    // return
    return buffer_projected_distance_.get(hash);
}

const Eigen::Vector3d& InteriorPoint::buffer_compute_projected_position() { return buffer_compute_projected_position(get_surface()); }
const double& InteriorPoint::buffer_compute_projected_distance() { return buffer_compute_projected_distance(get_surface()); }

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

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face->connect(shared_from_this());

    // update confirmed status
    if (inserted) update_confirmed_status();
}

void InteriorPoint::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect interior point with invalid surface.");

    // connect
    bool inserted = surface_ != surface;
    if (inserted) surface_ = surface;
    if (inserted) 
    {
        // update projected position
        projected_position_ = surface->compute_point_projective_position(get_origin(), get_original_position());
    }
    if (inserted) surface->connect(shared_from_this());    
}

void InteriorPoint::connect(const std::shared_ptr<InteriorPoint>& sibling_interior_point)
{
    // check input
    if (sibling_interior_point->is_expired()) throw std::runtime_error("Attempts to connect interior point with invalid sibling interior point.");

    // skip if try to connect to itself
    if (sibling_interior_point == shared_from_this()) return;

    // connect
    bool inserted = sibling_interior_points_.insert(sibling_interior_point).second;
    if (inserted) std::cout << "Connected interior point " << id_ << " with interiror point " << sibling_interior_point->get_id() << " as sibling."<< std::endl;
    if (inserted) sibling_interior_point->connect(shared_from_this());
    if (inserted)
    {
        for (const std::shared_ptr<InteriorPoint>& sibling_interior_point_ : sibling_interior_points_)
        {
            sibling_interior_point_->connect(sibling_interior_point);
        }
    }
}

void InteriorPoint::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());

    // udpate confirmed status
    if (erased) update_confirmed_status();

    // self destruct
    if (!deleting_ && faces_.empty()) storage_->delete_interior_point(shared_from_this());
}

void InteriorPoint::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surface_ == surface;
    if (erased) surface->disconnect(shared_from_this());
    if (erased) surface_ = nullptr;
    if (erased) projected_position_ = Eigen::Vector3d::Zero();

    // self destruct
    if (!deleting_ && erased && can_self_destruct_) storage_->delete_interior_point(shared_from_this());
}

void InteriorPoint::disconnect(const std::shared_ptr<InteriorPoint>& sibling_interior_point)
{
    // check input
    if (sibling_interior_point->is_expired()) return;

    // disconnect
    bool erased = sibling_interior_points_.erase(sibling_interior_point);
    if (erased) sibling_interior_point->disconnect(shared_from_this());
}

void InteriorPoint::delete_subscribers()
{
    // interior point subscribers
    std::vector<std::shared_ptr<Vertex>> interior_point_distance_subscribers_copy = interior_point_distance_subscribers_;
    for (const auto& interior_point_subscriber : interior_point_distance_subscribers_copy)
    {
        // try lock vertex and surface
        if(!omp_test_nest_lock(&interior_point_subscriber->vertex_lock)) continue;
        if (!omp_test_nest_lock(&interior_point_subscriber->get_surface()->lock)) 
        {
            // unlock vertex
            omp_unset_nest_lock(&interior_point_subscriber->vertex_lock);
            continue;
        }

        // delete
        delete_interior_point_distance_subscriber(interior_point_subscriber);

        // unlock vertex and surface
        omp_unset_nest_lock(&interior_point_subscriber->get_surface()->lock);
        omp_unset_nest_lock(&interior_point_subscriber->vertex_lock);
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
    interior_point_subscriber->upon_adding_publisher();
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

void InteriorPoint::update_confirmed_status()
{
    // update number of confirmed faces
    num_confirmed_faces = 0;
    for (const std::shared_ptr<Face>& face : faces_)
    {
        if (face->is_confirmed()) num_confirmed_faces++;
    }

    // update confirmed status
    if (num_confirmed_faces >= 1) is_confirmed_ = true;
    else is_confirmed_ = false;
}

bool InteriorPoint::is_confirmed() const
{
    return is_confirmed_;
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
