#include "MeshObject/Storage.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "utilities/covariance_math.hpp"

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius)
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
    radius_ = radius;

    num_deletes_ = 0;

    // log
    std::cout << "InteriorPoint " << id_ << " created.\n";
}

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
{
    std::shared_ptr<GenericPoint> generic_point = storage->add_generic_point(position, origin);
    initialize_(storage, generic_point);
    storage->delete_generic_point(generic_point);
}

void InteriorPoint::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point)
{
    initialize_(storage, generic_point->get_position(), generic_point->get_origin(), generic_point->get_radius());
    num_deletes_ = generic_point->get_num_deletes();
}

void InteriorPoint::delete_()
{
    // log
    std::cout << "Destroying InteriorPoint " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces = surfaces_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> sibling_interior_points = sibling_interior_points_;
    for (const auto& face : faces) disconnect(face);
    for (const auto& surface : surfaces) disconnect(surface);
    for (const auto& sibling_interior_point : sibling_interior_points) disconnect(sibling_interior_point);

    // compute radius
    if (storage_->has_penetrating_point())
    {
        // compute radius from storage
        double radius = (storage_->get_penetrating_point() - get_position()).norm();
        if (radius < radius_) radius_ = radius;

        // add to storage as penetrated point
        storage_->add_penetrated_point(shared_from_this());
    }
    else
    {
        // add to storage as generic point
        storage_->add_generic_point(shared_from_this());
    }

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
    if (surfaces_.empty()) throw std::runtime_error("InteriorPoint has no surface.");

    // Select the surface with the lowest projective std, return as reference
    double min_std = std::numeric_limits<double>::max();
    const std::shared_ptr<Surface>* selected_surface = nullptr;
    for (const std::shared_ptr<Surface>& surface : surfaces_) 
    {
        // get stats
        const std::vector<double>& stats = surface->get_projective_distance_stats();
        double std = compute_std(stats);
        if (std < min_std) 
        {
            min_std = std;
            selected_surface = &surface;
        }
    }

    return *selected_surface;
}


const std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& InteriorPoint::get_surfaces() const
{
    // if more than one surface, throw error
    if (surfaces_.size() > 1) throw std::runtime_error("Interior point connected to more than one surface.");

    return surfaces_;
}

const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& InteriorPoint::get_sibling_interior_points() const
{
    return sibling_interior_points_;
}

const double& InteriorPoint::get_radius() const
{
    return radius_;
}

void InteriorPoint::try_update_surface_projection(const std::shared_ptr<Surface> surface)
{
    // update if surface changes
    if (normal_used_ != surface->get_normal())
    {
        normal_used_ = surface->get_normal();
        projected_position_ = surface->compute_point_projective_position(get_origin(), get_position());
        projected_distance_ = surface->compute_point_projective_distance(get_origin(), get_position());
    }
}

void InteriorPoint::try_update_surface_projection()
{
    try_update_surface_projection(get_surface());
}

const Eigen::Vector3d& InteriorPoint::get_projected_position(const std::shared_ptr<Surface> surface)
{
    try_update_surface_projection(surface);
    return projected_position_;
}

const Eigen::Vector3d& InteriorPoint::get_projected_position()
{
    return get_projected_position(get_surface());
}

const double& InteriorPoint::get_projected_distance(const std::shared_ptr<Surface> surface)
{
    try_update_surface_projection(surface);
    return projected_distance_;
}

const double& InteriorPoint::get_projected_distance()
{
    return get_projected_distance(get_surface());
}

bool InteriorPoint::is_expired() const
{
    return is_expired_;
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
    bool inserted = surfaces_.insert(surface).second;
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
    bool erased = surfaces_.erase(surface);
    if (erased) surface->disconnect(shared_from_this());

    // self destruct
    if (!deleting_ && surfaces_.empty()) storage_->delete_interior_point(shared_from_this());
}

void InteriorPoint::disconnect(const std::shared_ptr<InteriorPoint>& sibling_interior_point)
{
    // check input
    if (sibling_interior_point->is_expired()) return;

    // disconnect
    bool erased = sibling_interior_points_.erase(sibling_interior_point);
    if (erased) sibling_interior_point->disconnect(shared_from_this());
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
    if (surfaces_.find(surface1) != surfaces_.end())
    {
        connect(surface2);
        disconnect(surface1);

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
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired InteriorPoints");
    return lhs->get_id() == rhs->get_id();
}
