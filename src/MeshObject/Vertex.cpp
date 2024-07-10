#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include <iostream>

#include "MeshObject/GenericPoint.hpp"
#include <set>
#include "utilities/covariance_math.hpp"

Settings Vertex::settings_;

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
    direction_ = (position_ - origin_).normalized();

    // num of deletes
    num_deletes_ = 0;

    // set reverse search radius based on input parameter
    set_reverse_radius_search_radius(radius);

    // update boundary state
    update_boundary_state();

    // log
    std::cout << "Vertex " << id_ << " created.\n";
}

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point)
{
    initialize_(storage, generic_point->get_position(), generic_point->get_origin(), generic_point->get_radius());
    num_deletes_ = generic_point->get_num_deletes();
}

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
{
    std::shared_ptr<GenericPoint> generic_point = storage->add_generic_point(position, origin);
    initialize_(storage, generic_point);
    storage->delete_generic_point(generic_point);
}

void Vertex::delete_()
{
    // log
    std::cout << "Destroying vertex " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges = edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces = surfaces_;
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& face : faces) disconnect(face);
    for (const auto& surface : surfaces) disconnect(surface);

    // remove from search tree
    if (is_searchable_)
    {
        storage_->remove_searchable_vertex(shared_from_this());
        is_searchable_ = false;
    }
    
    // update delete count
    num_deletes_++;

    // only create penetrated point / generic point if sibling is empty
    if (sibling_vertices_.empty())
    {
        if (storage_->has_penetrating_point())
        {
            // compute radius from storage
            double radius = (storage_->get_penetrating_point() - get_position()).norm();
            if (radius < reverse_search_radius_) reverse_search_radius_ = radius;

            // add to storage as penetrated point
            storage_->add_penetrated_point(shared_from_this());
        }
        else
        {
            // add to storage as generic point
            storage_->add_generic_point(shared_from_this());
        }
    }

    // disconnect from sibling vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> sibling_vertices = sibling_vertices_;
    for (const auto& sibling_vertex : sibling_vertices) disconnect(sibling_vertex);

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

void Vertex::try_update_surface_projection(const std::shared_ptr<Surface> surface)
{
    // update if surface changes
    if (mean_used_ != surface->get_mean())
    {
        mean_used_ = surface->get_mean();
        projected_position_ = surface->compute_point_projective_position(get_origin(), get_position());
        projected_distance_ = surface->compute_point_projective_distance(get_origin(), get_position());
    }
}

void Vertex::try_update_surface_projection()
{
    try_update_surface_projection(get_surface());
}

const Eigen::Vector3d& Vertex::get_projected_position(const std::shared_ptr<Surface> surface)
{
    try_update_surface_projection(surface);
    return projected_position_;
}

const Eigen::Vector3d& Vertex::get_projected_position()
{
    return get_projected_position(get_surface());
}

const double& Vertex::get_projected_distance(const std::shared_ptr<Surface> surface)
{
    try_update_surface_projection(surface);
    return projected_distance_;
}

const double& Vertex::get_projected_distance()
{
    return get_projected_distance(get_surface());
}

const Eigen::Vector3d& Vertex::get_origin() const 
{ 
    return origin_; 
}

const Eigen::Vector3d& Vertex::get_direction() const 
{ 
    return direction_; 
}

const std::shared_ptr<Surface>& Vertex::get_surface() const
{    
    if (surfaces_.empty()) throw std::runtime_error("Vertex has no surface.");

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

const std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& Vertex::get_surfaces() const 
{ 
    // if more than one surface, throw error
    if (surfaces_.size() > 1) throw std::runtime_error("Vertex connected to more than one surface.");

    return surfaces_; 
}

bool Vertex::has_surface() const
{
    return !surfaces_.empty();
}

const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Vertex::get_edges() const 
{ 
    return edges_; 
}

const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& Vertex::get_sibling_vertices() const 
{ 
    return sibling_vertices_; 
}

std::size_t Vertex::get_num_deletes() const
{
    return num_deletes_;
}

void Vertex::try_merge_surfaces()
{
    // check if there is only one surface
    if (surfaces_.size() == 1) return;

    bool merge_happened = false;

    // merge surfaces
    while (true)
    {
        bool merge_again = false;

        // generate pairs of surfaces
        std::set<std::pair<std::shared_ptr<Surface>, std::shared_ptr<Surface>>> surface_pairs;
        for (std::shared_ptr<Surface> surface1 : surfaces_) 
        {
            for (std::shared_ptr<Surface> surface2 : surfaces_) 
            {
                if (surface1 >= surface2) continue;
                surface_pairs.insert(std::make_pair(surface1, surface2));
            }
        }

        // merge surfaces
        for (const auto& pairs : surface_pairs) 
        {
            // skip if combined surface have large eigenvalue
            const std::shared_ptr<Surface>& surface1 = pairs.first;
            const std::shared_ptr<Surface>& surface2 = pairs.second;
            const Eigen::Matrix3d& cov1 = surface1->get_covariance();
            const Eigen::Matrix3d& cov2 = surface2->get_covariance();
            const Eigen::Vector3d& mean1 = surface1->get_mean();
            const Eigen::Vector3d& mean2 = surface2->get_mean();
            int size1 = surface1->get_total_point_size();
            int size2 = surface2->get_total_point_size();
            Eigen::Matrix3d covariance_matrix = merge_covariance(cov1, cov2, mean1, mean2, size1, size2);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance_matrix);
            
            double eigenvalue = solver.eigenvalues()[0];
            if (eigenvalue > Vertex::settings_.merged_eigenvalue_threshold) continue;

            // merge by changing surface1 into surface2
            std::cout << ">> Merging surface " << surface1->get_id() << " with " << surface1->get_total_point_size() << " points into surface " << surface2->get_id() << " with " << surface2->get_total_point_size() << " points." << std::endl;
            swap(surface1, surface2);
            merge_happened = true;
            std::cout << ">> resultant surface " << surface1->get_id() << " has " << surface1->get_total_point_size() << " points." << std::endl;
            std::cout << ">> resultant surface " << surface2->get_id() << " has " << surface2->get_total_point_size() << " points." << std::endl;

            // reset
            merge_again = true;
            break;
        }

        if (!merge_again) break;
    }

    if (merge_happened) 
    {
        std::cout << ">> Merging surfaces done." << std::endl;
    }
    else
    {
        std::cout << ">> No merging happened." << std::endl;
    }
}

const Eigen::Vector2d& Vertex::get_surface_coordinate(const std::shared_ptr<Surface> surface)
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
        Eigen::Vector3d projected_position = get_projected_position(surface);
        surface_coordinate_ = (projection_matrix.transpose() * projected_position).head<2>();
        eigenvectors_used_ = eigenvectors;
        return surface_coordinate_;
    }
}

const Eigen::Vector2d& Vertex::get_surface_coordinate()
{
    return get_surface_coordinate(get_surface());
}

bool Vertex::is_expired() const
{
    return is_expired_;
}

bool Vertex::is_boundary(const std::shared_ptr<Surface>& surface) const
{
    return is_boundary_map_.at(surface);
}

bool Vertex::is_boundary() const
{
    for (const auto& pair : is_boundary_map_)
    {
        if (pair.second) return true;
    }
    return false;
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

    // update confirmed status
    if (inserted) update_confirmed_status();
    if (inserted) update_singular_state();
}

void Vertex::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface->connect(shared_from_this());
    if (inserted) is_boundary_map_[surface] = false;
    if (inserted) update_boundary_state(surface);
    if (inserted) is_singular_map_[surface] = true;
    if (inserted) is_matched_surface_map_[surface] = false;
}

void Vertex::connect(const std::shared_ptr<Vertex>& sibling_vertex)
{
    // check input
    if (sibling_vertex->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid sibling vertex.");

    // skip if try to connect to itself
    if (sibling_vertex == shared_from_this()) return;

    // connect
    bool inserted = sibling_vertices_.insert(sibling_vertex).second;
    // if (inserted) std::cout << "Connected vertex " << id_ << " with vertex " << sibling_vertex->get_id() << " as sibling."<< std::endl;
    if (inserted) sibling_vertex->connect(shared_from_this());
    if (inserted)
    {
        for (const std::shared_ptr<Vertex>& sibling_vertex_ : sibling_vertices_)
        {
            sibling_vertex_->connect(sibling_vertex);
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

    // // check self destruct
    // if (!deleting_ && edges_.empty()) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Face>& face)
{
    // check pointer validity
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());

    // update confirmed status
    if (erased) update_confirmed_status();
    if (erased) update_singular_state();

    // // check self destruct
    // if (!deleting_ && faces_.empty()) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surfaces_.erase(surface);
    if (erased) surface->disconnect(shared_from_this());
    if (erased) is_boundary_map_.erase(surface);
    if (erased) is_singular_map_.erase(surface);
    if (erased) is_matched_surface_map_.erase(surface);

    // check self destruct
    if (!deleting_ && surfaces_.empty()) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Vertex>& sibling_vertex)
{
    // check input
    if (sibling_vertex->is_expired()) return;

    // disconnect
    bool erased = sibling_vertices_.erase(sibling_vertex);
    if (erased) sibling_vertex->disconnect(shared_from_this());
}

void Vertex::add_matched_surface(const std::shared_ptr<Surface>& surface)
{
    is_matched_surface_map_.at(surface) = true;
}

bool Vertex::is_matched_surface(const std::shared_ptr<Surface>& surface) const
{
    return is_matched_surface_map_.at(surface);
}

void Vertex::update_confirmed_status()
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

void Vertex::update_singular_state(const std::shared_ptr<Surface>& surface)
{
    // count number of faces in this surface
    int num_faces_in_surface = 0;
    for (const std::shared_ptr<Face>& face : faces_)
    {
        if (face->get_surfaces().find(surface) != face->get_surfaces().end()) num_faces_in_surface++;
    }

    // update singular state
    if (num_faces_in_surface == 0) is_singular_map_.at(surface) = true;
    else is_singular_map_.at(surface) = false;
}

void Vertex::update_singular_state()
{
    for (const std::shared_ptr<Surface>& surface : surfaces_)
    {
        update_singular_state(surface);
    }
}

bool Vertex::is_confirmed() const
{
    return is_confirmed_;
}

bool Vertex::is_singular(const std::shared_ptr<Surface>& surface) const
{
    return is_singular_map_.at(surface);
}

bool Vertex::is_singular() const
{
    // singular if all singular
    for (const auto& pair : is_singular_map_)
    {
        if (!pair.second) return false;
    }
    return true;
}

// swap surface1 with surface2
void Vertex::swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2)
{
    // if contains surfacce1
    bool contains_surface1 = surfaces_.find(surface1) != surfaces_.end();
    
    if (contains_surface1)
    {
        connect(surface2);
        disconnect(surface1);

        // cascade swap
        for (const std::shared_ptr<Edge>& edge : edges_)
        {
            edge->swap(surface1, surface2);
        }
        for (const std::shared_ptr<Face>& face : faces_)
        {
            face->swap(surface1, surface2);
        }
    }
}

void Vertex::update_boundary_state(const std::shared_ptr<Surface>& surface)
{
    if (deleting_) return;

    // skip if surface not yet in surfaces_ list
    if (surfaces_.find(surface) == surfaces_.end()) return;

    // becomes boundary when one of the connected edges is boundary, or when the point is alone
    is_boundary_map_.at(surface) = false;
    int num_edge_in_same_surface = 0;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        // skip if edge is not in the same surface
        if (edge->get_surfaces().find(surface) == edge->get_surfaces().end()) continue;

        num_edge_in_same_surface++;

        if (edge->is_boundary(surface))
        {
            is_boundary_map_.at(surface) = true;
            break;
        }
    }
    if (num_edge_in_same_surface == 0)
    {
        is_boundary_map_.at(surface) = true;
    }

    // update searchable state
    update_searchable_state();
}

void Vertex::update_boundary_state()
{
    for (const std::shared_ptr<Surface>& surface : surfaces_)
    {
        update_boundary_state(surface);
    }
}

void Vertex::update_searchable_state()
{
    // check if is boundary
    if (is_boundary() && !is_searchable_)
    {
        storage_->add_searchable_vertex(shared_from_this());
        is_searchable_ = true;
    }
    else if (!is_boundary() && is_searchable_)
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

void Vertex::reduce_reverse_radius_search_radius(double radius)
{
    if (radius < reverse_search_radius_) set_reverse_radius_search_radius(radius);
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
    if (!lhs && !rhs) return true; // true if both are nullptr
    if (!lhs || !rhs) return false; // false if either is nullptr
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs->get_id() == rhs->get_id();
}