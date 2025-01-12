#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Storage.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>
#include "MeshObject/InteriorPoint.hpp"
#include "utilities/covariance_math.hpp"
#include "MeshObject/UnionFind.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "utilities/utilities.hpp"

#include "utilities/gtsam_plane.hpp"

Settings Surface::settings_;

Surface::Surface()
{
    omp_init_nest_lock(&lock);
}

Surface::~Surface()
{
    omp_destroy_nest_lock(&lock);
}

void Surface::initialize_(const std::shared_ptr<Storage>& storage)
{
    // set expired
    is_expired_ = false;
    
    // check pointer validity
    if (storage->is_expired()) throw std::runtime_error("Attempts to create surface with invalid storage.");
    auto storage_valid = storage;

    // get id
    id_ = storage_valid->get_next_surface_id();

    // store
    storage_ = storage_valid;

    // initialize surface fitting
    mean_ = Eigen::Vector3d::Zero();
    covariance_ = Eigen::Matrix3d::Zero();
    eigenvectors_ = Eigen::Matrix3d::Identity();
    eigenvalues_ = Eigen::Vector3d::Zero();
    normal_ = Eigen::Vector3d(0, 0, 1);

    // initialize surface color
    set_random_color();

    // initialize edge_bvh
    edge_bvh_.set_surface(shared_from_this());
    
    // log
    if (settings_.log.initialize) std::cout << "Surface " << id_ << " created.\n";
}

void Surface::delete_()
{
    // log
    if (settings_.log.deletion) std::cout << "Destroying surface " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // make copies first since disconnect will modify the set
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices = vertices_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges = edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points = interior_points_;
    for (const auto& vertex : vertices) disconnect(vertex);
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& face : faces) disconnect(face);
    for (const auto& interior_point : interior_points) disconnect(interior_point);

    // log
    if (settings_.log.deletion) std::cout << "---------- surface " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& Surface::get_id() const
{
    return id_;
}

double Surface::compute_point_to_plane_distance(const Eigen::Vector3d& point) const
{
    if (get_total_point_size() < settings_.fit_plane_threshold) throw std::runtime_error("Surface is seed surface.");

    return (point - mean_).dot(normal_);
}

double Surface::compute_point_projective_distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& position) const
{
    // if perpendicular, return NaN
    Eigen::Vector3d rayDirection = (position - origin).normalized();
    double distance = std::abs( (mean_ - position).dot(normal_) / rayDirection.dot(normal_) );

    // negative distance if origin and position are on the opposite side of the plane
    const bool opposite_side = (origin - mean_).dot(normal_) * (position - mean_).dot(normal_) < 0;
    if (opposite_side) distance *= -1;

    // return
    return distance;
}

double Surface::compute_point_projective_distance(const std::shared_ptr<GenericPoint>& generic_point) const
{
    return compute_point_projective_distance(generic_point->get_origin(), generic_point->get_position());
}

double Surface::compute_point_projective_distance(const std::shared_ptr<Vertex>& vertex) const
{
    return compute_point_projective_distance(vertex->get_origin(), vertex->get_position());
}

double Surface::compute_point_projective_distance_with_improved_covariance(const Eigen::Vector3d& origin, const Eigen::Vector3d& position) const
{
    // set
    int size1 = get_total_point_size()-1;
    Eigen::Vector3d mean1 = mean_;
    Eigen::Matrix3d cov1 = covariance_;

    // point
    int size2 = 1;
    Eigen::Vector3d mean2 = position;
    Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();

    // set + point
    Eigen::Vector3d new_mean = merge_mean(mean1, mean2, size1, size2);
    Eigen::Matrix3d new_cov = merge_covariance(cov1, cov2, mean1, mean2, size1, size2);

    // plane estimate
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(new_cov);
    Eigen::Matrix3d new_eigenvectors = solver.eigenvectors();
    Eigen::Vector3d new_normal = new_eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal
    Eigen::Vector3d vector_towards_origin = origin - position;
    if (new_normal.dot(vector_towards_origin) < 0) new_normal *= -1; // normal should points towards the origin

    // compute distance
    Eigen::Vector3d rayDirection = (position - origin).normalized();
    double distance = (new_mean - position).dot(new_normal) / rayDirection.dot(new_normal);

    // return
    return distance;
}

Eigen::Vector3d Surface::compute_point_projective_position(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const
{
    // read lock
    std::shared_lock lock(rwlock_surface_fitting_);

    // if surface is seed surface, return original point as projected point
    if (get_total_point_size() < settings_.fit_plane_threshold) return point;

    // compute
    Eigen::Vector3d rayDirection = (point - origin).normalized();
    double distance = (mean_ - point).dot(normal_) / rayDirection.dot(normal_);
    Eigen::Vector3d intersection = point + distance * rayDirection;

    // return
    return intersection;
}

RelativePosition Surface::check_relative_position(double distance_travelled, const Eigen::Vector3d& origin, const Eigen::Vector3d& point, const Eigen::Vector3d& direction, double& projected_uncertainty)
{
    // lock surface fitting
    std::shared_lock lock(rwlock_surface_fitting_); // read lock

    // compute point to plane projective distance std
    const bool use_improved_covariance = true;
    if (use_improved_covariance)
    {
        // return no_relative_position if not enough points
        if (get_total_point_size() < settings_.fit_plane_threshold) return RelativePosition::NO_RELATIVE_POSITION;

        // compute d(range)/d(...)
        Eigen::Vector3d d_range_d_origin = - normal_ / normal_.dot(direction);
        Eigen::Vector3d d_range_d_mean = normal_ / normal_.dot(direction);
        Eigen::Vector3d d_range_d_direction = (mean_ - origin).dot(normal_) / std::pow(normal_.dot(direction), 2) * normal_;
        Eigen::Vector3d d_range_d_normal = ( (mean_ - origin) * normal_.dot(direction) - (mean_ - origin).dot(normal_) * direction ) / std::pow(normal_.dot(direction), 2);

        // compute covariance of ...
        // - settings
        double odometry_position_uncertainty_rate = settings_.odometry_position_uncertainty_rate;
        double odometry_angular_uncertainty_rate = settings_.odometry_angular_uncertainty_rate;
        double epsilon = 1e-6;
        // - uncertainties
        double odometry_position_uncertainty = std::abs(distance_travelled - get_smallest_distance_travelled()) * odometry_position_uncertainty_rate;
        double odometry_angular_uncertainty = std::abs(distance_travelled - get_smallest_distance_travelled()) * odometry_angular_uncertainty_rate / 180.0 * M_PI;
        double normal_angular_uncertainty = normal_uncertainty_;
        double plane_position_uncertainty = get_surface_position_std_in_normal_direction();
        // - computation
        Eigen::Matrix3d cov_mean = Eigen::Matrix3d::Identity() * std::pow(plane_position_uncertainty, 2);
        Eigen::Matrix3d cov_origin = Eigen::Matrix3d::Identity() * std::pow(odometry_position_uncertainty, 2);
        Eigen::Matrix3d cov_normal = generate_unit_vector_covariance(normal_, normal_angular_uncertainty, epsilon);
        Eigen::Matrix3d cov_direction = generate_unit_vector_covariance(direction, odometry_angular_uncertainty, epsilon);

        // compute variance of range due to ...
        double variance_range_origin = d_range_d_origin.transpose() * cov_origin * d_range_d_origin;
        double variance_range_mean = d_range_d_mean.transpose() * cov_mean * d_range_d_mean;
        double variance_range_direction = d_range_d_direction.transpose() * cov_direction * d_range_d_direction;
        double variance_range_normal = d_range_d_normal.transpose() * cov_normal * d_range_d_normal;

        // // compute std of range due to ...
        // double std_range_origin = std::sqrt(variance_range_origin);
        // double std_range_mean = std::sqrt(variance_range_mean);
        // double std_range_direction = std::sqrt(variance_range_direction);
        // double std_range_normal = std::sqrt(variance_range_normal);

        // compute combined std
        double combined_std = std::sqrt(variance_range_origin + variance_range_mean + variance_range_direction + variance_range_normal + settings_.range_precision * settings_.range_precision);
        double partial_combined_std = std::sqrt(variance_range_origin + variance_range_direction + settings_.range_precision * settings_.range_precision);

        // compute point to plane projective distance
        double projective_distance = compute_point_projective_distance(origin, point);

        // multiplier for confidence interval
        // double multiplier = 2.576; // for 99% confidence interval
        double multiplier = 1.96; // for 95% confidence interval

        // confidence interval values
        double threshold_in_front = multiplier * combined_std;
        double threshold_behind = - multiplier * combined_std;

        // check
        bool points_in_front_of_surface = projective_distance > threshold_in_front;
        bool points_behind_surface = projective_distance < threshold_behind;
        bool points_within_surface = !points_in_front_of_surface && !points_behind_surface;

        // return
        if (points_in_front_of_surface) return RelativePosition::IN_FRONT;
        else if (points_behind_surface) return RelativePosition::BEHIND;
        else if (points_within_surface)
        {
            projected_uncertainty = partial_combined_std;
            return RelativePosition::WITHIN;
        } 
        else throw std::runtime_error("Invalid relative position.");

        // // print all

        // double range = (point - origin).norm();
        // double angle = std::acos(direction.dot(normal_)) * 180.0 / M_PI - 90.0;

        // std::cout << "d_range_d_origin: \n" << d_range_d_origin << std::endl;
        // std::cout << "d_range_d_mean: \n" << d_range_d_mean << std::endl;
        // std::cout << "d_range_d_direction: \n" << d_range_d_direction << std::endl;
        // std::cout << "d_range_d_normal: \n" << d_range_d_normal << std::endl;

        // std::cout << "cov_origin: \n" << cov_origin << std::endl;
        // std::cout << "cov_mean: \n" << cov_mean << std::endl;
        // std::cout << "cov_direction: \n" << cov_direction << std::endl;
        // std::cout << "cov_normal: \n" << cov_normal << std::endl;

        // std::cout << "variance_range_origin: \n" << variance_range_origin << std::endl;
        // std::cout << "variance_range_mean: \n" << variance_range_mean << std::endl;
        // std::cout << "variance_range_direction: \n" << variance_range_direction << std::endl;
        // std::cout << "variance_range_normal: \n" << variance_range_normal << std::endl;

        // std::cout << "std_range_origin: \n" << std_range_origin << std::endl;
        // std::cout << "std_range_mean: \n" << std_range_mean << std::endl;
        // std::cout << "std_range_direction: \n" << std_range_direction << std::endl;
        // std::cout << "std_range_normal: \n" << std_range_normal << std::endl;

        // std::cout << range << " | " << projective_distance << "|" << angle << " | " << std_range_origin << " | " << std_range_mean << " | " << std_range_direction << " | " << std_range_normal << " | " << combined_std << std::endl;
    }

    // projective distance
    double projective_distance = compute_point_projective_distance(origin, point);

    // modified projective distance (taking into account accuracy in favor for points inside planes)
    double projective_distance_modified = sign(projective_distance) * std::max(0.0, std::fabs(projective_distance) - settings_.range_accuracy);

    // surface projective std
    double surface_position_std = get_surface_position_std_in_normal_direction();
    double surface_projective_std = surface_position_std / std::fabs(normal_.dot(direction));

    // point projective std
    double point_projective_std = settings_.range_precision;

    // combined projective std
    double combined_projective_std = std::sqrt(surface_projective_std * surface_projective_std + point_projective_std * point_projective_std);

    // multiplier for confidence interval
    double multiplier = 1.96; // for 95% confidence interval
    // double multiplier = settings_.envelope_size;

    // confidence interval values
    double threshold_in_front = multiplier * combined_projective_std;
    double threshold_behind = - multiplier * combined_projective_std;

    // check
    bool points_in_front_of_surface = projective_distance_modified > threshold_in_front;
    bool points_behind_surface = projective_distance_modified < threshold_behind;
    bool points_within_surface = !points_in_front_of_surface && !points_behind_surface;

    // return
    if (points_in_front_of_surface) return RelativePosition::IN_FRONT;
    else if (points_behind_surface) return RelativePosition::BEHIND;
    else if (points_within_surface) return RelativePosition::WITHIN;
    else throw std::runtime_error("Invalid relative position.");
}

RelativePosition Surface::check_relative_position(const std::shared_ptr<GenericPoint>& generic_point)
{
    return check_relative_position(generic_point->get_distance_travelled(), generic_point->get_origin(), generic_point->get_position(), generic_point->get_direction(), generic_point->get_projected_uncertainty());
}

RelativePosition Surface::check_relative_position(const std::shared_ptr<Vertex>& vertex)
{
    return check_relative_position(vertex->get_distance_travelled(), vertex->get_origin(), vertex->get_original_position(), vertex->get_direction(), vertex->get_projected_uncertainty());
}

// RelativePosition Surface::check_relative_position(const std::shared_ptr<Vertex>& vertex)
// {
//     // compute
//     // double projective_distance = vertex->buffer_compute_projected_distance(shared_from_this());
//     double projective_distance = compute_point_projective_distance(vertex);

//     double surface_position_std = get_surface_position_std_in_normal_direction();
//     double surface_projective_std = surface_position_std / std::fabs(normal_.dot(vertex->get_direction()));

//     // given projective_distance and surface_projective_std and range_precision and range_accuracy
//     double new_std = surface_projective_std + settings_.range_precision;

//     // distance is positive when in front of the surface
//     double threshold_in_front = 3.0 * new_std + settings_.range_accuracy;
//     double threshold_behind = - 3.0 * new_std - settings_.range_accuracy;
//     std::cout << "threshold_in_front: " << threshold_in_front << std::endl;
//     std::cout << "threshold_behind: " << threshold_behind << std::endl;

//     // check
//     bool points_in_front_of_surface = projective_distance > threshold_in_front;
//     bool points_behind_surface = projective_distance < threshold_behind;
//     bool points_within_surface = !points_in_front_of_surface && !points_behind_surface;

//     // return
//     if (points_in_front_of_surface) return RelativePosition::IN_FRONT;
//     else if (points_behind_surface) return RelativePosition::BEHIND;
//     else if (points_within_surface) return RelativePosition::WITHIN;
//     else throw std::runtime_error("Invalid relative position.");
// }

RelativePosition Surface::check_relative_position(const std::shared_ptr<InteriorPoint>& interior_point)
{
    return check_relative_position(interior_point->get_distance_travelled(), interior_point->get_origin(), interior_point->get_original_position(), interior_point->get_direction(), interior_point->get_projected_uncertainty());
}

// RelativePosition Surface::check_relative_position(const std::shared_ptr<InteriorPoint>& interior_point)
// {
//     // compute
//     double projective_distance = interior_point->buffer_compute_projected_distance(shared_from_this());

//     double surface_position_std = get_surface_position_std_in_normal_direction();
//     double surface_projective_std = surface_position_std / std::fabs(normal_.dot(interior_point->get_direction()));

//     bool points_in_front_of_surface = projective_distance > 3 * (settings_.range_precision + surface_projective_std);
//     bool points_behind_surface = projective_distance < - 3 * (settings_.range_precision + surface_projective_std);
//     bool points_within_surface = !points_in_front_of_surface && !points_behind_surface;

//     // return
//     if (points_in_front_of_surface) return RelativePosition::IN_FRONT;
//     else if (points_behind_surface) return RelativePosition::BEHIND;
//     else if (points_within_surface) return RelativePosition::WITHIN;
//     else throw std::runtime_error("Invalid relative position.");
// }

void Surface::merge_surface(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to merge surface with invalid surface.");

    // merge
    const auto& surface_valid = surface;
    for (const auto& vertex : surface_valid->vertices_) connect(vertex);

    // log
    if (settings_.log.merge_surface) std::cout << "Surface " << surface_valid->get_id() << " merged into surface " << id_ << std::endl;

    // delete
    storage_->delete_surface(surface);
}

const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& Surface::get_vertices() const
{
    return vertices_;
}

const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& Surface::get_interior_points() const
{
    return interior_points_;
}

const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Surface::get_edges() const
{
    return edges_;
}

const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Surface::get_faces() const
{
    return faces_;
}

const Eigen::Vector3d& Surface::get_mean() const
{
    return mean_;
}

const Eigen::Matrix3d& Surface::get_covariance() const
{
    return covariance_;
}

const Eigen::Matrix3d& Surface::get_eigenvectors() const
{
    return eigenvectors_;
}

const Eigen::Vector3d& Surface::get_eigenvalues() const
{
    return eigenvalues_;
}

const Eigen::Vector3d& Surface::get_normal() const
{
    return normal_;
}

std::size_t Surface::get_total_point_size() const
{
    return vertices_.size() + interior_points_.size();
}

double Surface::get_average_distance_travelled() const
{
    return sum_of_average_distance_travelled_ / get_total_point_size();
}

double Surface::get_smallest_distance_travelled() const
{
    return smallest_distance_travelled_;
}

const std::tuple<int, int, int>& Surface::get_color() const
{
    return color_;
}

const std::vector<double>& Surface::get_point_to_plane_distance_stats()
{
    // update if point count changes
    if (get_total_point_size() != previous_total_point_size_for_point_to_plane_) 
    {
        // clear
        stored_point_to_plane_distance_stats_.clear();

        // add
        for (const auto& vertex : vertices_)
        {
            double point2plane_distance = (vertex->get_original_position() - get_mean()).dot(get_normal());
            stored_point_to_plane_distance_stats_.push_back(point2plane_distance);
        }
        for (const auto& interior_point : interior_points_)
        {
            double point2plane_distance = (interior_point->get_original_position() - get_mean()).dot(get_normal());
            stored_point_to_plane_distance_stats_.push_back(point2plane_distance);
        }

        // update
        previous_total_point_size_for_point_to_plane_ = get_total_point_size();        
    }

    // return stored value
    return stored_point_to_plane_distance_stats_;
}

const std::vector<double>& Surface::get_projective_distance_stats()
{
    // update if point count changes
    if (get_total_point_size() != previous_total_point_size_for_projective_) 
    {
        // clear
        stored_projective_distance_stats_.clear();

        // add
        for (const auto& vertex : vertices_)
        {
            double projective_distance = vertex->compute_projected_distance();
            stored_projective_distance_stats_.push_back(projective_distance);
        }
        for (const auto& interior_point : interior_points_)
        {
            double projective_distance = interior_point->compute_projected_distance();
            stored_projective_distance_stats_.push_back(projective_distance);
        }

        // update
        previous_total_point_size_for_projective_ = get_total_point_size();        
    }

    // return stored value
    return stored_projective_distance_stats_;
}

double Surface::get_average_projective_distance()
{
    // compute
    double sum = 0;
    for (const double& distance : get_projective_distance_stats()) 
    {
        sum += std::fabs(distance);
    }

    // return
    return sum / get_total_point_size();
}

bool Surface::is_expired() const
{
    return is_expired_;
}

// decide wheather to remove this surface completely
bool Surface::is_abnormal()
{
    bool do_abnormal_check = false;
    if (!do_abnormal_check) return false;
    
    // not abnormal if low confidence surface
    if (get_total_point_size() < settings_.fit_plane_threshold) return false;

    // not abnormal if within range
    std::vector <double> projective_distance_stats = get_projective_distance_stats();
    // subtract accuracy from each distance
    std::vector<double> projective_distance_stats_modified = projective_distance_stats;
    for (double& distance : projective_distance_stats_modified) 
    {
        distance = sign(distance) * std::max(0.0, std::fabs(distance) - settings_.range_accuracy);
    }

    double new_projective_std = compute_std(projective_distance_stats_modified);

    // check if abnormal
    bool abnormal = new_projective_std > settings_.abnormal_size * settings_.range_precision;
    return abnormal;
}

bool Surface::can_merge(const std::shared_ptr<Surface>& surface) const
{
    // should use the abnormal surface test        
    // collect all points in this surface and all points in the other surface, 
    // compute the covariance matrix, 
    // decompose to find the normal, compute the projective distance.
    // compute stats of projective distance, 
    // if not abnormal, means can merge

    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to check mergeability with invalid surface.");

    // get mean and covariance
    Eigen::Vector3d mean1 = mean_;
    Eigen::Matrix3d cov1 = covariance_;
    int size1 = get_total_point_size();

    Eigen::Vector3d mean2 = surface->get_mean();
    Eigen::Matrix3d cov2 = surface->get_covariance();
    int size2 = surface->get_total_point_size();

    // combined mean and covariance
    Eigen::Vector3d new_mean = merge_mean(mean1, mean2, size1, size2);
    Eigen::Matrix3d new_cov = merge_covariance(cov1, cov2, mean1, mean2, size1, size2);

    // compute normal
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(new_cov);
    Eigen::Matrix3d new_eigenvectors = solver.eigenvectors();
    Eigen::Vector3d new_normal = new_eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal

    // compute projective distance stats for the combined surface
    std::vector<double> projective_distance_list;
    for (const auto& vertex : vertices_)                            projective_distance_list.push_back((new_mean - vertex->get_original_position()).dot(new_normal) / vertex->get_direction().dot(new_normal));
    for (const auto& interior_point : interior_points_)             projective_distance_list.push_back((new_mean - interior_point->get_original_position()).dot(new_normal) / interior_point->get_direction().dot(new_normal));
    for (const auto& vertex : surface->vertices_)                   projective_distance_list.push_back((new_mean - vertex->get_original_position()).dot(new_normal) / vertex->get_direction().dot(new_normal));
    for (const auto& interior_point : surface->interior_points_)    projective_distance_list.push_back((new_mean - interior_point->get_original_position()).dot(new_normal) / interior_point->get_direction().dot(new_normal));

    // subtract accuracy from each distance
    std::vector<double> projective_distance_list_modified = projective_distance_list;
    for (double& distance : projective_distance_list_modified) 
    {
        distance = sign(distance) * std::max(0.0, std::fabs(distance) - settings_.range_accuracy);
    }

    double new_projective_std = compute_std(projective_distance_list_modified);

    // check if mergable
    bool mergeable = new_projective_std <= settings_.range_precision;
    if (!mergeable)
    {
        if (settings_.log.can_merge) std::cout << "Surface " << id_ << " with " << get_total_point_size() << " points and surface " << surface->get_id() << " with " << surface->get_total_point_size() << " points are not mergable." << std::endl;
        if (settings_.log.can_merge) std::cout << "New projective std: " << new_projective_std << std::endl;
    }

    // return
    return mergeable;
}

void Surface::connect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid vertex.");

    {
        // write lock
        std::unique_lock lock(rwlock_vertices_); 

        // skip if already connected
        if (vertices_.find(vertex) != vertices_.end()) return;

        // insert
        vertices_.insert(vertex);
    }

    // reverse connection
    vertex->connect(shared_from_this());

    // update uncertainty
    if (get_total_point_size() <= settings_.fit_plane_threshold) 
    {
        vertex->weight_ = 1.0 / (settings_.range_precision * settings_.range_precision);
    }
    else
    {
        // within the check relative position, the uncertainty will be updated
        check_relative_position(vertex);
        vertex->weight_ = 1.0 / (vertex->get_projected_uncertainty() * vertex->get_projected_uncertainty());
    }

    add_point_to_surface_fitting(vertex->get_original_position(), vertex->get_origin(), vertex->get_distance_travelled(), vertex->weight_);
}

bool Surface::tree_intersect_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
{
    // read lock
    std::shared_lock lock(rwlock_edge_bvh_);

    return edge_bvh_.tree_intersect_edge(vertex0, vertex1);
}

bool Surface::connect_by_edges_and_faces(const std::shared_ptr<Vertex>& vertex, const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& all_nearby_vertices)
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(vertex->rwlock_lifecycle_);
    // no need to check expired since it is not yet connected to anything

    // try create edge with nearby vertices
    for (const auto& nearby_vertex : all_nearby_vertices)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(nearby_vertex->rwlock_lifecycle_);

        // check input
        if (nearby_vertex->is_expired()) continue;

        // skip if does not belong to the same surface
        if (nearby_vertex->get_surface() != shared_from_this()) continue;

        // skip if not boundary
        if (!nearby_vertex->is_boundary()) continue;

        // skip if same vertex
        if (nearby_vertex == vertex) continue;
        
        // skip if edge is longer than any of the radius of vertices
        const double distance = (vertex->get_position() - nearby_vertex->get_position()).norm();
        const double radius0 = vertex->get_radius(shared_from_this());
        const double radius1 = nearby_vertex->get_radius();
        if (!settings_.edge_is_short_enough(distance, radius0, radius1)) continue;

        // skip if edge intersects
        if (tree_intersect_edge(vertex, nearby_vertex)) continue;

        // create edge
        std::shared_ptr<Edge> new_edge = storage_->add_edge(vertex, nearby_vertex);

        // read lock edge
        std::shared_lock<std::shared_mutex> lock_edge(new_edge->rwlock_lifecycle_);

        // skip if edge is expired
        if (new_edge->is_expired()) continue;

        // connect
        connect(new_edge);
    }

    // false if no edge is created
    if (vertex->get_edges().empty()) return false;

    // try close holes
    vertex->try_close_holes_repeatedly();

    // else
    return true;
}

void Surface::compute_surface_position_std_in_normal_direction()
{
    // compute information weighted mean and std

    // get mean from all points and std from all points
    std::vector<double> mean_list;
    std::vector<double> information_list;
    for (auto& vertex : vertices_)
    {
        // compute conversion ratio
        double ratio = std::fabs(normal_.dot(vertex->get_direction()));

        // mean and informaiton
        double mean = ratio * vertex->compute_projected_distance();
        double std = ratio * settings_.range_precision;
        double information = 1.0 / (std * std);
        
        // store
        mean_list.push_back(mean);
        information_list.push_back(information);
    }
    for (auto& interior_point : interior_points_)
    {
        // compute conversion ratio
        double ratio = std::fabs(normal_.dot(interior_point->get_direction()));

        // mean and informaiton
        double mean = ratio * interior_point->compute_projected_distance();
        double std = ratio * settings_.range_precision;
        double information = 1.0 / (std * std);
        
        // store
        mean_list.push_back(mean);
        information_list.push_back(information);
    }

    // compute information weighted mean and std
    double weighted_mean = 0;
    double weighted_information = 0;
    double weighted_std = 0;
    for (std::size_t i = 0; i < mean_list.size(); i++)
    {
        weighted_mean += mean_list[i] * information_list[i];
        weighted_information += information_list[i];
    }
    weighted_mean /= weighted_information;
    weighted_std = 1.0 / std::sqrt(weighted_information);
    
    // // log
    // std::cout << "computed surface position with Bayesian mean: " << weighted_mean << ", Bayesian std: " << weighted_std << std::endl;
    
    // store
    previous_normal_distance_ = weighted_mean;
    previous_normal_std_ = weighted_std;
}

double Surface::get_surface_position_std_in_normal_direction()
{
    return previous_normal_std_;
}

void Surface::connect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid edge.");

    {
        // write lock
        std::unique_lock lock(rwlock_edges_);

        // skip if already connected
        if (edges_.find(edge) != edges_.end()) return;

        // insert
        edges_.insert(edge);
    }

    // add to searchable edge
    add_searchable_edge(edge);

    // reverse connection
    edge->connect(shared_from_this());
}

void Surface::connect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid face.");

    {
        // write lock
        std::unique_lock lock(rwlock_faces_);

        // skip if already connected
        if (faces_.find(face) != faces_.end()) return;

        // insert
        faces_.insert(face);
    }

    // reverse connection
    face->connect(shared_from_this());
}

void Surface::connect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid interior point.");

    {
        // write lock
        std::unique_lock lock(rwlock_interior_points_);

        // skip if already connected
        if (interior_points_.find(interior_point) != interior_points_.end()) return;

        // insert
        interior_points_.insert(interior_point);
    }
    
    // update uncertainty
    if (get_total_point_size() <= settings_.fit_plane_threshold) 
    {
        interior_point->weight_ = 1.0 / (settings_.range_precision * settings_.range_precision);
    }
    else
    {
        // within the check relative position, the uncertainty will be updated
        check_relative_position(interior_point);
        interior_point->weight_ = 1.0 / (interior_point->get_projected_uncertainty() * interior_point->get_projected_uncertainty());
    }

    add_point_to_surface_fitting(interior_point->get_original_position(), interior_point->get_origin(), interior_point->get_distance_travelled(), interior_point->weight_);
    
    // reverse connection
    interior_point->connect(shared_from_this());
}

void Surface::disconnect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) return;

    {
        // write lock
        std::unique_lock lock(rwlock_vertices_);

        // skip if not connected
        auto it = vertices_.find(vertex);
        if (it == vertices_.end()) return;

        // erase
        vertices_.erase(it);
    }

    remove_point_from_surface_fitting(vertex->get_original_position(), vertex->get_origin(), vertex->get_distance_travelled(), vertex->weight_);

    // reverse disconnect
    vertex->disconnect(shared_from_this());
}

void Surface::disconnect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) return;

    {
        // write lock
        std::unique_lock lock(rwlock_edges_);

        // skip if not connected
        auto it = edges_.find(edge);
        if (it == edges_.end()) return;

        // erase
        edges_.erase(it);
    }

    remove_searchable_edge(edge);

    // reverse disconnect
    edge->disconnect(shared_from_this());
}

void Surface::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    {
        // write lock
        std::unique_lock lock(rwlock_faces_);

        // skip if not connected
        auto it = faces_.find(face);
        if (it == faces_.end()) return;

        // erase
        faces_.erase(it);
    }

    // reverse disconnect
    face->disconnect(shared_from_this());
}

void Surface::disconnect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    {
        // write lock
        std::unique_lock lock(rwlock_interior_points_);

        // skip if not connected
        auto it = interior_points_.find(interior_point);
        if (it == interior_points_.end()) return;

        // erase
        interior_points_.erase(it);
    }

    remove_point_from_surface_fitting(interior_point->get_original_position(), interior_point->get_origin(), interior_point->get_distance_travelled(), interior_point->weight_);
    
    // reverse disconnect
    interior_point->disconnect(shared_from_this());
}

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> Surface::get_boundary_vertices()
{
    // initialize
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> boundary_vertices;

    // get boundary vertices
    for (const auto& vertex : vertices_)
    {
        if (vertex->is_boundary()) boundary_vertices.insert(vertex);
    }

    // return
    return boundary_vertices;
}

bool Surface::try_close_holes_repeatedly()
{
    // initialize
    bool changed = false;

    // try close holes repeatedly
    while (try_close_holes()) 
    {
        changed = true;
    }

    // return
    return changed;
}

bool Surface::try_close_holes()
{
    // initialize
    bool changed = false;

    // boundary vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> boundary_vertices = get_boundary_vertices();

    // for each boundary vertex
    for (const auto& vertex : boundary_vertices)
    {
        if (vertex->try_close_holes_repeatedly()) changed = true;
    }

    // return
    return changed;
}

void Surface::remove_non_manifold_edges()
{
    // make copy of edges
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_copy = edges_;

    for (const auto& edge : edges_copy)
    {
        if (edge->is_non_manifold()) disconnect(edge);
    }
}

void Surface::swap(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{
    // if contains vertex1
    if (vertices_.find(vertex1) != vertices_.end())
    {
        disconnect(vertex1);
        connect(vertex2);
    }
}

void Surface::pause_normal_std_update()
{
    update_normal_position_std_ = false;
}

void Surface::resume_normal_std_update()
{
    update_normal_position_std_ = true;
    compute_surface_position_std_in_normal_direction();
}

void Surface::add_searchable_edge(const std::shared_ptr<Edge>& edge)
{
    // write lock
    std::unique_lock lock(rwlock_edge_bvh_);

    edge_bvh_.tree_add_edge(edge);
}

void Surface::remove_searchable_edge(const std::shared_ptr<Edge>& edge)
{
    // write lock
    std::unique_lock lock(rwlock_edge_bvh_);
    
    edge_bvh_.tree_delete_edge(edge);
}

void Surface::print_info()
{
    std::cout << "============================= info for surface " << id_ << " =============================" << std::endl;
    // print surface info
    std::cout << "Surface " << id_ << " has " << vertices_.size() << " vertices, " << interior_points_.size() << " interior points, " << edges_.size() << " edges, " << faces_.size() << " faces." << std::endl;

    // print covariance etc.
    std::cout << "Mean: " << mean_.transpose() << std::endl;
    std::cout << "Covariance: " << covariance_ << std::endl;
    std::cout << "Normal: " << normal_.transpose() << std::endl;
    std::cout << "Eigenvalues: " << eigenvalues_.transpose() << std::endl;
    std::cout << "Eigenvectors: " << eigenvectors_ << std::endl;
    std::cout << "Color: " << std::get<0>(color_) << ", " << std::get<1>(color_) << ", " << std::get<2>(color_) << std::endl;
    std::cout << "Total point size: " << get_total_point_size() << std::endl;
    std::cout << "Average projective distance: " << get_average_projective_distance() << std::endl;
    std::cout << "Projective distance std: " << compute_std(get_projective_distance_stats()) << std::endl;
    std::cout << "Projective distance mean: " << compute_mean(get_projective_distance_stats()) << std::endl;

    std::cout << "======================================================================================" << std::endl;
}

void Surface::add_point_to_surface_fitting(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled, double weight)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_surface_fitting_);

    // surface
    double weight1 = weight_;
    Eigen::Vector3d mean1 = mean_;
    Eigen::Matrix3d cov1 = covariance_;

    // point
    double weight2 = weight;
    Eigen::Vector3d mean2 = position;
    Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();

    // set + point
    Eigen::Vector3d new_mean = weighted_merge_mean(mean1, mean2, weight1, weight2);
    Eigen::Matrix3d new_cov = weighted_merge_covariance(cov1, cov2, mean1, mean2, weight1, weight2);
    double new_weight = weight1 + weight2;

    // plane estimate
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(new_cov);
    Eigen::Matrix3d new_eigenvectors = solver.eigenvectors();
    Eigen::Vector3d new_eigenvalues = solver.eigenvalues();
    Eigen::Vector3d new_normal = new_eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal
    Eigen::Vector3d vector_towards_origin = origin - position;
    if (new_normal.dot(vector_towards_origin) < 0) new_normal *= -1; // normal should points towards the origin

    // store
    weight_ = new_weight;
    mean_ = new_mean;
    covariance_ = new_cov;
    eigenvectors_ = new_eigenvectors;
    eigenvalues_ = new_eigenvalues;
    normal_ = new_normal;

    sum_of_average_distance_travelled_ += distance_travelled;
    if (distance_travelled < smallest_distance_travelled_) smallest_distance_travelled_ = distance_travelled;

    // store characteristic length
    characteristic_length_ = std::max(characteristic_length_, (position -  mean_).norm());
    normal_uncertainty_ = settings_.range_precision / characteristic_length_;

    // update approximate uncertainty envelope
    // if approximate normal is the same as last time, incrementally update the uncertianty envelope
    // if approximate normal is not the same as last time, recompute the uncertainty envelope all together

    if (!update_normal_position_std_) return;

    // incrementally update
    double old_distance = previous_normal_distance_;
    double old_std = previous_normal_std_;
    double old_information = 1.0 / (old_std * old_std);

    double new_distance = compute_point_projective_distance(origin, position);
    double new_std = settings_.range_precision;
    double new_information = 1.0 / (new_std * new_std);

    double combined_distance = merge_information_weighted_mean(old_distance, new_distance, old_information, new_information);
    double combined_information = merge_information(old_information, new_information);
    double combined_std = 1.0 / std::sqrt(combined_information);

    previous_normal_distance_ = combined_distance;
    previous_normal_std_ = combined_std;
}

void Surface::remove_point_from_surface_fitting(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled, double weight)
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_surface_fitting_);
    
    // surface
    double combined_weight = weight_;
    Eigen::Vector3d combined_mean = mean_;
    Eigen::Matrix3d combined_cov = covariance_;

    // point
    double weight2 = weight;
    Eigen::Vector3d mean2 = position;
    Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();

    // set + point
    Eigen::Vector3d mean1 = weighted_remove_mean(combined_mean, mean2, combined_weight, weight2);
    Eigen::Matrix3d cov1 = weighted_remove_covariance(combined_cov, cov2, combined_mean, mean2, combined_weight, weight2);
    double weight1 = combined_weight - weight2;

    // plane estimate
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov1);
    Eigen::Matrix3d eigenvectors1 = solver.eigenvectors();
    Eigen::Vector3d eigenvalues1 = solver.eigenvalues();
    Eigen::Vector3d normal1 = eigenvectors1.col(0); // Assuming the smallest eigenvalue corresponds to the normal
    Eigen::Vector3d vector_towards_origin = origin - position;
    if (normal1.dot(vector_towards_origin) < 0) normal1 *= -1; // normal should points towards the origin

    // store
    weight_ = weight1;
    mean_ = mean1;
    covariance_ = cov1;
    eigenvectors_ = eigenvectors1;
    eigenvalues_ = eigenvalues1;
    normal_ = normal1;

    sum_of_average_distance_travelled_ -= distance_travelled;
    if (distance_travelled == smallest_distance_travelled_)
    {
        smallest_distance_travelled_ = std::numeric_limits<double>::max();
        for (const auto& vertex : vertices_)
        {
            if (vertex->get_distance_travelled() < smallest_distance_travelled_) smallest_distance_travelled_ = vertex->get_distance_travelled();
        }
        for (const auto& interior_point : interior_points_)
        {
            if (interior_point->get_distance_travelled() < smallest_distance_travelled_) smallest_distance_travelled_ = interior_point->get_distance_travelled();
        }
    }

    // update approximate uncertainty envelope
    // if approximate normal is the same as last time, incrementally update the uncertianty envelope
    // if approximate normal is not the same as last time, recompute the uncertainty envelope all together

    if (!update_normal_position_std_) return;

    // incrementally update
    double combined_distance = previous_normal_distance_;
    double combined_std = previous_normal_std_;
    double combined_information = 1.0 / (combined_std * combined_std);

    double new_distance = compute_point_projective_distance(origin, position);
    double new_std = settings_.range_precision;
    double new_information = 1.0 / (new_std * new_std);

    double old_distance = remove_information_weighted_mean(combined_distance, new_distance, combined_information, new_information);
    double old_information = remove_information(combined_information, new_information);
    double old_std = 1.0 / std::sqrt(old_information);

    previous_normal_distance_ = old_distance;
    previous_normal_std_ = old_std;
}

void Surface::set_random_color()
{
    // random color
    color_ = std::make_tuple(rand() % 256, rand() % 256, rand() % 256);
}

void Surface::optimize_surface_normal()
{
    // FIT PLANE
    double bearing_noise = 0.01;
    double range_noise = 0.01;

    // generate dataset
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> dataset;
    for (const auto& vertex : vertices_)
    {
        dataset.push_back(std::make_pair(vertex->get_origin(), vertex->get_original_position()));
    }
    for (const auto& interior_point : interior_points_)
    {
        dataset.push_back(std::make_pair(interior_point->get_origin(), interior_point->get_original_position()));
    }

    fit_plane_to_points(dataset, mean_, normal_, bearing_noise, range_noise);
}

bool Surface::remove_unmatched_points()
{
    // skip if less than 3 points   
    if (get_total_point_size() < 3) return false;

    // collect vertices to delete
    std::vector<std::shared_ptr<Vertex>> vertices_to_delete;
    for (const auto& vertex : vertices_)
    {
        // check input
        if (vertex->is_expired()) throw std::runtime_error("Surface contains expired vertex.");

        // if vertex is not within the surface, add to remove list
        if (check_relative_position(vertex) != RelativePosition::WITHIN) vertices_to_delete.push_back(vertex);
    }

    // collect interior points to delete
    std::vector<std::shared_ptr<InteriorPoint>> interior_points_to_delete;
    for (const auto& interior_point : interior_points_)
    {
        // check input
        if (interior_point->is_expired()) throw std::runtime_error("Surface contains expired interior point.");

        // if interior point is not within the surface, add to remove list
        if (check_relative_position(interior_point) != RelativePosition::WITHIN) interior_points_to_delete.push_back(interior_point);
    }

    // delete points
    for (const auto& vertex : vertices_to_delete)
    {
        if (vertex->is_expired()) continue;
        storage_->delete_vertex(vertex);
    }

    for (const auto& interior_point : interior_points_to_delete)
    {
        if (interior_point->is_expired()) continue;
        storage_->delete_interior_point(interior_point);
    }

    // return
    return !vertices_to_delete.empty();
}

void Surface::remove_singular_components()
{
    // get singular components
    std::vector<std::shared_ptr<Vertex>> singular_vertices;
    std::vector<std::shared_ptr<Edge>> singular_edges;
    for (const auto& vertex : vertices_) if (vertex->is_singular()) singular_vertices.push_back(vertex);
    for (const auto& edge : edges_) if (edge->is_singular()) singular_edges.push_back(edge);

    // delete singular components
    for (const auto& vertex : singular_vertices)
    {
        if (vertex->is_expired()) continue;
        storage_->delete_vertex(vertex);
    }
    for (const auto& edge : singular_edges)
    {
        if (edge->is_expired()) continue;
        storage_->delete_edge(edge);
    }
}

void Surface::split_surface_by_connected_components()
{
    // split surface
    UnionFind uf;
    uf.add_vertices(vertices_);
    uf.add_edges(edges_);
    std::vector<std::pair<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Vertex>>>> sorted_grouped_vertices = uf.compute_sorted_grouped_vertices();

    // create new surfaces from the second group onwards
    for (std::size_t i = 1; i < sorted_grouped_vertices.size(); i++)
    {
        // get root vertex
        const auto& root_vertex = sorted_grouped_vertices[i].first;

        // connect to new surface
        std::shared_ptr<Surface> new_surface = storage_->add_surface();
        root_vertex->swap(shared_from_this(), new_surface);
    }
}

bool operator<(const std::shared_ptr<Surface> &lhs, const std::shared_ptr<Surface> &rhs)
{
    // check pointer validity
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs)
{
    if (!lhs && !rhs) return true; // true if both are nullptr
    if (!lhs || !rhs) return false; // false if either is nullptr
    // check pointer validity
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs->get_id() == rhs->get_id();
}

bool operator>= (const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs)
{
    // check pointer validity
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs->get_id() >= rhs->get_id();
}

bool operator!= (const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs)
{
    // check pointer validity
    if (!lhs && !rhs) return false; // false if both are nullptr
    if (!lhs || !rhs) return true; // true if either is nullptr
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs->get_id() != rhs->get_id();
}