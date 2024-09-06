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
    omp_init_lock(&lock_);
}

Surface::~Surface()
{
    omp_destroy_lock(&lock_);
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

double Surface::compute_point_projective_distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& position) const
{
    // if perpendicular, return NaN
    Eigen::Vector3d rayDirection = (position - origin).normalized();
    double distance = (mean_ - position).dot(normal_) / rayDirection.dot(normal_);

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
    // compute
    Eigen::Vector3d rayDirection = (point - origin).normalized();
    double distance = (mean_ - point).dot(normal_) / rayDirection.dot(normal_);
    Eigen::Vector3d intersection = point + distance * rayDirection;

    // return
    return intersection;
}

RelativePosition Surface::check_relative_position(const Eigen::Vector3d& origin, const Eigen::Vector3d& point, const Eigen::Vector3d& direction)
{
    // compute
    double projective_distance = compute_point_projective_distance(origin, point);

    double surface_position_std = get_surface_position_std_in_normal_direction();
    double surface_projective_std = surface_position_std / std::fabs(normal_.dot(direction));

    // given projective_distance and surface_projective_std and range_precision and range_accuracy
    double new_std = std::sqrt(surface_projective_std * surface_projective_std + settings_.range_precision * settings_.range_precision);

    // distance is positive when in front of the surface
    // modify projective_distance
    projective_distance = sign(projective_distance) * std::max(0.0, std::fabs(projective_distance) - settings_.range_accuracy);

    double threshold_in_front = settings_.envelope_size * new_std;
    double threshold_behind = - settings_.envelope_size * new_std;

    // check
    bool points_in_front_of_surface = projective_distance > threshold_in_front;
    bool points_behind_surface = projective_distance < threshold_behind;
    bool points_within_surface = !points_in_front_of_surface && !points_behind_surface;

    // return
    if (points_in_front_of_surface) return RelativePosition::IN_FRONT;
    else if (points_behind_surface) return RelativePosition::BEHIND;
    else if (points_within_surface) return RelativePosition::WITHIN;
    else throw std::runtime_error("Invalid relative position.");
}

RelativePosition Surface::check_relative_position(const std::shared_ptr<GenericPoint>& generic_point)
{
    return check_relative_position(generic_point->get_origin(), generic_point->get_position(), generic_point->get_direction());
}

RelativePosition Surface::check_relative_position(const std::shared_ptr<Vertex>& vertex)
{
    return check_relative_position(vertex->get_origin(), vertex->get_position(), vertex->get_direction());
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
    return check_relative_position(interior_point->get_origin(), interior_point->get_position(), interior_point->get_direction());
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

std::size_t Surface::get_approximate_normal_hash()
{
    double factor = 100.0;
    Eigen::Vector3d approximate_normal = (normal_ * factor).array().round() / factor;
    approximate_normal = approximate_normal.normalized();

    std::size_t h1 = std::hash<double>{}(approximate_normal.x());
    std::size_t h2 = std::hash<double>{}(approximate_normal.y());
    std::size_t h3 = std::hash<double>{}(approximate_normal.z());
    std::size_t hash = h1 ^ (h2 << 1) ^ (h3 << 2); // Combining hashes

    return hash;
}

std::size_t Surface::get_total_point_size() const
{
    return vertices_.size() + interior_points_.size();
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
            double point2plane_distance = (vertex->get_position() - get_mean()).dot(get_normal());
            stored_point_to_plane_distance_stats_.push_back(point2plane_distance);
        }
        for (const auto& interior_point : interior_points_)
        {
            double point2plane_distance = (interior_point->get_position() - get_mean()).dot(get_normal());
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
            double projective_distance = vertex->buffer_compute_projected_distance(shared_from_this());
            stored_projective_distance_stats_.push_back(projective_distance);
        }
        for (const auto& interior_point : interior_points_)
        {
            double projective_distance = interior_point->buffer_compute_projected_distance(shared_from_this());
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
    for (const auto& vertex : vertices_)                            projective_distance_list.push_back((new_mean - vertex->get_position()).dot(new_normal) / vertex->get_direction().dot(new_normal));
    for (const auto& interior_point : interior_points_)             projective_distance_list.push_back((new_mean - interior_point->get_position()).dot(new_normal) / interior_point->get_direction().dot(new_normal));
    for (const auto& vertex : surface->vertices_)                   projective_distance_list.push_back((new_mean - vertex->get_position()).dot(new_normal) / vertex->get_direction().dot(new_normal));
    for (const auto& interior_point : surface->interior_points_)    projective_distance_list.push_back((new_mean - interior_point->get_position()).dot(new_normal) / interior_point->get_direction().dot(new_normal));

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

    // connect
    if (vertices_.insert(vertex).second)
    {
        vertex->connect(shared_from_this());
        add_point_to_surface_fitting(vertex->get_position(), vertex->get_origin());
    }
}

bool Surface::connect_by_edges_and_faces(const std::shared_ptr<Vertex>& vertex, const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& all_nearby_vertices)
{
    // initialize
    bool connected = false;

    // get nearby vertices in the same surface
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> nearby_vertices;
    for (const auto& nearby_vertex : all_nearby_vertices)
    {
        // check input
        if (nearby_vertex->is_expired()) continue;

        // skip if same vertex
        if (nearby_vertex == vertex) continue;

        // skip if does not belong to the same surface
        if (nearby_vertex->get_surface() != shared_from_this()) continue;

        // add to nearby vertices
        nearby_vertices.insert(nearby_vertex);
    }

    // update search radius of the vertex
    for (const auto& nearby_vertex : nearby_vertices)
    {
        double distance = (vertex->get_position() - nearby_vertex->get_position()).norm();
        vertex->reduce_reverse_radius_search_radius(distance + nearby_vertex->get_radius());
    }

    // create edges
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> used_vertices;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> new_edges;
    for (const auto& nearby_vertex : nearby_vertices)
    {
        // skip if edge is longer than any of the radius of vertices
        double distance = (vertex->get_position() - nearby_vertex->get_position()).norm();
        
        // skip if edge is longer than either vertices radius
        if (distance > vertex->get_radius(shared_from_this()) || distance > nearby_vertex->get_radius()) continue;

        // if edge intersects
        if (edge_bvh_.tree_intersect_edge(vertex, nearby_vertex)) 
        {   
            // log
            if (settings_.log.connect_by_edges_and_faces) std::cout << "Try to create edge between " << vertex->get_id() << " and " << nearby_vertex->get_id() << " but is intersected." << std::endl;

            // should not reduce search radius!!!
            // radius represents extends of flat surface -> edge intersection within the same plane is still flat surface!

            // // that means the nearby_vertex have too large of search radius
            // // so we should reduce it
            // double distance = (vertex->get_position() - nearby_vertex->get_position()).norm();
            // nearby_vertex->reduce_reverse_radius_search_radius(distance);
        }
        else // if edge does not intersect
        {
            // create edge
            std::shared_ptr<Edge> new_edge = storage_->add_edge(vertex, nearby_vertex);
            connect(new_edge);
            connect(vertex);
            used_vertices.insert(nearby_vertex);
            new_edges.insert(new_edge);

            connected = true;

            // check if the new edge have any sibling edges
            for (const auto& sibling_vertex : vertex->get_sibling_vertices())
            {
                for (const auto& edge : sibling_vertex->get_edges())
                {
                    if (edge->has_vertex(nearby_vertex))
                    {
                        new_edge->connect(edge);
                    }
                }
            }
        }
    }

    // create faces
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> new_faces;
    for (const std::shared_ptr<Vertex>& nearby_vertex0 : used_vertices)
    {
        for (const std::shared_ptr<Vertex>& nearby_vertex1 : used_vertices)
        {
            // skip if repeated
            if (nearby_vertex1 <= nearby_vertex0) continue;

            // skip if edge does not exist between the two vertices
            bool edge_exist = false;
            std::shared_ptr<Edge> existing_edge;
            for (const std::shared_ptr<Edge>& edge : nearby_vertex0->get_edges())
            {
                if (edge->get_surface() != shared_from_this()) continue;
                if (edge->has_vertex(nearby_vertex1))
                {
                    edge_exist = true;
                    existing_edge = edge;
                    break;
                }
            }
            if (!edge_exist) continue;

            // skip if edge is not boundary
            if (!existing_edge->is_boundary()) continue;

            // skip if face contains nearby vertices
            // get surface coordinate of the vertices
            Eigen::Vector2d surface_coordinate0 = vertex->get_surface_coordinate(shared_from_this());
            Eigen::Vector2d surface_coordinate1 = nearby_vertex0->get_surface_coordinate(shared_from_this());
            Eigen::Vector2d surface_coordinate2 = nearby_vertex1->get_surface_coordinate(shared_from_this());
            // check if the triangle formed by the vertices contains any other nearby_vertices
            bool triangle_contains_nearby_vertices = false;
            for (const std::shared_ptr<Vertex>& nearby_vertex : nearby_vertices)
            {
                if (nearby_vertex == nearby_vertex0 || nearby_vertex == nearby_vertex1) continue;
                Eigen::Vector2d surface_coordinate = nearby_vertex->get_surface_coordinate(shared_from_this());
                if (is_point_in_triangle(surface_coordinate0, surface_coordinate1, surface_coordinate2, surface_coordinate))
                {
                    triangle_contains_nearby_vertices = true;
                    break;
                }
            }
            if (triangle_contains_nearby_vertices) continue;

            // don't remove silver triangle, need better ways

            // // skip if face is too thin (silver triangle)
            // Eigen::Vector3d edge1 = nearby_vertex0->get_position() - vertex->get_position();
            // Eigen::Vector3d edge2 = nearby_vertex1->get_position() - vertex->get_position();
            // double angle = std::acos(edge1.normalized().dot(edge2.normalized())) * 180 / M_PI;
            // if (angle < settings_.min_face_angle) continue;
            // if (angle > (180 - 2.0*settings_.min_face_angle)) continue;


            // // check if face already exists
            // bool face_exist = false;
            // std::shared_ptr<Face> existing_face;
            // for (const std::shared_ptr<Face>& face : existing_edge->get_faces())
            // {
            //     if (face->has_vertex(vertex))
            //     {
            //         // log
            //         std::cout << "Face already exists between " << vertex->get_id() << " and " << nearby_vertex0->get_id() << " and " << nearby_vertex1->get_id() << std::endl;
                    
            //         face_exist = true;
            //         existing_face = face;
            //         break;
            //     }
            // }
            // if (face_exist)
            // {
            //     connect(existing_face);
            //     continue;
            // }

            // if face not already exists, create face
            std::shared_ptr<Face> new_face = storage_->add_face(shared_from_this(), vertex, nearby_vertex0, nearby_vertex1);
            new_faces.insert(new_face);

            // connnect new face to its sibling faces
            for (const auto& sibling_edge : existing_edge->get_sibling_edges())
            {
                for (const auto& sibling_face : sibling_edge->get_faces())
                {
                    if (sibling_face->has_vertex(vertex))
                    {
                        new_face->connect(sibling_face);
                    }
                }
            }
        }
    }

    // return
    return connected;
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
        double mean = ratio * vertex->buffer_compute_projected_distance(shared_from_this());
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
        double mean = ratio * interior_point->buffer_compute_projected_distance(shared_from_this());
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

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge->connect(shared_from_this());
}

void Surface::connect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face->connect(shared_from_this());
}

void Surface::connect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid interior point.");

    // connect
    bool inserted = interior_points_.insert(interior_point).second;
    if (inserted) interior_point->connect(shared_from_this());

    // update surface fitting
    if (inserted) add_point_to_surface_fitting(interior_point->get_position(), interior_point->get_origin());
}

void Surface::disconnect(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) return;

    // disconnect
    bool erased = vertices_.erase(vertex);
    if (erased) vertex->disconnect(shared_from_this());

    // remove from surface fitting
    if (erased) remove_point_from_surface_fitting(vertex->get_position(), vertex->get_origin());
}

void Surface::disconnect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) return;

    // disconnect
    bool erased = edges_.erase(edge);
    if (erased) edge->disconnect(shared_from_this());
}

void Surface::disconnect(const std::shared_ptr<Face>& face)
{
    // check input
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());
}

void Surface::disconnect(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    // disconnect
    bool erased = interior_points_.erase(interior_point);
    if (erased) interior_point->disconnect(shared_from_this());

    // remove from surface fitting
    if (erased) remove_point_from_surface_fitting(interior_point->get_position(), interior_point->get_origin());
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
    edge_bvh_.tree_add_edge(edge);
}

void Surface::remove_searchable_edge(const std::shared_ptr<Edge>& edge)
{
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

void Surface::add_point_to_surface_fitting(const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
{
    // surface
    int size1 = get_total_point_size()-1; // need to exclude the new point
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
    Eigen::Vector3d new_eigenvalues = solver.eigenvalues();
    Eigen::Vector3d new_normal = new_eigenvectors.col(0); // Assuming the smallest eigenvalue corresponds to the normal
    Eigen::Vector3d vector_towards_origin = origin - position;
    if (new_normal.dot(vector_towards_origin) < 0) new_normal *= -1; // normal should points towards the origin

    // store
    mean_ = new_mean;
    covariance_ = new_cov;
    eigenvectors_ = new_eigenvectors;
    eigenvalues_ = new_eigenvalues;
    normal_ = new_normal;

    // update approximate uncertainty envelope
    // if approximate normal is the same as last time, incrementally update the uncertianty envelope
    // if approximate normal is not the same as last time, recompute the uncertainty envelope all together

    if (!update_normal_position_std_) return;

    // approximate_normal hash
    std::size_t hash = get_approximate_normal_hash();

    if (hash != previous_approximate_normal_hash_)
    {
        // recompute
        compute_surface_position_std_in_normal_direction();
        previous_approximate_normal_hash_ = hash;
    }
    else
    {
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
}

void Surface::remove_point_from_surface_fitting(const Eigen::Vector3d& position, const Eigen::Vector3d& origin)
{
    // surface
    int combined_size = get_total_point_size()+1; // need to include the point just removed
    Eigen::Vector3d combined_mean = mean_;
    Eigen::Matrix3d combined_cov = covariance_;

    // point
    int size2 = 1;
    Eigen::Vector3d mean2 = position;
    Eigen::Matrix3d cov2 = Eigen::Matrix3d::Zero();

    // set + point
    Eigen::Vector3d mean1 = remove_mean(combined_mean, mean2, combined_size, size2);
    Eigen::Matrix3d cov1 = remove_covariance(combined_cov, cov2, combined_mean, mean2, combined_size, size2);

    // plane estimate
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov1);
    Eigen::Matrix3d eigenvectors1 = solver.eigenvectors();
    Eigen::Vector3d eigenvalues1 = solver.eigenvalues();
    Eigen::Vector3d normal1 = eigenvectors1.col(0); // Assuming the smallest eigenvalue corresponds to the normal
    Eigen::Vector3d vector_towards_origin = origin - position;
    if (normal1.dot(vector_towards_origin) < 0) normal1 *= -1; // normal should points towards the origin

    mean_ = mean1;
    covariance_ = cov1;
    eigenvectors_ = eigenvectors1;
    eigenvalues_ = eigenvalues1;
    normal_ = normal1;

    // update approximate uncertainty envelope
    // if approximate normal is the same as last time, incrementally update the uncertianty envelope
    // if approximate normal is not the same as last time, recompute the uncertainty envelope all together

    if (!update_normal_position_std_) return;

    // approximate_normal hash
    std::size_t hash = get_approximate_normal_hash();

    if (hash != previous_approximate_normal_hash_)
    {
        // recompute
        compute_surface_position_std_in_normal_direction();
        previous_approximate_normal_hash_ = hash;
    }
    else
    {
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
        dataset.push_back(std::make_pair(vertex->get_origin(), vertex->get_position()));
    }
    for (const auto& interior_point : interior_points_)
    {
        dataset.push_back(std::make_pair(interior_point->get_origin(), interior_point->get_position()));
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