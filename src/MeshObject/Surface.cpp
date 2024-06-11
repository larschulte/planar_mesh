#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/InteriorPoint.hpp"
#include "MeshObject/Storage.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>
#include "MeshObject/InteriorPoint.hpp"

void Surface::initialize_(std::shared_ptr<Storage> storage)
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
    
    // log
    std::cout << "Surface " << id_ << " created.\n";
}

void Surface::delete_()
{
    // log
    std::cout << "Destroying surface " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    while (!vertices_.empty())
    {
        disconnect(*vertices_.begin());
    }
    while (!edges_.empty())
    {
        disconnect(*edges_.begin());
    }
    while (!faces_.empty())
    {
        disconnect(*faces_.begin());
    }
    while (!interior_points_.empty())
    {
        disconnect(*interior_points_.begin());
    }

    // log
    std::cout << "---------- surface " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

const int& Surface::get_id() const
{
    return id_;
}

double Surface::compute_point_to_surface_distance(const Eigen::Vector3d& origin, const Eigen::Vector3d& position) const
{
    // if perpendicular, return NaN
    Eigen::Vector3d rayDirection = (position - origin).normalized();
    double distance = (mean_ - position).dot(normal_) / rayDirection.dot(normal_);

    // return
    return distance;
}

double Surface::compute_point_to_surface_distance_with_improved_covariance(const Eigen::Vector3d& origin, const Eigen::Vector3d& position) const
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
    Eigen::Vector3d new_mean = merge_means(mean1, mean2, size1, size2);
    Eigen::Matrix3d new_cov = merge_covariances(cov1, cov2, mean1, mean2, size1, size2);

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

Eigen::Vector3d Surface::compute_point_to_surface_position(const Eigen::Vector3d& origin, const Eigen::Vector3d& point) const
{
    // compute
    Eigen::Vector3d rayDirection = (point - origin).normalized();
    double distance = (mean_ - point).dot(normal_) / rayDirection.dot(normal_);
    Eigen::Vector3d intersection = point + distance * rayDirection;

    // return
    return intersection;
}

void Surface::merge_surface(std::shared_ptr<Surface> surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to merge surface with invalid surface.");

    // merge
    auto surface_valid = surface;
    for (auto vertex : surface_valid->vertices_) connect(vertex);
    for (auto edge : surface_valid->edges_) connect(edge);
    for (auto face : surface_valid->faces_) connect(face);
    for (auto interior_point : surface_valid->interior_points_) connect(interior_point);

    // merge EdgeBVH, since this is maintained per surface
    for (auto edge : surface_valid->edges_)
    {
        edge_bvh_.add_edge(edge);
    }

    // log
    std::cout << "Surface " << surface_valid->get_id() << " merged into surface " << id_ << std::endl;

    // delete
    storage_->delete_surface(surface);
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

const std::tuple<int, int, int>& Surface::get_color() const
{
    return color_;
}

bool Surface::is_expired() const
{
    return is_expired_;
}

void Surface::connect(std::shared_ptr<Vertex> vertex)
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid vertex.");

    // connect
    if (vertices_.insert(vertex).second)
    {
        vertex->connect(shared_from_this());
        add_point_to_surface_fitting(vertex->get_position(), vertex->get_origin());
    }

    // will need to create edges and faces

}

void Surface::connect(std::shared_ptr<Vertex> vertex, std::set<std::shared_ptr<Vertex>> all_nearby_vertices)
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid vertex.");

    // connect
    if (vertices_.insert(vertex).second)
    {
        vertex->connect(shared_from_this());
        add_point_to_surface_fitting(vertex->get_position(), vertex->get_origin());
    }

    // get nearby vertices in the same surface
    std::set<std::shared_ptr<Vertex>> nearby_vertices;
    for (auto nearby_vertex : all_nearby_vertices)
    {
        // check input
        if (nearby_vertex->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid nearby vertex.");

        // skip if same vertex
        if (nearby_vertex == vertex) continue;

        // skip if not in the same surface
        if (nearby_vertex->get_surface() != shared_from_this()) continue;

        // add to nearby vertices
        nearby_vertices.insert(nearby_vertex);
    }

    // create edges
    std::set<std::shared_ptr<Vertex>> used_vertices;
    for (auto nearby_vertex : nearby_vertices)
    {
        // skip if edge is intersected
        if (edge_bvh_.intersect_edges(vertex, nearby_vertex)) continue;

        // create edge
        std::shared_ptr<Edge> new_edge = storage_->add_edge(shared_from_this(), vertex, nearby_vertex);
        used_vertices.insert(nearby_vertex);
    }

    // create faces
    for (std::shared_ptr<Vertex> nearby_vertex0 : used_vertices)
    {
        for (std::shared_ptr<Vertex> nearby_vertex1 : used_vertices)
        {
            // skip if repeated
            if (nearby_vertex1 <= nearby_vertex0) continue;

            // skip if edge does not exist
            bool edge_exist = false;
            std::shared_ptr<Edge> existing_edge;
            for (auto edge : nearby_vertex0->get_edges())
            {
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

            // skip if face have intersections
            // // todo

            // create face
            std::shared_ptr<Face> new_face = storage_->add_face(vertex, nearby_vertex0, nearby_vertex1);
            connect(new_face);
        }
    }
}

void Surface::connect(std::shared_ptr<Edge> edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge->connect(shared_from_this());
}

void Surface::connect(std::shared_ptr<Face> face)
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid face.");

    // connect
    bool inserted = faces_.insert(face).second;
    if (inserted) face->connect(shared_from_this());
}

void Surface::connect(std::shared_ptr<InteriorPoint> interior_point)
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to connect surface with invalid interior point.");

    // connect
    bool inserted = interior_points_.insert(interior_point).second;
    if (inserted) interior_point->connect(shared_from_this());

    // update surface fitting
    if (inserted) add_point_to_surface_fitting(interior_point->get_position(), interior_point->get_origin());
}

void Surface::disconnect(std::shared_ptr<Vertex> vertex)
{
    // check input
    if (vertex->is_expired()) return;

    // disconnect
    bool erased = vertices_.erase(vertex);
    if (erased) vertex->disconnect(shared_from_this());
}

void Surface::disconnect(std::shared_ptr<Edge> edge)
{
    // check input
    if (edge->is_expired()) return;

    // disconnect
    bool erased = edges_.erase(edge);
    if (erased) edge->disconnect(shared_from_this());
}

void Surface::disconnect(std::shared_ptr<Face> face)
{
    // check input
    if (face->is_expired()) return;

    // disconnect
    bool erased = faces_.erase(face);
    if (erased) face->disconnect(shared_from_this());
}

void Surface::disconnect(std::shared_ptr<InteriorPoint> interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    // disconnect
    bool erased = interior_points_.erase(interior_point);
    if (erased) interior_point->disconnect(shared_from_this());
}

void Surface::add_searchable_edge(std::shared_ptr<Edge> edge)
{
    edge_bvh_.add_edge(edge);
}

void Surface::remove_searchable_edge(std::shared_ptr<Edge> edge)
{
    edge_bvh_.delete_edge(edge);
}

void Surface::add_point_to_surface_fitting(Eigen::Vector3d position, Eigen::Vector3d origin)
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
    Eigen::Vector3d new_mean = merge_means(mean1, mean2, size1, size2);
    Eigen::Matrix3d new_cov = merge_covariances(cov1, cov2, mean1, mean2, size1, size2);

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
}

void Surface::set_random_color()
{
    // random color
    color_ = std::make_tuple(rand() % 256, rand() % 256, rand() % 256);
}

bool operator<(const std::shared_ptr<Surface> &lhs, const std::shared_ptr<Surface> &rhs)
{
    // check pointer validity
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs->get_id() < rhs->get_id();
}

bool operator==(const std::shared_ptr<Surface>& lhs, const std::shared_ptr<Surface>& rhs)
{
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
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired surfaces");
    return lhs->get_id() != rhs->get_id();
}

Eigen::Vector3d merge_means(const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, int size1, int size2) 
{
    return (size1 * mean1 + size2 * mean2) / (size1 + size2);
}

Eigen::Matrix3d merge_covariances(const Eigen::Matrix3d& cov1, const Eigen::Matrix3d& cov2, 
                                const Eigen::Vector3d& mean1, const Eigen::Vector3d& mean2, 
                                int size1, int size2) 
{
    // Handle the edge case where one of the sizes is zero
    if (size1 == 0 && size2 == 0) throw std::invalid_argument("Both sizes are zero");
    if (size1 == 0) return cov2;
    if (size2 == 0) return cov1;

    Eigen::Vector3d combined_mean = merge_means(mean1, mean2, size1, size2);
    Eigen::Matrix3d mean_diff1 = (mean1 - combined_mean) * (mean1 - combined_mean).transpose();
    Eigen::Matrix3d mean_diff2 = (mean2 - combined_mean) * (mean2 - combined_mean).transpose();
    Eigen::Matrix3d combined_covariance = (size1 * cov1 + size2 * cov2 + size1 * mean_diff1 + size2 * mean_diff2) / (size1 + size2);

    return combined_covariance;
}