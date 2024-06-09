#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include <iostream>

#include "MeshObject/InteriorPoint.hpp"

void Face::initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vertex> vertex0, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2)
{
    // check pointer validity
    if (storage.expired()) throw std::runtime_error("Attempts to create face with invalid storage.");
    if (vertex0.expired()) throw std::runtime_error("Attempts to create face with invalid vertex0.");
    if (vertex1.expired()) throw std::runtime_error("Attempts to create face with invalid vertex1.");
    if (vertex2.expired()) throw std::runtime_error("Attempts to create face with invalid vertex2.");

    // store
    storage_ = storage;

    // get id
    id_ = storage_.lock()->get_next_face_id();


    // connect
    connect(vertex0);
    connect(vertex1);
    connect(vertex2);
    std::weak_ptr<Edge> edge0 = storage.lock()->get_edge(vertex0, vertex1);
    std::weak_ptr<Edge> edge1 = storage.lock()->get_edge(vertex1, vertex2);
    std::weak_ptr<Edge> edge2 = storage.lock()->get_edge(vertex2, vertex0);
    connect(edge0);
    connect(edge1);
    connect(edge2);

    // compute center
    Eigen::Vector3d pos0 = vertex0.lock()->get_position();
    Eigen::Vector3d pos1 = vertex1.lock()->get_position();
    Eigen::Vector3d pos2 = vertex2.lock()->get_position();
    center_ = (pos0 + pos1 + pos2) / 3;

    // log
    std::cout << "Face " << id_ << " created between vertex " << vertex0.lock()->get_id() << ", vertex " << vertex1.lock()->get_id() << " and vertex " << vertex2.lock()->get_id() << std::endl;
}

void Face::delete_()
{
    // log
    std::cout << "Destroying face " << id_ << std::endl;

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
    while (!surfaces_.empty())
    {
        disconnect(*surfaces_.begin());
    }
    while (!interior_points_.empty())
    {
        disconnect(*interior_points_.begin());
    }

    // log
    std::cout << "---------- face " << id_ << " destroyed" << std::endl;
}

int Face::get_id() const
{
    return id_;
}

Eigen::Vector3d Face::get_center() const
{
    return center_;
}

std::set<std::weak_ptr<Vertex>> Face::get_vertices() const
{
    return vertices_;
}

std::weak_ptr<Vertex> Face::get_vertex(int index) const
{
    if (index < 0 || index > 2) throw std::runtime_error("Invalid index for vertex.");
    auto it = vertices_.begin();
    for (int i = 0; i < index; i++) it++;
    return *it;
}

std::weak_ptr<Surface> Face::get_surface() const
{
    if (surfaces_.empty()) return std::weak_ptr<Surface>();
    return *surfaces_.begin();
}

bool Face::intersects_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction)
{
    // get first Vertex in vertices_
    auto it = vertices_.begin();
    std::weak_ptr vertex0 = *(it++);
    std::weak_ptr vertex1 = *(it++);
    std::weak_ptr vertex2 = *(it++);
    Eigen::Vector3d v0 = vertex0.lock()->get_position();
    Eigen::Vector3d v1 = vertex1.lock()->get_position();
    Eigen::Vector3d v2 = vertex2.lock()->get_position();


    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    Eigen::Vector3d pvec = direction.cross(edge2);
    double det = edge1.dot(pvec);
    if (std::fabs(det) < EPSILON) return false;

    double invDet = 1.0 / det;

    Eigen::Vector3d tvec = origin - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) return false;
    
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = direction.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) return false;

    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) return false;

    return true;
}   

Eigen::Vector3d Face::compute_intersection_point(const Eigen::Vector3d& origin, const Eigen::Vector3d& direction)
{
    // get first Vertex in vertices_
    auto it = vertices_.begin();
    std::weak_ptr vertex0 = *(it++);
    std::weak_ptr vertex1 = *(it++);
    std::weak_ptr vertex2 = *(it++);
    Eigen::Vector3d v0 = vertex0.lock()->get_position();
    Eigen::Vector3d v1 = vertex1.lock()->get_position();
    Eigen::Vector3d v2 = vertex2.lock()->get_position();

    const double EPSILON = 1e-8;
    Eigen::Vector3d edge1 = v1 - v0;
    Eigen::Vector3d edge2 = v2 - v0;
    
    Eigen::Vector3d pvec = direction.cross(edge2);
    double det = edge1.dot(pvec);
    if (std::fabs(det) < EPSILON) throw std::runtime_error("No intersection point found.");

    double invDet = 1.0 / det;

    Eigen::Vector3d tvec = origin - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0.0 || u > 1.0) throw std::runtime_error("No intersection point found.");
    
    Eigen::Vector3d qvec = tvec.cross(edge1);
    double v = direction.dot(qvec) * invDet;
    if (v < 0.0 || u + v > 1.0) throw std::runtime_error("No intersection point found.");

    double t = edge2.dot(qvec) * invDet;
    if (t < EPSILON) throw std::runtime_error("No intersection point found.");

    return origin + direction * t;
}


void Face::connect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) throw std::runtime_error("Attempts to connect face with invalid vertex.");

    // connect
    bool inserted = vertices_.insert(vertex).second;
    if (inserted) vertex.lock()->connect(shared_from_this());

    // check size
    if (vertices_.size() > 3) throw std::runtime_error("Face connected to more than 3 vertices.");
}

void Face::connect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) throw std::runtime_error("Attempts to connect face with invalid edge.");

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge.lock()->connect(shared_from_this());

    // check size
    if (edges_.size() > 3) throw std::runtime_error("Face connected to more than 3 edges.");
}

void Face::connect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) throw std::runtime_error("Attempts to connect face with invalid surface.");

    // connect
    bool inserted = surfaces_.insert(surface).second;
    if (inserted) surface.lock()->connect(shared_from_this());
}

void Face::connect(std::weak_ptr<InteriorPoint> interior_point)
{
    // check input
    if (interior_point.expired()) throw std::runtime_error("Attempts to connect face with invalid interior point.");

    // connect
    bool inserted = interior_points_.insert(interior_point).second;
    if (inserted) interior_point.lock()->connect(shared_from_this());
}

void Face::disconnect(std::weak_ptr<Vertex> vertex)
{
    // check input
    if (vertex.expired()) return;

    // delete
    bool erased = vertices_.erase(vertex);
    if (erased) vertex.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_.lock()->delete_face(shared_from_this());
}

void Face::disconnect(std::weak_ptr<Edge> edge)
{
    // check input
    if (edge.expired()) return;

    // delete
    bool erased = edges_.erase(edge);
    if (erased) edge.lock()->disconnect(shared_from_this());

    // self destruct
    if (!deleting_) storage_.lock()->delete_face(shared_from_this());
}

void Face::disconnect(std::weak_ptr<Surface> surface)
{
    // check input
    if (surface.expired()) return;

    // delete
    bool erased = surfaces_.erase(surface);
    if (erased) surface.lock()->disconnect(shared_from_this());
}

void Face::disconnect(std::weak_ptr<InteriorPoint> interior_point)
{
    // check input
    if (interior_point.expired()) return;

    // delete
    bool erased = interior_points_.erase(interior_point);
    if (erased) interior_point.lock()->disconnect(shared_from_this());
}

bool operator<(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs)
{
    if (lhs.expired() || rhs.expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}