#include "MeshObject/Storage.hpp"
#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"
#include "MeshObject/Surface.hpp"
#include "MeshObject/GenericPoint.hpp"
#include "MeshObject/InteriorPoint.hpp"

#include "MeshObject/RRSTree.hpp"
#include "MeshObject/TriangleBVH.hpp"

#include <queue>

#include <omp.h>
#include "utilities/queue_or_stack.hpp"

Settings Storage::settings_;

Storage::Storage()
{
    is_expired_ = false;

    // resize to num_threads
    smaller_queues_.resize(settings_.num_threads);
    smaller_repeated_queues_.resize(settings_.num_threads);
    smaller_abort_queues_.resize(settings_.num_threads);
    
    // initialize with queue or stack
    for (size_t i = 0; i < settings_.num_threads; ++i)
    {
        smaller_queues_[i] = queue_or_stack<std::shared_ptr<GenericPoint>>(settings_.use_queue);
        smaller_repeated_queues_[i] = queue_or_stack<std::shared_ptr<GenericPoint>>(settings_.use_queue);
        smaller_abort_queues_[i] = queue_or_stack<std::shared_ptr<GenericPoint>>(settings_.use_queue);
    }
}

Storage::~Storage()
{
    is_expired_ = true;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position) 
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), surface, position, origin);

    // lock
    std::lock_guard<std::mutex> lock(vertices_mutex_);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius)
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), surface, position, origin, radius);

    // lock
    std::lock_guard<std::mutex> lock(vertices_mutex_);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const std::shared_ptr<Surface>& surface, const std::shared_ptr<GenericPoint>& generic_point)
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), surface, generic_point);

    // lock
    std::lock_guard<std::mutex> lock(vertices_mutex_);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Edge>& Storage::add_edge(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{    
    // create
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), vertex1, vertex2);

    // lock
    std::lock_guard<std::mutex> lock(edges_mutex_);

    // store
    return *edges_.insert(edge).first;
}

const std::shared_ptr<Face>& Storage::add_face(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3) 
{
    // create
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), surface, vertex1, vertex2, vertex3);

    // lock 
    std::lock_guard<std::mutex> lock(faces_mutex_);

    // store
    return *faces_.insert(face).first;
}

const std::shared_ptr<Surface>& Storage::add_surface() 
{
    // create
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this());

    // lock
    std::lock_guard<std::mutex> lock(surfaces_mutex_);

    // store
    return *surfaces_.insert(surface).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), position, origin);

    // lock
    std::lock_guard<std::mutex> lock(genertic_points_mutex_);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const std::shared_ptr<Vertex>& vertex) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), vertex);

    // lock
    std::lock_guard<std::mutex> lock(genertic_points_mutex_);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const std::shared_ptr<InteriorPoint>& interiror_point) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), interiror_point);

    // lock
    std::lock_guard<std::mutex> lock(genertic_points_mutex_);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<InteriorPoint>& Storage::add_interior_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin) 
{
    // create
    std::shared_ptr<InteriorPoint> interior_point = std::make_shared<InteriorPoint>();
    interior_point->initialize_(shared_from_this(), position, origin);

    // lock
    std::lock_guard<std::mutex> lock(interior_points_mutex_);

    // store
    return *interior_points_.insert(interior_point).first;
}

const std::shared_ptr<InteriorPoint>& Storage::add_interior_point(const std::shared_ptr<GenericPoint>& generic_point) 
{
    // create
    std::shared_ptr<InteriorPoint> interior_point = std::make_shared<InteriorPoint>();
    interior_point->initialize_(shared_from_this(), generic_point);

    // lock
    std::lock_guard<std::mutex> lock(interior_points_mutex_);

    // store
    return *interior_points_.insert(interior_point).first;
}

// need to ensure the vertex/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vertex(const std::shared_ptr<Vertex>& vertex) 
{
    // check input
    if (vertex->is_expired()) throw std::runtime_error("Attempts to delete expired vertex.");

    {
        // lock
        std::lock_guard<std::mutex> lock(vertices_mutex_);
        
        // storage delete
        vertices_.erase(vertex);
    }
    
    // member delete
    vertex->delete_();    
}

void Storage::delete_edge(const std::shared_ptr<Edge>& edge) 
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to delete expired edge.");

    {
        // lock
        std::lock_guard<std::mutex> lock(edges_mutex_);

        // storage delete
        edges_.erase(edge);
    }
    
    // member delete
    edge->delete_();
}

void Storage::delete_face(const std::shared_ptr<Face>& face) 
{
    // check input
    if (face->is_expired()) return; // face might be already deleted due to reducion in radius

    {
        // lock 
        std::lock_guard<std::mutex> lock(faces_mutex_);

        // storage delete
        faces_.erase(face);
    }

    // member delete
    face->delete_();
}

void Storage::delete_surface(const std::shared_ptr<Surface>& surface) 
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to delete expired surface.");

    {    
        // lock
        std::lock_guard<std::mutex> lock(surfaces_mutex_);
        
        // storage delete
        surfaces_.erase(surface);
    }

    // member delete
    surface->delete_();
}

void Storage::delete_generic_point(const std::shared_ptr<GenericPoint>& genertic_point) 
{
    // check input
    if (genertic_point->is_expired()) throw std::runtime_error("Attempts to delete expired genertic point.");

    // make a copy of the generic point
    std::shared_ptr<GenericPoint> genertic_point_copy = genertic_point;

    {
        // lock
        std::lock_guard<std::mutex> lock(genertic_points_mutex_);

        // storage delete
        genertic_points_.erase(genertic_point);
    }

    // member delete
    genertic_point_copy->delete_();
}

void Storage::delete_interior_point(const std::shared_ptr<InteriorPoint>& interior_point) 
{
    // check input
    if (interior_point->is_expired()) throw std::runtime_error("Attempts to delete expired interior point.");

    {
        // lock
        std::lock_guard<std::mutex> lock(interior_points_mutex_);

        // storage delete
        interior_points_.erase(interior_point);
    }
    
    // member delete
    interior_point->delete_();
}

void Storage::add_to_main_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin) 
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), position, origin);

    main_queue_.push(queue_point);
}

void Storage::add_points_in_smaller_repeated_queues_to_main_queue()
{
    // for each smaller repeated queue
    for (queue_or_stack<std::shared_ptr<GenericPoint>>& smaller_repeated_queue : smaller_repeated_queues_)
    {
        while (!smaller_repeated_queue.empty())
        {
            main_queue_.push(smaller_repeated_queue.get());
            smaller_repeated_queue.pop();
        }
    }
}

void Storage::add_points_in_smaller_abort_queues_to_main_queue()
{
    bool any_non_empty = true;

    // Keep looping until all queues are empty
    while (any_non_empty)
    {
        any_non_empty = false; // Reset flag to check if any queue still has elements

        // Iterate over each smaller abort queue
        for (queue_or_stack<std::shared_ptr<GenericPoint>> &smaller_abort_queue : smaller_abort_queues_)
        {
            // If the queue is not empty, process one element
            if (!smaller_abort_queue.empty())
            {
                // Move the first element of the smaller queue to the main queue
                main_queue_.push(smaller_abort_queue.get());
                smaller_abort_queue.pop(); // Remove the processed element
                any_non_empty = true;      // Indicate that there are still elements to process
            }
        }
    }
}

void Storage::split_main_queue_into_smaller_queues()
{
    // Calculate base number of points per queue
    unsigned int total_points = main_queue_.size();
    unsigned int average = total_points / settings_.num_threads;  // Average points
    unsigned int remainder = total_points % settings_.num_threads;        // Extra points to distribute

    for (unsigned int i = 0; i < settings_.num_threads; ++i) {
        // Calculate the number of points for this thread
        unsigned int num_points = average + (i < remainder ? 1 : 0);

        // Move points from main_queue_ to smaller_queues_[i]
        for (unsigned int j = 0; j < num_points; ++j) {
            if (!main_queue_.empty()) {
                smaller_queues_[i].push(main_queue_.front());  // Move point to the smaller queue
                main_queue_.pop();  // Remove the point from main_queue_
            }
        }
    }
}


void Storage::add_to_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin) 
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), position, origin);

    smaller_queues_[omp_get_thread_num()].push(queue_point);
}

void Storage::add_to_queue(const std::shared_ptr<GenericPoint>& generic_point) 
{
    smaller_queues_[omp_get_thread_num()].push(generic_point);
}

void Storage::add_to_queue(const std::shared_ptr<InteriorPoint>& interior_point) 
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), interior_point);

    if (queue_point->get_num_deletes() <= settings_.num_of_delete_before_put_to_repeated_queue)
    {
        smaller_queues_[omp_get_thread_num()].push(queue_point);
    }
    else
    {
        // if queue_point number of delete exceeds 5, 
        // reset number of delete and add to repeated_queue_
        queue_point->reset_num_deletes();

        smaller_repeated_queues_[omp_get_thread_num()].push(queue_point);
    }
}

void Storage::add_to_queue(const std::shared_ptr<Vertex>& vertex) 
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), vertex);

    if (queue_point->get_num_deletes() <= settings_.num_of_delete_before_put_to_repeated_queue)
    {
        smaller_queues_[omp_get_thread_num()].push(queue_point);
    }
    else
    {
        // if queue_point number of delete exceeds 5, 
        // reset number of delete and add to repeated_queue_
        queue_point->reset_num_deletes();

        smaller_repeated_queues_[omp_get_thread_num()].push(queue_point);
    }
}

void Storage::add_to_abort_queue(const std::shared_ptr<GenericPoint>& generic_point) 
{
    smaller_abort_queues_[omp_get_thread_num()].push(generic_point);
}

std::shared_ptr<GenericPoint> Storage::pop_from_queue()
{
    if (smaller_queues_[omp_get_thread_num()].empty()) return nullptr;

    std::shared_ptr<GenericPoint> queue_point = smaller_queues_[omp_get_thread_num()].get();
    smaller_queues_[omp_get_thread_num()].pop();
    return queue_point;
}

unsigned int Storage::get_queue_size()
{
    return smaller_queues_[omp_get_thread_num()].size();
}

unsigned int Storage::get_repeated_queue_size()
{    
    unsigned int size = 0;
    for (const queue_or_stack<std::shared_ptr<GenericPoint>>& repeated_queue : smaller_repeated_queues_)
    {
        size += repeated_queue.size();
    }
    return size;
}

unsigned int Storage::get_abort_queue_size()
{
    unsigned int size = 0;
    for (const queue_or_stack<std::shared_ptr<GenericPoint>>& abort_queue : smaller_abort_queues_)
    {
        size += abort_queue.size();
    }
    return size;
}

void Storage::add_searchable_vertex(const std::shared_ptr<Vertex>& vertex)
{
    // add to rrs_tree
    rrs_tree_.tree_add_vertex(vertex);
}

void Storage::remove_searchable_vertex(const std::shared_ptr<Vertex>& vertex)
{
    // remove from rrs_tree
    rrs_tree_.tree_delete_vertex(vertex);
}

void Storage::add_searchable_face(const std::shared_ptr<Face>& face)
{
    // add to triangle_bvh
    triangle_bvh_.tree_add_face(face);
}

void Storage::remove_searchable_face(const std::shared_ptr<Face>& face)
{
    // remove from triangle_bvh
    triangle_bvh_.tree_delete_face(face);
}

// get id
int Storage::get_next_vertex_id() { return next_vertex_id_++; }
int Storage::get_next_edge_id() { return next_edge_id_++; }
int Storage::get_next_face_id() { return next_face_id_++; }
int Storage::get_next_surface_id() { return next_surface_id_++; }
int Storage::get_next_generic_point_id() { return next_genertic_point_id_++; }
int Storage::get_next_interior_point_id() { return next_interior_point_id_++; }

bool Storage::can_reverse_radius_search() 
{ 
    return rrs_tree_.can_reverse_radius_search(); 
}

RRSReturnType Storage::reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& result) 
{
    return rrs_tree_.tree_reverse_radius_search(generic_point, result);
}

RRSReturnType Storage::reverse_radius_search_find_node(const Eigen::Vector3d& point, std::shared_ptr<RRSNode>& return_node)
{
    return rrs_tree_.tree_find_leaf_node(point, return_node);
}

BVHReturnType Storage::face_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& result) 
{
    return triangle_bvh_.tree_intersection_search(generic_point, result);
}

BVHReturnType Storage::face_intersection_search_find_node(const Eigen::Vector3d& origin, const Eigen::Vector3d& point, std::shared_ptr<Node>& return_node)
{
    return triangle_bvh_.tree_find_leaf_node(origin, point, return_node);
}

const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& Storage::get_vertices() const
{
    return vertices_;
}

const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Storage::get_edges() const
{
    return edges_;
}

const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Storage::get_faces() const
{
    return faces_;
}

const std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& Storage::get_surfaces() const
{
    return surfaces_;
}

const std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash>& Storage::get_generic_points() const
{
    return genertic_points_;
}

const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& Storage::get_interior_points() const
{
    return interior_points_;
}

std::vector<std::shared_ptr<Vertex>> Storage::get_rrs_vertices()
{
    return rrs_tree_.compute_vertices_list();
}

std::map<std::shared_ptr<Vertex>, int> Storage::get_vertex_to_cloud_indices_map() const
{
    // initialize
    std::map<std::shared_ptr<Vertex>, int> vertex_to_cloud_indices_map;

    // fill
    int id = 0;
    for (const auto& vertex : vertices_)
    {
        vertex_to_cloud_indices_map[vertex] = id;
        id++;
    }

    // return
    return vertex_to_cloud_indices_map;
} 

bool Storage::is_expired() const
{
    return is_expired_;
}

// get edge
const std::shared_ptr<Edge>& Storage::get_edge(std::shared_ptr<Vertex> vertex1, std::shared_ptr<Vertex> vertex2) const
{
    // check input
    if (vertex1->is_expired() || vertex2->is_expired()) throw std::runtime_error("Attempts to get edge with invalid vertex.");

    // search
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->has_vertex(vertex1) && edge->has_vertex(vertex2)) return edge;
    }

    // not found
    throw std::runtime_error("Edge not found.");
}

void Storage::print_rrs() const
{
    rrs_tree_.tree_print();
}

void Storage::print_bvh() const
{
    triangle_bvh_.tree_print();
}

void Storage::check_tree_rebuild()
{
    rrs_tree_.check_rebuild();
    triangle_bvh_.check_rebuild();
}

void Storage::rebuild_tree()
{
    rrs_tree_.rebuild();
    triangle_bvh_.rebuild();
}

unsigned int Storage::get_bvh_size() const
{
    return triangle_bvh_.get_size();
}
