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
    smaller_add_searchable_vertices_queue_.resize(settings_.num_threads);
    smaller_set_of_vertices_to_update_rrs_tree.resize(settings_.num_threads);
    smaller_set_of_faces_to_update_rrs_tree.resize(settings_.num_threads);
    smaller_set_of_edges_to_update_edgeBVH_tree.resize(settings_.num_threads);
    thread_vertices_to_be_deleted_.resize(settings_.num_threads);
    thread_edges_to_be_deleted_.resize(settings_.num_threads);
    thread_faces_to_be_deleted_.resize(settings_.num_threads);
    thread_interior_points_to_be_deleted_.resize(settings_.num_threads);
    thread_vertices_that_have_deleted_publishers_.resize(settings_.num_threads);
    thread_vertices_that_have_added_publishers_.resize(settings_.num_threads);
    thread_vertices_that_have_changed_box_.resize(settings_.num_threads);
    thread_nodes_to_be_deleted_.resize(settings_.num_threads);
    
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

const std::shared_ptr<Vertex>& Storage::add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position, double distance_travelled) 
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), surface, position, origin, distance_travelled);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_vertices_);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius, double distance_travelled) 
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), surface, position, origin, radius, distance_travelled);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_vertices_);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Vertex>& Storage::add_vertex(const std::shared_ptr<Surface>& surface, const std::shared_ptr<GenericPoint>& generic_point)
{
    // create
    std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>();
    vertex->initialize_(shared_from_this(), surface, generic_point);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_vertices_);

    // store
    return *vertices_.insert(vertex).first;
}

const std::shared_ptr<Edge>& Storage::add_edge(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2)
{    
    // create
    std::shared_ptr<Edge> edge = std::make_shared<Edge>();
    edge->initialize_(shared_from_this(), surface, vertex1, vertex2);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_edges_);

    // store
    return *edges_.insert(edge).first;
}

const std::shared_ptr<Face>& Storage::add_face(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3) 
{
    // create
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), surface, vertex1, vertex2, vertex3);

    // write lock 
    std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

    // store
    return *faces_.insert(face).first;
}

const std::shared_ptr<Face>& Storage::add_face(
    const std::shared_ptr<Surface>& surface, 
    const std::shared_ptr<Vertex>& vertex1, 
    const std::shared_ptr<Vertex>& vertex2, 
    const std::shared_ptr<Vertex>& vertex3, 
    const std::shared_ptr<Edge>& edge1, 
    const std::shared_ptr<Edge>& edge2, 
    const std::shared_ptr<Edge>& edge3) 
{
    // create
    std::shared_ptr<Face> face = std::make_shared<Face>();
    face->initialize_(shared_from_this(), surface, vertex1, vertex2, vertex3, edge1, edge2, edge3);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

    // store
    return *faces_.insert(face).first;
} 

const std::shared_ptr<Surface>& Storage::add_surface() 
{
    // create
    std::shared_ptr<Surface> surface = std::make_shared<Surface>();
    surface->initialize_(shared_from_this());

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_surfaces_);
    
    // store
    return *surfaces_.insert(surface).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), position, origin, distance_travelled);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_genertic_points_);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const std::shared_ptr<Vertex>& vertex) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), vertex);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_genertic_points_);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<GenericPoint>& Storage::add_generic_point(const std::shared_ptr<InteriorPoint>& interiror_point) 
{
    // create
    std::shared_ptr<GenericPoint> genertic_point = std::make_shared<GenericPoint>();
    genertic_point->initialize_(shared_from_this(), interiror_point);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_genertic_points_);

    // store
    return *genertic_points_.insert(genertic_point).first;
}

const std::shared_ptr<InteriorPoint>& Storage::add_interior_point(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Face>& face, const std::shared_ptr<GenericPoint>& generic_point) 
{
    // create
    std::shared_ptr<InteriorPoint> interior_point = std::make_shared<InteriorPoint>();
    interior_point->initialize_(shared_from_this(), surface, face, generic_point);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_interior_points_);

    // store
    return *interior_points_.insert(interior_point).first;
}

const std::shared_ptr<InteriorPoint>& Storage::add_interior_point(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex, const std::shared_ptr<GenericPoint>& generic_point) 
{
    // create
    std::shared_ptr<InteriorPoint> interior_point = std::make_shared<InteriorPoint>();
    interior_point->initialize_(shared_from_this(), surface, vertex, generic_point);

    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_interior_points_);

    // store
    return *interior_points_.insert(interior_point).first;
}

// need to ensure the vertex/edge/face are only stored using shared_ptr here and nowhere else
void Storage::delete_vertex(const std::shared_ptr<Vertex> vertex) 
{
    // member delete
    vertex->delete_();    

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_vertices_);
        
        // storage delete
        vertices_.erase(vertex);
    }
}

void Storage::delete_edge(const std::shared_ptr<Edge> edge) 
{
    // member delete
    edge->delete_();

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_edges_);

        // storage delete
        edges_.erase(edge);
    }    
}

void Storage::delete_face(const std::shared_ptr<Face> face) 
{
    // member delete
    face->delete_();

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

        // storage delete
        faces_.erase(face);
    }
}

void Storage::delete_surface(const std::shared_ptr<Surface> surface) 
{
    // member delete
    surface->delete_();

    {    
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_surfaces_);
        
        // storage delete
        surfaces_.erase(surface);
    }    
}

void Storage::delete_generic_point(const std::shared_ptr<GenericPoint> genertic_point) 
{
    // member delete
    genertic_point->delete_();
    
    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_genertic_points_);

        // storage delete
        genertic_points_.erase(genertic_point);
    }
}

void Storage::delete_interior_point(const std::shared_ptr<InteriorPoint> interior_point) 
{
    // member delete
    interior_point->delete_();

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_points_);

        // storage delete
        interior_points_.erase(interior_point);
    }    
}

void Storage::add_to_main_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled)
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), position, origin, distance_travelled);

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

void Storage::print_main_queue_stats()
{
    // print individual point
    const bool print_individual_point = false;
    const bool print_surface_related = true;

    if (print_individual_point)
    {
        unsigned int size = main_queue_.size();
        for (unsigned int i = 0; i < size; i++)
        {
            std::shared_ptr<GenericPoint> generic_point = main_queue_.front();
            main_queue_.pop();
            main_queue_.push(generic_point);

            std::stringstream ss;

            for (const auto& pair : generic_point->contented_surfaces)
            {
                ss << "surface " << pair.first->get_id() << " contented for " << pair.second << " times | ";
            }

            ss << std::endl;
            std::cout << ss.str();
        }
    }

    if (print_surface_related)
    {
        std::unordered_map<std::shared_ptr<Surface>, std::vector<std::shared_ptr<GenericPoint>>, MeshObjectHash> surface_to_generic_points;

        unsigned int num_point_with_no_surface = 0;

        unsigned int size = main_queue_.size();
        for (unsigned int i = 0; i < size; i++)
        {
            std::shared_ptr<GenericPoint> generic_point = main_queue_.front();
            main_queue_.pop();
            main_queue_.push(generic_point);

            std::stringstream ss;

            // skip if no surface
            if (generic_point->contented_surfaces.empty()) 
            {
                num_point_with_no_surface++;
                continue;
            }

            // find the surface with the largest count
            std::pair<std::shared_ptr<Surface>, unsigned int> pair = std::make_pair(nullptr, 0);
            for (const auto& surface : generic_point->contented_surfaces)
            {
                if (surface.second > pair.second)
                {
                    pair = surface;
                }
            }

            // store in the map
            surface_to_generic_points[pair.first].push_back(generic_point);
        }

        // print 
        for (const auto& pair : surface_to_generic_points)
        {
            std::cout << pair.second.size() << "            points failed to lock surface         " << pair.first->get_id() << std::endl;
        }

        std::cout << "Num of point with no surface " << num_point_with_no_surface << std::endl;
    }
}

void Storage::split_main_queue_into_smaller_queues_by_angle(Eigen::Vector3d origin)
{
    // Define angle thresholds for splitting
    double angle_step = 360.0 / settings_.num_threads; // Divide the sphere into equal segments

    // Initialize the smaller queues
    for (auto& queue : smaller_queues_) {
        while (!queue.empty()) {
            queue.pop(); // Ensure the smaller queues are empty before starting
        }
    }

    // Process each point in the main queue
    while (!main_queue_.empty()) {
        // Get the point from the main queue
        std::shared_ptr<GenericPoint> point = main_queue_.front();
        main_queue_.pop();

        // Compute the angle of the point
        Eigen::Vector3d position = point->get_position();
        Eigen::Vector3d direction = position - origin;
        double angle = std::atan2(direction[1], direction[0]) * 180.0 / M_PI;

        // Determine the target queue based on the angle
        int queue_index = static_cast<int>(angle / angle_step) % settings_.num_threads;

        // Add the point to the appropriate smaller queue
        smaller_queues_[queue_index].push(point);
    }
}

void Storage::split_main_queue_into_smaller_queues_by_contention()
{
    // group points by contented surface
    std::unordered_map<std::shared_ptr<Surface>, std::vector<std::shared_ptr<GenericPoint>>, MeshObjectHash> surface_to_generic_points;
    {
        unsigned int size = main_queue_.size();
        for (unsigned int i = 0; i < size; i++)
        {
            // get point from main queue
            std::shared_ptr<GenericPoint> generic_point = main_queue_.front();
            main_queue_.pop();

            // throw if no surface
            if (generic_point->contented_surfaces.empty()) throw std::runtime_error("No surface for generic point.");

            // find the surface with the largest count
            std::pair<std::shared_ptr<Surface>, unsigned int> pair = std::make_pair(nullptr, 0);
            for (const auto& surface : generic_point->contented_surfaces)
            {
                if (surface.second > pair.second)
                {
                    pair = surface;
                }
            }

            // clear contented surface
            generic_point->contented_surfaces.clear();

            // store in the map
            surface_to_generic_points[pair.first].push_back(generic_point);
        }
    }

    // sort by group size, largest first
    std::vector<std::vector<std::shared_ptr<GenericPoint>>> sorted_surface_to_generic_points;
    {
        for (const auto &pair : surface_to_generic_points)
        {
            sorted_surface_to_generic_points.push_back(pair.second);
        }
        std::sort(sorted_surface_to_generic_points.begin(), sorted_surface_to_generic_points.end(), 
            [](const std::vector<std::shared_ptr<GenericPoint>> &a, const std::vector<std::shared_ptr<GenericPoint>> &b)
            { 
                return a.size() > b.size(); 
            });
    }
    
    // allocte the sorted groups to smaller lists
    std::vector<std::vector<std::shared_ptr<GenericPoint>>> smallerLists(settings_.num_threads);
    {
        for (const auto &group : sorted_surface_to_generic_points)
        {
            // Find the smallest list to distribute the group
            unsigned int minListSize = smallerLists[0].size();
            int minIndex = 0;
            for (unsigned int i = 0; i < settings_.num_threads; ++i)
            {
                if (smallerLists[i].size() < minListSize)
                {
                    minListSize = smallerLists[i].size();
                    minIndex = i;
                }
            }

            // Add the group to the smallest list
            smallerLists[minIndex].insert(smallerLists[minIndex].end(), group.begin(), group.end());
        }
    }
    
    // add list to queue
    for (unsigned int i = 0; i < smallerLists.size(); ++i)
    {
        for (const auto &point : smallerLists[i])
        {
            smaller_queues_[i].push(point);
        }
    }
}

void Storage::add_to_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled) 
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), position, origin, distance_travelled);

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

    // if (queue_point->get_num_deletes() <= settings_.num_of_delete_before_put_to_repeated_queue)
    // {
    //     smaller_queues_[omp_get_thread_num()].push(queue_point);
    // }
    // else
    // {
    //     // if queue_point number of delete exceeds 5, 
    //     // reset number of delete and add to repeated_queue_
    //     queue_point->reset_num_deletes();

        smaller_repeated_queues_[omp_get_thread_num()].push(queue_point);
    // }
}

void Storage::add_to_queue(const std::shared_ptr<Vertex>& vertex) 
{
    std::shared_ptr<GenericPoint> queue_point = std::make_shared<GenericPoint>();
    queue_point->initialize_(shared_from_this(), vertex);

    // if (queue_point->get_num_deletes() <= settings_.num_of_delete_before_put_to_repeated_queue)
    // {
    //     smaller_queues_[omp_get_thread_num()].push(queue_point);
    // }
    // else
    // {
    //     // if queue_point number of delete exceeds 5, 
    //     // reset number of delete and add to repeated_queue_
    //     queue_point->reset_num_deletes();

        smaller_repeated_queues_[omp_get_thread_num()].push(queue_point);
    // }
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

unsigned int Storage::get_main_queue_size()
{
    return main_queue_.size();
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

void Storage::clear_queues()
{
    // clear main queue
    while (!main_queue_.empty()) main_queue_.pop();

    for (queue_or_stack<std::shared_ptr<GenericPoint>>& queue : smaller_queues_)
    {
        while (!queue.empty()) queue.pop();
    }

    for (queue_or_stack<std::shared_ptr<GenericPoint>>& repeated_queue : smaller_repeated_queues_)
    {
        while (!repeated_queue.empty()) repeated_queue.pop();
    }

    for (queue_or_stack<std::shared_ptr<GenericPoint>>& abort_queue : smaller_abort_queues_)
    {
        while (!abort_queue.empty()) abort_queue.pop();
    }
}

void Storage::cleanup_surfaces()
{
    // skip if settings is -1
    if (settings_.cleanup_seed_surface_after_ith_cloud == -1) return;

    // make copy of surfaces
    std::vector<std::shared_ptr<Surface>> surfaces_copy;
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(rwlock_surfaces_);

        // copy
        surfaces_copy = std::vector<std::shared_ptr<Surface>>(surfaces_.begin(), surfaces_.end());
    }

    // for each surface
    for (const std::shared_ptr<Surface>& surface : surfaces_copy)
    {
        // skip if surface is expired
        if (surface->is_expired()) continue;

        // skip if modulus is not equal to thread number
        if (surface->get_id() % settings_.num_threads != omp_get_thread_num()) continue;

        // skip if surface recently get updated
        // if (distance_travelled_ - surface->get_max_distance_travelled() < settings_.cleanup_seed_surface_after_distance_travelled) continue;
        if (ith_cloud_ - surface->get_ith_cloud() < settings_.cleanup_seed_surface_after_ith_cloud) continue;

        // // delete if surface is seed
        // if (surface->is_seed()) 
        // {
            delete_surface(surface);
            continue;
        // }
    }

    update_vertices_that_have_added_publishers();
    delete_to_be_deleted_repeatedly();
    update_vertices_that_have_changed_box();
    add_or_remove_vertices_from_rrs_tree();
    add_or_remove_faces_from_bvh_tree();
    add_or_remove_edges_from_edgeBVH_tree();
}

void Storage::remove_non_manifold_edges()
{
    // make copy of edges
    std::vector<std::shared_ptr<Edge>> edges_copy(edges_.begin(), edges_.end());

    // for each edge
    for (const std::shared_ptr<Edge>& edge : edges_copy)
    {
        // delete if edge is non-manifold
        if (edge->is_non_manifold()) delete_edge(edge);
    }
}

void Storage::remove_non_manifold_faces()
{
    // initialize
    std::vector<std::shared_ptr<Face>> non_manifold_faces;
    
    // collect all non manifold faces
    for (const std::shared_ptr<Face>& face : faces_)
    {
        if (face->is_non_manifold()) non_manifold_faces.push_back(face);
    }

    // delete all non manifold faces
    for (const std::shared_ptr<Face>& face : non_manifold_faces)
    {
        // skip if face is expired
        if (face->is_expired()) continue;

        // delete face
        delete_face(face);
    }
}

void Storage::remove_non_manifold_vertices()
{
    // make copy of vertices
    std::vector<std::shared_ptr<Vertex>> vertices_copy(vertices_.begin(), vertices_.end());

    // for each vertex
    for (const std::shared_ptr<Vertex>& vertex : vertices_copy)
    {
        // delete if vertex is non-manifold
        if (vertex->is_non_manifold()) delete_vertex(vertex);
    }
}

void Storage::update_radius()
{
    // make copy of vertices
    std::vector<std::shared_ptr<Vertex>> vertices_copy(vertices_.begin(), vertices_.end());

    for (const std::shared_ptr<Vertex>& vertex : vertices_copy)
    {
        // skip if vertex is expired
        if (vertex->is_expired()) continue;

        // try update radius
        vertex->try_update_radius();
    }
}

void Storage::add_searchable_vertex(const std::shared_ptr<Vertex>& vertex)
{
    rrs_tree_.tree_add_vertex(vertex);
}

void Storage::remove_searchable_vertex(const std::shared_ptr<Vertex>& vertex)
{
    // get node
    std::shared_ptr<RRSNode> node = vertex->node;

    // skip if node does not exist
    if (node == nullptr) return;

    // add node to thread_nodes_to_be_deleted_
    thread_nodes_to_be_deleted_[omp_get_thread_num()].insert(node);

    // remove from rrs_tree
    rrs_tree_.tree_delete_vertex(vertex);
}

void Storage::tree_update_vertex_box(const std::shared_ptr<Vertex>& vertex)
{
    rrs_tree_.tree_update_vertex_box(vertex);
}

void Storage::add_points_in_add_searchable_vertex_queue()
{
    const unsigned int num_points = smaller_add_searchable_vertices_queue_[omp_get_thread_num()].size();

    for (unsigned int i = 0; i < num_points; i++)
    {
        std::shared_ptr<Vertex> vertex = smaller_add_searchable_vertices_queue_[omp_get_thread_num()].front();
        smaller_add_searchable_vertices_queue_[omp_get_thread_num()].pop();

        // add to rrs_tree only add if the vertex is searchable
        if (!vertex->is_searchable()) continue; 

        rrs_tree_.tree_add_vertex(vertex);
    }
}

void Storage::add_or_remove_vertices_from_rrs_tree()
{
    // use while not empty erase idiom
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& vertices_to_update_rrs_tree = smaller_set_of_vertices_to_update_rrs_tree[omp_get_thread_num()];

    // update
    while (!vertices_to_update_rrs_tree.empty())
    {
        // get vertex
        std::shared_ptr<Vertex> vertex = *vertices_to_update_rrs_tree.begin();
        vertices_to_update_rrs_tree.erase(vertex);

        // read lock the vertex
        std::shared_lock<std::shared_mutex> lock(vertex->rwlock_lifecycle_);

        // get flags
        const bool is_expired = vertex->is_expired();
        const bool is_searchable = vertex->is_searchable();

        // make updates according to flags
        if (is_expired || is_searchable)
        {
            rrs_tree_.tree_delete_vertex(vertex);
        }
        else if (!is_expired && !is_searchable)
        {
            rrs_tree_.tree_add_vertex(vertex);
        }
    }
}

void Storage::add_or_remove_faces_from_bvh_tree()
{
    // copy
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_to_update_rrs_tree = smaller_set_of_faces_to_update_rrs_tree[omp_get_thread_num()];

    // clear
    smaller_set_of_faces_to_update_rrs_tree[omp_get_thread_num()].clear();

    // update
    for (const std::shared_ptr<Face>& face : faces_to_update_rrs_tree)
    {
        // read lock the face
        std::shared_lock<std::shared_mutex> lock(face->rwlock_lifecycle_);

        // get flags
        const bool is_expired = face->is_expired();

        // make updates according to flags
        if (is_expired)
        {
            triangle_bvh_.tree_delete_face(face);
        }
        else
        {
            triangle_bvh_.tree_add_face(face);
        }
    }
}

void Storage::add_or_remove_edges_from_edgeBVH_tree()
{
    for (const std::pair<std::shared_ptr<Edge>, std::shared_ptr<Surface>>& edge_surface_pair : smaller_set_of_edges_to_update_edgeBVH_tree[omp_get_thread_num()])
    {
        // at this point in time, is the edge expired?
        bool is_expired;
        {
            // read lock the edge
            std::shared_lock<std::shared_mutex> lock(edge_surface_pair.first->rwlock_lifecycle_);

            // check if edge is expired
            is_expired = edge_surface_pair.first->is_expired();
        }

        // if expired, remove from tree
        if (is_expired)
        {
            edge_surface_pair.second->remove_searchable_edge(edge_surface_pair.first);
        }
        // if not expired, remove from tree
        else
        {
            edge_surface_pair.second->add_searchable_edge(edge_surface_pair.first);
        }
    }

    // clear
    smaller_set_of_edges_to_update_edgeBVH_tree[omp_get_thread_num()].clear();
}

void Storage::add_vertex_to_be_deleted(const std::shared_ptr<Vertex>& vertex)
{
    thread_vertices_to_be_deleted_[omp_get_thread_num()].insert(vertex);
}

void Storage::add_edge_to_be_deleted(const std::shared_ptr<Edge>& edge)
{
    thread_edges_to_be_deleted_[omp_get_thread_num()].insert(edge);
}

void Storage::add_face_to_be_deleted(const std::shared_ptr<Face>& face)
{
    thread_faces_to_be_deleted_[omp_get_thread_num()].insert(face);
}

void Storage::add_interior_point_to_be_deleted(const std::shared_ptr<InteriorPoint>& interior_point)
{
    thread_interior_points_to_be_deleted_[omp_get_thread_num()].insert(interior_point);
}

void Storage::add_surface_to_be_split(const std::shared_ptr<Surface>& surface)
{
    surfaces_to_be_split_.insert(surface);
}

void Storage::clear_surfaces_to_be_split()
{
    surfaces_to_be_split_.clear();
}

void Storage::delete_to_be_deleted_repeatedly()
{
    bool repeat = true;
    while (repeat)
    {
        // delete
        delete_to_be_deleted();

        // update vertices that have deleted publishers
        update_vertices_that_have_deleted_publishers();

        // check if there is any new vertex to be deleted
        repeat = thread_vertices_to_be_deleted_[omp_get_thread_num()].size() > 0
            || thread_edges_to_be_deleted_[omp_get_thread_num()].size() > 0
            || thread_faces_to_be_deleted_[omp_get_thread_num()].size() > 0
            || thread_interior_points_to_be_deleted_[omp_get_thread_num()].size() > 0;
    }
}

void Storage::delete_to_be_deleted()
{
    // use iterators to avoid invalidation issues
    int thread_id = omp_get_thread_num();

    auto& vertex_set = thread_vertices_to_be_deleted_[thread_id];
    auto& edge_set = thread_edges_to_be_deleted_[thread_id];
    auto& face_set = thread_faces_to_be_deleted_[thread_id];
    auto& interior_point_set = thread_interior_points_to_be_deleted_[thread_id];

    while (!vertex_set.empty()) {
        std::shared_ptr<Vertex> vertex = *vertex_set.begin();
        vertex_set.erase(vertex);
        delete_vertex(vertex);
    }

    while (!edge_set.empty()) {
        std::shared_ptr<Edge> edge = *edge_set.begin();
        edge_set.erase(edge);
        delete_edge(edge);
    }

    while (!face_set.empty()) {
        std::shared_ptr<Face> face = *face_set.begin();
        face_set.erase(face);
        delete_face(face);
    }

    while (!interior_point_set.empty()) {
        std::shared_ptr<InteriorPoint> ip = *interior_point_set.begin();
        interior_point_set.erase(ip);
        delete_interior_point(ip);
    }

    // // move content out instead of copying
    // int thread_id = omp_get_thread_num();
    // std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_to_be_deleted = std::move(thread_vertices_to_be_deleted_[thread_id]);
    // std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_to_be_deleted = std::move(thread_edges_to_be_deleted_[thread_id]);
    // std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_to_be_deleted = std::move(thread_faces_to_be_deleted_[thread_id]);
    // std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points_to_be_deleted = std::move(thread_interior_points_to_be_deleted_[thread_id]);

    // // containers are already cleared due to move semantics, no need to call clear()

    // // delete
    // for (const std::shared_ptr<Vertex>& vertex : vertices_to_be_deleted) delete_vertex(vertex);
    // for (const std::shared_ptr<Edge>& edge : edges_to_be_deleted) delete_edge(edge);
    // for (const std::shared_ptr<Face>& face : faces_to_be_deleted) delete_face(face);
    // for (const std::shared_ptr<InteriorPoint>& interior_point : interior_points_to_be_deleted) delete_interior_point(interior_point);


    
    // // make copy
    // std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_to_be_deleted = thread_vertices_to_be_deleted_[omp_get_thread_num()];
    // std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_to_be_deleted = thread_edges_to_be_deleted_[omp_get_thread_num()];
    // std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_to_be_deleted = thread_faces_to_be_deleted_[omp_get_thread_num()];
    // std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points_to_be_deleted = thread_interior_points_to_be_deleted_[omp_get_thread_num()];

    // // clear
    // thread_vertices_to_be_deleted_[omp_get_thread_num()].clear();
    // thread_edges_to_be_deleted_[omp_get_thread_num()].clear();
    // thread_faces_to_be_deleted_[omp_get_thread_num()].clear();
    // thread_interior_points_to_be_deleted_[omp_get_thread_num()].clear();

    // // delete
    // for (const std::shared_ptr<Vertex>& vertex : vertices_to_be_deleted) delete_vertex(vertex);
    // for (const std::shared_ptr<Edge>& edge : edges_to_be_deleted) delete_edge(edge);
    // for (const std::shared_ptr<Face>& face : faces_to_be_deleted) delete_face(face);
    // for (const std::shared_ptr<InteriorPoint>& interior_point : interior_points_to_be_deleted) delete_interior_point(interior_point);
}

void Storage::split_surfaces_per_thread()
{
    for (const std::shared_ptr<Surface>& surface : surfaces_to_be_split_)
    {
        // skip if surface is expired
        if (surface->is_expired()) continue;

        // skip if modulus is not equal to thread number
        if (surface->get_id() % settings_.num_threads != omp_get_thread_num()) continue;

        // split surface
        surface->split_surface_by_connected_components();
    }
}

void Storage::split_surfaces()
{
    for (const std::shared_ptr<Surface>& surface : surfaces_)
    {
        // skip if surface is expired
        if (surface->is_expired()) continue;

        // skip if surface does not have recently removed edges
        // [todo]

        // split surface
        surface->split_surface_by_connected_components();
    }
}

void Storage::add_vertex_that_have_deleted_publishers(const std::shared_ptr<Vertex>& vertex)
{
    thread_vertices_that_have_deleted_publishers_[omp_get_thread_num()].insert(vertex);
}

void Storage::add_vertex_that_have_added_publishers(const std::shared_ptr<Vertex>& vertex)
{
    thread_vertices_that_have_added_publishers_[omp_get_thread_num()].insert(vertex);
}

void Storage::add_vertex_that_have_changed_box(const std::shared_ptr<Vertex>& vertex)
{
    thread_vertices_that_have_changed_box_[omp_get_thread_num()].insert(vertex);
}

void Storage::update_vertices_that_have_deleted_publishers()
{
    // use while not empty erase idiom
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& vertices_that_have_deleted_publishers = thread_vertices_to_be_deleted_[omp_get_thread_num()];

    // update
    while (!vertices_that_have_deleted_publishers.empty())
    {
        // get vertex
        std::shared_ptr<Vertex> vertex = *vertices_that_have_deleted_publishers.begin();
        vertices_that_have_deleted_publishers.erase(vertex);

        // read lock
        std::shared_lock<std::shared_mutex> lock(vertex->rwlock_lifecycle_);

        // update
        vertex->upon_deleting_publisher();
    }
}

void Storage::update_vertices_that_have_added_publishers()
{
    // use while not empty erase idiom
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& vertices_that_have_added_publishers = thread_vertices_that_have_added_publishers_[omp_get_thread_num()];

    // update
    while (!vertices_that_have_added_publishers.empty())
    {
        // get vertex
        std::shared_ptr<Vertex> vertex = *vertices_that_have_added_publishers.begin();
        vertices_that_have_added_publishers.erase(vertex);

        // read lock
        std::shared_lock<std::shared_mutex> lock(vertex->rwlock_lifecycle_);

        // update
        vertex->upon_adding_publisher();
    }
}

void Storage::update_vertices_that_have_changed_box()
{
    // make copy
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_that_have_changed_box = thread_vertices_that_have_changed_box_[omp_get_thread_num()];

    // clear
    thread_vertices_that_have_changed_box_[omp_get_thread_num()].clear();

    // update 
    for (const std::shared_ptr<Vertex>& vertex : vertices_that_have_changed_box)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(vertex->rwlock_lifecycle_);

        // update
        rrs_tree_.tree_update_vertex_box(vertex);
    }
}

void Storage::add_to_set_of_vertices_to_update_rrs_tree(const std::shared_ptr<Vertex>& vertex)
{
    // add to affected vertices set
    smaller_set_of_vertices_to_update_rrs_tree[omp_get_thread_num()].insert(vertex);
}

void Storage::add_to_set_of_faces_to_update_bvh_tree(const std::shared_ptr<Face>& face)
{
    // add to affected vertices set
    smaller_set_of_faces_to_update_rrs_tree[omp_get_thread_num()].insert(face);
}

void Storage::add_to_set_of_edge_to_update_edgeBVH_tree(const std::shared_ptr<Edge>& edge, const std::shared_ptr<Surface>& surface)
{
    // add to affected vertices set
    smaller_set_of_edges_to_update_edgeBVH_tree[omp_get_thread_num()].emplace_back(edge, surface);
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

RRSReturnType Storage::reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& result) 
{
    return rrs_tree_.tree_reverse_radius_search(generic_point, result);
}

BVHReturnType Storage::face_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& result) 
{
    return triangle_bvh_.tree_intersection_search(generic_point, result);
}

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> Storage::get_boundary_vertices() const
{
    // initialize
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> boundary_vertices;

    // get vertices
    const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& vertices = get_vertices();

    for (const std::shared_ptr<Vertex>& vertex : vertices)
    {
        if (vertex->is_boundary()) boundary_vertices.insert(vertex);
    }    

    // return
    return boundary_vertices;
}

std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> Storage::get_boundary_edges() const
{
    // initialize
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> boundary_edges;
    {
        // lock
        std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

        // get edges
        for (const std::shared_ptr<Edge>& edge : edges_)
        {
            // read lock
            std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

            if (edge->is_boundary()) boundary_edges.insert(edge);
        }    
    }
    
    // return
    return boundary_edges;
}

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> Storage::get_vertices() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_vertices_);

    // return
    return vertices_;
}

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& Storage::get_vertices_ref()
{
    // return
    return vertices_;
}

std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> Storage::get_edges() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

    // return
    return edges_;
}

std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> Storage::get_faces() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_faces_);

    // return
    return faces_;
}

std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Storage::get_faces_ref()
{
    // return
    return faces_;
}

std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> Storage::get_surfaces() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_surfaces_);

    // return
    return surfaces_;
}

std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> Storage::get_generic_points() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_genertic_points_);

    // return
    return genertic_points_;
}

std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> Storage::get_interior_points() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_interior_points_);

    // return
    return interior_points_;
}

unsigned int Storage::get_boundary_vertices_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_vertices_);

    // return
    unsigned int size = 0;
    for (const std::shared_ptr<Vertex>& vertex : vertices_)
    {
        if (vertex->is_boundary()) size++;
    }

    return size;
}

unsigned int Storage::get_boundary_edges_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

    // return
    unsigned int size = 0;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->is_boundary()) size++;
    }

    return size;
}

unsigned int Storage::get_vertices_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_vertices_);

    // return
    return vertices_.size();
}

unsigned int Storage::get_edges_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

    // return
    return edges_.size();
}

unsigned int Storage::get_faces_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_faces_);

    // return
    return faces_.size();
}

unsigned int Storage::get_surfaces_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_surfaces_);

    // return
    return surfaces_.size();
}

unsigned int Storage::get_generic_points_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_genertic_points_);

    // return
    return genertic_points_.size();
}

unsigned int Storage::get_interior_points_size() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_interior_points_);

    // return
    return interior_points_.size();
}

std::vector<std::shared_ptr<Vertex>> Storage::get_rrs_vertices()
{
    return rrs_tree_.compute_vertices_list();
}

bool Storage::is_expired() const
{
    return is_expired_;
}

unsigned int Storage::get_rrs_size() const
{
    return rrs_tree_.get_size();
}

unsigned int Storage::get_rrs_node_size() const
{
    return rrs_tree_.get_node_size();
}

unsigned int Storage::get_bvh_size() const
{
    return triangle_bvh_.get_size();
}

void Storage::print_rrs() const
{
    rrs_tree_.tree_print();
}

void Storage::print_bvh() const
{
    triangle_bvh_.tree_print();
}

void Storage::set_distance_travelled(double distance_travelled)
{ 
    distance_travelled_ = distance_travelled; 
}

void Storage::set_ith_cloud(unsigned int ith_cloud)
{ 
    ith_cloud_ = ith_cloud;
}