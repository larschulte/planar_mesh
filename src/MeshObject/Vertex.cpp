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
#include "MeshObject/RRSTree.hpp"

Settings Vertex::settings_;

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius, double distance_travelled)
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
    distance_travelled_ = distance_travelled;
    direction_ = (position_ - origin_).normalized();

    // num of deletes
    num_deletes_ = 0;

    // connect
    connect(surface);

    // check if update search tree
    check_if_update_search_tree();

    // set reverse search radius based on input parameter
    try_update_radius();
    try_update_node_box();

    // log
    if (settings_.log.initialize) std::cout << "Vertex " << id_ << " created.\n";
}

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<GenericPoint>& generic_point)
{
    initialize_(storage, surface, generic_point->get_position(), generic_point->get_origin(), generic_point->get_radius(), generic_point->get_distance_travelled());
    num_deletes_ = generic_point->get_num_deletes();
    projected_uncertainty_ = generic_point->get_projected_uncertainty();
}

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled)
{
    std::shared_ptr<GenericPoint> generic_point = storage->add_generic_point(position, origin, distance_travelled);
    initialize_(storage, surface, generic_point);
    storage->delete_generic_point(generic_point);
}

// [todo]
// i initially wanted each vertex to be strongly bounded with node, so that i can use node's lock to present vertex's lock
// however, since now each vertex may not get a node after process point
// or since each vertex may not disconnect from node after process point

// it shouldn't hurt to have a lock for each vertex
// this prevents other thread from trying to read the vertex while it is being deleted

// when adding vertex to search tree, we can postpone and add it to the tree after process point
// when removing vertex from search tree, there are two causes
    // when deleting vertex
        // thus the surface will be nullptr
        // some check need to be done before accessing its surface during search
    // when vertex is not boundary
        // the surface will not be null ptr

void Vertex::delete_()
{
    // write lock
    std::unique_lock<std::shared_mutex> lock(rwlock_lifecycle_);

    // add to update rrs tree vertex set
    storage_->add_to_set_of_vertices_to_update_rrs_tree(shared_from_this());

    // log
    if (settings_.log.deletion) std::cout << "Destroying vertex " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;
    
    // // cascade radius reduction
    // cascade_radius_reduction_to_connected_vertices();

    // disconnect
    delete_publishers();
    delete_subscribers();
    std::vector<std::shared_ptr<Edge>> edges = edges_;
    std::vector<std::shared_ptr<Face>> faces = faces_;
    // std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices_that_affect_radius = neighboring_vertices_that_affect_radius_;
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& face : faces) disconnect(face);
    // for (const auto& neighboring_vertex : neighboring_vertices_that_affect_radius) disconnect_neighboring_vertex(neighboring_vertex);

    std::shared_ptr<Surface> surface = surface_; // make copy to prevent cyclic reference and nullptr access
    if (surface) disconnect(surface);

    // update delete count
    num_deletes_++;

    // only create penetrated point / generic point if sibling is empty
    if (can_create_generic_point_)
    {
        storage_->add_to_queue(shared_from_this());
    }

    // log
    if (settings_.log.deletion) std::cout << "---------- vertex " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

void Vertex::temp_initialize(const Eigen::Vector3d& position, unsigned int id)
{
    // set expired
    is_expired_ = false;

    // set id
    id_ = id;

    // set position
    position_ = position;

    // set radius
    try_update_radius();
    try_update_node_box();
}

const int& Vertex::get_id() const 
{ 
    return id_; 
}

const Eigen::Vector3d& Vertex::get_original_position() const 
{ 
    return position_; 
}

const Eigen::Vector3d& Vertex::get_position() const 
{ 
    if (projected_position_.isZero()) throw std::runtime_error("Vertex projected position is not set.");
    
    return projected_position_;
}

const Eigen::Vector3d& Vertex::buffer_compute_projected_position(const std::shared_ptr<Surface> surface)
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

const double& Vertex::buffer_compute_projected_distance(const std::shared_ptr<Surface> surface)
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

const Eigen::Vector3d& Vertex::buffer_compute_projected_position() { return buffer_compute_projected_position(get_surface()); }
const double& Vertex::buffer_compute_projected_distance() { return buffer_compute_projected_distance(get_surface()); }

const Eigen::Vector3d& Vertex::get_origin() const 
{ 
    return origin_; 
}

const double& Vertex::get_distance_travelled() const 
{ 
    return distance_travelled_; 
}

const Eigen::Vector3d& Vertex::get_direction() const 
{ 
    return direction_; 
}

const std::shared_ptr<Surface>& Vertex::get_surface() const
{    
    return surface_;
}

const std::shared_ptr<Surface>& Vertex::get_surface_check() const
{    
    return surface_;
}

bool Vertex::has_surface() const
{
    return surface_ != nullptr;
}

const std::vector<std::shared_ptr<Edge>>& Vertex::get_edges() const 
{ 
    return edges_; 
}

const std::vector<std::shared_ptr<Face>>& Vertex::get_faces() const 
{ 
    return faces_; 
}

std::size_t Vertex::get_num_deletes() const
{
    return num_deletes_;
}

double Vertex::get_current_surface_uncertainty() const
{
    return current_surface_uncertainty_;
}

double& Vertex::get_projected_uncertainty()
{
    return projected_uncertainty_;
}

const std::shared_ptr<Edge>& Vertex::get_edge(const std::shared_ptr<Vertex>& vertex) const
{
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->has_vertex(vertex)) return edge;
    }

    throw std::runtime_error("Edge not found.");
}

// void Vertex::try_merge_surfaces()
// {
    // // check if there is only one surface
    // if (surfaces_.size() == 1) return;

    // bool merge_happened = false;

    // // merge surfaces
    // while (true)
    // {
    //     bool merge_again = false;

    //     // generate pairs of surfaces
    //     std::set<std::pair<std::shared_ptr<Surface>, std::shared_ptr<Surface>>> surface_pairs;
    //     for (std::shared_ptr<Surface> surface1 : surfaces_) 
    //     {
    //         for (std::shared_ptr<Surface> surface2 : surfaces_) 
    //         {
    //             if (surface1 >= surface2) continue;
    //             surface_pairs.insert(std::make_pair(surface1, surface2));
    //         }
    //     }

    //     // merge surfaces
    //     for (const auto& pairs : surface_pairs) 
    //     {
    //         // skip if combined surface have large eigenvalue
    //         const std::shared_ptr<Surface>& surface1 = pairs.first;
    //         const std::shared_ptr<Surface>& surface2 = pairs.second;
    //         const Eigen::Matrix3d& cov1 = surface1->get_covariance();
    //         const Eigen::Matrix3d& cov2 = surface2->get_covariance();
    //         const Eigen::Vector3d& mean1 = surface1->get_mean();
    //         const Eigen::Vector3d& mean2 = surface2->get_mean();
    //         int size1 = surface1->get_total_point_size();
    //         int size2 = surface2->get_total_point_size();
    //         Eigen::Matrix3d covariance_matrix = merge_covariance(cov1, cov2, mean1, mean2, size1, size2);
    //         Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance_matrix);
            
    //         double eigenvalue = solver.eigenvalues()[0];
    //         if (eigenvalue > Vertex::settings_.merged_eigenvalue_threshold) continue;

    //         // merge by changing surface1 into surface2
    //         std::cout << ">> Merging surface " << surface1->get_id() << " with " << surface1->get_total_point_size() << " points into surface " << surface2->get_id() << " with " << surface2->get_total_point_size() << " points." << std::endl;
    //         swap(surface1, surface2);
    //         merge_happened = true;
    //         std::cout << ">> resultant surface " << surface1->get_id() << " has " << surface1->get_total_point_size() << " points." << std::endl;
    //         std::cout << ">> resultant surface " << surface2->get_id() << " has " << surface2->get_total_point_size() << " points." << std::endl;

    //         // reset
    //         merge_again = true;
    //         break;
    //     }

    //     if (!merge_again) break;
    // }

    // if (merge_happened) 
    // {
    //     std::cout << ">> Merging surfaces done." << std::endl;
    // }
    // else
    // {
    //     std::cout << ">> No merging happened." << std::endl;
    // }
// }

const Eigen::Vector2d& Vertex::get_surface_coordinate(const std::shared_ptr<Surface>& surface)
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

const Eigen::Vector2d& Vertex::get_surface_coordinate()
{
    return get_surface_coordinate(get_surface());
}

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> Vertex::compute_connected_vertices()
{
    // iterate through all edges and get connected vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_vertices;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        // skip if edge is deleting
        if (edge->is_deleting()) continue;
        
        connected_vertices.insert(edge->get_vertex(0));
        connected_vertices.insert(edge->get_vertex(1));
    }

    // remove itself
    connected_vertices.erase(shared_from_this());

    // return
    return connected_vertices;
}

std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> Vertex::compute_connected_interior_points()
{
    // iterate through all faces and get connected interior points
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> connected_interior_points;
    for (const std::shared_ptr<Face>& face : faces_)
    {
        connected_interior_points.insert(face->get_interior_points().begin(), face->get_interior_points().end());
    }

    // return
    return connected_interior_points;
}

bool Vertex::is_expired() const
{
    return is_expired_;
}

bool Vertex::is_boundary() const
{
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

    // becomes boundary when one of the connected edges is boundary, or when the point is alone
    bool is_boundary_flag = false;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->is_boundary())
        {
            is_boundary_flag = true;
            break;
        }
    }
    if (edges_.empty()) is_boundary_flag = true;

    return is_boundary_flag;
}

bool Vertex::is_searchable() const
{
    return node != nullptr;
}

bool Vertex::is_deleting() const
{
    return deleting_;
}

void Vertex::connect(const std::shared_ptr<Edge>& edge)
{
    // check input
    if (edge->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid edge.");

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_edges_);

        // skip if already connected
        for (const std::shared_ptr<Edge>& edge_ : edges_) if (edge_ == edge) return;    

        // connect
        edges_.push_back(edge);

    }
    
    check_if_update_search_tree();

    // reverse connect
    edge->connect(shared_from_this());
}

void Vertex::connect(const std::shared_ptr<Face>& face) 
{
    // check input
    if (face->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid face.");

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

        // skip if already connected
        for (const std::shared_ptr<Face>& face_ : faces_) if (face_ == face) return;

        // connect
        faces_.push_back(face);
    }
    
    update_singular_state();

    // reverse connect
    face->connect(shared_from_this());
}

void Vertex::connect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) throw std::runtime_error("Attempts to connect vertex with invalid surface.");

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_surface_);

        // skip if already connected
        if (surface_ == surface) return;

        // connect
        surface_ = surface;
    }
        
    // update projected position
    projected_position_ = surface->compute_point_projective_position(get_origin(), get_original_position());

    // reverse connect
    surface->connect(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Edge>& edge) 
{
    // check input
    if (edge->is_expired()) return;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_edges_);

        // skip if not connected
        auto it = std::find(edges_.begin(), edges_.end(), edge);
        if (it == edges_.end()) return;

        // disconnect
        edges_.erase(it);
        
    }

    // reverse disconnect
    edge->disconnect(shared_from_this());

    // update boundary state
    check_if_update_search_tree();

    // check self destruct
    if (!deleting_ && edges_.empty() && can_self_destruct_) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Face>& face)
{
    // check pointer validity
    if (face->is_expired()) return;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_faces_);

        // skip if not connected
        auto it = std::find(faces_.begin(), faces_.end(), face);
        if (it == faces_.end()) return;

        // disconnect
        faces_.erase(it);
    }

    // reverse disconnect
    face->disconnect(shared_from_this());

    // do not self destruct when have no face
    // check self destruct
    if (!deleting_ && faces_.empty() && can_self_destruct_) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Surface>& surface)
{
    // check input
    if (surface->is_expired()) return;

    {
        // write lock
        std::unique_lock<std::shared_mutex> lock(rwlock_surface_);

        // skip if not connected
        if (surface_ != surface) return;

        // disconnect
        surface_ = nullptr;
    }

    // update projected position
    projected_position_ = Eigen::Vector3d::Zero();

    // reverse disconnect
    surface->disconnect(shared_from_this());

    // check self destruct
    if (!deleting_ && can_self_destruct_) storage_->delete_vertex(shared_from_this());
}

std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> Vertex::get_connected_boundary_edges() const
{
    // initialize
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> connected_boundary_edges;

    // get boundary edges
    for (const auto& edge : get_edges())
    {
        if (edge->is_boundary()) connected_boundary_edges.insert(edge);
    }

    // return
    return connected_boundary_edges;
}

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> Vertex::get_connected_boundary_vertices()
{
    // get connected boundary edges
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> connected_boundary_edges = get_connected_boundary_edges();

    // get connected boundary vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_boundary_vertices;
    for (const auto& edge : connected_boundary_edges)
    {
        connected_boundary_vertices.insert(edge->get_vertex(0));
        connected_boundary_vertices.insert(edge->get_vertex(1));
    }
    connected_boundary_vertices.erase(shared_from_this()); // remove self

    // return
    return connected_boundary_vertices;
}

bool Vertex::check_connected_by_edge(const std::shared_ptr<Vertex>& vertex)
{
    // check if connected by edge
    for (const auto& edge : edges_)
    {
        if (edge->get_vertex(0) == vertex || edge->get_vertex(1) == vertex) return true;
    }

    // return
    return false;
}

bool Vertex::check_connected_by_face(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1)
{
    // check if connected by face
    for (const auto& face : faces_)
    {
        if (face->has_vertex(vertex0) && face->has_vertex(vertex1)) return true;
    }

    // return
    return false;
}

bool Vertex::try_close_holes_repeatedly()
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

bool Vertex::try_close_holes_between_self_and(std::shared_ptr<Vertex>& vertex0, std::shared_ptr<Vertex>& vertex1) 
{
    // initialize flag
    bool changed = false;

    // skip if edge is not short enough
    const double edge_length = (vertex0->get_position() - vertex1->get_position()).norm();
    const double radius0 = vertex0->get_radius();
    const double radius1 = vertex1->get_radius();
    if (!settings_.edge_is_short_enough(edge_length, radius0, radius1)) return changed;

    // skip if edge intersects
    if (get_surface()->tree_intersect_edge(vertex0, vertex1)) return changed;

    // create edge if edge does not exist
    const bool edge_exist = vertex0->check_connected_by_edge(vertex1);
    if (!edge_exist)
    {
        // create edge
        std::shared_ptr<Edge> new_edge = storage_->add_edge(vertex0, vertex1);
        get_surface()->connect(new_edge);

        // update flag
        changed = true;
    }

    // skip if edge between 0 and 1 is not boundary
    if (!vertex0->get_edge(vertex1)->is_boundary()) return changed;

    // create face if face does not exist
    const bool face_exist = check_connected_by_face(vertex0, vertex1);
    if (!face_exist)
    {
        // create face
        std::shared_ptr<Face> new_face = storage_->add_face(get_surface(), shared_from_this(), vertex0, vertex1);
        get_surface()->connect(new_face);

        // un-add the face if face is non-manifold
        if (new_face->is_non_manifold())
        {
            // un-add face
            new_face->un_add_face();

            // update flag
            changed = false;
        }
        else
        {
            // update flag
            changed = true;
        }
    }

    // return changed
    return changed;
}

bool Vertex::try_close_holes()
{
    // initialize
    bool changed = false;

    // initialize
    bool repeat_loop = true;
    while (repeat_loop)
    {
        repeat_loop = false;

        // get connected boundary vertices
        std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_boundary_vertices = get_connected_boundary_vertices();

        // skip if less than two connected boundary vertex (this includes when the vertex itself is not boundary)
        if (connected_boundary_vertices.size() < 2)
        {
            break;
        }

        // compute connected boundary vertices with angle
        std::vector<std::pair<std::shared_ptr<Vertex>, double>> connected_boundary_vertices_with_angle;
        for (const auto& vertex : connected_boundary_vertices)
        {
            // use projected position
            Eigen::Vector2d direction = vertex->get_surface_coordinate() - get_surface_coordinate();
            double angle = std::atan2(direction.y(), direction.x());

            // add to list
            connected_boundary_vertices_with_angle.push_back(std::make_pair(vertex, angle));
        }

        // sort by angle
        std::sort(connected_boundary_vertices_with_angle.begin(), connected_boundary_vertices_with_angle.end(), [](const std::pair<std::shared_ptr<Vertex>, double>& a, const std::pair<std::shared_ptr<Vertex>, double>& b) { return a.second < b.second; });

        // create pairs of connected boundary vertices
        std::vector<std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>>> connected_boundary_vertex_pairs;
        for (int i = 0; i < connected_boundary_vertices_with_angle.size(); i++)
        {
            connected_boundary_vertex_pairs.push_back(std::make_pair(connected_boundary_vertices_with_angle[i].first, connected_boundary_vertices_with_angle[(i + 1) % connected_boundary_vertices_with_angle.size()].first));
        }

        // try close holes between the two vertices
        for (auto& [vertex0, vertex1] : connected_boundary_vertex_pairs)
        {
            if (try_close_holes_between_self_and(vertex0, vertex1)) 
            {
                changed = true;
                repeat_loop = true;
                break;
            }
        }
    }

    return changed;
}

void Vertex::remove_all_edges()
{
    // make copy of edges
    std::vector<std::shared_ptr<Edge>> edges_copy = edges_;

    // remove all edges
    for (const auto& edge : edges_copy)
    {
        disconnect(edge);
    }
}

void Vertex::set_can_self_destruct(bool can_self_destruct)
{
    can_self_destruct_ = can_self_destruct;
}

bool Vertex::is_non_manifold() const
{
    // non manifold if 
    // 1. connected by one boundary edge
    // 2. connected by more than 2 boundary edges
    return get_connected_boundary_edges().size() > 2 || get_connected_boundary_edges().size() == 1;
}

void Vertex::delete_publishers()
{
    // vertex point publisher
    std::vector<std::pair<std::shared_ptr<Vertex>, double>> vertex_point_distance_publishers_copy = vertex_point_distance_publishers_;
    for (const auto& [vertex_point_publisher, distance] : vertex_point_distance_publishers_copy)
    {
        delete_vertex_point_distance_publisher(vertex_point_publisher);
    }

    // interior point publisher
    std::vector<std::pair<std::shared_ptr<InteriorPoint>, double>> interior_point_distance_publishers_copy = interior_point_distance_publishers_;
    for (const auto& [interior_point_publisher, distance] : interior_point_distance_publishers_copy)
    {
        delete_interior_point_distance_publisher(interior_point_publisher);
    }
}

void Vertex::delete_subscribers()
{
    // vertex point subscribers
    std::vector<std::shared_ptr<Vertex>> vertex_point_distance_subscribers_copy = vertex_point_distance_subscribers_;
    for (const auto& vertex_point_subscriber : vertex_point_distance_subscribers_copy)
    {
        // delete
        delete_vertex_point_distance_subscriber(vertex_point_subscriber);
    }
}

void Vertex::upon_adding_publisher()
{
    // previous and current radius
    const double previous_radius = get_radius();
    try_update_radius();
    const double current_radius = get_radius();
    
    // only try break edge and update node box if radius is reduced
    if (current_radius < previous_radius) 
    {
        try_break_edges();
        if (!is_expired()) try_update_node_box();
    }    
}

void Vertex::upon_deleting_publisher()
{
    // skip if deleting
    if (is_deleting()) return;

    // previous and current radius
    const double previous_radius = get_radius();
    try_update_radius();
    const double current_radius = get_radius();

    // only try close holes and update node box if radius is increased
    if (current_radius > previous_radius) 
    {
        try_close_holes_repeatedly();
        try_update_node_box();
    }
}

void Vertex::add_vertex_point_distance_publisher(const std::shared_ptr<Vertex> vertex_point_publisher)
{
    // check input
    if (vertex_point_publisher->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_vertex_point_distance_publishers_);

        // skip if already exist
        for (const auto& pair : vertex_point_distance_publishers_) if (pair.first == vertex_point_publisher) return; // Already exists

        // compute distance
        const double distance = (get_position() - vertex_point_publisher->get_position()).norm(); 

        // add publisher
        vertex_point_distance_publishers_.emplace_back(vertex_point_publisher, distance);
    }

    // add self to publisher vertex as subscriber
    vertex_point_publisher->add_vertex_point_distance_subscriber(shared_from_this());
}

void Vertex::delete_vertex_point_distance_publisher(const std::shared_ptr<Vertex> vertex_point_publisher)
{
    // check input
    if (vertex_point_publisher->is_expired()) return;
    
    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_vertex_point_distance_publishers_);

        // skip if not exist
        auto it = std::find_if(vertex_point_distance_publishers_.begin(), vertex_point_distance_publishers_.end(), [&](const std::pair<std::shared_ptr<Vertex>, double>& pair) {return pair.first == vertex_point_publisher;});
        if (it == vertex_point_distance_publishers_.end()) return;

        // delete publisher
        vertex_point_distance_publishers_.erase(it);
    }
    
    // upon deleting publisher
    upon_deleting_publisher();

    // delete self from publisher vertex as subscriber
    vertex_point_publisher->delete_vertex_point_distance_subscriber(shared_from_this());
}

void Vertex::add_vertex_point_distance_subscriber(const std::shared_ptr<Vertex> vertex_point_subscriber)
{
    // check input
    if (vertex_point_subscriber->is_expired()) return;
    
    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_vertex_point_distance_subscribers_);

        // skip if already exist
        for (const auto& vertex_point_subscriber_ : vertex_point_distance_subscribers_) if (vertex_point_subscriber_ == vertex_point_subscriber) return; // Already exists

        // add subscriber
        vertex_point_distance_subscribers_.push_back(vertex_point_subscriber);
    }

    // add self to subscriber vertex as publisher
    vertex_point_subscriber->add_vertex_point_distance_publisher(shared_from_this());
    vertex_point_subscriber->upon_adding_publisher();
}

void Vertex::delete_vertex_point_distance_subscriber(const std::shared_ptr<Vertex> vertex_point_subscriber)
{
    // check input
    if (vertex_point_subscriber->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_vertex_point_distance_subscribers_);

        // skip if not exist
        auto it = std::find(vertex_point_distance_subscribers_.begin(), vertex_point_distance_subscribers_.end(), vertex_point_subscriber);
        if (it == vertex_point_distance_subscribers_.end()) return; // skip if not exist

        // delete subscriber
        vertex_point_distance_subscribers_.erase(it);
    }
    
    // delete self from subscriber vertex as publisher
    vertex_point_subscriber->delete_vertex_point_distance_publisher(shared_from_this());
}

void Vertex::add_interior_point_distance_publisher(const std::shared_ptr<InteriorPoint> interior_point_publisher)
{
    // check input
    if (interior_point_publisher->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_publishers_);

        // skip if already exist
        for (const auto& pair : interior_point_distance_publishers_) if (pair.first == interior_point_publisher) return; // Already exists
    
        // compute distance
        const double distance = (get_position() - interior_point_publisher->get_position()).norm(); 

        // add publisher
        interior_point_distance_publishers_.emplace_back(interior_point_publisher, distance);
    }

    // add self to publisher vertex as subscriber
    interior_point_publisher->add_interior_point_distance_subscriber(shared_from_this());
}

void Vertex::delete_interior_point_distance_publisher(const std::shared_ptr<InteriorPoint> interior_point_publisher)
{
    // check input
    if (interior_point_publisher->is_expired()) return;

    {
        // lock
        std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_publishers_);

        // skip if not exist
        auto it = std::find_if(interior_point_distance_publishers_.begin(), interior_point_distance_publishers_.end(), [&](const std::pair<std::shared_ptr<InteriorPoint>, double>& pair) {return pair.first == interior_point_publisher;});
        if (it == interior_point_distance_publishers_.end()) return;

        // delete publisher
        interior_point_distance_publishers_.erase(it);
    }
    
    // upon deleting publisher
    upon_deleting_publisher();

    // delete self from publisher vertex as subscriber
    interior_point_publisher->delete_interior_point_distance_subscriber(shared_from_this());
}

double Vertex::compute_radius()
{
    // reset to default
    double new_radius = settings_.radius_value;

    // remove expired entries
    for (auto it = vertex_point_distance_publishers_.begin(); it != vertex_point_distance_publishers_.end();)
    {
        if (it->first->is_expired())
        {
            it = vertex_point_distance_publishers_.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // remove expired entries
    for (auto it = interior_point_distance_publishers_.begin(); it != interior_point_distance_publishers_.end();)
    {
        if (it->first->is_expired())
        {
            it = interior_point_distance_publishers_.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // reduce value
    for (const auto& [neighboring_vertex, distance] : vertex_point_distance_publishers_)
    {
        const double extra_radius = settings_.extra_radius;
        if (distance + extra_radius < new_radius) new_radius = distance + extra_radius;
    }

    // reduce value
    for (const auto& [interior_point, distance] : interior_point_distance_publishers_)
    {
        if (distance < new_radius) new_radius = distance;
    }

    // return radius
    return new_radius;
}

void Vertex::try_update_radius()
{
    reverse_search_radius_ = compute_radius();
}

void Vertex::try_break_edges()
{
    // collect edges to delete
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_to_delete;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->is_expired()) continue; // could turn expired if below deletes an edge which then deletes a face
        if (edge->is_deleting()) continue; // could be deleting

        const double edge_length = edge->get_length();
        const double radius0 = edge->get_vertex(0)->get_radius();
        const double radius1 = edge->get_vertex(1)->get_radius();
        if (!settings_.edge_is_short_enough(edge_length, radius0, radius1))
        {
            edges_to_delete.insert(edge);
        }
    }

    // delete edges
    for (const std::shared_ptr<Edge>& edge : edges_to_delete)
    {
        // skip if expired
        if (is_expired()) return;

        // skip if edge is expired
        if (edge->is_expired()) continue; 

        storage_->delete_edge(edge);
    }
}

void Vertex::try_update_node_box()
{
    // previous radius
    const double previous_rrs_half_size = (max_.x() - min_.x()) / 2.0;

    // update min and max
    const double current_rrs_half_size = settings_.compute_rrs_half_size(reverse_search_radius_);
    const Eigen::Vector3d half_size_vecotr = Eigen::Vector3d(current_rrs_half_size, current_rrs_half_size, current_rrs_half_size);
    min_ = get_position() - half_size_vecotr;
    max_ = get_position() + half_size_vecotr;

    if (node)
    {
        // update node
        if (current_rrs_half_size > previous_rrs_half_size)
        {
            node->box_ = RRSBoundingBox(min_, max_);
            node->recursive_expand_parent_box();
        }
        else if (current_rrs_half_size < previous_rrs_half_size)
        {
            node->box_ = RRSBoundingBox(min_, max_);
            node->recursive_shrink_parent_box();
        }
    }
}

void Vertex::update_singular_state()
{
    // count number of faces in this surface
    int num_faces_in_surface = faces_.size();

    // update singular state
    if (num_faces_in_surface == 0) is_singular_ = true;
    else is_singular_ = false;
}

bool Vertex::is_singular() const
{
    return is_singular_;
}

// swap surface1 with surface2
void Vertex::swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2)
{
    // if contains surface1
    if (surface_ == surface1)
    {
        // std::cout << "Swapping vertex " << id_ << " surface " << surface1->get_id() << " with surface " << surface2->get_id() << std::endl;

        can_self_destruct_ = false;
        disconnect(surface1);
        connect(surface2);
        can_self_destruct_ = true;

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

// replace this vertex with the input vertex
void Vertex::absorbs(const std::shared_ptr<Vertex>& input_vertex)
{
    // change input_vertex's surface to this vertex's surface
    std::shared_ptr<Surface> current_surface = input_vertex->get_surface();
    std::shared_ptr<Surface> new_surface = surface_;
    input_vertex->swap(current_surface, new_surface);

    // make input_vertex's edges and faces connect to this vertex
    std::vector<std::shared_ptr<Face>> faces_copy = input_vertex->get_faces();
    std::vector<std::shared_ptr<Edge>> edges_copy = input_vertex->get_edges();
    for (const std::shared_ptr<Face>& face : faces_copy) face->swap(input_vertex, shared_from_this());            
    for (const std::shared_ptr<Edge>& edge : edges_copy) edge->swap(input_vertex, shared_from_this());
}

void Vertex::check_if_update_search_tree()
{
    // check if is boundary
    if (is_boundary() && !is_searchable())
    {
        storage_->add_to_set_of_vertices_to_update_rrs_tree(shared_from_this());
    }
    else if (!is_boundary() && is_searchable())
    {
        storage_->add_to_set_of_vertices_to_update_rrs_tree(shared_from_this());
    }
}

void Vertex::print_info()
{
    std::cout << "Vertex " << id_ << " at " << position_.transpose() << std::endl;
    std::cout << "Boundary state: " << is_boundary() << std::endl;
    std::cout << "Singular state: " << is_singular_ << std::endl;
    std::cout << "Searchable state: " << is_searchable() << std::endl;
    std::cout << "Expired: " << is_expired_ << std::endl;
}

void Vertex::can_create_generic_point(bool state)
{
    can_create_generic_point_ = state;
}

const Eigen::Vector3d& Vertex::get_min() const
{
    return min_;
}

const Eigen::Vector3d& Vertex::get_max() const
{
    return max_;
}

const double& Vertex::get_radius() const
{
    return reverse_search_radius_;
}

const double& Vertex::get_radius(const std::shared_ptr<Surface>& surface) const
{
    // only used when being connected to a surface by edge and faces
    return reverse_search_radius_;
}

bool Vertex::contains(const Eigen::Vector3d& point) const
{
    return (point - get_position()).norm() < settings_.compute_rrs_half_size(reverse_search_radius_);
}

bool Vertex::approx_contains(const Eigen::Vector3d& point) const
{
    // by comparing to max and min
    return (point.x() > min_.x() && point.x() < max_.x() &&
            point.y() > min_.y() && point.y() < max_.y() &&
            point.z() > min_.z() && point.z() < max_.z());
}

bool Vertex::approx_contains(const std::shared_ptr<GenericPoint>& generic_point) const
{
    return approx_contains(generic_point->get_position());
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
    return lhs->get_id() == rhs->get_id();
}