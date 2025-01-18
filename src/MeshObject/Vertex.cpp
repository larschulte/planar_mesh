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

    // connect to surface
    {
        // connect to surface
        surface_ = surface;
        surface_->connect(shared_from_this());

        // projected position
        projected_position_ = surface_->compute_point_projective_position(origin_, position_);
    }

    // check if update search tree
    check_if_update_search_tree();

    // set reverse search radius based on input parameter
    try_update_radius();

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

    // set deletion flag
    deleting_ = true;
    
    // log
    if (settings_.log.deletion) std::cout << "Destroying vertex " << id_ << std::endl;
    
    // publishers (disconnect)
    {
        // lock, copy and clear
        std::vector<std::pair<std::shared_ptr<Vertex>, double>> vertex_point_distance_publishers_copy;
        {
            std::unique_lock<std::shared_mutex> lock(rwlock_vertex_point_distance_publishers_);
            vertex_point_distance_publishers_copy = vertex_point_distance_publishers_;
            vertex_point_distance_publishers_.clear();
        }

        // disconnect
        for (const auto& [vertex_point_publisher, distance] : vertex_point_distance_publishers_copy) vertex_point_publisher->delete_vertex_point_distance_subscriber(shared_from_this());
    }

    // publishers (disconnect)
    {
        // lock, copy and clear
        std::vector<std::pair<std::shared_ptr<InteriorPoint>, double>> interior_point_distance_publishers_copy;
        {
            std::unique_lock<std::shared_mutex> lock(rwlock_interior_point_distance_publishers_);
            interior_point_distance_publishers_copy = interior_point_distance_publishers_;
            interior_point_distance_publishers_.clear();
        }

        // disconnect
        for (const auto& [interior_point_publisher, distance] : interior_point_distance_publishers_copy) interior_point_publisher->delete_interior_point_distance_subscriber(shared_from_this());
    }

    // subscribers (disconnect)
    {
        // lock, copy and clear
        std::vector<std::shared_ptr<Vertex>> vertex_point_distance_subscribers_copy;
        {
            std::unique_lock<std::shared_mutex> lock(rwlock_vertex_point_distance_subscribers_);
            vertex_point_distance_subscribers_copy = vertex_point_distance_subscribers_;
            vertex_point_distance_subscribers_.clear();
        }

        // disconnect
        for (const auto& vertex_point_subscriber : vertex_point_distance_subscribers_copy) vertex_point_subscriber->delete_vertex_point_distance_publisher(shared_from_this());

        // add vertex to list to update
        for (const auto& vertex_point_subscriber : vertex_point_distance_subscribers_copy) storage_->add_vertex_that_have_deleted_publishers(vertex_point_subscriber);
    }

    // surface (disconnect)
    {
        // ask to be removed from surface
        surface_->disconnect(shared_from_this());

        // clear surface
        surface_ = nullptr;
    }

    // edges (delete)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock_edges(rwlock_edges_);

        // delete edge
        for (const auto& edge : edges_) storage_->add_edge_to_be_deleted(edge);

        // clear
        edges_.clear();
    }
    
    // faces (delete)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock_faces(rwlock_faces_);

        // delete face
        for (const auto& face : faces_) storage_->add_face_to_be_deleted(face);

        // clear
        faces_.clear();
    }

    // create generic point
    {
        // update delete count
        num_deletes_++;

        // only create penetrated point / generic point if sibling is empty
        if (can_create_generic_point_)
        {
            storage_->add_to_queue(shared_from_this());
        }
    }

    // update tree
    {
        // add to update rrs tree vertex set
        storage_->add_to_set_of_vertices_to_update_rrs_tree(shared_from_this());
    }

    // log
    if (settings_.log.deletion) std::cout << "---------- vertex " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;
}

Vertex::~Vertex()
{
    // std::cout << "Vertex " << id_ << " deleted." << std::endl;
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

Eigen::Vector3d Vertex::compute_projected_position() 
{ 
    return get_surface()->compute_point_projective_position(get_origin(), get_original_position());
}

double Vertex::compute_projected_distance() 
{ 
    return get_surface()->compute_point_projective_distance(get_origin(), get_original_position());
}

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

std::vector<std::shared_ptr<Edge>> Vertex::get_edges() const 
{ 
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

    return edges_; 
}

std::vector<std::shared_ptr<Face>> Vertex::get_faces() const 
{ 
    // read lock
    std::shared_lock<std::shared_mutex> lock(rwlock_faces_);

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

std::shared_ptr<Edge> Vertex::get_edge(const std::shared_ptr<Vertex>& vertex) const
{
    // copy edge list
    std::vector<std::shared_ptr<Edge>> edges;
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

        // copy
        edges = edges_;
    }
    
    for (const std::shared_ptr<Edge>& edge : edges)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

        // skip if expired
        if (edge->is_expired()) continue; 

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
    // make copy of edges
    std::vector<std::shared_ptr<Edge>> edges;
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

        // copy
        edges = edges_;
    }

    // iterate through all edges and get connected vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_vertices;
    for (const std::shared_ptr<Edge>& edge : edges)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

        // skip if edge is deleting
        if (edge->is_expired()) continue;
        
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

    // update boundary state
    check_if_update_search_tree();

    // check self destruct
    if (!deleting_ && edges_.empty() && can_self_destruct_ && !connecting_to_edges_and_faces_) storage_->add_vertex_to_be_deleted(shared_from_this());
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
    // make copy of edges 
    std::vector<std::shared_ptr<Edge>> edges;
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

        // copy
        edges = edges_;
    }

    // check if connected by edge
    for (const auto& edge : edges)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

        // skip if expired
        if (edge->is_expired()) continue;

        // check if connected
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

bool Vertex::try_close_holes_between_self_and(std::shared_ptr<Vertex>& vertex0, std::shared_ptr<Vertex>& vertex1, std::shared_ptr<Edge>& edge0, std::shared_ptr<Edge>& edge1)
{
    // read lock
    std::shared_lock<std::shared_mutex> vertex0_lock(vertex0->rwlock_lifecycle_, std::defer_lock);
    std::shared_lock<std::shared_mutex> vertex1_lock(vertex1->rwlock_lifecycle_, std::defer_lock);
    std::lock(vertex0_lock, vertex1_lock);

    // skip if expired
    if (vertex0->is_expired() || vertex1->is_expired()) return false;

    // try lock the edge connecting the two vertices
    std::shared_lock<std::shared_mutex> edge0_lock(edge0->rwlock_lifecycle_, std::defer_lock);
    std::shared_lock<std::shared_mutex> edge1_lock(edge1->rwlock_lifecycle_, std::defer_lock);
    std::lock(edge0_lock, edge1_lock);

    // skip if expired
    if (edge0->is_expired() || edge1->is_expired()) return false;
    

    // skip if edge is not short enough
    const double edge_length = (vertex0->get_position() - vertex1->get_position()).norm();
    const double radius0 = vertex0->get_radius();
    const double radius1 = vertex1->get_radius();
    if (!settings_.edge_is_short_enough(edge_length, radius0, radius1)) return false;

    // skip if edge intersects
    if (get_surface()->tree_intersect_edge(vertex0, vertex1)) return false;

    // get inter edge and its lock
    std::shared_ptr<Edge> inter_edge = nullptr;
    std::shared_lock<std::shared_mutex> inter_edge_lock;

    // make copy of edges 
    std::vector<std::shared_ptr<Edge>> vertex_0_edges = vertex0->get_edges();

    // get inter edge if exists
    for (const auto& vertex_0_edge : vertex_0_edges)
    {
        // read lock
        inter_edge_lock = std::shared_lock<std::shared_mutex>(vertex_0_edge->rwlock_lifecycle_);

        // skip if expired
        if (vertex_0_edge->is_expired()) continue;

        std::shared_ptr<Vertex> vertex_0_edge_vertex0 = vertex_0_edge->get_vertex(0);
        std::shared_ptr<Vertex> vertex_0_edge_vertex1 = vertex_0_edge->get_vertex(1);

        // skip if nullptr
        if (!vertex_0_edge_vertex0 || !vertex_0_edge_vertex1) continue;

        // skip if none is vertex1
        if (vertex_0_edge_vertex0 != vertex1 && vertex_0_edge_vertex1 != vertex1) continue;

        // get inter edge
        inter_edge = vertex_0_edge;
    }

    if (!inter_edge)
    {
        // create edge
        inter_edge = storage_->add_edge(surface_, vertex0, vertex1);

        // read lock on the inter edge
        inter_edge_lock = std::shared_lock<std::shared_mutex>(inter_edge->rwlock_lifecycle_);

        // skip if inter_edge is expired
        if (inter_edge->is_expired()) return false;

        // connect
        get_surface()->connect(inter_edge);
    }

    // skip if inter_edge is not boundary
    if (!inter_edge->is_boundary()) return false;

    // make copy of faces
    std::vector<std::shared_ptr<Face>> faces;
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(rwlock_faces_);

        // copy
        faces = faces_;
    }

    // skip if face already exists
    for (const auto& face : faces)
    {
        // read lock face
        std::shared_lock<std::shared_mutex> lock(face->rwlock_lifecycle_);

        // skip if exists
        if (face->has_vertex(vertex0) && face->has_vertex(vertex1)) return false;
    }

    // skip if inter_edge is expired
    inter_edge_lock = std::shared_lock<std::shared_mutex>(inter_edge->rwlock_lifecycle_);
    if (inter_edge->is_expired()) return false;
    
    // create face
    std::shared_ptr<Face> new_face = storage_->add_face(get_surface(), shared_from_this(), vertex0, vertex1, edge0, edge1, inter_edge);

    // read lock face
    std::shared_lock<std::shared_mutex> lock(new_face->rwlock_lifecycle_);

    // skip if expired
    if (new_face->is_expired()) return false;

    // connect
    get_surface()->connect(new_face);

    // // un-add the face if face is non-manifold
    // if (new_face->is_non_manifold())
    // {
    //     // un-add face
    //     new_face->un_add_face();

    //     // update flag
    //     changed = false;
    // }
    // else
    // {
    //     // update flag
    //     changed = true;
    // }

    // return changed
    return true;
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

        // copy edges
        std::vector<std::shared_ptr<Edge>> edges;
        {
            // read lock
            std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

            // copy
            edges = edges_;
        }

        // get boundary edges
        std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> connected_boundary_edges;
        for (const auto& edge : edges)
        {
            // read lock
            std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

            // skip if expired
            if (edge->is_expired()) continue;

            // skip if not boundary
            if (!edge->is_boundary()) continue;

            // add to list
            connected_boundary_edges.insert(edge);
        }

        // get boundary vertices
        std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_boundary_vertices;
        // boundary vertices to edge map
        std::unordered_map<std::shared_ptr<Vertex>, std::shared_ptr<Edge>, MeshObjectHash> boundary_vertices_to_edge_map;
        for (const auto& edge : connected_boundary_edges)
        {
            // read lock
            std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

            // skip if expired
            if (edge->is_expired()) continue;

            // vertex 
            std::shared_ptr<Vertex> vertex0 = edge->get_vertex(0);
            std::shared_ptr<Vertex> vertex1 = edge->get_vertex(1);

            // skip if nullptr
            if (!vertex0 || !vertex1) continue;

            if (vertex0 != shared_from_this())
            {
                connected_boundary_vertices.insert(vertex0);
                boundary_vertices_to_edge_map[vertex0] = edge;
            }
            
            if (vertex1 != shared_from_this())
            {
                connected_boundary_vertices.insert(vertex1);
                boundary_vertices_to_edge_map[vertex1] = edge;
            }
        }

        // skip if less than two connected boundary vertex (this includes when the vertex itself is not boundary)
        if (connected_boundary_vertices.size() < 2)
        {
            break;
        }

        // compute connected boundary vertices with angle
        std::vector<std::pair<std::shared_ptr<Vertex>, double>> connected_boundary_vertices_with_angle;
        for (const auto& vertex : connected_boundary_vertices)
        {
            // read lock
            std::shared_lock<std::shared_mutex> lock(vertex->rwlock_lifecycle_);

            // skip if expired
            if (vertex->is_expired()) continue;

            // use projected position
            Eigen::Vector2d direction = vertex->get_surface_coordinate() - get_surface_coordinate();
            double angle = std::atan2(direction.y(), direction.x());

            // add to list
            connected_boundary_vertices_with_angle.push_back(std::make_pair(vertex, angle));
        }

        // sort by angle
        std::sort(connected_boundary_vertices_with_angle.begin(), connected_boundary_vertices_with_angle.end(), 
            [](const std::pair<std::shared_ptr<Vertex>, double>& a, const std::pair<std::shared_ptr<Vertex>, double>& b) 
            { 
                return a.second < b.second; 
            });

        // create pairs of connected boundary vertices
        std::vector<std::pair<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>>> connected_boundary_vertex_pairs;
        for (int i = 0; i < connected_boundary_vertices_with_angle.size(); i++)
        {
            connected_boundary_vertex_pairs.push_back(std::make_pair(connected_boundary_vertices_with_angle[i].first, connected_boundary_vertices_with_angle[(i + 1) % connected_boundary_vertices_with_angle.size()].first));
        }

        // try close holes between the two vertices
        for (auto& [vertex0, vertex1] : connected_boundary_vertex_pairs)
        {
            if (try_close_holes_between_self_and(vertex0, vertex1, boundary_vertices_to_edge_map[vertex0], boundary_vertices_to_edge_map[vertex1]))
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

void Vertex::set_connecting_to_edges_and_faces(bool connecting_to_edges_and_faces)
{
    connecting_to_edges_and_faces_ = connecting_to_edges_and_faces;
}

bool Vertex::is_non_manifold() const
{
    // non manifold if 
    // 1. connected by one boundary edge
    // 2. connected by more than 2 boundary edges
    return get_connected_boundary_edges().size() > 2 || get_connected_boundary_edges().size() == 1;
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
        storage_->add_vertex_that_have_changed_box(shared_from_this());
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
        storage_->add_vertex_that_have_changed_box(shared_from_this());
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
}

double Vertex::compute_radius()
{
    // reset to default
    double new_radius = settings_.radius_value;

    // read lock the lists
    std::shared_lock<std::shared_mutex> lock_vertex_point_distance_publishers(rwlock_vertex_point_distance_publishers_);
    std::shared_lock<std::shared_mutex> lock_interior_point_distance_publishers(rwlock_interior_point_distance_publishers_);

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
    // update radius
    reverse_search_radius_ = compute_radius();

    // update max and min
    const double current_rrs_half_size = settings_.compute_rrs_half_size(reverse_search_radius_);
    const Eigen::Vector3d half_size_vecotr = Eigen::Vector3d(current_rrs_half_size, current_rrs_half_size, current_rrs_half_size);
    min_ = get_position() - half_size_vecotr;
    max_ = get_position() + half_size_vecotr;
}

void Vertex::try_break_edges()
{
    // copy of edges
    std::vector<std::shared_ptr<Edge>> edges_copy;
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(rwlock_edges_);

        edges_copy = edges_;
    }

    // collect edges to delete
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_to_delete;
    for (const std::shared_ptr<Edge>& edge : edges_copy)
    {
        // read lock
        std::shared_lock<std::shared_mutex> lock(edge->rwlock_lifecycle_);

        if (edge->is_expired()) continue; // could turn expired if below deletes an edge which then deletes a face

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

        storage_->add_edge_to_be_deleted(edge);
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