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
    // lock vertex
    while (!omp_test_nest_lock(&vertex_lock)) 
    {
        std::cout << "waiting to lock vertex " << id_ << std::endl;
    }

    // add to update rrs tree vertex set
    storage_->add_to_set_of_vertices_to_update_rrs_tree(shared_from_this());

    // log
    if (settings_.log.deletion) std::cout << "Destroying vertex " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;
    
    // // cascade radius reduction
    // cascade_radius_reduction_to_connected_vertices();

    // disconnect
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges = edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    // std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> neighboring_vertices_that_affect_radius = neighboring_vertices_that_affect_radius_;
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& face : faces) disconnect(face);
    // for (const auto& neighboring_vertex : neighboring_vertices_that_affect_radius) disconnect_neighboring_vertex(neighboring_vertex);
    delete_self_from_nearby_vertices();
    delete_self_from_penetrated_vertices();
    delete_self_from_penetrating_vertex_points();
    if (surface_) disconnect(surface_);

    // update delete count
    num_deletes_++;

    // only create penetrated point / generic point if sibling is empty
    if (sibling_vertices_.empty() && can_create_generic_point_)
    {
        storage_->add_to_queue(shared_from_this());
    }

    // disconnect from sibling vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> sibling_vertices = sibling_vertices_;
    for (const auto& sibling_vertex : sibling_vertices) disconnect(sibling_vertex);

    // log
    if (settings_.log.deletion) std::cout << "---------- vertex " << id_ << " destroyed" << std::endl;

    // set expired
    is_expired_ = true;

    // release vertex lock
    omp_unset_nest_lock(&vertex_lock);
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
    // check if have node lock
    if (!omp_test_nest_lock(&node->omp_lock)) throw std::runtime_error("Can't lock node in BVH.");
    // release
    omp_unset_nest_lock(&node->omp_lock);
    return surface_;
}

bool Vertex::has_surface() const
{
    return surface_ != nullptr;
}

const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& Vertex::get_edges() const 
{ 
    return edges_; 
}

const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& Vertex::get_faces() const 
{ 
    return faces_; 
}

const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& Vertex::get_sibling_vertices() const 
{ 
    return sibling_vertices_; 
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

std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> Vertex::get_penetrating_vertices() const
{
    // initialize
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> penetrating_vertices;
    
    // add
    for (const auto& vertex : distances_to_ray_of_penetrating_vertex_points_)
    {
        penetrating_vertices.insert(vertex.first);
    }

    // return
    return penetrating_vertices;
}

std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> Vertex::get_penetrating_interior_points() const
{
    // initialize
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> penetrating_interior_points;
    
    // add
    for (const auto& interior_point : distances_to_ray_of_penetrating_interior_points_)
    {
        penetrating_interior_points.insert(interior_point.first);
    }

    // return
    return penetrating_interior_points;
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

    // connect
    bool inserted = edges_.insert(edge).second;
    if (inserted) edge->connect(shared_from_this());

    // update boundary state
    if (inserted) check_if_update_search_tree();
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
    bool inserted = surface_ != surface;
    if (inserted) surface_ = surface;
    if (inserted) 
    {
        // update projected position
        projected_position_ = surface->compute_point_projective_position(get_origin(), get_original_position());
    }
    if (inserted) surface->connect(shared_from_this());
    if (inserted) is_singular_ = true;
    if (inserted) update_singular_state();
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
    check_if_update_search_tree();

    // check self destruct
    if (!deleting_ && edges_.empty() && can_self_destruct_) storage_->delete_vertex(shared_from_this());
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

    // do not self destruct when have no face
    // check self destruct
    if (!deleting_ && faces_.empty() && can_self_destruct_) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Surface>& surface)
{
    // lock node
    std::shared_ptr<RRSNode> node_copy = node ? node : std::make_shared<RRSNode>(); // lock if node exists
    while (!omp_test_nest_lock(&node_copy->omp_lock)) 
    {
        std::cout << "disconnect vertex waiting " << id_ << std::endl;
    }

    // check input
    if (surface->is_expired()) return;
    
    // disconnect
    bool erased = surface_ == surface;
    if (erased) surface->disconnect(shared_from_this());
    if (erased) is_singular_ = false;
    if (erased) surface_ = nullptr;
    if (erased) projected_position_ = Eigen::Vector3d::Zero();

    // check self destruct
    if (!deleting_ && erased && can_self_destruct_) storage_->delete_vertex(shared_from_this());

    // release lock
    omp_unset_nest_lock(&node_copy->omp_lock);
}

void Vertex::disconnect(const std::shared_ptr<Vertex>& sibling_vertex)
{
    // check input
    if (sibling_vertex->is_expired()) return;

    // disconnect
    bool erased = sibling_vertices_.erase(sibling_vertex);
    if (erased) sibling_vertex->disconnect(shared_from_this());
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

        // un-add the face if face is non-manifold or penetrated
        if (new_face->is_non_manifold() || new_face->is_penetrated())
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
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_copy = edges_;

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

void Vertex::add_nearby_vertex(const std::shared_ptr<Vertex>& rrs_vertex)
{
    // check input 
    if (rrs_vertex->is_expired()) return;

    // compute distance
    const double distance = (get_position() - rrs_vertex->get_position()).norm(); 

    // add to map
    distances_to_nearby_vertices_[rrs_vertex] = distance;
}

void Vertex::delete_nearby_vertex(const std::shared_ptr<Vertex>& rrs_vertex)
{
    distances_to_nearby_vertices_.erase(rrs_vertex);
}

void Vertex::add_self_to_nearby_vertices()
{
    // make copy as list may change
    std::unordered_map<std::shared_ptr<Vertex>, double, MeshObjectHash> distance_to_neighboring_rrs_vertices_copy = distances_to_nearby_vertices_;

    // update
    for (const auto& [neighboring_vertex, distance] : distance_to_neighboring_rrs_vertices_copy)
    {
        // skip if expired
        if (neighboring_vertex->is_expired()) continue;

        // add to neighboring vertex
        neighboring_vertex->add_nearby_vertex(shared_from_this());

        // try update
        neighboring_vertex->try_update_radius();
        neighboring_vertex->try_break_edges();

        // skip if expired
        if (neighboring_vertex->is_expired()) continue;

        // try update
        neighboring_vertex->try_update_node_box();
    }
}

void Vertex::delete_self_from_nearby_vertices()
{
    // make copy as list may change
    std::unordered_map<std::shared_ptr<Vertex>, double, MeshObjectHash> distance_to_neighboring_rrs_vertices_copy = distances_to_nearby_vertices_;

    // update
    for (const auto& [neighboring_vertex, distance] : distance_to_neighboring_rrs_vertices_copy)
    {
        // try lock the vertex
        if (!omp_test_nest_lock(&neighboring_vertex->vertex_lock)) continue;

        // skip if expired
        if (neighboring_vertex->is_expired()) 
        {
            // release lock
            omp_unset_nest_lock(&neighboring_vertex->vertex_lock);
            continue;
        }

        // try lock the surface
        std::shared_ptr<Surface> surface_copy = neighboring_vertex->get_surface();
        if (!omp_test_nest_lock(&surface_copy->lock)) 
        {
            // release vertex lock
            omp_unset_nest_lock(&neighboring_vertex->vertex_lock);
            continue;
        }

        // delete from neighboring vertex
        neighboring_vertex->delete_nearby_vertex(shared_from_this());        

        // try update
        neighboring_vertex->try_update_radius();
        neighboring_vertex->try_break_edges();
        
        // skip if expired
        if (neighboring_vertex->is_expired()) 
        {
            // release lock
            omp_unset_nest_lock(&neighboring_vertex->vertex_lock);
            omp_unset_nest_lock(&surface_copy->lock);
            continue;
        }

        // try close holes
        neighboring_vertex->try_close_holes_repeatedly();

        // try update
        neighboring_vertex->try_update_node_box();

        // release lock
        omp_unset_nest_lock(&neighboring_vertex->vertex_lock);
        omp_unset_nest_lock(&surface_copy->lock);
    }
}

void Vertex::add_penetrating_interior_point(const std::shared_ptr<InteriorPoint>& interior_point)
{
    // check input
    if (interior_point->is_expired()) return;

    // add self to interior point as well
    interior_point->add_penetrated_vertex(shared_from_this());

    // compute projected position of interior point
    const Eigen::Vector3d projected_position = get_surface()->compute_point_projective_position(interior_point->get_origin(), interior_point->get_original_position());

    // compute distance
    const double distance = (projected_position - get_position()).norm();

    // add to list
    distances_to_ray_of_penetrating_interior_points_[interior_point] = distance;
}

void Vertex::delete_penetrating_interior_point(const std::shared_ptr<InteriorPoint>& interior_point)
{
    distances_to_ray_of_penetrating_interior_points_.erase(interior_point);
}

void Vertex::add_penetrating_vertex_point(const std::shared_ptr<Vertex>& vertex_point)
{
    // check input
    if (vertex_point->is_expired()) return;

    // add self to interior point as well
    vertex_point->add_penetrated_vertex(shared_from_this());

    // compute projected position of vertex point
    const Eigen::Vector3d projected_position = get_surface()->compute_point_projective_position(vertex_point->get_origin(), vertex_point->get_original_position());

    // compute distance
    const double distance = (projected_position - get_position()).norm();

    // add to list
    distances_to_ray_of_penetrating_vertex_points_[vertex_point] = distance;
}


void Vertex::delete_penetrating_vertex_point(const std::shared_ptr<Vertex>& vertex_point)
{
    distances_to_ray_of_penetrating_vertex_points_.erase(vertex_point);
}

void Vertex::delete_self_from_penetrating_vertex_points()
{
    // make copy as list may change
    std::unordered_map<std::shared_ptr<Vertex>, double, MeshObjectHash> distances_to_ray_of_penetrating_vertex_points_copy = distances_to_ray_of_penetrating_vertex_points_;

    // update
    for (const auto& [vertex, distance] : distances_to_ray_of_penetrating_vertex_points_copy)
    {
        // try lock the vertex
        if (!omp_test_nest_lock(&vertex->vertex_lock)) continue;

        // skip if expired
        if (vertex->is_expired()) 
        {
            // release lock
            omp_unset_nest_lock(&vertex->vertex_lock);
            continue;
        }

        // try lock the surface
        std::shared_ptr<Surface> surface_copy = vertex->get_surface();
        if (!omp_test_nest_lock(&surface_copy->lock)) 
        {
            // release vertex lock
            omp_unset_nest_lock(&vertex->vertex_lock);
            continue;
        }

        // delete from neighboring vertex
        vertex->delete_penetrated_vertex(shared_from_this());

        // try update
        vertex->try_update_radius();
        vertex->try_break_edges();
        
        // skip if expired
        if (vertex->is_expired()) 
        {
            // release lock
            omp_unset_nest_lock(&vertex->vertex_lock);
            omp_unset_nest_lock(&surface_copy->lock);
            continue;
        }

        // try close holes
        vertex->try_close_holes_repeatedly();

        // try update
        vertex->try_update_node_box();

        // release lock
        omp_unset_nest_lock(&vertex->vertex_lock);
        omp_unset_nest_lock(&surface_copy->lock);
    }
}

void Vertex::add_penetrated_vertex(const std::shared_ptr<Vertex>& vertex)
{
    // check input
    if (vertex->is_expired()) return;

    // if enough points, compute point to plane distance
    double distance = settings_.radius_value;
    if (vertex->get_surface()->get_total_point_size() >= settings_.fit_plane_threshold)
    {
        distance = vertex->get_surface()->compute_point_to_plane_distance(get_position());
    }

    // add to list
    distances_to_plane_of_penetrated_vertex_points_[vertex] = distance;
}

void Vertex::delete_penetrated_vertex(const std::shared_ptr<Vertex>& vertex)
{
    distances_to_plane_of_penetrated_vertex_points_.erase(vertex);
}

void Vertex::add_self_to_penetrated_vertices()
{
    // make copy as list may change
    std::unordered_map<std::shared_ptr<Vertex>, double, MeshObjectHash> distances_to_plane_of_penetrated_vertex_points_copy = distances_to_plane_of_penetrated_vertex_points_;

    // update
    for (const auto& [vertex, distance] : distances_to_plane_of_penetrated_vertex_points_copy)
    {
        // skip if expired
        if (vertex->is_expired()) continue;

        // add to neighboring vertex
        vertex->add_penetrating_vertex_point(shared_from_this());

        // try update
        vertex->try_update_radius();
        vertex->try_break_edges();

        // skip if expired
        if (vertex->is_expired()) continue;

        // try update
        vertex->try_update_node_box();
    }
}

void Vertex::delete_self_from_penetrated_vertices()
{
    // make copy as list may change
    std::unordered_map<std::shared_ptr<Vertex>, double, MeshObjectHash> distances_to_plane_of_penetrated_vertex_points_copy = distances_to_plane_of_penetrated_vertex_points_;

    // update
    for (const auto& [vertex, distance] : distances_to_plane_of_penetrated_vertex_points_copy)
    {
        // try lock the vertex
        if (!omp_test_nest_lock(&vertex->vertex_lock)) continue;

        // skip if expired
        if (vertex->is_expired()) 
        {
            // release lock
            omp_unset_nest_lock(&vertex->vertex_lock);
            continue;
        }

        // try lock the surface
        std::shared_ptr<Surface> surface_copy = vertex->get_surface();
        if (!omp_test_nest_lock(&surface_copy->lock)) 
        {
            // release vertex lock
            omp_unset_nest_lock(&vertex->vertex_lock);
            continue;
        }

        // delete from neighboring vertex
        vertex->delete_penetrating_vertex_point(shared_from_this());

        // try update
        vertex->try_update_radius();
        vertex->try_break_edges();
        
        // skip if expired
        if (vertex->is_expired()) 
        {
            // release lock
            omp_unset_nest_lock(&vertex->vertex_lock);
            omp_unset_nest_lock(&surface_copy->lock);
            continue;
        }

        // try close holes
        vertex->try_close_holes_repeatedly();

        // try update
        vertex->try_update_node_box();

        // release lock
        omp_unset_nest_lock(&vertex->vertex_lock);
        omp_unset_nest_lock(&surface_copy->lock);
    }
}

void Vertex::cascade_radius_reduction_to_connected_vertices()
{
    // get connected vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_vertices = compute_connected_vertices();

    // update
    for (const std::shared_ptr<Vertex>& connected_vertex : connected_vertices)
    {
        // skip if expired
        if (connected_vertex->is_expired()) continue;

        // add nearby vertices to connected vertex
        for (const auto& [nearby_vertex, distance] : distances_to_nearby_vertices_)
        {
            connected_vertex->add_nearby_vertex(nearby_vertex);
        }

        // add penetrating interior points to connected vertex
        for (const auto& [interior_point, distance] : distances_to_ray_of_penetrating_interior_points_)
        {
            connected_vertex->add_penetrating_interior_point(interior_point);
        }

        // try update
        connected_vertex->try_update_radius();
        connected_vertex->try_break_edges();
        
        // skip if expired
        if (connected_vertex->is_expired()) continue;

        // try update
        connected_vertex->try_update_node_box();
    }
}

double Vertex::compute_radius()
{
    // reset to default
    double new_radius = settings_.radius_value;

    // remove expired entries
    for (auto it = distances_to_nearby_vertices_.begin(); it != distances_to_nearby_vertices_.end();)
    {
        if (it->first->is_expired())
        {
            it = distances_to_nearby_vertices_.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // remove expired entries
    for (auto it = distances_to_ray_of_penetrating_interior_points_.begin(); it != distances_to_ray_of_penetrating_interior_points_.end();)
    {
        if (it->first->is_expired())
        {
            it = distances_to_ray_of_penetrating_interior_points_.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // remove expired entries
    for (auto it = distances_to_ray_of_penetrating_vertex_points_.begin(); it != distances_to_ray_of_penetrating_vertex_points_.end();)
    {
        if (it->first->is_expired())
        {
            it = distances_to_ray_of_penetrating_vertex_points_.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // remove expired entries
    for (auto it = distances_to_plane_of_penetrated_vertex_points_.begin(); it != distances_to_plane_of_penetrated_vertex_points_.end();)
    {
        if (it->first->is_expired())
        {
            it = distances_to_plane_of_penetrated_vertex_points_.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // reduce value
    for (const auto& [neighboring_vertex, distance] : distances_to_nearby_vertices_)
    {
        const double extra_radius = 0.2f;
        if (distance + extra_radius < new_radius) new_radius = distance + extra_radius;
    }

    // reduce value
    for (const auto& [interior_point, distance] : distances_to_ray_of_penetrating_interior_points_)
    {
        if (distance < new_radius) new_radius = distance;
    }

    // reduce value
    for (const auto& [vertex_point, distance] : distances_to_ray_of_penetrating_vertex_points_)
    {
        if (distance < new_radius) new_radius = distance;
    }

    // reduce value
    for (const auto& [vertex_point, distance] : distances_to_plane_of_penetrated_vertex_points_)
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
    // if for an edge, any vertices have smaller radius than the length of the edge, delete the edge
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_copy = edges_;
    for (const std::shared_ptr<Edge>& edge : edges_copy)
    {
        if (edge->is_expired()) continue; // could turn expired if below deletes an edge which then deletes a face
        if (edge->is_deleting()) continue; // could be deleting

        const double edge_length = edge->get_length();
        const double radius0 = edge->get_vertex(0)->get_radius();
        const double radius1 = edge->get_vertex(1)->get_radius();
        if (!settings_.edge_is_short_enough(edge_length, radius0, radius1))
        {
            storage_->delete_edge(edge);
        }
    }

    // skip if expired
    if (is_expired()) return;

    // delete face if penetrated
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_copy = faces_;
    for (const std::shared_ptr<Face>& face : faces_copy)
    {
        if (face->is_expired()) continue; // could turn expired if below deletes a face
        if (face->is_deleting()) continue; // could be deleting

        if (face->is_penetrated())
        {
            storage_->delete_face(face);
        }
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
        // skip if can't lock node
        if (!omp_test_nest_lock(&node->omp_lock)) return;

        // update node
        if (current_rrs_half_size > previous_rrs_half_size)
        {
            node->box = RRSBoundingBox(min_, max_);
            node->recursive_expand_parent_box();
        }
        else if (current_rrs_half_size < previous_rrs_half_size)
        {
            node->box = RRSBoundingBox(min_, max_);
            node->recursive_shrink_parent_box();
        }

        // release lock
        omp_unset_nest_lock(&node->omp_lock);
    }
}

void Vertex::review_surfaces()
{
    // review have the following cases
    // - high confidence and not match -> delete
    // - high confidence and match -> need sibling info
    // - low confidence -> need sibling info

    // what about sibling interior points??
    
    // skip if already under review
    if (under_review_) return;
    under_review_ = true;

    if (settings_.log.review_surfaces) std::cout << "reviewing vertex " << id_ << std::endl;

    // delete if surface is high confidence and mismatched
    std::shared_ptr<Surface> surface = get_surface();

    // disconnect surface from this vertex when reviewing
    can_self_destruct_ = false;
    disconnect(surface);
    can_self_destruct_ = true;
    

    bool low_confidence = surface->get_total_point_size() < settings_.fit_plane_threshold;
    if (!low_confidence)
    {
        // // mismatch if observed from behind
        // Eigen::Vector3d normal = surface->get_normal();
        // Eigen::Vector3d direction = get_direction();
        // bool observed_from_behind = normal.dot(direction) > 0;
        // mismatch if not within surface
       bool not_within_surface = surface->check_relative_position(shared_from_this()) != RelativePosition::WITHIN;
        
        // delete if surface is high confidence and mismatched
        // bool mismatch = observed_from_behind || not_within_surface;
        bool mismatch = not_within_surface;
        if (mismatch)
        {
            // find connected vertices
            std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_vertices = compute_connected_vertices();

            // find connected interior points
            std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> connected_interior_points = compute_connected_interior_points();

            storage_->delete_vertex(shared_from_this());
            under_review_ = false;
            return;
        }
    }

    // if reached here, means either low confidence or high confidence but matched, record the surface uncertainty measure
    current_surface_uncertainty_ = (surface->get_total_point_size() < settings_.fit_plane_threshold) ? std::numeric_limits<double>::max() : surface->get_surface_position_std_in_normal_direction();

    // connect surface back
    connect(surface);

    // ask siblings to review themselves
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> sibling_vertices_copy = sibling_vertices_;
    if (settings_.log.review_surfaces) std::cout << ">> Reviewing sibling vertices" << std::endl;
    for (const std::shared_ptr<Vertex>& sibling : sibling_vertices_copy)
    {
        // skip if expired
        if (sibling->is_expired()) continue; // some sibling vertex may be expired during previous review 

        // skip if sibling is already under review
        if (sibling->is_under_review()) continue;

        // review
        sibling->can_create_generic_point(false);
        sibling->review_surfaces();
        sibling->can_create_generic_point(true);
    }
    if (settings_.log.review_surfaces) std::cout << ">> Finished reviewing sibling vertices" << std::endl;

    // skip if current one is expired during sibling review
    if (is_expired()) return;

    // given correct uncertainty envelope computation, can assume no surface will overlap each other
    // a point will only be connected to a surface if there is edge connecting to the surface
    
    // Currently
    // if a point is connected to multiple high confidence matched surface, will delete the current positional uncertainty if not lowest.

    // Proposed
    // if point is in high uncertainty surface, check if can be merged into low uncertainty surface, if not, delete the point
    // if can merged, merge into the low uncertainty surface

    // Procedure
    // 1. collected a list of all surfaces and positional uncertainty
    // 2. sort so surface with smallest positional uncertainty is first
    // 3. check if the current surface can merge with the first uncertainty surface, merge into it, break
    // 4. if can't merge, check next surface that have smaller positional uncertainty
    // 5. if there is a surface with smaller positional uncertainty, but the current larger positoinal uncertianty surface can't merge into any, delete this vertex

    // get list of siblings and their surface uncertainty
    std::vector<std::pair<std::shared_ptr<Vertex>, double>> sibling_surface_uncertainty_list;
    for (const std::shared_ptr<Vertex>& sibling : sibling_vertices_copy)
    {
        // skip if expired
        if (sibling->is_expired()) continue;

        // record
        sibling_surface_uncertainty_list.push_back(std::make_pair(sibling, sibling->get_current_surface_uncertainty()));
    }

    // sort by surface uncertainty
    std::sort(sibling_surface_uncertainty_list.begin(), sibling_surface_uncertainty_list.end(), 
        [](const std::pair<std::shared_ptr<Vertex>, double>& a, const std::pair<std::shared_ptr<Vertex>, double>& b) { return a.second < b.second; });

    // start from the smallest uncertainty, check if the current surface can merge into the sibling surface
    bool exists_smaller_uncertainty = false;
    bool merge_happened = false;
    for (const auto& sibling_surface_uncertainty_pair : sibling_surface_uncertainty_list)
    {
        // extract
        const std::shared_ptr<Vertex> sibling_vertex = sibling_surface_uncertainty_pair.first;
        const std::shared_ptr<Surface> sibling_surface = sibling_surface_uncertainty_pair.first->get_surface();
        double sibling_surface_uncertainty = sibling_surface_uncertainty_pair.second;

        // skip if sibling surface have higher positional uncertainty than current one
        if (sibling_surface_uncertainty > current_surface_uncertainty_) continue;
        exists_smaller_uncertainty = true;
        
        // check if can merge
        can_self_destruct_ = false;
        disconnect(surface); // remove dupliate point before checking if can merge
        bool can_merge = surface->can_merge(sibling_surface);
        connect(surface);
        can_self_destruct_ = true;

        // merging
        if (can_merge) 
        {
            if (settings_.log.review_surfaces) std::cout << ">> Merging between current vertex " << id_ << " with surface " << surface->get_id() << " and sibling vertex " << sibling_vertex->get_id() << " with surface " << sibling_surface->get_id() << std::endl;

            // flag
            merge_happened = true;

            // merge (the smaller one is absorbed by the larger one)
            if (sibling_surface->get_total_point_size() <= surface_->get_total_point_size())
            {
                absorbs(sibling_vertex);
                storage_->delete_vertex(sibling_vertex);
            }
            else
            {
                sibling_vertex->absorbs(shared_from_this());
                storage_->delete_vertex(shared_from_this());
            }

            // break
            break;
        }
    }

    if (exists_smaller_uncertainty && !merge_happened)
    {
        storage_->delete_vertex(shared_from_this());
    }

    // return
    under_review_ = false;
    return;
}

bool Vertex::is_under_review() const
{
    return under_review_;
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

void Vertex::update_singular_state()
{
    // count number of faces in this surface
    int num_faces_in_surface = faces_.size();

    // update singular state
    if (num_faces_in_surface == 0) is_singular_ = true;
    else is_singular_ = false;
}

bool Vertex::is_confirmed() const
{
    return is_confirmed_;
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
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_copy = input_vertex->get_faces();
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_copy = input_vertex->get_edges();
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
    std::cout << "Connected to " << edges_.size() << " edges, " << faces_.size() << " faces, 1 surface, " << sibling_vertices_.size() << " sibling vertices." << std::endl;
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
    // if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired vertices");
    return lhs->get_id() == rhs->get_id();
}