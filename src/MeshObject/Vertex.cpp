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

    // add to rrs tree
    storage_->add_searchable_vertex(shared_from_this());

    // update boundary state
    update_boundary_state();

    // log
    if (settings_.log.initialize) std::cout << "Vertex " << id_ << " created.\n";
}

void Vertex::initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point)
{
    initialize_(storage, generic_point->get_position(), generic_point->get_origin(), generic_point->get_radius());
    previous_surface_ = generic_point->get_previous_surface();
    previous_radius_ = generic_point->get_previous_radius();
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
    if (settings_.log.deletion) std::cout << "Destroying vertex " << id_ << std::endl;

    // set deletion flag
    deleting_ = true;

    // disconnect
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges = edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces = faces_;
    for (const auto& edge : edges) disconnect(edge);
    for (const auto& face : faces) disconnect(face);
    if (surface_) disconnect(surface_);

    // remove from rrs tree
    storage_->remove_searchable_vertex(shared_from_this());
    is_searchable_ = false;
    
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
}

void Vertex::temp_initialize(const Eigen::Vector3d& position, unsigned int id)
{
    // set expired
    is_expired_ = false;

    // set id
    id_ = id;

    // set position
    position_ = position;
}

const int& Vertex::get_id() const 
{ 
    return id_; 
}

const Eigen::Vector3d& Vertex::get_position() const 
{ 
    return position_; 
}

const Eigen::Vector3d& Vertex::buffer_compute_projected_position(const std::shared_ptr<Surface> surface)
{
    // do cartersian rounding now, swtich to Locality Sensitive Hashing later

    // compute hash
    std::size_t hash = surface->get_approximate_normal_hash();

    // add to cache if not exist
    if (!buffer_projected_position_.exists(hash)) 
    {
        const Eigen::Vector3d computedResult = surface->compute_point_projective_position(get_origin(), get_position());
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
        const double computedResult = surface->compute_point_projective_distance(get_origin(), get_position());
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
        Eigen::Vector3d projected_position = surface->compute_point_projective_position(get_origin(), get_position());
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
    return is_boundary_;
}

bool Vertex::is_searchable() const
{
    return is_searchable_;
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
    bool inserted = surface_ != surface;
    if (inserted) surface_ = surface;
    if (inserted) 
    {
        // if new surface is the same as the previous surface, set the radius to the updated previous radius
        if (surface == previous_surface_) 
        {
            // set radius to the previous radius
            reduce_reverse_radius_search_radius(previous_radius_);
        }

        previous_surface_ = nullptr;
        previous_radius_ = 0;
    }
    if (inserted) surface->connect(shared_from_this());
    if (inserted) is_boundary_ = false;
    if (inserted) update_boundary_state();
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
    update_boundary_state();

    // check self destruct
    if (!deleting_ && edges_.empty()) storage_->delete_vertex(shared_from_this());
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
    // check input
    if (surface->is_expired()) return;

    // disconnect
    bool erased = surface_ == surface;
    if (erased) surface->disconnect(shared_from_this());
    if (erased) is_boundary_ = false;
    if (erased) is_singular_ = false;
    if (erased) surface_ = nullptr;

    // check self destruct
    if (!deleting_ && erased && can_self_destruct_) storage_->delete_vertex(shared_from_this());
}

void Vertex::disconnect(const std::shared_ptr<Vertex>& sibling_vertex)
{
    // check input
    if (sibling_vertex->is_expired()) return;

    // disconnect
    bool erased = sibling_vertices_.erase(sibling_vertex);
    if (erased) sibling_vertex->disconnect(shared_from_this());
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

            // reduce the search radius of the connected vertices
            for (std::shared_ptr<Vertex> vertex : connected_vertices)
            {
                // distance
                double distance = (vertex->get_position() - get_position()).norm();
                
                // reduce the search radius of the searched vertex
                if (settings_.log.review_surfaces) std::cout << ">>   reducing search radius of vertex " << vertex->get_id() << std::endl;
                vertex->reduce_reverse_radius_search_radius(distance);
            }

            // reduce the search radius of the connected interior points
            for (std::shared_ptr<InteriorPoint> interior_point : connected_interior_points)
            {
                // distance
                double distance = (interior_point->get_position() - get_position()).norm();
                
                // reduce the search radius of the searched interior point
                if (settings_.log.review_surfaces) std::cout << ">>   reducing search radius of interior point " << interior_point->get_id() << std::endl;
                interior_point->reduce_reverse_radius_search_radius(distance);
            }

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

void Vertex::update_boundary_state()
{
    if (deleting_) return;

    // becomes boundary when one of the connected edges is boundary, or when the point is alone
    is_boundary_ = false;
    for (const std::shared_ptr<Edge>& edge : edges_)
    {
        if (edge->is_boundary())
        {
            is_boundary_ = true;
            break;
        }
    }
    if (edges_.empty()) is_boundary_ = true;

    // update searchable state
    update_searchable_state();
}

void Vertex::update_searchable_state()
{
    // check if is boundary
    if (is_boundary() && !is_searchable_)
    {
        is_searchable_ = true;
    }
    else if (!is_boundary() && is_searchable_)
    {
        is_searchable_ = false;
    }
}

void Vertex::print_info()
{
    std::cout << "Vertex " << id_ << " at " << position_.transpose() << std::endl;
    std::cout << "Connected to " << edges_.size() << " edges, " << faces_.size() << " faces, 1 surface, " << sibling_vertices_.size() << " sibling vertices." << std::endl;
    std::cout << "Boundary state: " << is_boundary_ << std::endl;
    std::cout << "Singular state: " << is_singular_ << std::endl;
    std::cout << "Searchable state: " << is_searchable_ << std::endl;
    std::cout << "Expired: " << is_expired_ << std::endl;
}

void Vertex::can_create_generic_point(bool state)
{
    can_create_generic_point_ = state;
}

void Vertex::set_reverse_radius_search_radius(double radius)
{
    // set radius
    double previous_radius = reverse_search_radius_;
    reverse_search_radius_ = radius;

    // update min and max
    min_ = position_ - Eigen::Vector3d(radius, radius, radius);
    max_ = position_ + Eigen::Vector3d(radius, radius, radius);

    // should update search tree if expand radius (if new radius is larger than old one, delete then re-add the searchable vertex if it is searchable)
    if (reverse_search_radius_ > previous_radius && is_searchable_)
    {
        storage_->remove_searchable_vertex(shared_from_this());
        storage_->add_searchable_vertex(shared_from_this());
    }
}

void Vertex::reduce_reverse_radius_search_radius(double radius)
{
    if (radius >= reverse_search_radius_) return;

    // update radius
    set_reverse_radius_search_radius(radius);

    // cascade to connected vertices
    // get connected vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> connected_vertices = compute_connected_vertices();
    for (const std::shared_ptr<Vertex>& vertex : connected_vertices)
    {
        // distance
        double distance = (vertex->get_position() - get_position()).norm();
        
        // reduce the search radius of the searched vertex
        vertex->reduce_reverse_radius_search_radius(distance+radius);
    }

    // cascade to interior point of connected faces
    // get connected interior points
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> connected_interior_points = compute_connected_interior_points();
    for (const std::shared_ptr<InteriorPoint>& interior_point : connected_interior_points)
    {
        // distance
        double distance = (interior_point->get_position() - get_position()).norm();
        
        // reduce the search radius of the searched interior point
        interior_point->reduce_reverse_radius_search_radius(distance+radius);
    }


    // if for an edge, any vertices have smaller radius than the length of the edge, delete the edge
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_copy = edges_;
    for (const std::shared_ptr<Edge>& edge : edges_copy)
    {
        if (edge->is_expired()) continue; // could turn expired if below deletes an edge which then deletes a face
        if (edge->is_deleting()) continue; // could be deleting

        if (edge->get_length() > edge->get_vertex(0)->get_radius() || edge->get_length() > edge->get_vertex(1)->get_radius())
        {
            storage_->delete_edge(edge);
        }
    }
}

void Vertex::reduce_previous_radius(double radius)
{
    if (radius >= previous_radius_) return;

    // update radius
    previous_radius_ = radius;
}

Eigen::Vector3d Vertex::get_min() const
{
    return min_;
}

Eigen::Vector3d Vertex::get_max() const
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
    if (surface == previous_surface_)
    {
        return previous_radius_;
    }
    else
    {
        return reverse_search_radius_;
    }
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
    if (lhs->is_expired() || rhs->is_expired()) throw std::runtime_error("Comparing expired vertices");
    return lhs->get_id() == rhs->get_id();
}