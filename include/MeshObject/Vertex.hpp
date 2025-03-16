#pragma once

#include <memory>
#include <Eigen/Dense>
#include <unordered_set>

#include "MeshObject/MeshObject.hpp"
#include "MeshObject/Settings.hpp"
#include <map>

#include "Cache/FIFOCache.hpp"

#include "MeshObject/RRSTree.hpp"

// Forward declarations
class Edge;
class Face;
class Storage;
class Surface;
class GenericPoint;
class InteriorPoint;

class Vertex : public std::enable_shared_from_this<Vertex>, public MeshObject
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius, double distance_travelled);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<GenericPoint>& generic_point);
    void delete_();

public:

    ~Vertex();

    // read write locks for shared resources
    mutable std::shared_mutex rwlock_vertex_point_distance_publishers_;
    mutable std::shared_mutex rwlock_interior_point_distance_publishers_;
    mutable std::shared_mutex rwlock_vertex_point_distance_subscribers_;

    mutable std::shared_mutex rwlock_edges_;
    mutable std::shared_mutex rwlock_faces_;

    mutable std::shared_mutex rwlock_lifecycle_;

    void temp_initialize(const Eigen::Vector3d& position, unsigned int id);
    std::shared_ptr<RRSNode> node;

    const int& get_id() const;
    const Eigen::Vector3d& get_original_position() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const double& get_distance_travelled() const;
    const Eigen::Vector3d& get_direction() const;
    const std::shared_ptr<Surface>& get_surface() const;
    std::vector<std::shared_ptr<Edge>> get_edges() const;
    std::vector<std::shared_ptr<Face>> get_faces() const;
    std::size_t get_num_deletes() const;
    double get_current_surface_uncertainty() const;

    double& get_projected_uncertainty();

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& get_vertex_point_distance_publishers();
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& get_interior_point_distance_publishers();

    std::shared_ptr<Edge> get_edge(const std::shared_ptr<Vertex>& vertex) const;

    // void try_merge_surfaces();

    Eigen::Vector3d compute_projected_position();
    double compute_projected_distance();
    const Eigen::Vector2d& get_surface_coordinate(const std::shared_ptr<Surface>& surface);
    const Eigen::Vector2d& get_surface_coordinate();

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> compute_connected_vertices();
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> compute_connected_interior_points();

    bool is_expired() const;
    bool is_boundary() const;
    bool is_searchable() const;
    bool is_deleting() const;

    void connect(const std::shared_ptr<Edge>& edge);
    void connect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Edge>& edge);
    void disconnect(const std::shared_ptr<Face>& face);

    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> get_connected_boundary_edges() const;
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> get_connected_boundary_vertices();
    bool check_connected_by_edge(const std::shared_ptr<Vertex>& vertex);
    bool check_connected_by_face(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1);
    bool try_close_holes_repeatedly();
    bool try_close_holes_between_self_and(std::shared_ptr<Vertex>& vertex0, std::shared_ptr<Vertex>& vertex1, std::shared_ptr<Edge>& edge0, std::shared_ptr<Edge>& edge1);
    bool try_close_holes();
    void remove_all_edges();

    // can self destruct flag
    void set_can_self_destruct(bool can_self_destruct);
    void set_connecting_to_edges_and_faces(bool connecting_to_edges_and_faces);

    // non manifold
    bool is_non_manifold() const;

    void upon_adding_publisher();
    void upon_deleting_publisher();

    // vertex point distance publisher and subscriber
    void add_vertex_point_distance_publisher(const std::shared_ptr<Vertex> vertex_point_publisher);
    void delete_vertex_point_distance_publisher(const std::shared_ptr<Vertex> vertex_point_publisher);
    void add_vertex_point_distance_subscriber(const std::shared_ptr<Vertex> vertex_point_subscriber);
    void delete_vertex_point_distance_subscriber(const std::shared_ptr<Vertex> vertex_point_subscriber);

    // interior point distance publisher
    void add_interior_point_distance_publisher(const std::shared_ptr<InteriorPoint> interior_point_publisher);
    void delete_interior_point_distance_publisher(const std::shared_ptr<InteriorPoint> interior_point_publisher);

    // update radius
    double compute_radius();
    void try_update_radius();
    void try_break_edges();

    void update_singular_state();
    bool is_singular() const;
    
    void check_if_update_search_tree();

    void print_info();

    void set_do_not_add_back_due_to_not_connected(bool do_not_add_back_due_to_not_connected);
    void set_do_not_add_back_due_to_seed_surface(bool do_not_add_back_due_to_seed_surface);

public: // for reverse radius search
    const Eigen::Vector3d& get_min() const;
    const Eigen::Vector3d& get_max() const;
    const double& get_radius() const;
    const double& get_radius(const std::shared_ptr<Surface>& surface) const;
    bool contains(const Eigen::Vector3d& point) const;
    bool approx_contains(const Eigen::Vector3d& point) const;
    bool approx_contains(const std::shared_ptr<GenericPoint>& generic_point) const;

public: // for interactive viewer
    unsigned int index_in_cloud_;

private: // for reverse radius search
    double reverse_search_radius_;
    Eigen::Vector3d min_;
    Eigen::Vector3d max_;

private:
    static Settings settings_;

    bool deleting_ = false;
    bool is_expired_ = true;
    bool is_singular_;
    bool can_self_destruct_ = true;
    bool do_not_add_back_due_to_not_connected_ = false;
    bool do_not_add_back_due_to_seed_surface_ = false;
    double current_surface_uncertainty_;

    std::size_t num_deletes_;

    int id_;
    std::shared_ptr<Storage> storage_;

    bool connecting_to_edges_and_faces_ = false;

    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<Face>> faces_;
    std::shared_ptr<Surface> surface_;

    Eigen::Matrix3d eigenvectors_used_;
    Eigen::Vector2d surface_coordinate_;

    FIFOCache<std::size_t, Eigen::Vector3d> buffer_projected_position_{3};
    FIFOCache<std::size_t, double> buffer_projected_distance_{3};
    

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    double distance_travelled_;
    Eigen::Vector3d direction_;

    double projected_uncertainty_;

    Eigen::Vector3d projected_position_ = Eigen::Vector3d::Zero();

    // a publisher is one that reduces radius of other vertices
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertex_point_distance_publishers_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_point_distance_publishers_;

    // a subscriber is one that is reduced by the current vertex
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertex_point_distance_subscribers_;

public:
    double weight_;
};

bool operator<(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator<=(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator==(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);