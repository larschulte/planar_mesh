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
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point);
    void delete_();

public:
    std::shared_ptr<RRSNode> node;

    const int& get_id() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const Eigen::Vector3d& get_direction() const;
    const std::shared_ptr<Surface>& get_surface() const;
    bool has_surface() const;
    const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& get_edges() const;
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_faces() const;
    const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& get_sibling_vertices() const;
    std::size_t get_num_deletes() const;
    double get_current_surface_uncertainty() const;

    // void try_merge_surfaces();

    const Eigen::Vector3d& buffer_compute_projected_position(const std::shared_ptr<Surface> surface);
    const Eigen::Vector3d& buffer_compute_projected_position();
    const double& buffer_compute_projected_distance(const std::shared_ptr<Surface> surface);
    const double& buffer_compute_projected_distance();
    const Eigen::Vector2d& get_surface_coordinate(const std::shared_ptr<Surface>& surface);
    const Eigen::Vector2d& get_surface_coordinate();

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> compute_connected_vertices();
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> compute_connected_interior_points();

    bool is_expired() const;
    bool is_boundary() const;
    bool is_searchable() const;

    void connect(const std::shared_ptr<Edge>& edge);
    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void connect(const std::shared_ptr<Vertex>& sibling_vertex);
    void disconnect(const std::shared_ptr<Edge>& edge);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<Vertex>& sibling_vertex);

    void review_surfaces();
    bool is_under_review() const;

    void update_confirmed_status();
    void update_singular_state();
    bool is_confirmed() const;
    bool is_singular() const;
    
    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);
    void absorbs(const std::shared_ptr<Vertex>& input_vertex);

    void update_boundary_state();
    void update_searchable_state();

    void print_info();

    void can_create_generic_point(bool can_create);

public: // for reverse radius search
    void set_reverse_radius_search_radius(double radius);
    void reduce_reverse_radius_search_radius(double radius);
    void reduce_previous_radius(double radius);
    Eigen::Vector3d get_min() const;
    Eigen::Vector3d get_max() const;
    const double& get_radius() const;
    const double& get_radius(const std::shared_ptr<Surface>& surface) const;
    bool contains(const Eigen::Vector3d& point) const;
    bool approx_contains(const Eigen::Vector3d& point) const;
    bool approx_contains(const std::shared_ptr<GenericPoint>& generic_point) const;

private: // for reverse radius search
    double reverse_search_radius_;
    Eigen::Vector3d min_;
    Eigen::Vector3d max_;

private:
    static Settings settings_;

    bool deleting_ = false;
    bool under_review_ = false;
    bool is_boundary_;
    bool is_searchable_ = false;
    bool is_expired_ = true;
    bool is_singular_;
    bool can_self_destruct_ = true;
    bool can_create_generic_point_ = true;
    double current_surface_uncertainty_;

    std::size_t num_deletes_;

    std::size_t num_confirmed_faces = 0;
    bool is_confirmed_ = false;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::shared_ptr<Surface> surface_;
    std::shared_ptr<Surface> previous_surface_;
    double previous_radius_;

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> sibling_vertices_;

    Eigen::Matrix3d eigenvectors_used_;
    Eigen::Vector2d surface_coordinate_;

    FIFOCache<std::size_t, Eigen::Vector3d> buffer_projected_position_{3};
    FIFOCache<std::size_t, double> buffer_projected_distance_{3};
    

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    Eigen::Vector3d direction_;
};

bool operator<(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator<=(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator==(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);