#pragma once

#include <memory>
#include <Eigen/Dense>
#include <unordered_set>

#include "MeshObject/MeshObject.hpp"
// #include "Cache/FIFOCache.hpp"
#include "MeshObject/Settings.hpp"

// Forward declarations
class Storage;
class Face;
class Surface;
class GenericPoint;

class InteriorPoint : public std::enable_shared_from_this<InteriorPoint>, public MeshObject 
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<Face>& face, const std::shared_ptr<GenericPoint>& generic_point);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex, const std::shared_ptr<GenericPoint>& generic_point);
    void delete_(); 

public:

    ~InteriorPoint();

    // read write locks
    mutable std::shared_mutex rwlock_interior_point_distance_subscribers_;

    mutable std::shared_mutex rwlock_lifecycle_;

    const int& get_id() const;
    const Eigen::Vector3d& get_original_position() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const double& get_distance_travelled() const;
    const Eigen::Vector3d& get_direction() const;
    const std::shared_ptr<Surface>& get_surface() const;
    const double& get_radius() const;
    bool is_expired() const;

    double& get_projected_uncertainty();

    // get surface coordinate
    const Eigen::Vector2d& get_surface_coordinate(const std::shared_ptr<Surface>& surface);
    const Eigen::Vector2d& get_surface_coordinate();

    std::size_t get_num_deletes() const;

    Eigen::Vector3d compute_projected_position();
    double compute_projected_distance();
    
    void add_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber);
    void delete_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber);

    void set_reverse_radius_search_radius(double radius);

    void set_do_not_add_back_due_to_seed_surface(bool do_not_add_back_due_to_seed_surface);

private:
    static Settings settings_;

    bool deleting_ = false;
    bool is_expired_ = true;

    std::size_t num_deletes_;

    bool can_self_destruct_ = true;
    bool do_not_add_back_due_to_seed_surface_ = false;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::shared_ptr<Face> face_;
    std::shared_ptr<Vertex> vertex_;
    std::shared_ptr<Surface> surface_;

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    double distance_travelled_;
    Eigen::Vector3d direction_;
    double radius_;

    double projected_uncertainty_;

    Eigen::Matrix3d eigenvectors_used_;
    Eigen::Vector2d surface_coordinate_;

    // FIFOCache<std::size_t, Eigen::Vector3d> buffer_projected_position_{3};
    // FIFOCache<std::size_t, double> buffer_projected_distance_{3};

    Eigen::Vector3d projected_position_ = Eigen::Vector3d::Zero();

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> interior_point_distance_subscribers_;

public:
    double weight_;
};

bool operator<(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);
bool operator==(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);