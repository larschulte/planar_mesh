#pragma once

#include <memory>
#include <Eigen/Dense>
#include <unordered_set>

#include "MeshObject/MeshObject.hpp"
#include "Cache/FIFOCache.hpp"
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
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius, double distance_travelled);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point);
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled);
    void delete_(); 

public:

    // read write locks
    mutable std::shared_mutex rwlock_interior_point_distance_subscribers_;

    const int& get_id() const;
    const Eigen::Vector3d& get_original_position() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const double& get_distance_travelled() const;
    const Eigen::Vector3d& get_direction() const;
    const std::shared_ptr<Surface>& get_surface() const;
    const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& get_sibling_interior_points() const;
    const double& get_radius() const;
    bool is_expired() const;

    double& get_projected_uncertainty();

    std::size_t get_num_deletes() const;

    const Eigen::Vector3d& buffer_compute_projected_position(const std::shared_ptr<Surface> surface);
    const Eigen::Vector3d& buffer_compute_projected_position();
    const double& buffer_compute_projected_distance(const std::shared_ptr<Surface> surface);
    const double& buffer_compute_projected_distance();

    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void connect(const std::shared_ptr<InteriorPoint>& sibling_interior_point);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<InteriorPoint>& sibling_interior_point);

    void delete_subscribers();
    
    void add_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber);
    void delete_interior_point_distance_subscriber(const std::shared_ptr<Vertex> interior_point_subscriber);

    void set_reverse_radius_search_radius(double radius);

    void update_confirmed_status();
    bool is_confirmed() const;
    
    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);

private:
    static Settings settings_;

    bool deleting_ = false;
    bool is_expired_ = true;

    std::size_t num_deletes_;

    std::size_t num_confirmed_faces = 0;
    bool is_confirmed_ = false;

    bool can_self_destruct_ = true;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::shared_ptr<Surface> surface_;

    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> sibling_interior_points_;

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    double distance_travelled_;
    Eigen::Vector3d direction_;
    double radius_;

    double projected_uncertainty_;

    FIFOCache<std::size_t, Eigen::Vector3d> buffer_projected_position_{3};
    FIFOCache<std::size_t, double> buffer_projected_distance_{3};

    Eigen::Vector3d projected_position_ = Eigen::Vector3d::Zero();

    std::vector<std::shared_ptr<Vertex>> interior_point_distance_subscribers_;

public:
    double weight_;
};

bool operator<(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);
bool operator==(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);