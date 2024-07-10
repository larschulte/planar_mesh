#pragma once

#include <memory>
#include <Eigen/Dense>
#include <unordered_set>

#include "MeshObject/MeshObject.hpp"

// Forward declarations
class Storage;
class Face;
class Surface;
class GenericPoint;

class InteriorPoint : public std::enable_shared_from_this<InteriorPoint>, public MeshObject 
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<GenericPoint>& generic_point);
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    void delete_(); 

public:
    const int& get_id() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const Eigen::Vector3d& get_direction() const;
    const std::shared_ptr<Surface>& get_surface() const;
    const std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& get_surfaces() const;
    const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& get_sibling_interior_points() const;
    const double& get_radius() const;
    bool is_expired() const;

    std::size_t get_num_deletes() const;

    void try_update_surface_projection(const std::shared_ptr<Surface> surface);
    void try_update_surface_projection();
    const Eigen::Vector3d& get_projected_position(const std::shared_ptr<Surface> surface);
    const Eigen::Vector3d& get_projected_position();
    const double& get_projected_distance(const std::shared_ptr<Surface> surface);
    const double& get_projected_distance();

    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void connect(const std::shared_ptr<InteriorPoint>& sibling_interior_point);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<InteriorPoint>& sibling_interior_point);

    void update_confirmed_status();
    bool is_confirmed() const;
    
    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);

private:
    bool deleting_ = false;
    bool is_expired_ = true;

    std::size_t num_deletes_;

    std::size_t num_confirmed_faces = 0;
    bool is_confirmed_ = false;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_;

    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> sibling_interior_points_;

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    Eigen::Vector3d direction_;
    double radius_;

    Eigen::Vector3d normal_used_;
    Eigen::Vector3d projected_position_;
    double projected_distance_;
};

bool operator<(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);
bool operator==(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);