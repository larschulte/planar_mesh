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
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Face>& face, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Face>& face, const std::shared_ptr<GenericPoint>& generic_point);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Face>& face, const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    void delete_(); 

public:
    const int& get_id() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const std::shared_ptr<Surface>& get_surface() const;
    const double& get_radius() const;
    bool is_expired() const;

    std::size_t get_num_deletes() const;

    void try_update_surface_projection();
    const Eigen::Vector3d& get_projected_position();
    const double& get_projected_distance();

    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);

    void swap(const std::shared_ptr<Surface>& surface1, const std::shared_ptr<Surface>& surface2);

private:
    bool deleting_ = false;
    bool is_expired_ = true;

    std::size_t num_deletes_;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    
    std::shared_ptr<Surface> surface_ = nullptr;

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    double radius_;

    Eigen::Vector3d normal_used_;
    Eigen::Vector3d projected_position_;
    double projected_distance_;
};

bool operator<(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);
bool operator==(const std::shared_ptr<InteriorPoint>& lhs, const std::shared_ptr<InteriorPoint>& rhs);