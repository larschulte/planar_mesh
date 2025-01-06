#pragma once

#include <memory>
#include <Eigen/Dense>
#include <unordered_set>

#include "MeshObject/MeshObject.hpp"
#include "MeshObject/Settings.hpp"
#include "MeshObject/Surface.hpp"
// Forward declarations
class Storage;

class GenericPoint : public std::enable_shared_from_this<GenericPoint>, public MeshObject
{
public:
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled);

protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Vertex>& vertex);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<InteriorPoint>& interior_point);
    void delete_(); 

public:
    const int& get_id() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const double& get_distance_travelled() const;
    const Eigen::Vector3d& get_direction() const;
    const Eigen::Vector3d& get_inv_direction() const;
    const double& get_radius() const;
    bool is_expired() const;

    double& get_projected_uncertainty();

    std::size_t get_num_deletes() const;
    void reset_num_deletes();

    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> intersected_surfaces;
    std::unordered_map<std::shared_ptr<Surface>, unsigned int, MeshObjectHash> contented_surfaces;
    unsigned int get_contention_count() const;
    
private:
    static Settings settings_;

    bool deleting_ = false;
    bool is_expired_ = true;

    std::size_t num_deletes_;

    int id_;
    std::shared_ptr<Storage> storage_;
    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    double distance_travelled_;
    Eigen::Vector3d direction_;
    Eigen::Vector3d inv_direction_;
    double radius_;

    double projected_uncertainty_;
};

bool operator<(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator==(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);