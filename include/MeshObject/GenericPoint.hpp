#pragma once

#include <memory>
#include <Eigen/Dense>
#include <unordered_set>

#include "MeshObject/MeshObject.hpp"
#include "MeshObject/Settings.hpp"
// Forward declarations
class Storage;

class GenericPoint : public std::enable_shared_from_this<GenericPoint>, public MeshObject
{
public:
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin);

protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Vertex>& vertex);
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<InteriorPoint>& interior_point);
    void delete_(); 

public:
    const int& get_id() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const double& get_radius() const;
    bool is_expired() const;

    std::size_t get_num_deletes() const;
    
private:
    static Settings settings_;

    bool deleting_ = false;
    bool is_expired_ = true;

    std::size_t num_deletes_;

    int id_;
    std::shared_ptr<Storage> storage_;
    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
    double radius_;
};

bool operator<(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator==(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);