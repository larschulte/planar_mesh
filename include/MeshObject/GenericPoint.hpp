#pragma once

#include <memory>
#include <Eigen/Dense>
#include <set>

// Forward declarations
class Storage;

class GenericPoint : public std::enable_shared_from_this<GenericPoint> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d position, Eigen::Vector3d origin);
    void delete_(); 

public:
    int get_id() const;
    Eigen::Vector3d get_position() const;
    Eigen::Vector3d get_origin() const;

private:
    bool deleting_ = false;

    std::weak_ptr<Storage> storage_;
    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
};

bool operator<(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs);
bool operator==(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs);