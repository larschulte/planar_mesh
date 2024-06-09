#pragma once

#include <memory>
#include <Eigen/Dense>
#include <set>

// Forward declarations
class Storage;
class Face;
class Surface;

class InteriorPoint : public std::enable_shared_from_this<InteriorPoint> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Face> face, Eigen::Vector3d position, Eigen::Vector3d origin);
    void delete_(); 

public:
    int get_id() const;
    Eigen::Vector3d get_position() const;
    Eigen::Vector3d get_origin() const;

    void connect(std::weak_ptr<Face> face);
    void connect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<Face> face);
    void disconnect(std::weak_ptr<Surface> surface);

private:
    bool deleting_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Face>> faces_;
    std::set<std::weak_ptr<Surface>> surfaces_;

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
};

bool operator<(const std::weak_ptr<InteriorPoint>& lhs, const std::weak_ptr<InteriorPoint>& rhs);
bool operator==(const std::weak_ptr<InteriorPoint>& lhs, const std::weak_ptr<InteriorPoint>& rhs);