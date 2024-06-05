#pragma once

#include <memory>
#include <Eigen/Dense>
#include <set>

// Forward declarations
class Edge;
class Face;
class Storage;
class Surface;

class Vertex : public std::enable_shared_from_this<Vertex> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, Eigen::Vector3d pos);
    void delete_();

public:
    int get_id() const;
    Eigen::Vector3d get_pos() const;

    void connect(std::weak_ptr<Edge> edge);
    void connect(std::weak_ptr<Face> face);
    void connect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<Edge> edge);
    void disconnect(std::weak_ptr<Face> face);
    void disconnect(std::weak_ptr<Surface> surface);

private:
    bool deleting_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Edge>> edges_;
    std::set<std::weak_ptr<Face>> faces_;
    std::set<std::weak_ptr<Surface>> surfaces_;

    Eigen::Vector3d pos_;
};

bool operator<(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs);
bool operator==(const std::weak_ptr<Vertex>& lhs, const std::weak_ptr<Vertex>& rhs);