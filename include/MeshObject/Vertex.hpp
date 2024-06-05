#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

// Forward declarations
class Edge;
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
    void connect_edge(std::weak_ptr<Edge> edge);
    void disconnect_edge(std::weak_ptr<Edge> edge);
    void connect_surface(std::weak_ptr<Surface> surface);
    void disconnect_surface(std::weak_ptr<Surface> surface);

private:
    int id_;
    std::weak_ptr<Storage> storage_;
    std::vector<std::weak_ptr<Edge>> edges_;
    std::vector<std::weak_ptr<Surface>> surfaces_;
    Eigen::Vector3d pos_;
};