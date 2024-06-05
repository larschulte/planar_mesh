#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

// forward declarations
class Vertex;
class Edge;
class Face;
class Storage;

class Surface : public std::enable_shared_from_this<Surface> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage);
    void delete_();

public:
    int get_id() const;

    void connect(std::weak_ptr<Vertex> vertex);
    void connect(std::weak_ptr<Edge> edge);
    void connect(std::weak_ptr<Face> face);

    void disconnect(std::weak_ptr<Vertex> vertex);
    void disconnect(std::weak_ptr<Edge> edge);
    void disconnect(std::weak_ptr<Face> face);
    
private:
    bool deleting_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Vertex>> vertices_;
    std::set<std::weak_ptr<Edge>> edges_;
    std::set<std::weak_ptr<Face>> faces_;

    Eigen::Vector3d mean_;
    Eigen::Matrix3d covariance_;
    Eigen::Matrix3d eigenvectors_;
    Eigen::Vector3d eigenvalues_;
    Eigen::Vector3d normal_;
};

bool operator<(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs);
bool operator==(const std::weak_ptr<Surface>& lhs, const std::weak_ptr<Surface>& rhs);