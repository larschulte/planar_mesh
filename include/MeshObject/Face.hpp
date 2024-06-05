#pragma once

#include <memory>
#include <vector>

// Forward declarations
class Edge;
class Storage;
class Surface;

class Face : public std::enable_shared_from_this<Face> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3);
    void delete_();

public:
    int get_id() const;

    void connect(std::weak_ptr<Vertex> vertex);
    void connect(std::weak_ptr<Edge> edge);
    void connect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<Vertex> vertex);
    void disconnect(std::weak_ptr<Edge> edge);
    void disconnect(std::weak_ptr<Surface> surface);
    

private:
    bool deleting_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Vertex>> vertices_;
    std::set<std::weak_ptr<Edge>> edges_;
    std::set<std::weak_ptr<Surface>> surfaces_;
};

bool operator<(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs);
bool operator==(const std::weak_ptr<Face>& lhs, const std::weak_ptr<Face>& rhs);