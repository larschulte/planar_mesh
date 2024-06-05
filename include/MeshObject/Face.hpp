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
    void cascade_delete_from_edge(std::weak_ptr<Edge> edge);
    void connect_surface(std::weak_ptr<Surface> surface);
    void disconnect_surface(std::weak_ptr<Surface> surface);

private:
    int id_;
    std::weak_ptr<Edge> edge1_;
    std::weak_ptr<Edge> edge2_;
    std::weak_ptr<Edge> edge3_;
    std::vector<std::weak_ptr<Surface>> surfaces_;
    std::weak_ptr<Storage> storage_;
};