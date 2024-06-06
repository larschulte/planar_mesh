#pragma once

#include <memory>
#include <set>

// Forward declarations
class Vertex;
class Face;
class Storage;
class Surface;

class Edge : public std::enable_shared_from_this<Edge> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Vertex> vertex1, std::weak_ptr<Vertex> vertex2);
    void delete_();

public:
    int get_id() const;

    void connect(std::weak_ptr<Vertex> vertex);
    void connect(std::weak_ptr<Face> face);
    void connect(std::weak_ptr<Surface> surface);
    void disconnect(std::weak_ptr<Vertex> vertex);
    void disconnect(std::weak_ptr<Face> face);
    void disconnect(std::weak_ptr<Surface> surface);

    bool has_vertex(std::weak_ptr<Vertex> vertex) const;
    bool is_boundary() const;
    void update_boundary_state();

private:
    bool deleting_ = false;
    bool is_boundary_ = false;

    int id_;
    std::weak_ptr<Storage> storage_;

    std::set<std::weak_ptr<Vertex>> vertices_;
    std::set<std::weak_ptr<Face>> faces_;
    std::set<std::weak_ptr<Surface>> surfaces_;
};

bool operator<(const std::weak_ptr<Edge>& lhs, const std::weak_ptr<Edge>& rhs);
bool operator==(const std::weak_ptr<Edge>& lhs, const std::weak_ptr<Edge>& rhs);