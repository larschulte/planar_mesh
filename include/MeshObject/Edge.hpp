#pragma once

#include <memory>
#include <set>
#include <map>

// Forward declarations
class Vertex;
class Face;
class Storage;
class Surface;

class Edge : public std::enable_shared_from_this<Edge> 
{
protected:
    friend class Storage;
    void initialize_(const std::shared_ptr<Storage>& storage, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    void delete_();

public:
    const int& get_id() const;
    const std::shared_ptr<Vertex>& get_vertex(int index) const;
    const std::shared_ptr<Surface>& get_surface() const;
    const Eigen::Vector3d& get_center() const;
    const Eigen::Vector3d& get_max() const;
    const Eigen::Vector3d& get_min() const;
    bool is_expired() const;

    void connect(const std::shared_ptr<Vertex>& vertex);
    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<Vertex>& vertex);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);

    bool has_vertex(const std::shared_ptr<Vertex>& vertex) const;
    bool is_boundary() const;
    void update_boundary_state();
    void update_searchable_state();
    void update_searchable_state(const std::shared_ptr<Surface>& surface);

    bool intersects_edge(const std::shared_ptr<Vertex>& vertex0, const std::shared_ptr<Vertex>& vertex1);

private:
    bool deleting_ = false;
    bool is_boundary_ = false;
    bool is_expired_ = true;

    Eigen::Vector3d center_;
    Eigen::Vector3d max_;
    Eigen::Vector3d min_;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::set<std::shared_ptr<Vertex>> vertices_;
    std::set<std::shared_ptr<Face>> faces_;
    std::set<std::shared_ptr<Surface>> surfaces_;
    std::map<std::shared_ptr<Surface>, bool> is_searchable_in_surface_;
};

bool operator<(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs);
bool operator==(const std::shared_ptr<Edge>& lhs, const std::shared_ptr<Edge>& rhs);

bool segments_intersect(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);