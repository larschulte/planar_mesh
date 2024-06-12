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
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    void initialize_(const std::shared_ptr<Storage>& storage, const Eigen::Vector3d& position, const Eigen::Vector3d& origin, const double& radius);
    void delete_();

public:
    const int& get_id() const;
    const Eigen::Vector3d& get_position() const;
    const Eigen::Vector3d& get_origin() const;
    const std::shared_ptr<Surface>& get_surface() const;
    const std::set<std::shared_ptr<Edge>>& get_edges() const;

    void try_update_surface_projection();
    const Eigen::Vector3d& get_projected_position();
    const double& get_projected_distance();

    Eigen::Vector2d get_surface_coordinate();
    bool is_expired() const;
    bool is_to_be_deleted() const;

    void connect(const std::shared_ptr<Edge>& edge);
    void connect(const std::shared_ptr<Face>& face);
    void connect(const std::shared_ptr<Surface>& surface);
    void disconnect(const std::shared_ptr<Edge>& edge);
    void disconnect(const std::shared_ptr<Face>& face);
    void disconnect(const std::shared_ptr<Surface>& surface);

    void update_boundary_state();

public: // for reverse radius search
    void set_reverse_radius_search_radius(double radius);
    Eigen::Vector3d get_min() const;
    Eigen::Vector3d get_max() const;
    double get_radius() const;
    bool contains(const Eigen::Vector3d& point) const;
    bool approx_contains(const Eigen::Vector3d& point) const;

private: // for reverse radius search
    double reverse_search_radius_;
    Eigen::Vector3d min_;
    Eigen::Vector3d max_;

private:
    bool deleting_ = false;
    bool is_boundary_ = false;
    bool is_searchable_ = false;
    bool is_expired_ = true;
    bool is_to_be_deleted_ = false;

    int id_;
    std::shared_ptr<Storage> storage_;

    std::set<std::shared_ptr<Edge>> edges_;
    std::set<std::shared_ptr<Face>> faces_;
    std::set<std::shared_ptr<Surface>> surfaces_;

    Eigen::Matrix3d eigenvectors_used_;
    Eigen::Vector2d surface_coordinate_;
    Eigen::Vector3d normal_used_;
    Eigen::Vector3d projected_position_;
    double projected_distance_;
    

    Eigen::Vector3d position_;
    Eigen::Vector3d origin_;
};

bool operator<(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator<=(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);
bool operator==(const std::shared_ptr<Vertex>& lhs, const std::shared_ptr<Vertex>& rhs);