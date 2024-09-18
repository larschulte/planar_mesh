#pragma once

#include <memory>
#include <unordered_set>
#include <Eigen/Dense> 

#include "MeshObject/RRSTree.hpp"
#include "MeshObject/TriangleBVH.hpp"

#include "utilities/queue_or_stack.hpp"

#include "MeshObject/Settings.hpp"

#include <omp.h>
#include <mutex>

// forward declarations
class Vertex;
class Edge;
class Face;
class Surface;
class GenericPoint;
class InteriorPoint;

class Storage : public std::enable_shared_from_this<Storage> 
{
public: // to user
    Storage();
    ~Storage();

    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position);
    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius);
    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<Surface>& surface, const std::shared_ptr<GenericPoint>& generic_point);
    const std::shared_ptr<Edge>& add_edge(const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    const std::shared_ptr<Face>& add_face(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3);
    const std::shared_ptr<Surface>& add_surface();
    const std::shared_ptr<GenericPoint>& add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    const std::shared_ptr<GenericPoint>& add_generic_point(const std::shared_ptr<Vertex>& vertex);
    const std::shared_ptr<GenericPoint>& add_generic_point(const std::shared_ptr<InteriorPoint>& interior_point);
    const std::shared_ptr<InteriorPoint>& add_interior_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    const std::shared_ptr<InteriorPoint>& add_interior_point(const std::shared_ptr<GenericPoint>& generic_point);
    
    void delete_vertex(const std::shared_ptr<Vertex>& vertex);
    void delete_edge(const std::shared_ptr<Edge>& edge);
    void delete_face(const std::shared_ptr<Face>& face);
    void delete_surface(const std::shared_ptr<Surface>& surface);
    void delete_generic_point(const std::shared_ptr<GenericPoint>& genertic_point);
    void delete_interior_point(const std::shared_ptr<InteriorPoint>& interior_point);

    void add_to_main_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    void split_main_queue_into_smaller_queues();

    void add_points_in_smaller_repeated_queues_to_main_queue();
    void add_points_in_smaller_abort_queues_to_main_queue();

    void add_to_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin);
    void add_to_queue(const std::shared_ptr<GenericPoint>& generic_point);
    void add_to_queue(const std::shared_ptr<InteriorPoint>& interior_point);
    void add_to_queue(const std::shared_ptr<Vertex>& vertex);
    void add_to_abort_queue(const std::shared_ptr<GenericPoint>& generic_point);
    std::shared_ptr<GenericPoint> pop_from_queue();
    unsigned int get_queue_size();
    unsigned int get_repeated_queue_size();
    unsigned int get_abort_queue_size();

    bool can_reverse_radius_search();
    RRSReturnType reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& result);
    RRSReturnType reverse_radius_search_find_node(const Eigen::Vector3d& point, std::shared_ptr<RRSNode>& return_node);
    BVHReturnType face_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& result);
    BVHReturnType face_intersection_search_find_node(const Eigen::Vector3d& origin, const Eigen::Vector3d& point, std::shared_ptr<Node>& return_node);

    const std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>& get_vertices() const;
    const std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>& get_edges() const;
    const std::shared_ptr<Edge>& get_edge(std::shared_ptr<Vertex> vertex1, std::shared_ptr<Vertex> vertex2) const;
    const std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>& get_faces() const;
    const std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash>& get_surfaces() const;
    const std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash>& get_generic_points() const;
    const std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>& get_interior_points() const;

    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices();
    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map() const;
    bool is_expired() const;

    void print_rrs() const;
    void print_bvh() const;
    void check_tree_rebuild();
    void rebuild_tree();

    unsigned int get_bvh_size() const;

private: // to Vertex and Face class
    friend class Vertex;
    friend class Face;
    void add_searchable_vertex(const std::shared_ptr<Vertex>& vertex);
    void remove_searchable_vertex(const std::shared_ptr<Vertex>& vertex);

    void add_searchable_face(const std::shared_ptr<Face>& face);
    void remove_searchable_face(const std::shared_ptr<Face>& face);

public: // to MeshObject class
    int get_next_vertex_id();
    int get_next_edge_id();
    int get_next_face_id();
    int get_next_surface_id();
    int get_next_generic_point_id();
    int get_next_interior_point_id();

private:
    static Settings settings_;

    bool is_expired_ = true;

    RRSTree rrs_tree_;
    TriangleBVH triangle_bvh_;

    std::queue<std::shared_ptr<GenericPoint>> main_queue_;
    std::queue<std::shared_ptr<GenericPoint>> main_repeated_queue_;
    std::vector<queue_or_stack<std::shared_ptr<GenericPoint>>> smaller_queues_;
    std::vector<queue_or_stack<std::shared_ptr<GenericPoint>>> smaller_repeated_queues_;
    std::vector<queue_or_stack<std::shared_ptr<GenericPoint>>> smaller_abort_queues_;
    unsigned int num_delete_before_put_to_repeated_queue_ = 2;

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points_;
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> genertic_points_;

    std::mutex vertices_mutex_;
    std::mutex edges_mutex_;
    std::mutex faces_mutex_;
    std::mutex surfaces_mutex_;
    std::mutex interior_points_mutex_;
    std::mutex genertic_points_mutex_;

    std::atomic<int> next_vertex_id_{0};
    std::atomic<int> next_edge_id_{0};
    std::atomic<int> next_face_id_{0};
    std::atomic<int> next_surface_id_{0};
    std::atomic<int> next_genertic_point_id_{0};
    std::atomic<int> next_interior_point_id_{0};
};