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

    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position, double distance_travelled);
    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<Surface>& surface, const Eigen::Vector3d& origin, const Eigen::Vector3d& position, const double& radius, double distance_travelled);
    const std::shared_ptr<Vertex>& add_vertex(const std::shared_ptr<Surface>& surface, const std::shared_ptr<GenericPoint>& generic_point);
    const std::shared_ptr<Edge>& add_edge(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2);
    const std::shared_ptr<Face>& add_face(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Vertex>& vertex1, const std::shared_ptr<Vertex>& vertex2, const std::shared_ptr<Vertex>& vertex3);
    const std::shared_ptr<Face>& add_face(
        const std::shared_ptr<Surface>& surface, 
        const std::shared_ptr<Vertex>& vertex1, 
        const std::shared_ptr<Vertex>& vertex2, 
        const std::shared_ptr<Vertex>& vertex3,
        const std::shared_ptr<Edge>& edge1,
        const std::shared_ptr<Edge>& edge2,
        const std::shared_ptr<Edge>& edge3);
    const std::shared_ptr<Surface>& add_surface();
    const std::shared_ptr<GenericPoint>& add_generic_point(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled);
    const std::shared_ptr<GenericPoint>& add_generic_point(const std::shared_ptr<Vertex>& vertex);
    const std::shared_ptr<GenericPoint>& add_generic_point(const std::shared_ptr<InteriorPoint>& interior_point);
    const std::shared_ptr<InteriorPoint>& add_interior_point(const std::shared_ptr<Surface>& surface, const std::shared_ptr<Face>& face, const std::shared_ptr<GenericPoint>& generic_point);
    
    void delete_vertex(const std::shared_ptr<Vertex> vertex);
    void delete_edge(const std::shared_ptr<Edge> edge);
    void delete_face(const std::shared_ptr<Face> face);
    void delete_surface(const std::shared_ptr<Surface> surface);
    void delete_generic_point(const std::shared_ptr<GenericPoint> genertic_point);
    void delete_interior_point(const std::shared_ptr<InteriorPoint> interior_point);

    void add_to_main_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled);
    void split_main_queue_into_smaller_queues();
    void split_main_queue_into_smaller_queues_by_angle(Eigen::Vector3d origin);
    void split_main_queue_into_smaller_queues_by_contention();
    void print_main_queue_stats();

    void add_points_in_smaller_repeated_queues_to_main_queue();
    void add_points_in_smaller_abort_queues_to_main_queue();

    void add_to_queue(const Eigen::Vector3d& position, const Eigen::Vector3d& origin, double distance_travelled);
    void add_to_queue(const std::shared_ptr<GenericPoint>& generic_point);
    void add_to_queue(const std::shared_ptr<InteriorPoint>& interior_point);
    void add_to_queue(const std::shared_ptr<Vertex>& vertex);
    void add_to_abort_queue(const std::shared_ptr<GenericPoint>& generic_point);
    std::shared_ptr<GenericPoint> pop_from_queue();
    unsigned int get_main_queue_size();
    unsigned int get_queue_size();
    unsigned int get_repeated_queue_size();
    unsigned int get_abort_queue_size();
    void clear_queues();

    void cleanup_surfaces();
    void remove_non_manifold_edges();
    void remove_non_manifold_vertices();
    void remove_non_manifold_faces();
    void update_radius();

    void add_points_in_add_searchable_vertex_queue();

    void add_or_remove_vertices_from_rrs_tree();
    void add_or_remove_faces_from_bvh_tree();
    void add_or_remove_edges_from_edgeBVH_tree();

    void add_vertex_to_be_deleted(const std::shared_ptr<Vertex>& vertex);
    void add_edge_to_be_deleted(const std::shared_ptr<Edge>& edge);
    void add_face_to_be_deleted(const std::shared_ptr<Face>& face);
    void add_interior_point_to_be_deleted(const std::shared_ptr<InteriorPoint>& interior_point);

    void delete_to_be_deleted_repeatedly();
    void delete_to_be_deleted();

    void add_vertex_that_have_deleted_publishers(const std::shared_ptr<Vertex>& vertex);
    void add_vertex_that_have_added_publishers(const std::shared_ptr<Vertex>& vertex);
    void update_vertices_that_have_deleted_publishers();
    void update_vertices_that_have_added_publishers();

    void add_to_set_of_vertices_to_update_rrs_tree(const std::shared_ptr<Vertex>& vertex);
    void add_to_set_of_faces_to_update_bvh_tree(const std::shared_ptr<Face>& face);
    void add_to_set_of_edge_to_update_edgeBVH_tree(const std::shared_ptr<Edge>& edge, const std::shared_ptr<Surface>& surface);

    RRSReturnType reverse_radius_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Vertex>>& result);
    BVHReturnType face_intersection_search(const std::shared_ptr<GenericPoint>& generic_point, std::vector<std::shared_ptr<Face>>& result);

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> get_boundary_vertices() const;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> get_boundary_edges() const;
    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> get_vertices() const;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> get_edges() const;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> get_faces() const;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> get_surfaces() const;
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> get_generic_points() const;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> get_interior_points() const;
    unsigned int get_boundary_vertices_size() const;
    unsigned int get_boundary_edges_size() const;
    unsigned int get_vertices_size() const;
    unsigned int get_edges_size() const;
    unsigned int get_faces_size() const;
    unsigned int get_surfaces_size() const;
    unsigned int get_generic_points_size() const;
    unsigned int get_interior_points_size() const;

    std::vector<std::shared_ptr<Vertex>> get_rrs_vertices();
    std::map<std::shared_ptr<Vertex>, int> get_vertex_to_cloud_indices_map() const;
    bool is_expired() const;

    unsigned int get_rrs_size() const;
    unsigned int get_bvh_size() const;

    void print_rrs() const;
    void print_bvh() const;

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

    std::vector<std::queue<std::shared_ptr<Vertex>>> smaller_add_searchable_vertices_queue_;
    std::vector<std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>> smaller_set_of_vertices_to_update_rrs_tree;
    std::vector<std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>> smaller_set_of_faces_to_update_rrs_tree;
    std::vector<std::vector<std::pair<std::shared_ptr<Edge>, std::shared_ptr<Surface>>> > smaller_set_of_edges_to_update_edgeBVH_tree;

    std::vector<std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>> thread_vertices_to_be_deleted_;
    std::vector<std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash>> thread_edges_to_be_deleted_;
    std::vector<std::unordered_set<std::shared_ptr<Face>, MeshObjectHash>> thread_faces_to_be_deleted_;
    std::vector<std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash>> thread_interior_points_to_be_deleted_;

    std::vector<std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>> thread_vertices_that_have_deleted_publishers_;
    std::vector<std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash>> thread_vertices_that_have_added_publishers_;

    std::unordered_set<std::shared_ptr<Vertex>, MeshObjectHash> vertices_;
    std::unordered_set<std::shared_ptr<Edge>, MeshObjectHash> edges_;
    std::unordered_set<std::shared_ptr<Face>, MeshObjectHash> faces_;
    std::unordered_set<std::shared_ptr<Surface>, MeshObjectHash> surfaces_;
    std::unordered_set<std::shared_ptr<InteriorPoint>, MeshObjectHash> interior_points_;
    std::unordered_set<std::shared_ptr<GenericPoint>, MeshObjectHash> genertic_points_;

    // read write lock
    mutable std::shared_mutex rwlock_vertices_;
    mutable std::shared_mutex rwlock_edges_;
    mutable std::shared_mutex rwlock_faces_;
    mutable std::shared_mutex rwlock_surfaces_;
    mutable std::shared_mutex rwlock_interior_points_;
    mutable std::shared_mutex rwlock_genertic_points_;

    std::atomic<int> next_vertex_id_{0};
    std::atomic<int> next_edge_id_{0};
    std::atomic<int> next_face_id_{0};
    std::atomic<int> next_surface_id_{0};
    std::atomic<int> next_genertic_point_id_{0};
    std::atomic<int> next_interior_point_id_{0};
};