#pragma once

#include <memory>
#include <vector>

#include "MeshObject/Vert.hpp"
#include "MeshObject/Edge.hpp"
#include "MeshObject/Face.hpp"

class Storage : public std::enable_shared_from_this<Storage> 
{
public:
    void add_vert(Eigen::Vector3d pos);
    void add_edge(std::weak_ptr<Vert> vert1, std::weak_ptr<Vert> vert2);
    void add_face(std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3);

    void delete_vert(std::weak_ptr<Vert> vert);
    void delete_edge(std::weak_ptr<Edge> edge);
    void delete_face(std::weak_ptr<Face> face);

    int get_next_vert_id();
    int get_next_edge_id();
    int get_next_face_id();

    std::vector<std::shared_ptr<Vert>> verts_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::vector<std::shared_ptr<Face>> faces_;

private:
    int next_vert_id_ = 0;
    int next_edge_id_ = 0;
    int next_face_id_ = 0;
};