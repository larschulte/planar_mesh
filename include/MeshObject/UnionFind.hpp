#pragma once

#include <unordered_set>
#include "MeshObject/Vertex.hpp"

class Vertex;
class Edge;

class UnionFind 
{
public:
    void add_vertices(const std::unordered_set<std::weak_ptr<Vertex>, MeshObjectHash>& vertices);
    void add_edges(const std::unordered_set<std::weak_ptr<Edge>, MeshObjectHash>& edges);
    std::vector<std::pair<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Vertex>>>> compute_sorted_grouped_vertices();
    void print_sorted_grouped_vertices();

private:
    void make_share_root(const std::shared_ptr<Vertex>& u, const std::shared_ptr<Vertex>& v);
    const std::shared_ptr<Vertex>& find_root(const std::shared_ptr<Vertex>& u); 

    std::unordered_map<std::shared_ptr<Vertex>, std::shared_ptr<Vertex>> parent;  // Parent array
    std::unordered_map<std::shared_ptr<Vertex>, int> rank;    // Rank array
};