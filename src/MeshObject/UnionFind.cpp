#include "MeshObject/UnionFind.hpp"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <unordered_set>
#include <algorithm>

#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"

void UnionFind::add_vertices(const std::unordered_set<std::weak_ptr<Vertex>, MeshObjectHash>& vertices) 
{
    for (const auto& vertex : vertices) 
    {
        parent[vertex.lock()] = vertex.lock();
        rank[vertex.lock()] = 0;
    }
}

void UnionFind::add_edges(const std::unordered_set<std::weak_ptr<Edge>, MeshObjectHash>& edges) 
{
    for (const auto& edge : edges) 
    {
        // edge locked
        std::shared_ptr<Edge> edge_locked = edge.lock();

        // skip if nullptr
        if (!edge_locked) continue;

        std::weak_ptr<Vertex> u = edge_locked->get_vertex_weak_ptr(0);
        std::weak_ptr<Vertex> v = edge_locked->get_vertex_weak_ptr(1);

        std::shared_ptr<Vertex> u_locked = u.lock();
        std::shared_ptr<Vertex> v_locked = v.lock();

        // skip if nullptr
        if (!u_locked || !v_locked) continue;

        make_share_root(u_locked, v_locked);
    }
}

// make u and v share the same root
void UnionFind::make_share_root(const std::shared_ptr<Vertex>& u, const std::shared_ptr<Vertex>& v) 
{
    const std::shared_ptr<Vertex>& root_u = find_root(u);
    const std::shared_ptr<Vertex>& root_v = find_root(v);

    if (root_u != root_v) 
    {
        // Attach smaller tree under root of larger tree
        if (rank.at(root_u) > rank.at(root_v)) 
        {
            parent.at(root_v) = root_u;
        } 
        else if (rank.at(root_u) < rank.at(root_v)) 
        {
            parent.at(root_u) = root_v;
        } 
        else 
        {
            parent.at(root_v) = root_u;
            rank.at(root_u)++;
        }
    }
}

// find root vertex of vertex u 
const std::shared_ptr<Vertex>& UnionFind::find_root(const std::shared_ptr<Vertex>& u) 
{
    if (parent.at(u) != u) 
    {
        parent.at(u) = find_root(parent.at(u));  // Path compression
    }
    return parent.at(u);
}

std::vector<std::pair<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Vertex>>>> UnionFind::compute_sorted_grouped_vertices() 
{
    // group
    std::unordered_map<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Vertex>>> groups;
    for (const auto& pair : parent) 
    {
        groups[find_root(pair.second)].push_back(pair.first);
    }

    // sorted group
    std::vector<std::pair<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Vertex>>>> sorted_groups(groups.begin(), groups.end());
    std::sort(sorted_groups.begin(), sorted_groups.end(), 
        [](const auto& a, const auto& b)
        {
            return a.second.size() > b.second.size();
        });

    // return
    return sorted_groups;
}

void UnionFind::print_sorted_grouped_vertices()
{
    std::vector<std::pair<std::shared_ptr<Vertex>, std::vector<std::shared_ptr<Vertex>>>> sorted_groups = compute_sorted_grouped_vertices();

    std::cout << "Total number of groups: " << sorted_groups.size() << std::endl;

    for (const auto& pair : sorted_groups) 
    {
        std::cout << "root: " << pair.first->get_id() <<  ", size: " << pair.second.size() << std::endl;
        for (const auto& vertex : pair.second) 
        {
            std::cout << vertex->get_id() << " ";
        }
        std::cout << std::endl;
    }
}