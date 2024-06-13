#include "MeshObject/UnionFind.hpp"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <set>
#include <algorithm>

#include "MeshObject/Vertex.hpp"
#include "MeshObject/Edge.hpp"

void UnionFind::add_vertices(const std::set<std::shared_ptr<Vertex>>& vertices) 
{
    for (const auto& vertex : vertices) 
    {
        parent[vertex] = vertex;
        rank[vertex] = 0;
    }
}

void UnionFind::add_edges(const std::set<std::shared_ptr<Edge>>& edges) 
{
    for (const auto& edge : edges) 
    {
        const std::shared_ptr<Vertex>& u = edge->get_vertex(0);
        const std::shared_ptr<Vertex>& v = edge->get_vertex(1);
        make_share_root(u, v);
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
        groups[pair.second].push_back(pair.first);
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