#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>

// Forward declarations
class Edge;
class Storage;

class Face : public std::enable_shared_from_this<Face> 
{
protected:
    friend class Storage;
    void initialize_(std::weak_ptr<Storage> storage, std::weak_ptr<Edge> edge1, std::weak_ptr<Edge> edge2, std::weak_ptr<Edge> edge3);
    void delete_();

public:
    void cascade_delete_from_edge(std::weak_ptr<Edge> edge);

private:
    int id_;
    std::weak_ptr<Edge> edge1_;
    std::weak_ptr<Edge> edge2_;
    std::weak_ptr<Edge> edge3_;
    std::weak_ptr<Storage> storage_;
};