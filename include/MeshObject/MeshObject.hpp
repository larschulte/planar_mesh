#pragma once
#include <memory>
#include <type_traits>

class MeshObject
{
public:
    virtual const int& get_id() const = 0;
};

struct MeshObjectHash {
    template<typename T>
    std::size_t operator()(const std::shared_ptr<T>& v) const 
    {
        return std::hash<int>()(v->get_id());
    }
};