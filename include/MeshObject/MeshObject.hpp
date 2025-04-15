#pragma once
#include <memory>
#include <type_traits>

class MeshObject
{
public:
    virtual const int& get_id() const = 0;
    virtual bool is_expired() const = 0;
};

struct MeshObjectHash 
{
    template<typename T>
    std::size_t operator()(const std::shared_ptr<T>& v) const 
    {
        return std::hash<int>()(v->get_id());
    }

    template<typename T>
    std::size_t operator()(const std::weak_ptr<T>& v) const 
    {
        return std::hash<int>()(v.lock()->get_id());
    }
};

struct MeshObjectCompare 
{
    template<typename T>
    bool operator()(const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const {
        return lhs < rhs; // or any other comparison logic
    }

    template<typename T>
    bool operator()(const std::weak_ptr<T>& lhs, const std::weak_ptr<T>& rhs) const {
        return lhs.lock() < rhs.lock(); // or any other comparison logic
    }
};

bool operator<(const std::weak_ptr<MeshObject>& lhs, const std::weak_ptr<MeshObject>& rhs);
bool operator==(const std::weak_ptr<MeshObject>& lhs, const std::weak_ptr<MeshObject>& rhs);