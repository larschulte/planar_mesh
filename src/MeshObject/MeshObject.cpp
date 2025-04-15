#include "MeshObject/MeshObject.hpp"
#include <stdexcept>

bool operator<(const std::weak_ptr<MeshObject>& lhs, const std::weak_ptr<MeshObject>& rhs)
{
    // when updating the third point's boundary state (due to first point deleted), the first point is already deleted. 
    if (lhs.lock()->is_expired() || rhs.lock()->is_expired()) throw std::runtime_error("Comparing expired edges");
    return lhs.lock()->get_id() < rhs.lock()->get_id();
}

bool operator==(const std::weak_ptr<MeshObject>& lhs, const std::weak_ptr<MeshObject>& rhs)
{
    if (!lhs.lock() && !rhs.lock()) return true; // true if both are nullptr
    if (!lhs.lock() || !rhs.lock()) return false; // false if either is nullptr
    return lhs.lock()->get_id() == rhs.lock()->get_id();
}