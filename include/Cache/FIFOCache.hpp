#pragma once

#include <unordered_map>
#include <queue>
#include <utility>
#include <stdexcept>

template<typename Key, typename Value>
class FIFOCache {
public:
    FIFOCache(std::size_t max_size);

    void put(const Key& key, const Value& value);
    bool exists(const Key& key) const;
    const Value& get(const Key& key);

private:
    std::size_t max_size;
    std::queue<Key> cache_queue;
    std::unordered_map<Key, Value> cache_map;
};

#include "FIFOCache.tpp" // Include the implementation
