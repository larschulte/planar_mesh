#pragma once

#include <unordered_map>
#include <list>
#include <utility>

template<typename Key, typename Value>
class LRUCache 
{
public:
    LRUCache(std::size_t max_size);

    void put(const Key& key, const Value& value);
    bool exists(const Key& key);
    const Value& get(const Key& key);

private:
    size_t max_size;
    std::list<std::pair<Key, Value>> cache_list;
    std::unordered_map<Key, typename std::list<std::pair<Key, Value>>::iterator> cache_map;
};

#include "LRUCache.tpp" // Include the implementation