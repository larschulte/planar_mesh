#pragma once

#include "LRUCache.hpp"

template<typename Key, typename Value>
LRUCache<Key, Value>::LRUCache(std::size_t max_size) : max_size(max_size) {}

template<typename Key, typename Value>
void LRUCache<Key, Value>::put(const Key& key, const Value& value) 
{
    auto it = cache_map.find(key);
    if (it != cache_map.end()) 
    {
        cache_list.erase(it->second);
        cache_map.erase(it);
    }

    cache_list.push_front(std::make_pair(key, value));
    cache_map[key] = cache_list.begin();

    if (cache_map.size() > max_size) 
    {
        auto last = cache_list.end();
        last--;
        cache_map.erase(last->first);
        cache_list.pop_back();

        // // log
        // std::cout << "Cache of size " << max_size << " is full, removing least recent element" << std::endl;
    }
}

template <typename Key, typename Value>
bool LRUCache<Key, Value>::exists(const Key& key) 
{
    // check if key exists
    bool exists = cache_map.find(key) != cache_map.end();

    // log
    if (exists) std::cout << "found in cache" << std::endl;

    // return
    return exists;
}

template<typename Key, typename Value>
const Value& LRUCache<Key, Value>::get(const Key& key) 
{
    auto it = cache_map.find(key);
    if (it == cache_map.end()) 
    {
        throw std::runtime_error("Key not found in cache");
    } 
    else 
    {
        cache_list.splice(cache_list.begin(), cache_list, it->second);
        return it->second->second;
    }
}