#include "FIFOCache.hpp"
#include <iostream>

template<typename Key, typename Value>
FIFOCache<Key, Value>::FIFOCache(std::size_t max_size) : max_size(max_size) {}

template<typename Key, typename Value>
void FIFOCache<Key, Value>::put(const Key& key, const Value& value) 
{
    if (cache_map.size() >= max_size) 
    {
        Key oldest_key = cache_queue.front();
        cache_queue.pop();
        cache_map.erase(oldest_key);

        // // log
        // std::cout << "Cache of size " << max_size << " is full, removing oldest element" << std::endl;
    }

    cache_queue.push(key);
    cache_map[key] = value;
}

template<typename Key, typename Value>
bool FIFOCache<Key, Value>::exists(const Key& key) const 
{
    // check if key exists
    bool exists = cache_map.find(key) != cache_map.end();

    // // log
    // if (exists) std::cout << "found in cache" << std::endl;

    // return
    return exists;
}

template<typename Key, typename Value>
const Value& FIFOCache<Key, Value>::get(const Key& key) 
{
    auto it = cache_map.find(key);
    if (it == cache_map.end()) 
    {
        throw std::runtime_error("Key not found in cache");
    } 
    else 
    {
        return it->second;
    }
}
