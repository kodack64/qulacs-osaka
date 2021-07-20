#pragma once

#include <omp.h>

#include "state.hpp"
#include "type.hpp"

template <class T, class V>
class DllExport lockedHash {
    std::unordered_map<T, V> theMap;
    omp_lock_t theLock;

public:
    lockedHash(omp_sync_hint_t hint) {
        omp_init_lock_with_hint(&theLock, hint);
    }
    void insert(T key, V value) {
        omp_set_lock(&theLock);  // Claim the lock
        theMap.insert({key, value});
        omp_unset_lock(&theLock);  // Release the lock
    }
    V lookup(T key) {
        omp_set_lock(&theLock);  // Claim the lock
        auto result = theMap.find(key);
        omp_unset_lock(&theLock);  // Release the lock
        return result == theMap.end() ? 0 : result->second;
    }
};