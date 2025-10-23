/**
 * @file ParallelHashMap.h
 * @brief A thread-safe hash map supporting insertion and lookup operations.
 * @details The parallel hash map is built on top of a fixed-sized hash map
 *    object and features OpenMP concurrency structures. The underlying
 *    fixed-sized hash map handles collisions with chaining.
 * @date June 6, 2015
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef __PARALLEL_HASH_MAP__
#define __PARALLEL_HASH_MAP__
#include <iostream>
#include <stdexcept>
#include <functional>
#include <omp.h>

namespace antmoc
{

/**
 * @class FixedHashMap ParallelHashMap.h "src/ParallelHashMap.h"
 * @brief A fixed-size hash map supporting insertion and lookup operations.
 * @details The FixedHashMap class supports insertion and lookup operations
 *    but not deletion as deletion is not needed in the OpenMOC application.
 *    This hash table uses chaining for collisions and does not incorporate
 *    concurrency objects except for tracking the number of entries in the
 *    table for which an atomic increment is used. This hash table is not
 *    thread safe but is used as a building block for the ParallelHashMap
 *    class. This table guarantees O(1) insertions and lookups on average.
 */
template <class K, class V>
class FixedHashMap {
  struct node {
    node(K k_in, V v_in) : key(k_in), value(v_in), next(nullptr) {}
    K key;
    V value;
    node *next;
  };

  private:
    size_t _M;      /* table size */
    size_t _N;      /* number of elements present in table */
    node ** _buckets;   /* buckets of values stored in nodes */

  public:

    FixedHashMap(size_t M = 64);
    virtual ~FixedHashMap();
    bool contains(K& key);
    V& at(K& key);
    void insert(K key, V value);
    long insert_and_get_count(K key, V value);
    size_t size();
    bool empty();
    size_t bucket_count();
    K* keys();
    V* values();
    void clear();
    void print_buckets();
};


/**
 * @class ParallelHashMap ParallelHashMap.h "src/ParallelHashMap.h"
 * @brief A thread-safe hash map supporting insertion and lookup operations.
 * @details The ParallelHashMap class is built on top of the FixedHashMap
 *    class, supporting insertion and lookup operations but not deletion as
 *    deletion is not needed in the OpenMOC application. This hash table uses
 *    chaining for collisions, as defined in FixedHashMap. It offers lock
 *    free lookups in O(1) time on average and fine-grained locking for
 *    insertions in O(1) time on average as well. Resizing is conducted
 *    periodically during inserts, although the starting table size can be
 *    chosen to limit the number of resizing operations.
 */
template <class K, class V>
class ParallelHashMap {

  /* padded pointer to hash table to avoid false sharing */
  struct paddedPointer {
    volatile long pad_L1;
    volatile long pad_L2;
    volatile long pad_L3;
    volatile long pad_L4;
    volatile long pad_L5;
    volatile long pad_L6;
    volatile long pad_L7;
    FixedHashMap<K,V> volatile* value;
    volatile long pad_R1;
    volatile long pad_R2;
    volatile long pad_R3;
    volatile long pad_R4;
    volatile long pad_R5;
    volatile long pad_R6;
    volatile long pad_R7;
    volatile long pad_R8;
  };

  private:
    FixedHashMap<K,V> *_table;
    paddedPointer *_announce;
    size_t _num_threads;
    size_t _N;
    omp_lock_t * _locks;
    size_t _num_locks;
    bool _fixed_size;
    void resize();

  public:
    ParallelHashMap(size_t M = 64, size_t L = 64);
    virtual ~ParallelHashMap();
    bool contains(K& key);
    V at(K& key);
    void update(K& key, V value);
    void insert(K key, V value);
    long insert_and_get_count(K key, V value);
    size_t size();
    bool empty();
    size_t bucket_count();
    size_t num_locks();
    K* keys();
    V* values();
    void clear();
    void setFixedSize();
    void print_buckets();
    void realloc(size_t M);
};


} /* namespace antmoc */

#endif
