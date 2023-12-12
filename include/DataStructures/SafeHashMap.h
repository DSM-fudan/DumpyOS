//
// Created by pengwang5 on 2022/12/12.
//

#ifndef FADAS_SAFEHASHMAP_H
#define FADAS_SAFEHASHMAP_H
#include <cstdint>
#include <functional>
#include <iostream>
#include <mutex>
#include <shared_mutex>

template <typename K, typename V> class HashNode
{
public:
    HashNode()
    {
    }
    HashNode(K key_, V value_) : key(key_), value(value_)
    {
    }
    ~HashNode()
    {
        next = nullptr;
    }

    HashNode(const HashNode&) = delete;
    HashNode(HashNode&&)      = delete;
    HashNode& operator=(const HashNode&) = delete;
    HashNode& operator=(HashNode&&) = delete;


    const K &getKey() const
    {
        return key;
    }
    void setValue(V value_)
    {
        value = value_;
    }
    const V &getValue() const
    {
        return value;
    }

    HashNode *next = nullptr; // Pointer to the next node in the same bucket
private:
    K key;   // the hash key
    V value; // the value corresponding to the key
};

template <typename K, typename V> class HashBucket
{
public:
    HashBucket()
    {
    }

    ~HashBucket() // delete the bucket
    {
        clear();
    }

    // Function to find an entry in the bucket matching the key
    // If key is found, the corresponding value is copied into the parameter "value" and function returns true.
    // If key is not found, function returns false
    bool find(const K &key, V &value) const
    {
        // A shared mutex is used to enable mutiple concurrent reads
        std::shared_lock lock(mutex_);
        HashNode<K, V> *node = head;

        while (node != nullptr)
        {
            if (node->getKey() == key)
            {
                value = node->getValue();
                return true;
            }
            node = node->next;
        }
        return false;
    }

    // Function to insert into the bucket
    // If key already exists, update the value, else insert a new node in the bucket with the <key, value> pair
    void insert(const K &key, const V &value)
    {
        // Exclusive lock to enable single write in the bucket
        std::unique_lock lock(mutex_);
        HashNode<K, V> *prev = nullptr;
        HashNode<K, V> *node = head;

        while (node != nullptr && node->getKey() != key)
        {
            prev = node;
            node = node->next;
        }

        if (nullptr == node) // New entry, create a node and add to bucket
        {
            if (nullptr == head)
            {
                head = new HashNode<K, V>(key, value);
            }
            else
            {
                prev->next = new HashNode<K, V>(key, value);
            }
        }
        else
        {
            node->setValue(value); // Key found in bucket, update the value
        }
    }

    // Function to remove an entry from the bucket, if found
    void erase(const K &key)
    {
        // Exclusive lock to enable single write in the bucket
        std::unique_lock lock(mutex_);
        HashNode<K, V> *prev = nullptr;
        HashNode<K, V> *node = head;

        while (node != nullptr && node->getKey() != key)
        {
            prev = node;
            node = node->next;
        }

        if (nullptr == node) // Key not found, nothing to be done
        {
            return;
        }
        else // Remove the node from the bucket
        {
            if (head == node)
            {
                head = node->next;
            }
            else
            {
                prev->next = node->next;
            }
            delete node; // Free up the memory
        }
    }

    // Function to clear the bucket
    void clear()
    {
        // Exclusive lock to enable single write in the bucket
        std::unique_lock lock(mutex_);
        HashNode<K, V> *prev = nullptr;
        HashNode<K, V> *node = head;
        while (node != nullptr)
        {
            prev = node;
            node = node->next;
            delete prev;
        }
        head = nullptr;
    }

private:
    HashNode<K, V> *head = nullptr;         // The head node of the bucket
    mutable std::shared_timed_mutex mutex_; // The mutex for this bucket
};

constexpr size_t HASH_SIZE_DEFAULT = 1031;      // A prime number as hash size gives a better distribution of values in buckets

// The class representing the hash map.
// It is expected for user defined types, the hash function will be provided.
// By default, the std::hash function will be used
// If the hash size is not provided, then a default size of 1031 will be used
// The hash table itself consists of an array of hash buckets.
// Each hash bucket is implemented as singly linked list with the head as a dummy node created
// during the creation of the bucket. All the hash buckets are created during the construction of the map.
// Locks are taken per bucket, hence multiple threads can write simultaneously in different buckets in the hash map
template <typename K, typename V, typename F = std::hash<K>> class SafeHashMap
{
public:
    SafeHashMap(size_t hashSize_ = HASH_SIZE_DEFAULT) : hashSize(hashSize_)
    {
        hashTable = new HashBucket<K, V>[hashSize]; // create the hash table as an array of hash buckets
    }

    ~SafeHashMap()
    {
        delete[] hashTable;
    }
    // Copy and Move of the HashMap are not supported at this moment
    SafeHashMap(const SafeHashMap&) = delete;
    SafeHashMap(SafeHashMap&&)      = delete;
    SafeHashMap& operator=(const SafeHashMap&) = delete;
    SafeHashMap& operator=(SafeHashMap&&) = delete;

    // Function to find an entry in the hash map matching the key.
    // If key is found, the corresponding value is copied into the parameter "value" and function returns true.
    // If key is not found, function returns false.
    bool find(const K &key, V &value) const
    {
        size_t hashValue = hashFn(key) % hashSize;
        return hashTable[hashValue].find(key, value);
    }

    // Function to insert into the hash map.
    // If key already exists, update the value, else insert a new node in the bucket with the <key, value> pair.
    void insert(const K &key, const V &value)
    {
        size_t hashValue = hashFn(key) % hashSize;
        hashTable[hashValue].insert(key, value);
    }

    // Function to remove an entry from the bucket, if found
    void erase(const K &key)
    {
        size_t hashValue = hashFn(key) % hashSize;
        hashTable[hashValue].erase(key);
    }

    // Function to clean up the hasp map, i.e., remove all entries from it
    void clear()
    {
        for (size_t i = 0; i < hashSize; i++)
        {
            (hashTable[i]).clear();
        }
    }

private:
    HashBucket<K, V> *hashTable;
    F hashFn;
    const size_t hashSize;
};

#endif //FADAS_SAFEHASHMAP_H
