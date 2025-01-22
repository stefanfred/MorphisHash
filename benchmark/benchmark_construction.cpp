#include <chrono>
#include <iostream>
#include <bytehamster/util/XorShift64.h>
#include <tlx/cmdline_parser.hpp>
#include "BenchmarkData.h"


#include "MorphisHash.h"
#include "MorphisHashFlat.h"

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory")

bool rotate = false;
bool morphisHash = false;
bool morphisHashflat = false;
size_t numObjects = 1e6;
size_t numQueries = 1e6;
size_t leafSize = 20;
size_t relativeWidth = 4;
size_t bucketSize = 2000;
size_t threads = 1;

template<typename HashFunc>
void construct() {
    auto time = std::chrono::system_clock::now();
    long seed = std::chrono::duration_cast<std::chrono::milliseconds>(time.time_since_epoch()).count();
    bytehamster::util::XorShift64 prng(seed);
#define STRING_KEYS
#ifdef STRING_KEYS
    std::vector<std::string> keys = generateInputData(numObjects);
#else
    std::cout<<"Generating input data (Seed: "<<seed<<")"<<std::endl;
    std::vector<sux::function::hash128_t> keys;
    for (size_t i = 0; i < numObjects; i++) {
        keys.push_back(sux::function::hash128_t(prng(), prng()));
    }
#endif

    std::cout << "Constructing" << std::endl;
    sleep(1);
    auto beginConstruction = std::chrono::high_resolution_clock::now();
    HashFunc hashFunc(keys, bucketSize, threads);
    unsigned long constructionDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginConstruction).count();

    std::cout << "Testing" << std::endl;
    std::vector<bool> taken(keys.size(), false);
    for (size_t i = 0; i < keys.size(); i++) {
        size_t hash = hashFunc(keys.at(i));
        if (taken[hash]) {
            std::cerr << "Collision by key " << i << "!" << std::endl;
            exit(1);
        } else if (hash > numObjects) {
            std::cerr << "Out of range!" << std::endl;
            exit(1);
        }
        taken[hash] = true;
    }

    std::cout << "Preparing query plan" << std::endl;
#ifdef STRING_KEYS
    std::vector<std::string> queryPlan;
#else
    std::vector<sux::function::hash128_t> queryPlan;
#endif
    queryPlan.reserve(numQueries);
    for (size_t i = 0; i < numQueries; i++) {
        queryPlan.push_back(keys[prng(numObjects)]);
    }

    std::cout << "Querying" << std::endl;
    sleep(1);
    auto beginQueries = std::chrono::high_resolution_clock::now();
    for (const auto &key: queryPlan) {
        size_t retrieved = hashFunc(key);
        DO_NOT_OPTIMIZE(retrieved);
    }
    auto queryDurationMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - beginQueries).count();

    hashFunc.printBits();

#ifdef SIMD
    std::string method = "SIMD";
#else
    std::string method = "plain";
#endif
    if (rotate) {
        method += "Rotate";
    } else if (morphisHash) {
        method += "morphisHash";
    } else if (morphisHashflat) {
        method += "morphisHashFlat";
    }
    std::cout << "RESULT"
              << " method=" << method
              << " l=" << leafSize
              << " b=" << bucketSize
              << " N=" << numObjects
              << " w=" << relativeWidth
              << " threads=" << threads
              << " numQueries=" << numQueries
              << " queryTimeMilliseconds=" << queryDurationMs
              << " constructionTimeMilliseconds=" << constructionDurationMs
              << " bitsPerElement=" << (double) hashFunc.getBits() / numObjects
              << std::endl;
}

template<template<size_t, size_t> class HashFunc, size_t I, size_t WS>
void dispatchWidth() {
    if (WS == relativeWidth) {
        construct<HashFunc<I, WS>>();
    } else if constexpr (WS == 0) {
        std::cerr << "The relativeWidth " << relativeWidth << " was not compiled into this binary." << std::endl;
    } else {
        dispatchWidth<HashFunc, I, WS - 1>();
    }
}

template<template<size_t, size_t> class HashFunc, size_t I>
void dispatchLeafSize() {
    if constexpr (I <= 1) {
        std::cerr << "The leafSize " << leafSize << " was not compiled into this binary." << std::endl;
    } else if (I == leafSize) {
        dispatchWidth<HashFunc, I, std::min(morphishash::MAX_DIFF, I)>();
    } else {
        dispatchLeafSize<HashFunc, I - 1>();
    }
}

int main(int argc, const char *const *argv) {
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numObjects", numObjects, "Number of objects to construct with");
    cmd.add_bytes('q', "numQueries", numQueries, "Number of queries to measure");
    cmd.add_bytes('l', "leafSize", leafSize, "Leaf size to construct");
    cmd.add_bytes('w', "relativeWidth", relativeWidth, "Bits less than leaf size for fixed width retrieval.");
    cmd.add_bytes('b', "bucketSize", bucketSize, "Bucket size to construct");
    cmd.add_bool('r', "rotate", rotate, "Apply rotation fitting");
    cmd.add_bool('m', "morphisHash", morphisHash, "morphisHash");
    cmd.add_bool('f', "morphisHashFlat", morphisHashflat, "morphisHash flat");
    cmd.add_bytes('t', "threads", threads, "Number of threads");

    if (!cmd.process(argc, argv)) {
        return 1;
    }

    if (morphisHash) {
        dispatchLeafSize<morphishash::MorphisHash, morphishash::MAX_LEAF_SIZE>();
    } else if (morphisHashflat) {
        dispatchLeafSize<morphishash::MorphisHashFlat, morphishash::MAX_LEAF_SIZE>();
    }
    return 0;
}
