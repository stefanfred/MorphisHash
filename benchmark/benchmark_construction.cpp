#include <chrono>
#include <iostream>
#include <XorShift64.h>
#include <tlx/cmdline_parser.hpp>
#include "BenchmarkData.h"

#ifdef SIMD
#include "SIMDShockHash.hpp"
template <size_t l>
using ShockHash = shockhash::SIMDShockHash<l, false>;
template <size_t l>
using ShockHashRotate = shockhash::SIMDShockHash<l, true>;
#else
//#define STATS
//#define MORESTATS
#include "MorphisHash.h"

template<size_t l>
using ShockHash = morphishash::MorphisHash<l, false>;
template<size_t l>
using ShockHashRotate = morphishash::MorphisHash<l, true>;
#endif

#include "MorphisHash2.h"
#include "MorphisHash2Flat.h"

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory")

bool rotate = false;
bool shockhash2 = false;
bool shockhash2flat = false;
bool useBurr = false;
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
    util::XorShift64 prng(seed);
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
    } else if (shockhash2) {
        method += (useBurr ? "2Burr" : "2Fixed");
    } else if (shockhash2flat) {
        method += (useBurr ? "2flatBurr" : "2flatFixed");
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

template<template<size_t, bool, size_t> class RecSplit, size_t I, size_t WS>
void dispatchWidth() {
    if (useBurr) {
        construct<RecSplit<I, true, 0>>();
    } else if (WS == relativeWidth) {
        construct<RecSplit<I, false, WS>>();
    } else if constexpr (WS == 0) {
        std::cerr << "The relativeWidth " << relativeWidth << " was not compiled into this binary." << std::endl;
    } else {
        dispatchWidth<RecSplit, I, WS - 1>();
    }
}

template<template<size_t, bool, size_t> class RecSplit, size_t I>
void dispatchLeafSize() {
    if constexpr (I <= 1) {
        std::cerr << "The leafSize " << leafSize << " was not compiled into this binary." << std::endl;
    } else if (I == leafSize) {
        dispatchWidth<RecSplit, I, std::min(morphishash::MAX_DIFF, I)>();
    } else {
        dispatchLeafSize<RecSplit, I - 1>();
    }
}

int main(int argc, const char *const *argv) {
    tlx::CmdlineParser cmd;
    cmd.add_bytes('n', "numObjects", numObjects, "Number of objects to construct with");
    cmd.add_bytes('q', "numQueries", numQueries, "Number of queries to measure");
    cmd.add_bytes('l', "leafSize", leafSize, "Leaf size to construct");
    cmd.add_bytes('w', "relativeWidth", relativeWidth,
                  "Bits less than leaf size for fixed width retrieval. Use leafSize for BuRR mode.");
    cmd.add_bytes('b', "bucketSize", bucketSize, "Bucket size to construct");
    cmd.add_bool('r', "rotate", rotate, "Apply rotation fitting");
    cmd.add_bool('2', "shockhash2", shockhash2, "ShockHash2");
    cmd.add_bool('f', "shockhash2flat", shockhash2flat, "ShockHash2 flat");
    cmd.add_bool('u', "useburr", useBurr, "ShockHash2 burr mode");
    cmd.add_bytes('t', "threads", threads, "Number of threads");

    if (!cmd.process(argc, argv)) {
        return 1;
    }

    if (shockhash2) {
        dispatchLeafSize<morphishash::MorphisHash2, morphishash::MAX_LEAF_SIZE2>();
    } else if (shockhash2flat) {
        dispatchLeafSize<morphishash::MorphisHash2Flat, morphishash::MAX_LEAF_SIZE2>();
    }
    //construct<shockhash::ShockHash2Flat<70, false, 4>>();
    return 0;
}
