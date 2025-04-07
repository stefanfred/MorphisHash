#include <chrono>
#include <iostream>
#include <bytehamster/util/XorShift64.h>
#include <tlx/cmdline_parser.hpp>
#include "BenchmarkData.h"


#include "MorphisHash.h"
#include "MorphisHashFlat.h"

#define DO_NOT_OPTIMIZE(value) asm volatile("" : : "r,m"(value) : "memory")

size_t numObjects = 1e7;
size_t numQueries = 1e6;
size_t bucketSize = 2000;
size_t threads = 8;


template<size_t l, size_t w, size_t o>
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
    morphishash::MorphisHashFlat<l, w, o> hashFunc(keys, bucketSize, threads);
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

    std::cout << "RESULT"
              << " l=" << l
               << " w=" << w
               << " o=" << o
              << " b=" << bucketSize
              << " N=" << numObjects
              << " threads=" << threads
              << " numQueries=" << numQueries
              << " queryTimeMilliseconds=" << queryDurationMs
              << " constructionTimeMilliseconds=" << constructionDurationMs
              << " bitsPerElement=" << (double) hashFunc.getBits() / numObjects
              << std::endl;
}

int main() {
    construct<80,2,3>();
}