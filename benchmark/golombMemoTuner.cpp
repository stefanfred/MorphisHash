#include <vector>
#include <iostream>
#include "XorShift64.h"
#include "MorphisHash2-internal.h"
#include "MorphisHash2.h"

template<template<size_t, size_t> class T, size_t leafSize, size_t widthDiff>
void dispatchLeafSize() {
    if constexpr (leafSize > 1) {
        dispatchLeafSize<T, leafSize - 1, widthDiff>();
    }
    if constexpr (leafSize == 1) {
        std::cout << " 0, 0," << std::flush;
        return;
    }

    constexpr uint64_t width = leafSize >= widthDiff ? leafSize - widthDiff : 0;

    std::vector<uint64_t> leaf(leafSize);
    util::XorShift64 prng;
    std::vector<__uint128_t> seeds;
    uint64_t bitsexp = uint64_t(morphishash::bij_memoMorphis[widthDiff][leafSize]) - width;
    uint64_t bitsexpmax = uint64_t(morphishash::bij_memoMorphis[morphishash::MAX_DIFF][morphishash::MAX_LEAF_SIZE2]) - uint64_t(morphishash::MAX_LEAF_SIZE2) + uint64_t(morphishash::MAX_DIFF);
    uint64_t iterations = std::min(uint64_t(50) << (bitsexpmax - bitsexp), uint64_t(5000));
    seeds.reserve(iterations);
    for (size_t i = 0; i < iterations; i++) {
        for (size_t k = 0; k < leafSize; k++) {
            leaf[k] = prng();
        }
        seeds.push_back((T<leafSize, width>::findSeed(leaf)).first);
    }

    size_t spaceBest = std::numeric_limits<size_t>::max();
    size_t lowerBest = 0;
    for (size_t lower = 0; lower < 128; lower++) {
        size_t spaceBits = 0;
        for (__uint128_t seed: seeds) {
            spaceBits += (seed >> lower) + 1 + lower;
        }
        if (spaceBits < spaceBest) {
            spaceBest = spaceBits;
            lowerBest = lower;
        }
    }

    lowerBest += width;

    std::cout << (lowerBest < 10 ? " " : "") << lowerBest << ", " << std::flush;
}

template<template<size_t, size_t> class T, size_t widthDiff>
void dispatchWidth() {
    if constexpr (widthDiff <= morphishash::MAX_DIFF) {
        std::cout << "{" << std::flush;
        dispatchLeafSize<T, morphishash::MAX_LEAF_SIZE2, widthDiff>();
        std::cout << "}, " << std::flush;
        dispatchWidth<T, widthDiff + 1>();
    }
}


template<size_t leafSize, size_t width>
using ShockHash2 = std::conditional_t<(leafSize >= 10),
        morphishash::BijectionsShockHash2<leafSize, morphishash::QuadSplitCandidateFinderBuckets, true, false, width>,
        morphishash::BijectionsShockHash2<leafSize, morphishash::BasicSeedCandidateFinder::Finder, true, false, width>>;

int main() {
    dispatchWidth<ShockHash2, 0>();
}
