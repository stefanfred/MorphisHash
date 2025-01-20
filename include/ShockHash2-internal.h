#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <variant>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <TinyBinaryCuckooHashTable.h>
#include <UnionFind.h>
#include <experimental/simd>
#include <type_traits>
#include "PairingFunction.h"
#include "CuckooUnionFind.h"

namespace morphishash {
    namespace stdx = std::experimental;

    template<typename T>
    using simd_t = stdx::native_simd<T>;
//using simd_t = stdx::fixed_size_simd<T, 1>; // To disable SIMD

    template<size_t leafSize>
    static constexpr uint64_t MASK_HALF = ((leafSize + 1) / 2 >= 64) ? ~0ul : (1ul << ((leafSize + 1) / 2)) - 1;

    template<size_t leafSize, bool isolatedVertexFilter = false>
    struct SeedCache {
        [[no_unique_address]]
        std::conditional_t<isolatedVertexFilter, __uint128_t, std::monostate> isolatedVertices;
        uint64_t seed;
        uint8_t hashes[leafSize];
    };

    template<size_t leafSize>
    inline void calculateIsolatedVertices(SeedCache<leafSize, true> &seedCache) {
        seedCache.isolatedVertices = __uint128_t(0);
        uint64_t hitCount[leafSize] = {0};
        for (size_t i = 0; i < leafSize; i++) {
            hitCount[seedCache.hashes[i]]++;
        }
        for (size_t i = 0; i < leafSize; i++) {
            if (hitCount[seedCache.hashes[i]] == 1) {
                seedCache.isolatedVertices |= __uint128_t(1) << i;
            }
        }

        /*__uint128_t hitCount1 = 0;
        __uint128_t hitCountMoreThan1 = 0;
        for (size_t i = 0; i < leafSize; i++) {
            // hitCount[seedCache.hashes[i]]++;
            __uint128_t pos = __uint128_t(1) << seedCache.hashes[i];
            hitCountMoreThan1 |= (hitCount1 & pos);
            hitCount1 |= pos;
        }
        __uint128_t isolated = 0;
        for (size_t i = 0; i < leafSize; i++) {
            //if (hitCount[seedCache.hashes[i]] == 1) {
            //    seedCache.isolatedVertices |= __uint128_t(1) << i;
            //}
            __uint128_t pos = __uint128_t(1) << seedCache.hashes[i];
            if ((((~hitCount1) | hitCountMoreThan1) & pos) == 0) {
                isolated |= __uint128_t(1) << i;
            }
        }
        seedCache.isolatedVertices = isolated;*/
    }

    class BasicSeedCandidateFinder {
    public:
        template<size_t leafSize, bool isolatedVertexFilter = false>
        class Finder {
        private:
            size_t currentSeed = 0;
        public:
            std::array<uint64_t, leafSize> keys = {0};

            explicit Finder(const std::vector<uint64_t> &inKeys) {
                for (size_t i = 0; i < leafSize; ++i) {
                    keys[i] = inKeys[i];
                }
            }

            inline SeedCache<leafSize, isolatedVertexFilter> next() {
                SeedCache<leafSize, isolatedVertexFilter> seedCache; // NOLINT(cppcoreguidelines-pro-type-member-init)
                while (true) {
                    uint64_t taken = 0;
                    for (size_t i = 0; i < leafSize; i++) {
                        uint64_t hash = ::util::remix(keys.at(i) + currentSeed);
                        seedCache.hashes[i] = ::util::fastrange64(hash, (leafSize + 1) / 2);
                        // std::cout << keys.at(i) << " " << uint64_t(seedCache.hashes[i]) << " " << currentSeed << std::endl;
                        taken |= 1ul << seedCache.hashes[i];
                    }
                    if (taken == MASK_HALF<leafSize>) {
                        // Found a new seed candidate
                        seedCache.seed = currentSeed;
                        if constexpr (isolatedVertexFilter) {
                            calculateIsolatedVertices(seedCache);
                        }
                        currentSeed++;
                        return seedCache;
                    }
                    currentSeed++;
                }
            }

            static std::string name() {
                return "Basic";
            }
        };

    public:
        static size_t hash(uint64_t key, uint64_t seed, size_t leafSize) {
            return ::util::fastrange64(::util::remix(key + seed), (leafSize + 1) / 2);
        }
    };

    class RotatingSeedCandidateFinder {
    public:
        template<size_t leafSize, bool isolatedVertexFilter = false>
        class Finder {
        private:
            uint64_t sizeSetA = 0;
            size_t currentSeed = -1;
            size_t currentRotation = (leafSize + 1) / 2;
            uint64_t takenA = 0;
            uint64_t takenB = 0;
            SeedCache<leafSize, isolatedVertexFilter> seedCache; // NOLINT(cppcoreguidelines-pro-type-member-init)
        public:
            std::array<uint64_t, leafSize> keys = {0};

            explicit Finder(const std::vector<uint64_t> &keysIn) {
                for (size_t i = 0; i < leafSize; i++) {
                    if ((keysIn[i] & 1) == 0) {
                        keys.at(sizeSetA) = keysIn[i];
                        sizeSetA++;
                    } else {
                        keys.at(leafSize - i + sizeSetA - 1) = keysIn[i];
                    }
                }
            }

            static constexpr uint64_t rotate(uint64_t val, uint32_t x) {
                constexpr uint64_t l = (leafSize + 1) / 2;
                return ((val << x) | (val >> (l - x))) & ((1ul << l) - 1);
            }

            inline SeedCache<leafSize, isolatedVertexFilter> next() {
                while (true) {
                    while (currentRotation < (leafSize + 1) / 2) {
                        if ((takenA | rotate(takenB, currentRotation)) == MASK_HALF<leafSize>) {
                            // Found a new seed candidate
                            SeedCache<leafSize, isolatedVertexFilter> rotated = seedCache;
                            rotated.seed = currentSeed * ((leafSize + 1) / 2) + currentRotation;
                            for (size_t i = sizeSetA; i < leafSize; i++) {
                                rotated.hashes[i] = (rotated.hashes[i] + currentRotation) % ((leafSize + 1) / 2);
                            }
                            if constexpr (isolatedVertexFilter) {
                                calculateIsolatedVertices(rotated);
                            }
                            currentRotation++;
                            return rotated;
                        }
                        currentRotation++;
                    }
                    currentRotation = 0;
                    currentSeed++;

                    takenA = 0;
                    takenB = 0;
                    for (size_t i = 0; i < sizeSetA; i++) {
                        uint64_t hash = ::util::remix(keys[i] + currentSeed);
                        seedCache.hashes[i] = ::util::fastrange64(hash, (leafSize + 1) / 2);
                        takenA |= 1ul << seedCache.hashes[i];
                    }
                    for (size_t i = sizeSetA; i < leafSize; i++) {
                        uint64_t hash = ::util::remix(keys[i] + currentSeed);
                        seedCache.hashes[i] = ::util::fastrange64(hash, (leafSize + 1) / 2);
                        takenB |= 1ul << seedCache.hashes[i];
                    }
                }
            }


            static std::string name() {
                return "Rotate";
            }
        };

    public:
        static size_t hash(uint64_t key, uint64_t seed, size_t leafSize) {
            size_t hashSeed = seed / ((leafSize + 1) / 2);
            size_t rotation = seed % ((leafSize + 1) / 2);
            size_t baseHash = ::util::fastrange64(::util::remix(key + hashSeed), (leafSize + 1) / 2);
            if ((key & 1) == 0) {
                return baseHash;
            } else {
                return (baseHash + rotation) % ((leafSize + 1) / 2);
            }
        }
    };

    template<size_t leafSize>
    class CandidateList {
    private:
        static constexpr size_t NUM_SENTINELS = simd_t<uint64_t>::size() + 1;
        std::vector<uint64_t> candidates;
    public:
        explicit CandidateList(size_t expectedNumSeeds) {
            candidates.reserve(expectedNumSeeds);
            for (size_t i = 0; i < NUM_SENTINELS; i++) {
                candidates.emplace_back(MASK_HALF<leafSize>);
            }
        }

        inline void add(size_t seed, uint64_t mask) {
            assert(seed == candidates.size() - NUM_SENTINELS);
            (void) seed;
            *(candidates.end() - NUM_SENTINELS) = mask;
            candidates.emplace_back(MASK_HALF<leafSize>);
        }

        struct IteratorType {
            const CandidateList &candidateList;
            size_t currentIdx = -1;
            size_t size = 0;
            uint64_t filterMask;
            bool isEnd;

            IteratorType(const CandidateList &candidateList, uint64_t filterMask, bool isEnd)
                    : candidateList(candidateList), filterMask(filterMask), isEnd(isEnd) {
                size = candidateList.candidates.size();
                if (!isEnd) {
                    operator++(); // Initialized with -1, go to first actual item
                }
            }

            inline bool operator!=(IteratorType rhs) const {
                return isEnd != rhs.isEnd;
            }

            inline std::pair<size_t, uint64_t> operator*() const {
                uint64_t mask = candidateList.candidates.at(currentIdx);
                return std::make_pair(currentIdx, mask);
            }

            inline IteratorType &operator++() {
                ++currentIdx;
                simd_t<uint64_t> filterMaskSimd = filterMask;
                simd_t<uint64_t> maskHalfSimd = MASK_HALF<leafSize>;
                while (true) {
                    simd_t<uint64_t> read(&candidateList.candidates[currentIdx], stdx::element_aligned);
                    simd_t<uint64_t>::mask_type comparisonResult = (read | filterMaskSimd) == maskHalfSimd;
                    if (stdx::any_of(comparisonResult)) {
                        currentIdx += stdx::find_first_set(comparisonResult);
                        break;
                    }
                    currentIdx += simd_t<uint64_t>::size();
                }
                if (currentIdx < size - NUM_SENTINELS) {
                    // Last is sentinel, return only if not sentinel
                    return *this;
                }
                isEnd = true;
                return *this;
            }
        };

        struct FilteredListType {
            const CandidateList &candidateList;
            uint64_t mask;

            FilteredListType(CandidateList &candidateList, uint64_t mask)
                    : candidateList(candidateList), mask(mask) {
            }

            inline IteratorType begin() const {
                return IteratorType(candidateList, mask, false);
            }

            inline IteratorType end() const {
                return IteratorType(candidateList, mask, true);
            }
        };

        FilteredListType filter(uint64_t mask) {
            return FilteredListType(*this, mask);
        }

        static std::string name() {
            return "List";
        }
    };

    template<size_t leafSize>
    class CandidateBuckets {
    private:
        static constexpr uint64_t BUCKET_MASK = 0b11111;
        std::array<std::vector<uint64_t>, BUCKET_MASK + 1> candidateMasks;
        std::array<std::vector<uint64_t>, BUCKET_MASK + 1> candidateSeeds;
    public:
        explicit CandidateBuckets(size_t expectedNumSeeds) {
            size_t toReserve = expectedNumSeeds / candidateSeeds.size();
            for (size_t i = 0; i < candidateSeeds.size(); i++) {
                candidateSeeds[i].reserve(toReserve);
                candidateMasks[i].reserve(toReserve);
                candidateMasks[i].push_back(MASK_HALF<leafSize>); // Sentinel
            }
        }

        inline void add(size_t seed, uint64_t mask) {
            candidateMasks[mask & BUCKET_MASK].back() = mask;
            candidateMasks[mask & BUCKET_MASK].emplace_back(MASK_HALF<leafSize>);
            candidateSeeds[mask & BUCKET_MASK].emplace_back(seed);
        }

        struct IteratorType {
            const CandidateBuckets &candidateTree;
            size_t currentBucket = 0;
            size_t currentIdx = -1;
            uint64_t filterMask;
            uint64_t filterMaskRestrictedToBucketMask;
            bool isEnd;

            inline void nextRelevantBucket() {
                while ((filterMaskRestrictedToBucketMask | currentBucket) < BUCKET_MASK) {
                    currentBucket++;
                }
                if ((filterMaskRestrictedToBucketMask | currentBucket) != BUCKET_MASK && currentBucket <= BUCKET_MASK) {
                    currentBucket++;
                }
            }

            IteratorType(const CandidateBuckets &candidateTree, uint64_t filterMask, bool isEnd)
                    : candidateTree(candidateTree), filterMask(filterMask),
                      filterMaskRestrictedToBucketMask(filterMask & BUCKET_MASK), isEnd(isEnd) {
                if (!isEnd) {
                    nextRelevantBucket();
                    operator++(); // Initialized with -1, go to first actual item
                }
            }

            inline bool operator!=(IteratorType rhs) const {
                return isEnd != rhs.isEnd;
            }

            inline std::pair<size_t, uint64_t> operator*() const {
                return std::make_pair(candidateTree.candidateSeeds[currentBucket][currentIdx],
                                      candidateTree.candidateMasks[currentBucket][currentIdx]);
            }

            inline IteratorType &operator++() {
                ++currentIdx;
                while (currentBucket <= BUCKET_MASK) {
                    const std::vector<uint64_t> &bucket = candidateTree.candidateMasks[currentBucket];
                    while (true) {
                        if ((bucket[currentIdx] | filterMask) == MASK_HALF<leafSize>) {
                            break;
                        }
                        ++currentIdx;
                    }
                    if (currentIdx < bucket.size() - 1) {
                        // Last is sentinel, return only if not sentinel
                        return *this;
                    }
                    currentBucket++;
                    nextRelevantBucket();
                    currentIdx = 0;
                }
                isEnd = true;
                return *this;
            }
        };

        struct FilteredListType {
            const CandidateBuckets &candidateTree;
            uint64_t mask;

            FilteredListType(CandidateBuckets &candidateTree, uint64_t mask)
                    : candidateTree(candidateTree), mask(mask) {
            }

            inline IteratorType begin() const {
                return IteratorType(candidateTree, mask, false);
            }

            inline IteratorType end() const {
                return IteratorType(candidateTree, mask, true);
            }
        };

        FilteredListType filter(uint64_t mask) {
            return FilteredListType(*this, mask);
        }

        static std::string name() {
            return "Buckets" + std::to_string(BUCKET_MASK + 1);
        }
    };

    class QuadSplitCandidateFinder {
    private:

    public:
        template<template<size_t> typename CandidateList, size_t leafSize, bool isolatedVertexFilter = false>
        class Finder {
        public:
            static constexpr double E_HALF = 1.359140914;
            uint64_t sizeSetA = 0;
            size_t currentSeed = -1;
            CandidateList<leafSize> candidatesA;
            CandidateList<leafSize> candidatesB;
            std::vector<SeedCache<leafSize, isolatedVertexFilter>> extractedCandidates;
        public:
            std::array<uint64_t, leafSize> keys = {0};

            explicit Finder(const std::vector<uint64_t> &keysIn)
                    : candidatesA(std::pow(E_HALF, keysIn.size() / 2)),
                      candidatesB(std::pow(E_HALF, keysIn.size() / 2)) {
                for (size_t i = 0; i < leafSize; i++) {
                    if ((keysIn[i] & 1) == 0) {
                        keys.at(sizeSetA) = keysIn[i];
                        sizeSetA++;
                    } else {
                        keys.at(leafSize - i + sizeSetA - 1) = keysIn[i];
                    }
                }
            }


            inline SeedCache<leafSize, isolatedVertexFilter> next() {
                while (extractedCandidates.empty()) {
                    prepareNextSeed();
                }
                SeedCache<leafSize, isolatedVertexFilter> seedCache = extractedCandidates.back(); // Smallest seed is in the back
                extractedCandidates.pop_back();
                return seedCache;
            }

            inline SeedCache<leafSize, isolatedVertexFilter> makeCache(size_t seedA, size_t seedB) {
                SeedCache<leafSize, isolatedVertexFilter> cache = {};
                cache.seed = pairElegant(seedA, seedB);
                for (size_t i = 0; i < sizeSetA; i++) {
                    uint64_t hash = ::util::remix(keys[i] + seedA);
                    cache.hashes[i] = ::util::fastrange64(hash, (leafSize + 1) / 2);
                }
                for (size_t i = sizeSetA; i < leafSize; i++) {
                    uint64_t hash = ::util::remix(keys[i] + seedB);
                    cache.hashes[i] = ::util::fastrange64(hash, (leafSize + 1) / 2);
                }
                if constexpr (isolatedVertexFilter) {
                    calculateIsolatedVertices(cache);
                }
                return cache;
            }

            void prepareNextSeed() {
                currentSeed++;

                uint64_t takenA = 0;
                uint64_t takenB = 0;
                for (size_t i = 0; i < sizeSetA; i++) {
                    uint64_t hash = ::util::remix(keys[i] + currentSeed);
                    takenA |= 1ul << ::util::fastrange64(hash, (leafSize + 1) / 2);
                }
                for (size_t i = sizeSetA; i < leafSize; i++) {
                    uint64_t hash = ::util::remix(keys[i] + currentSeed);
                    takenB |= 1ul << ::util::fastrange64(hash, (leafSize + 1) / 2);
                }

                for (const auto [candidateSeed, candidateMask]: candidatesA.filter(takenB)) {
                    if (isFloatAccurateElegant(candidateSeed, currentSeed)) {
                        extractedCandidates.push_back(makeCache(candidateSeed, currentSeed));
                    } else {
                        std::cout << "Skipped seed because of floating point inaccuracy" << std::endl;
                    }
                }
                candidatesA.add(currentSeed,
                                takenA); // add after iterating, so we don't test the same seed combination twice
                candidatesB.add(currentSeed, takenB);
                for (const auto [candidateSeed, candidateMask]: candidatesB.filter(takenA)) {
                    if (isFloatAccurateElegant(currentSeed, candidateSeed)) {
                        extractedCandidates.push_back(makeCache(currentSeed, candidateSeed));
                    } else {
                        std::cout << "Skipped seed because of floating point inaccuracy" << std::endl;
                    }
                }
                if (extractedCandidates.size() > 1) {
                    // Smallest seed is in the back
                    std::sort(extractedCandidates.begin(), extractedCandidates.end(),
                              [](const SeedCache<leafSize, isolatedVertexFilter> &a,
                                 const SeedCache<leafSize, isolatedVertexFilter> &b) {
                                  return a.seed > b.seed;
                              });
                }
            }

            static std::string name() {
                return std::string("QuadSplit")
                       + CandidateList<leafSize>::name();
            }
        };

    public:
        static size_t hash(uint64_t key, uint64_t seed, size_t leafSize) {
            auto [seedA, seedB] = unpairElegant(seed);
            uint64_t seedToUse = ((key & 1) == 0) ? seedA : seedB;

            return ::util::fastrange64(::util::remix(key + seedToUse), (leafSize + 1) / 2);
        }
    };

    template<size_t leafSize, bool isolatedVertexFilter>
    using QuadSplitCandidateFinderList = QuadSplitCandidateFinder::Finder<CandidateList, leafSize, isolatedVertexFilter>;

    template<size_t leafSize, bool isolatedVertexFilter>
    using QuadSplitCandidateFinderBuckets = QuadSplitCandidateFinder::Finder<CandidateBuckets, leafSize, isolatedVertexFilter>;


    static constexpr size_t SockHash2SeedFinderLeafSizeThreshold = 10;


    static inline size_t queryCandidate(size_t seed, uint64_t key, size_t leafSize) {
        if (leafSize < SockHash2SeedFinderLeafSizeThreshold) {
            return BasicSeedCandidateFinder::hash(key, seed, leafSize);
        } else {
            return QuadSplitCandidateFinder::hash(key, seed, leafSize);
        }
    }

    static inline size_t queryHash(size_t seed, uint64_t key, uint64_t retrieved, size_t leafSize) {
        auto [seed1, seed2] = unpairTriangular(seed);

        size_t usedSeed;
        size_t result;
        if (retrieved == 0) {
            usedSeed = seed1;
            result = leafSize / 2;
        } else {
            usedSeed = seed2;
            result = 0;
        }
        result += queryCandidate(usedSeed, key, leafSize);

        assert(result <= leafSize);
        return result;
    }

    static void verify(size_t seed, const std::vector<uint64_t> &keys, size_t leafSize,
                       std::vector<std::pair<uint64_t, uint8_t>> &retrieval) {
        std::vector<bool> taken(leafSize, false);
        for (uint64_t key: keys) {
            size_t retrieved = ~0u;
            for (auto &retr: retrieval) {
                if (retr.first == key) {
                    retrieved = retr.second;
                }
            }
            if (retrieved == ~0u) {
                throw std::logic_error("Not in retrieval");
            }
            size_t hashValue = queryHash(seed, key, retrieved, leafSize);
            // size_t hashValue = hash(seed, key, 0);
            if (taken[hashValue]) {
                throw std::logic_error("Collision");
            }
            taken[hashValue] = true;
        }
    }

    static inline void constructRetrieval(const std::vector<uint64_t> &keys, size_t seed,
                                          std::vector<std::pair<uint64_t, uint8_t>> &retrieval, size_t leafSize) {

        auto [seed1, seed2] = unpairTriangular(seed >> 12);
        morphishash::TinyBinaryCuckooHashTable table(leafSize);
        for (uint64_t key: keys) {
            table.prepare(morphishash::HashedKey(key));
        }
        table.clearPlacement();
        for (size_t k = 0; k < leafSize; k++) {
            morphishash::TinyBinaryCuckooHashTable::CandidateCells candidateCells;
            candidateCells.cell1 = queryCandidate(seed1, table.heap[k].hash.mhc, leafSize) + leafSize / 2;
            candidateCells.cell2 = queryCandidate(seed2, table.heap[k].hash.mhc, leafSize);
            if (!table.insert(&table.heap[k], candidateCells)) {
                throw std::logic_error("Should be possible to construct");
            }
        }
        for (size_t k = 0; k < leafSize; k++) {
            size_t candidate2 = queryCandidate(seed2, table.heap[k].hash.mhc, leafSize);
            if (table.cells[candidate2] == &table.heap[k]) {
                retrieval.emplace_back(table.heap[k].hash.mhc, 1);
            } else {
                retrieval.emplace_back(table.heap[k].hash.mhc, 0);
            }
        }
    }

    static int parity(uint64_t val) {
        return __builtin_parityll(val);
    }

    static int parity(__uint128_t val) {
        return parity(uint64_t(val >> 64) ^ uint64_t(val));
    }

    static uint64_t inline sh2remix64(uint64_t z) {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
    }

    static __uint128_t inline sh2remix128(uint64_t z) {
        return (__uint128_t(sh2remix64(~z)) << 64) | sh2remix64(z);
    }

    template<bool bitMode128>
    static std::conditional<bitMode128, __uint128_t, uint64_t>::type inline sh2remix(uint64_t z) {
        if constexpr (bitMode128) {
            return sh2remix128(z);
        } else {
            return sh2remix64(z);
        }
    }


    static constexpr size_t MAX_LEAF_SIZE2 = 106;
    static constexpr size_t MAX_RETRIEVAL_WIDTH = MAX_LEAF_SIZE2;
    static constexpr size_t MAX_DIFF = 7;

/**
 * ShockHash2 base case.
 * Note that while this can be used with uneven leaf sizes, it achieves suboptimal space and time.
 */
    template<size_t leafSize, template<size_t, bool> typename SeedCandidateFinder, bool isolatedVertexFilter = false, bool useBurr = false, size_t matrix_width = leafSize>
    class BijectionsShockHash2 {

    private:
        using CandidateFinder = SeedCandidateFinder<leafSize, isolatedVertexFilter>;

        static constexpr bool bitMode128 = matrix_width >= size_t(64);
        using matrixRow = std::conditional<bitMode128, __uint128_t, uint64_t>::type;
        static constexpr matrixRow row_mask = (matrixRow(1) << matrix_width) - 1;


        template<typename t>
        static inline t check_bit(t number, size_t n) {
            return number & (t(1) << n);
        }

        template<typename t>
        static inline t set_bit(t number, size_t n, bool x) {
            return (number & ~(t(1) << n)) | (t(x) << n);
        }


    public:

        static inline std::pair<uint64_t, __uint128_t> findSeed(const std::vector<uint64_t> &inKeys) {
            assert(inKeys.size() == leafSize);
            if constexpr (not useBurr and matrix_width == 0) {
                size_t seed = 0;
                std::bitset<leafSize> taken;
                while (true) {
                    taken.reset();
                    bool works = true;
                    for (uint64_t k: inKeys) {
                        uint64_t pos = sh2remix64(k ^ seed) % leafSize;
                        if (taken[pos]) {
                            works = false;
                            seed++;
                            break;
                        } else {
                            taken[pos] = true;
                        }
                    }
                    if (works) {
                        return {seed, 0};
                    }
                }
            }
            CandidateFinder seedCandidateFinder(inKeys);
            std::array<uint64_t, leafSize> keys = seedCandidateFinder.keys;
            std::vector<SeedCache<leafSize, isolatedVertexFilter>> seedsCandidates;
            seedsCandidates.push_back(seedCandidateFinder.next());
            CuckooUnionFind unionFind(leafSize);
            while (true) {
                const SeedCache<leafSize, isolatedVertexFilter> newCandidate = seedCandidateFinder.next();
                SeedCache<leafSize, isolatedVertexFilter> newCandidateShifted = newCandidate;
                for (size_t i = 0; i < leafSize; i++) {
                    newCandidateShifted.hashes[i] += leafSize / 2;
                }
                for (const SeedCache<leafSize, isolatedVertexFilter> &other: seedsCandidates) {
                    if constexpr (isolatedVertexFilter) {
                        if ((other.isolatedVertices & newCandidate.isolatedVertices) != 0ul) {
                            continue;
                        }
                    }
                    size_t li = 0;
                    unionFind.clear();
                    for (; li < leafSize; li++) {
                        size_t end1 = newCandidateShifted.hashes[li];
                        size_t end2 = other.hashes[li];
                        /*std::cout << keys[li] << " " << end1 << " " << end2 << " "
                                  << (QuadSplitCandidateFinder::hash(keys[li], newCandidateShifted.seed, leafSize) +leafSize / 2) <<" "
                                  << QuadSplitCandidateFinder::hash(keys[li], other.seed, leafSize)
                                  << std::endl;*/
                        if (!unionFind.unionIsStillPseudoforest(end1, end2)) {
                            break;
                        }
                    }
                    if (li != leafSize) {
                        continue;
                    }

                    // Found working seed!
                    uint64_t seed1 = newCandidate.seed;
                    uint64_t seed2 = other.seed;
                    uint64_t fullSeed = pairTriangular(seed1, seed2);

                    if (!isFloatAccurateTriangular(seed1, seed2)) {
                        std::cout << "Skipped seed because of floating point inaccuracy" << std::endl;
                        continue;
                    }

                    if constexpr (useBurr) {
                        return {fullSeed, 0};
                    } else {
                        std::array<matrixRow, leafSize> matrix{};
                        std::bitset<leafSize> sol{};
                        sol.set();

                        // insert keys
                        for (size_t i = 0; i < leafSize; i++) {
                            //matrixRow hash = keys[i] & row_mask;
                            matrixRow hash = sh2remix<bitMode128>(keys[i] ^ fullSeed) & row_mask;
                            auto addCand = [&](size_t end, bool orientation) {
                                matrix[end] ^= hash;
                                if (!orientation)
                                    sol[end].flip();
                            };
                            addCand(newCandidateShifted.hashes[i], false);
                            addCand(other.hashes[i], true);
                        }


                        // gauss
                        std::bitset<leafSize> usedrow{};
                        for (size_t coloumn = 0; coloumn < matrix_width; ++coloumn) {
                            // find usable row
                            size_t pivotrowindex = 0;
                            for (; pivotrowindex < leafSize; ++pivotrowindex) {
                                if (not usedrow[pivotrowindex] and check_bit(matrix[pivotrowindex], coloumn)) {
                                    break;
                                }
                            }
                            if (pivotrowindex == leafSize) {
                                // zero coloumn
                                continue;
                            }
                            usedrow[pivotrowindex] = true;
                            matrixRow pivotrow = matrix[pivotrowindex];
                            bool pivotsol = sol[pivotrowindex];

                            // add to all other
                            for (size_t rowindex = 0; rowindex < leafSize; ++rowindex) {
                                if (rowindex != pivotrowindex && check_bit(matrix[rowindex], coloumn)) {
                                    matrix[rowindex] ^= pivotrow;
                                    sol[rowindex] = sol[rowindex] != pivotsol;
                                }
                            }
                        }

                        // check solvable
                        matrixRow res = 0;
                        bool fail = false;
                        for (size_t rowindex = 0; rowindex < leafSize and not fail; ++rowindex) {
                            matrixRow row = matrix[rowindex];
                            bool s = sol[rowindex];
                            if (row == 0) {
                                if (s) {
                                    // contradiction
                                    fail = true;
                                    break;
                                } else {
                                    // redundant
                                    continue;
                                }
                            }
                            int pivotIndex = std::countr_zero(row);
                            res = set_bit(res, pivotIndex, s);
                        }

                        if (fail) {
                            continue;
                        }


                        /*std::bitset<leafSize> occ;
                        for (size_t i = 0; i < leafSize; ++i) {
                            size_t pos = queryHash(fullSeed, keys[i],
                                                   parity(res & row_mask & sh2remix<bitMode128>(keys[i] ^ fullSeed)),
                                                   leafSize);
                            if (occ[pos]) {
                                std::cout << " FAIL " << std::endl;
                                exit(123);
                            }
                        }*/
                        //std::cout << " VALID " << std::endl;
/*
                        for (size_t i = 0; i < leafSize; ++i) {
                            std::cout << keys[i] << " " << queryHash(fullSeed, keys[i], 0, leafSize) << " "
                                      << queryHash(fullSeed, keys[i], 1, leafSize) << " "
                                      << queryHash(fullSeed, keys[i],
                                                   parity(res & row_mask & sh2remix(keys[i] ^ fullSeed)), leafSize)
                                      << std::endl;*/
                            //std::cout << keys[i] <<" "<< (parity(keys[i] & uint64_t(res & row_mask)) ? std::to_string(newCandidateShifted.hashes[i]) : std::to_string(other.hashes[i]))<<" "<<std::to_string(newCandidateShifted.hashes[i])<<" "<<std::to_string(other.hashes[i])<< " "<<seed1<<" "<<seed2<<  std::endl;
                            //std::cout<<(CandidateFinder::Finder::hash(keys[i], seed1))<<" "<<CandidateFinder::hash(keys[i], seed2)<<std::endl;
                            /*std::cout << keys[i] << " " << " " << fullSeed << " " << uint64_t(res & row_mask) << " "
                                      << leafSize
                                      << std::endl;
                        }*/
                        //std::cout << fullSeed << " " << uint64_t(res & row_mask) << std::endl;

                        //return fullSeed;
                        return {fullSeed, res};
                    }
                }
                seedsCandidates.push_back(newCandidate);
            }
            return {0, 0};
        }


        /*inline double calculateBijection(const std::vector<uint64_t> &keys) {
            assert(matrix_width == 0);
            size_t seed = findSeed(keys);
            // Begin: Validity check
#ifndef NDEBUG
            std::vector<std::pair<uint64_t, uint8_t>> retrieval;
            constructRetrieval(keys, seed, retrieval, leafSize);
            verify(seed, keys, leafSize, retrieval);
#endif
            // End: Validity check
            return seed;
        }*/

        static std::string name() {
            return std::string("ShockHash2")
                   + (isolatedVertexFilter ? "Filter" : "")
                   + CandidateFinder::name();
        }
    };
} // namespace shockhash
