#pragma once

#include <vector>
#include <SimpleRibbon.h>
#include <EliasFano.h>
#include <MurmurHash64.h>
#include <tlx/math/integer_log2.hpp>
#include <MorphisHash2-internal.h>
#include <MorphisHash2.h>
#include <sdsl/int_vector.hpp>
#include "MorphisHash2FlatBase.h"

namespace morphishash {

    template<size_t THRESHOLD_RANGE>
    constexpr std::array<uint32_t, THRESHOLD_RANGE> _fill_mapping() {
        const uint64_t ONE_THIRD = std::numeric_limits<uint32_t>::max() / 3;
        std::array<uint32_t, THRESHOLD_RANGE> array;
        if (THRESHOLD_RANGE == 1) {
            return array;
        } else if (THRESHOLD_RANGE == 2) {
            array.at(0) = 0;
            array.at(1) = std::numeric_limits<uint32_t>::max();
            return array;
        }
        array.at(0) = 0; // Last resort
        array.at(1) = ONE_THIRD; // Safeguard, so much bumping should never happen in practice
        size_t interpolationRange = THRESHOLD_RANGE - 3;
        for (size_t i = 0; i < interpolationRange; i++) {
            array.at(2 + i) = 2 * ONE_THIRD + ONE_THIRD * i / interpolationRange;
        }
        array.at(THRESHOLD_RANGE - 1) = std::numeric_limits<uint32_t>::max(); // Keep all
        return array;
    }

    template<size_t k, bool useBurr, size_t RETRIEVAL_DIFF>
    class MorphisHash2Flat {
        using BaseCase = BijectionsShockHash2<k, morphishash::QuadSplitCandidateFinderList, true, useBurr,
                k - RETRIEVAL_DIFF>;
        static constexpr double OVERLOAD_FACTOR = 0.9;
        static constexpr size_t THRESHOLD_BITS = tlx::integer_log2_floor(k) - 1;
        static constexpr size_t THRESHOLD_RANGE = 1ul << THRESHOLD_BITS;
        static constexpr std::array<uint32_t, THRESHOLD_RANGE> THRESHOLD_MAPPING = _fill_mapping<THRESHOLD_RANGE>();
        // static constexpr size_t SEED_BITS = std::ceil(0.442 * k - 0.2 + double(RETRIEVAL_DIFF)*0.65);
        static constexpr size_t SEED_BITS =
                bij_memoMorphis[RETRIEVAL_DIFF][k] + 4 - (RETRIEVAL_DIFF <= k ? k - RETRIEVAL_DIFF : 0);
        static constexpr size_t MAX_SEED = 1ul << SEED_BITS;
        static constexpr size_t SEED_FALLBACK_INDICATOR = 0;
        sdsl::int_vector<0> thresholdsAndSeeds;
        std::map<size_t, size_t> seedsFallback;
        std::vector<size_t> layerBases;

        RiceBitVector<> solutions;

        MorphisHash2<k, useBurr, RETRIEVAL_DIFF> fallbackPhf;
        size_t N;
        size_t nbuckets;
        pasta::BitVector freePositionsBv;
        pasta::FlatRankSelect <pasta::OptimizedFor::ONE_QUERIES> *freePositionsRankSelect = nullptr;
        using Ribbon = SimpleRibbon<1, (k > 24) ? 128 : 64>;
        Ribbon ribbon;
        size_t layers = 2;
    public:
        explicit MorphisHash2Flat(const std::vector<std::string> &keys, size_t ignore, size_t ignore2)
                : MorphisHash2Flat(keys) {
            (void) ignore;
            (void) ignore2;
        }

        explicit MorphisHash2Flat(const std::vector<std::string> &keys) {
            N = keys.size();
            nbuckets = N / k;
            size_t keysInEndBucket = N - nbuckets * k;
            //std::cout<<N<<" "<<nbuckets<<" "<<keysInEndBucket<<std::endl;
            size_t bucketsThisLayer = std::max(1ul, (size_t) std::ceil(OVERLOAD_FACTOR * nbuckets));
            std::vector<size_t> freePositions;
            std::vector<KeyInfo> hashes;
            hashes.reserve(keys.size());
            for (const std::string &key: keys) {
                uint64_t mhc = ::util::MurmurHash64(key);
                uint32_t bucket = ::util::fastrange32(mhc & 0xffffffff, bucketsThisLayer);
                uint32_t threshold = mhc >> 32;
                hashes.emplace_back(mhc, bucket, threshold);
            }
            std::vector<KeyInfo> allHashes = hashes;
            layerBases.push_back(0);
            thresholdsAndSeeds.bit_resize((THRESHOLD_BITS + SEED_BITS) * nbuckets);
            for (size_t layer = 0; layer < 2; layer++) {
                size_t layerBase = layerBases.back();
                if (layer != 0) {
                    bucketsThisLayer = OVERLOAD_FACTOR * (hashes.size() / k);
                    bucketsThisLayer = std::min(bucketsThisLayer, nbuckets - layerBase);
                    if (bucketsThisLayer == 0) {
                        layers = 1;
                        break;
                    }
                    // Rehash
                    for (auto &hash: hashes) {
                        hash.mhc = ::util::remix(hash.mhc);
                        hash.bucket = ::util::fastrange32(hash.mhc & 0xffffffff, bucketsThisLayer);
                        hash.threshold = hash.mhc >> 32;
                    }
                }
                layerBases.push_back(layerBase + bucketsThisLayer);
                sort_keyInfo(hashes);
                std::vector<KeyInfo> bumpedKeys;
                size_t bucketStart = 0;
                size_t previousBucket = 0;
                for (size_t i = 0; i < hashes.size(); i++) {
                    size_t bucket = hashes.at(i).bucket;
                    while (bucket != previousBucket) {
                        flushBucket(layerBase, bucketStart, i, previousBucket, hashes, bumpedKeys, freePositions);
                        previousBucket++;
                        bucketStart = i;
                    }
                }
                // Last bucket
                while (previousBucket < bucketsThisLayer) {
                    flushBucket(layerBase, bucketStart, hashes.size(), previousBucket, hashes, bumpedKeys,
                                freePositions);
                    previousBucket++;
                    bucketStart = hashes.size();
                }
                hashes = bumpedKeys;
            }

            std::vector<std::string> fallbackPhfContent;
            for (auto &hash: hashes) {
                fallbackPhfContent.push_back(std::to_string(hash.mhc));
            }
            //std::cout<<fallbackPhfContent.size()<<std::endl;
            fallbackPhf = MorphisHash2<k, useBurr, RETRIEVAL_DIFF>(fallbackPhfContent, 2000, 1);
            size_t additionalFreePositions = hashes.size() - freePositions.size();
            size_t nbucketsHandled = layerBases.back();
            {
                size_t i = 0;
                for (; i < additionalFreePositions - keysInEndBucket; i++) {
                    freePositions.push_back(nbucketsHandled + i / k);
                }
                for (; i < additionalFreePositions; i++) {
                    freePositions.push_back(nbuckets + i);
                }
            }
            freePositionsBv.resize(freePositions.size() + freePositions.back() + 1, false);
            for (size_t i = 0; i < freePositions.size(); i++) {
                freePositionsBv[i + freePositions.at(i)] = true;
            }
            freePositionsRankSelect = new pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES>(freePositionsBv);

            // Construct ShockHash within buckets
            std::vector<std::vector<uint64_t>> bucketContents(nbuckets);
            std::vector<uint64_t> lastBucket;
            for (KeyInfo key: allHashes) {
                size_t bucket = evaluateKPerfect(key.mhc).first;
                if (bucket >= nbuckets) {
                    lastBucket.push_back(key.mhc);
                } else {
                    bucketContents.at(bucket).push_back(key.mhc);
                }
            }
            std::vector<std::pair<uint64_t, uint8_t>> ribbonInput;
            RiceBitVector<>::Builder solBuilder;
            for (size_t i = 0; i < nbuckets; i++) {
                //std::cout<<bucketContents.at(i).size()<<std::endl;
                std::pair<uint64_t, __uint128_t> seed = BaseCase::findSeed(bucketContents.at(i));
                if (useBurr)
                    constructRetrieval(bucketContents.at(i), seed.first, ribbonInput, k);
                else
                    solBuilder.appendFixed128(seed.second, k - RETRIEVAL_DIFF);
                if (seed.first >= MAX_SEED || seed.first == SEED_FALLBACK_INDICATOR) {
                    seedsFallback.insert(std::make_pair(i, seed.first));
                    seed.first = SEED_FALLBACK_INDICATOR;
                }
                setSeed(i, seed.first);
            }
            // Construct last bucket
            if (useBurr)
                ribbon = Ribbon(ribbonInput);
            else
                solutions = solBuilder.build();
        }

        inline void setThreshold(size_t bucket, size_t value) {
            thresholdsAndSeeds.set_int(bucket * (THRESHOLD_BITS + SEED_BITS), value, THRESHOLD_BITS);
        }

        inline void setSeed(size_t bucket, size_t value) {
            thresholdsAndSeeds.set_int(bucket * (THRESHOLD_BITS + SEED_BITS) + THRESHOLD_BITS, value, SEED_BITS);
        }

        inline std::pair<size_t, size_t> getThresholdAndSeed(size_t bucket) {
            uint64_t thresholdAndSeed = thresholdsAndSeeds.get_int(
                    bucket * (THRESHOLD_BITS + SEED_BITS), SEED_BITS + THRESHOLD_BITS);
            size_t seed = thresholdAndSeed >> THRESHOLD_BITS;
            size_t threshold = thresholdAndSeed & (THRESHOLD_RANGE - 1);
            return std::make_pair(threshold, seed);
        }

        uint32_t compact_threshold(uint32_t threshold) {
            // Binary search would be better here, but this doesn't matter much for performance anyway
            for (size_t i = 0; i < THRESHOLD_RANGE; i++) {
                if (threshold <= THRESHOLD_MAPPING[i]) {
                    return i;
                }
            }
            return THRESHOLD_MAPPING.back();
        }

        void flushBucket(size_t layerBase, size_t bucketStart, size_t i, size_t bucketIdx,
                         std::vector<KeyInfo> &hashes, std::vector<KeyInfo> &bumpedKeys,
                         std::vector<size_t> &freePositions) {
            size_t bucketSize = i - bucketStart;
            if (bucketSize <= k) {
                size_t threshold = THRESHOLD_RANGE - 1;
                setThreshold(layerBase + bucketIdx, threshold);
                for (size_t b = bucketSize; b < k; b++) {
                    freePositions.push_back(layerBase + bucketIdx);
                }
            } else {
                size_t lastThreshold = compact_threshold(hashes.at(bucketStart + k - 1).threshold);
                size_t firstBumpedThreshold = compact_threshold(hashes.at(bucketStart + k).threshold);
                size_t threshold = lastThreshold;
                if (firstBumpedThreshold == lastThreshold) {
                    // Needs to bump more
                    threshold--;
                }
                setThreshold(layerBase + bucketIdx, threshold);
                uint32_t uncompressedThreshold = THRESHOLD_MAPPING[threshold];
                for (size_t l = 0; l < bucketSize; l++) {
                    if (hashes.at(bucketStart + l).threshold > uncompressedThreshold) {
                        bumpedKeys.push_back(hashes.at(bucketStart + l));
                        if (l < k) {
                            freePositions.push_back(layerBase + bucketIdx);
                        }
                    }
                }
            }
        }

        /** Estimate for the space usage of this structure, in bits */
        [[nodiscard]] size_t getBits() {
            return 8 * sizeof(*this)
                   + fallbackPhf.getBits()
                   + (freePositionsBv.size() + 8 * freePositionsRankSelect->space_usage())
                   + (useBurr ? 8 * ribbon.sizeBytes() : solutions.getBits())
                   + thresholdsAndSeeds.bit_size()
                   + 64 * seedsFallback.size();
        }

        void printBits() {
            std::cout << "Thresholds: " << 1.0f * THRESHOLD_BITS / k << std::endl;
            std::cout << "Fallback PHF keys: " << freePositionsBv.size() - N / k << std::endl;
            std::cout << "PHF: " << 1.0f * fallbackPhf.getBits() / N << std::endl;
            std::cout << "Fano: " << 1.0f * (freePositionsBv.size() + 8 * freePositionsRankSelect->space_usage()) / N
                      << std::endl;
            std::cout << "Base case seeds: " << 1.0f * SEED_BITS / k << std::endl;
            std::cout << "Base case seeds overflow: " << 1.0f * seedsFallback.size() * 64 / N << std::endl;
            std::cout << "Retrieval: " << (useBurr ? 8 * ribbon.sizeBytes() : solutions.getBits()) / N << std::endl;
        }

        size_t operator()(const std::string &key) {
            return operator()(::util::MurmurHash64(key));
        }

        size_t operator()(const hash128_t &hash) {
            return operator()(hash.second);
        }

        __attribute_noinline__ size_t operator()(uint64_t hash) {
            auto [bucket, seed] = evaluateKPerfect(hash);
            if (bucket >= nbuckets) {
                return bucket; // N that are not multiples of k
            }
            if (seed == SEED_FALLBACK_INDICATOR) {
                seed = seedsFallback.at(bucket);
            }
            size_t baseCase;
            if (useBurr) {
                size_t retrieved = ribbon.retrieve(hash);
                baseCase = queryHash(seed, hash, retrieved, k);
            } else {
                size_t width = k > RETRIEVAL_DIFF ? (k - RETRIEVAL_DIFF) : 0;
                if (width == 0) {
                    baseCase = sh2remix64(hash ^ seed) % k;
                } else {
                    __uint128_t remixed = sh2remix128(hash ^ seed);
                    RiceBitVector<>::Reader r = solutions.reader();
                    r.toFixedPos(width, bucket);
                    __uint128_t sol = r.readFixed128(width);
                    __uint128_t row_mask = (__uint128_t(1) << (width)) - 1;
                    uint64_t retrieved = parity(sol & row_mask & remixed);
                    baseCase = queryHash(seed, hash, retrieved, k);
                }
            }
            return bucket * k + baseCase;
        }

        /** Returns bucket and seed */
        inline std::pair<size_t, size_t> evaluateKPerfect(uint64_t mhc) {
            for (size_t layer = 0; layer < layers; layer++) {
                if (layer != 0) {
                    mhc = ::util::remix(mhc);
                }
                size_t base = layerBases.at(layer);
                size_t layerSize = layerBases.at(layer + 1) - base;
                uint32_t bucket = ::util::fastrange32(mhc & 0xffffffff, layerSize);
                uint32_t threshold = mhc >> 32;
                auto [storedThreshold, storedSeed] = getThresholdAndSeed(base + bucket);
                if (threshold <= THRESHOLD_MAPPING[storedThreshold]) {
                    return std::make_pair(base + bucket, storedSeed);
                }
            }
            size_t phf = fallbackPhf(std::to_string(mhc));
            size_t bucket = freePositionsRankSelect->select1(phf + 1) - phf;
            if (bucket >= nbuckets) { // Last half-filled bucket
                return std::make_pair(bucket - nbuckets + k * nbuckets, 1);
            }
            auto [_, storedSeed] = getThresholdAndSeed(bucket);
            return std::make_pair(bucket, storedSeed);
        }
    };
} // Namespace shockhash
