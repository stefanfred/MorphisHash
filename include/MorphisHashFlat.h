#pragma once

#include <vector>
#include <map>
#include <bytehamster/util/EliasFano.h>
#include <bytehamster/util/MurmurHash64.h>
#include <tlx/math/integer_log2.hpp>
#include <MorphisHash-internal.h>
#include <MorphisHash.h>
#include <bytehamster/util/IntVector.h>
#include "MorphisHashFlatBase.h"

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

    template<size_t k, size_t RETRIEVAL_DIFF, size_t EXTRA_SEED_BITS = 3>
    class MorphisHashFlat {
        using BaseCase = BijectionsMorphisHash<k, morphishash::QuadSplitCandidateFinderList, true,
                k - RETRIEVAL_DIFF>;
        static constexpr double OVERLOAD_FACTOR = 0.9;
        static constexpr size_t THRESHOLD_BITS = tlx::integer_log2_floor(k) - 1;
        static constexpr size_t THRESHOLD_RANGE = 1ul << THRESHOLD_BITS;
        static constexpr std::array<uint32_t, THRESHOLD_RANGE> THRESHOLD_MAPPING = _fill_mapping<THRESHOLD_RANGE>();
        // static constexpr size_t SEED_BITS = std::ceil(0.442 * k - 0.2 + double(RETRIEVAL_DIFF)*0.65);
        static constexpr size_t WIDTH = (RETRIEVAL_DIFF <= k ? k - RETRIEVAL_DIFF : 0);
        static constexpr size_t SEED_BITS = bij_memoMorphis[RETRIEVAL_DIFF][k] + EXTRA_SEED_BITS - WIDTH;
        static constexpr size_t MAX_SEED = 1ul << SEED_BITS;
        static constexpr size_t SEED_FALLBACK_INDICATOR = 0;
        static constexpr size_t BUCKET_DATA_BITS = THRESHOLD_BITS + SEED_BITS + WIDTH;
        //bytehamster::util::IntVector<THRESHOLD_BITS + SEED_BITS> thresholdsAndSeeds;
        std::map<size_t, size_t> seedsFallback;
        std::vector<size_t> layerBases;

        RiceBitVector<> bucketData;

        MorphisHash<k, RETRIEVAL_DIFF> fallbackPhf;
        size_t N;
        size_t nbuckets;
        pasta::BitVector freePositionsBv;
        pasta::FlatRankSelect <pasta::OptimizedFor::ONE_QUERIES> *freePositionsRankSelect = nullptr;
        size_t layers = 2;
    public:
        explicit MorphisHashFlat(const std::vector<std::string> &keys, size_t ignore, size_t ignore2)
                : MorphisHashFlat(keys) {
            (void) ignore;
            (void) ignore2;
        }

        explicit MorphisHashFlat(const std::vector<std::string> &keys) {
            N = keys.size();
            nbuckets = N / k;
            size_t keysInEndBucket = N - nbuckets * k;
            size_t
                    bucketsThisLayer = std::max(1ul, (size_t)
            std::ceil(OVERLOAD_FACTOR * nbuckets));
            std::vector<size_t> freePositions;
            std::vector<KeyInfo> hashes;
            hashes.reserve(keys.size());
            for (const std::string &key: keys) {
                uint64_t mhc = ::bytehamster::util::MurmurHash64(key);
                uint32_t bucket = ::bytehamster::util::fastrange32(mhc & 0xffffffff, bucketsThisLayer);
                uint32_t threshold = mhc >> 32;
                hashes.emplace_back(mhc, bucket, threshold);
            }
            std::vector<KeyInfo> allHashes = hashes;
            layerBases.push_back(0);

            bytehamster::util::IntVector<THRESHOLD_BITS> thresholds(nbuckets);
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
                        hash.mhc = ::bytehamster::util::remix(hash.mhc);
                        hash.bucket = ::bytehamster::util::fastrange32(hash.mhc & 0xffffffff, bucketsThisLayer);
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
                        flushBucket(layerBase, bucketStart, i, previousBucket, hashes, bumpedKeys, freePositions, thresholds);
                        previousBucket++;
                        bucketStart = i;
                    }
                }
                // Last bucket
                while (previousBucket < bucketsThisLayer) {
                    flushBucket(layerBase, bucketStart, hashes.size(), previousBucket, hashes, bumpedKeys,
                                freePositions, thresholds);
                    previousBucket++;
                    bucketStart = hashes.size();
                }
                hashes = bumpedKeys;
            }

            std::vector<std::string> fallbackPhfContent;
            for (auto &hash: hashes) {
                fallbackPhfContent.push_back(std::to_string(hash.mhc));
            }

            fallbackPhf = MorphisHash<k, RETRIEVAL_DIFF>(fallbackPhfContent, 2000, 1);
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

            // Construct MorphisHash within buckets
            std::vector<std::vector<uint64_t>> bucketContents(nbuckets);
            std::vector<uint64_t> lastBucket;
            for (KeyInfo key: allHashes) {
                size_t bucket = get<0>(evaluateKPerfect(key.mhc, thresholds));
                if (bucket >= nbuckets) {
                    lastBucket.push_back(key.mhc);
                } else {
                    bucketContents.at(bucket).push_back(key.mhc);
                }
            }
            RiceBitVector<>::Builder bucketDataBuilder;
            for (size_t i = 0; i < nbuckets; i++) {
                std::pair<uint64_t, __uint128_t> seed = BaseCase::findSeed(bucketContents.at(i));
                if (seed.first >= MAX_SEED || seed.first == SEED_FALLBACK_INDICATOR) {
                    seedsFallback.insert(std::make_pair(i, seed.first));
                    seed.first = SEED_FALLBACK_INDICATOR;
                }
                bucketDataBuilder.appendFixed(thresholds.at(i), THRESHOLD_BITS);
                bucketDataBuilder.appendFixed(seed.first, SEED_BITS);
                bucketDataBuilder.appendFixed128(seed.second, WIDTH);
            }

            bucketData = bucketDataBuilder.build();
        }

        inline std::tuple<size_t, size_t, __uint128_t> getThresholdAndSeedAndRetrieval(size_t bucket) {
            RiceBitVector<>::Reader r = bucketData.reader();
            r.toFixedPos(BUCKET_DATA_BITS, bucket);

            size_t threshold = r.readFixed(THRESHOLD_BITS);
            size_t seed = r.readFixed(SEED_BITS);
            __uint128_t sol = r.readFixed128(WIDTH);
            return {threshold, seed, sol};
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
                         std::vector<size_t> &freePositions, bytehamster::util::IntVector<THRESHOLD_BITS> &thresholds) {
            size_t bucketSize = i - bucketStart;
            if (bucketSize <= k) {
                size_t threshold = THRESHOLD_RANGE - 1;
                thresholds.set(layerBase + bucketIdx, threshold);
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
                thresholds.set(layerBase + bucketIdx, threshold);
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
                   + bucketData.getBits()
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
            std::cout << "Retrieval: " << bucketData.getBits() / N << std::endl;
        }

        size_t operator()(const std::string &key) {
            return operator()(::bytehamster::util::MurmurHash64(key));
        }

        size_t operator()(const hash128_t &hash) {
            return operator()(hash.second);
        }

        __attribute_noinline__ size_t operator()(uint64_t hash) {
            auto [bucket, seed, sol] = evaluateKPerfect(hash, {});
            if (bucket >= nbuckets) {
                return bucket; // N that are not multiples of k
            }
            if (seed == SEED_FALLBACK_INDICATOR) {
                seed = seedsFallback.at(bucket);
            }
            size_t baseCase;
            if constexpr (WIDTH == 0) {
                baseCase = sh2remix64(hash ^ seed) % k;
            } else {
                __uint128_t remixed = sh2remix128(hash ^ seed);
                constexpr __uint128_t row_mask = (__uint128_t(1) << (WIDTH)) - 1;
                uint64_t retrieved = parity(sol & row_mask & remixed);
                baseCase = queryHash(seed, hash, retrieved, k);
            }
            return bucket * k + baseCase;
        }

        /** Returns bucket and seed */
        inline std::tuple<size_t, size_t, __uint128_t> evaluateKPerfect(uint64_t mhc, const std::optional<const bytehamster::util::IntVector<THRESHOLD_BITS>> &constructThresholds) {
            for (size_t layer = 0; layer < layers; layer++) {
                if (layer != 0) {
                    mhc = ::bytehamster::util::remix(mhc);
                }
                size_t base = layerBases.at(layer);
                size_t layerSize = layerBases.at(layer + 1) - base;
                uint32_t bucket = ::bytehamster::util::fastrange32(mhc & 0xffffffff, layerSize);
                uint32_t threshold = mhc >> 32;
                if(constructThresholds.has_value()) {
                    if (threshold <= THRESHOLD_MAPPING[constructThresholds->at(base+bucket)]) {
                        return {base+bucket, 0, 0};
                    }
                } else {
                    auto [storedThreshold, storedSeed, sol] = getThresholdAndSeedAndRetrieval(base + bucket);
                    if (threshold <= THRESHOLD_MAPPING[storedThreshold]) {
                        return {base+bucket, storedSeed, sol};
                    }
                }
            }
            size_t phf = fallbackPhf(std::to_string(mhc));
            size_t bucket = freePositionsRankSelect->select1(phf + 1) - phf;
            if (bucket >= nbuckets) { // Last half-filled bucket
                return {bucket - nbuckets + k * nbuckets, 1,0};
            }
            if(constructThresholds.has_value()) {
                return {bucket,0,0};
            } else {
                auto [_, storedSeed, sol] = getThresholdAndSeedAndRetrieval(bucket);
                return {bucket, storedSeed, sol};
            }
        }
    };
} // Namespace morphishash
