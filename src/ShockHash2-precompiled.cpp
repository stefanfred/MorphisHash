#include "ShockHash2-precompiled.h"
#include "ShockHash2-internal.h"

namespace shockhash {

    template<size_t I, size_t WS, bool burr>
    __uint128_t construct(std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                     std::vector<uint64_t> &leafKeys) {
        using SH = std::conditional_t<I < SockHash2SeedFinderLeafSizeThreshold,
                BijectionsShockHash2<I, BasicSeedCandidateFinder::Finder, true, burr, I - WS>,
                BijectionsShockHash2<I, QuadSplitCandidateFinderList, true, burr, I - WS>>;
        __uint128_t x = SH::findSeed(leafKeys);
        if(burr) {
            constructRetrieval(leafKeys, x, ribbonInput, I);
        }
        return x;
    }

    template<size_t I, size_t WS>
    __uint128_t dispatchWidth(size_t width,
                         std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                         std::vector<uint64_t> &leafKeys, bool useBurr) {
        if (useBurr) {
            return construct<I, 0, true>(ribbonInput, leafKeys);
        } else if(width == 0) {
            return construct<I, I, false>(ribbonInput, leafKeys);
        } else if constexpr (WS >= I || WS > shockhash::MAX_DIFF || I - WS > shockhash::MAX_RETRIEVAL_WIDTH) {
            exit(1);
        } else if (I - WS == width) {
            return construct<I, WS, false>(ribbonInput, leafKeys);
        } else {
            return dispatchWidth<I, WS + 1>(width, ribbonInput, leafKeys, useBurr);
        }
    }

    template<size_t I>
    __uint128_t dispatchLeafSize(size_t leafSize, size_t width,
                            std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                            std::vector<uint64_t> &leafKeys, bool useBurr) {
        if constexpr (I <= 1) {
            std::cerr << "The leafSize " << leafSize << " was not compiled into this binary." << std::endl;
            return 0;
        } else if (I == leafSize) {
            return dispatchWidth<I, 4>(width, ribbonInput, leafKeys, useBurr);
        } else {
            return dispatchLeafSize<I - 1>(leafSize, width, ribbonInput, leafKeys, useBurr);
        }
    }


    __uint128_t shockhash2construct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys,
                               std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput, bool useBurr) {
        return dispatchLeafSize<shockhash::MAX_LEAF_SIZE2>(leafSize, width, ribbonInput, leafKeys, useBurr);
    }

}
