#include "ShockHash2-precompiled.h"
#include "ShockHash2-internal.h"

namespace morphishash {

    template<size_t I, size_t W, bool burr>
    std::pair<uint64_t, __uint128_t> construct(std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                                               std::vector<uint64_t> &leafKeys) {
        using SH = std::conditional_t<I < SockHash2SeedFinderLeafSizeThreshold,
                BijectionsShockHash2<I, BasicSeedCandidateFinder::Finder, true, burr, W>,
                BijectionsShockHash2<I, QuadSplitCandidateFinderList, true, burr, W>>;
        std::pair<uint64_t, __uint128_t> x = SH::findSeed(leafKeys);
        /*if(burr) {
            constructRetrieval(leafKeys, x, ribbonInput, I);
        }*/
        return x;
    }

    template<size_t I, size_t W>
    std::pair<uint64_t, __uint128_t> dispatchWidth(size_t width,
                                                   std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                                                   std::vector<uint64_t> &leafKeys, bool useBurr) {
        //std::cout << W << " " << I << " " << width << std::endl;
        /*if (useBurr) {
            return construct<I, 0, true>(ribbonInput, leafKeys);
        } else*/
        if (W == width) {
            return construct<I, W, false>(ribbonInput, leafKeys);
        }
        if constexpr (W > 0 && I - W < MAX_DIFF) {
            return dispatchWidth<I, W - 1>(width, ribbonInput, leafKeys, useBurr);
        }
    }

    template<size_t I>
    std::pair<uint64_t, __uint128_t> dispatchLeafSize(size_t leafSize, size_t width,
                                                      std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                                                      std::vector<uint64_t> &leafKeys, bool useBurr) {
        if constexpr (I <= 1) {
            std::cerr << "The leafSize " << leafSize << " was not compiled into this binary." << std::endl;
            exit(1);
        } else if (I == leafSize) {
            return dispatchWidth<I, I>(width, ribbonInput, leafKeys, useBurr);
        } else {
            return dispatchLeafSize<I - 1>(leafSize, width, ribbonInput, leafKeys, useBurr);
        }
    }


    std::pair<uint64_t, __uint128_t> shockhash2construct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys,
                                                         std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                                                         bool useBurr) {
        return dispatchLeafSize<morphishash::MAX_LEAF_SIZE2>(leafSize, width, ribbonInput, leafKeys, useBurr);
    }

}
