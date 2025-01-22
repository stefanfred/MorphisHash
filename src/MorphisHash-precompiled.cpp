#include "MorphisHash-precompiled.h"
#include "MorphisHash-internal.h"

namespace morphishash {

    template<size_t I, size_t W>
    std::pair<uint64_t, __uint128_t> construct(std::vector<uint64_t> &leafKeys) {
        using MH = std::conditional_t<I < SockHash2SeedFinderLeafSizeThreshold,
                BijectionsMorphisHash<I, BasicSeedCandidateFinder::Finder, true, W>,
                BijectionsMorphisHash<I, QuadSplitCandidateFinderList, true, W>>;
        std::pair<uint64_t, __uint128_t> x = MH::findSeed(leafKeys);
        return x;
    }

    template<size_t I, size_t W>
    std::pair<uint64_t, __uint128_t> dispatchWidth(size_t width, std::vector<uint64_t> &leafKeys) {
        if (W == width) {
            return construct<I, W>(leafKeys);
        }
        if constexpr (W > 0 && I - W < MAX_DIFF) {
            return dispatchWidth<I, W - 1>(width, leafKeys);
        }
        std::cerr<<"not included in compilation"<<std::endl;
        exit(1);
    }

    template<size_t I>
    std::pair<uint64_t, __uint128_t> dispatchLeafSize(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys) {
        if constexpr (I <= 1) {
            std::cerr << "The leafSize " << leafSize << " was not compiled into this binary." << std::endl;
            exit(1);
        } else if (I == leafSize) {
            return dispatchWidth<I, I>(width, leafKeys);
        } else {
            return dispatchLeafSize<I - 1>(leafSize, width, leafKeys);
        }
    }


    std::pair<uint64_t, __uint128_t> morphisHashconstruct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys) {
        return dispatchLeafSize<morphishash::MAX_LEAF_SIZE>(leafSize, width, leafKeys);
    }

}
