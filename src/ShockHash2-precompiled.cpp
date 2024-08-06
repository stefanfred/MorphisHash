#include "ShockHash2-precompiled.h"
#include "ShockHash2-internal.h"

namespace shockhash {

    template<size_t I, size_t WS, bool burr>
    size_t construct(std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                     std::vector<uint64_t> &leafKeys) {
        using SH = std::conditional_t<false,
                BijectionsShockHash2<I, QuadSplitCandidateFinderList, true, I - WS>,
                BijectionsShockHash2<I, BasicSeedCandidateFinder::Finder, true, I - WS>>;
        size_t x = SH::findSeed(leafKeys);
        if(burr) {
            constructRetrieval(leafKeys, x, ribbonInput, I);
        }
        return x;
    }

    template<size_t I, size_t WS>
    size_t dispatchWidth(size_t width,
                         std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                         std::vector<uint64_t> &leafKeys) {
        if constexpr (WS >= I || WS > shockhash::MAX_DIFF || I - WS > shockhash::MAX_RETRIEVAL_WIDTH) {
            std::cout << "Using Burr" << std::endl;
            return construct<I, I, true>(ribbonInput, leafKeys);
        } else if (I - WS == width) {
            return construct<I, WS, false>(ribbonInput, leafKeys);
        } else {
            return dispatchWidth<I, WS + 1>(width, ribbonInput, leafKeys);
        }
    }

    template<size_t I>
    size_t dispatchLeafSize(size_t leafSize, size_t width,
                            std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput,
                            std::vector<uint64_t> &leafKeys) {
        if constexpr (I <= 1) {
            std::cerr << "The leafSize " << leafSize << " was not compiled into this binary." << std::endl;
            return 0;
        } else if (I == leafSize) {
            return dispatchWidth<I, 0>(width, ribbonInput, leafKeys);
        } else {
            return dispatchLeafSize<I - 2>(leafSize, width, ribbonInput, leafKeys);
        }
    }

    size_t shockhash2construct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys,
                               std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput) {
        return dispatchLeafSize<shockhash::MAX_LEAF_SIZE2>(leafSize, width, ribbonInput, leafKeys);
    }

}
