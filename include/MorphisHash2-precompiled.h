#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace morphishash {
    std::pair<uint64_t, __uint128_t>  shockhash2construct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys,
                               std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput, bool burr);
}