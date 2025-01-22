#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace morphishash {
    std::pair<uint64_t, __uint128_t>  morphisHashconstruct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys);
}