#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace shockhash {
    size_t shockhash2construct(size_t leafSize, size_t width, std::vector<uint64_t> &leafKeys,
                               std::vector<std::pair<uint64_t, uint8_t>> &ribbonInput, bool burr);
}