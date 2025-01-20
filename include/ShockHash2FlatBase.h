#pragma once
#include <cstdint>

namespace morphishash {
struct KeyInfo {
    uint64_t mhc;
    uint32_t bucket;
    uint32_t threshold;
};
}
