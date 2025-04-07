#pragma once

#include <vector>
#include <cstdint>

namespace morphishash {
    class BitVec {
    private:
        std::vector<uint64_t> data;
    public:
        BitVec() : BitVec(0) {
        }

        explicit BitVec(const size_t bits) : data((bits + 63) / 64) {
        }

        template<uint64_t elementWidth>
        [[nodiscard]] std::conditional<elementWidth <= 64, uint64_t, __uint128_t>::type at(const size_t bitIdx) const {
            if constexpr (elementWidth == 0)
                return 0;

            if constexpr (elementWidth > 64) {
                static_assert(elementWidth <= 128);
                return __uint128_t(at<64>(bitIdx + elementWidth - 64)) |
                       (__uint128_t(at<elementWidth - 64>(bitIdx)) << 64);
            } else {
                static constexpr uint64_t MASK = elementWidth == 64 ? ~0ul : ((1ul << elementWidth) - 1);
                const size_t offset = bitIdx % 64;
                const uint64_t *word = &data[bitIdx / 64];
                uint64_t w1 = (*word) >> offset;
                if (offset + elementWidth > 64) {
                    return (w1 | (*(word + 1) << (64 - offset))) & MASK;
                } else {
                    return w1 & MASK;
                }
            }
        }

        template<uint64_t elementWidth>
        void set(const size_t bitIdx, std::conditional<elementWidth <= 64, uint64_t, __uint128_t>::type value) {
            if constexpr (elementWidth == 0)
                return;

            if constexpr (elementWidth > 64) {
                static_assert(elementWidth <= 128);
                set<elementWidth - 64>(bitIdx, uint64_t(value >> 64));
                set<64>(bitIdx + elementWidth - 64, uint64_t(value));
            } else {
                static constexpr uint64_t MASK = elementWidth == 64 ? ~0ul : ((1ul << elementWidth) - 1);
                uint64_t *word = &data[bitIdx / 64];
                const size_t offset = bitIdx % 64;
                value &= MASK;
                *word &= ~(MASK << offset);
                *word |= (value << offset);
                if (offset + elementWidth >= 64) {
                    *(word + 1) &= ~(MASK >> (64 - offset));
                    *(word + 1) |= (value >> (64 - offset));
                }
            }

        }

        [[nodiscard]] size_t dataSizeBytes() const {
            return data.size() * sizeof(uint64_t);
        }
    };
}