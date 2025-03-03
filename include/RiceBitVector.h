/*
 * Sux: Succinct data structures
 *
 * Copyright (C) 2019-2020 Emmanuel Esposito and Sebastiano Vigna
 *
 *  This library is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by the Free
 *  Software Foundation; either version 3 of the License, or (at your option)
 *  any later version.
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * Under Section 7 of GPL version 3, you are granted additional permissions
 * described in the GCC Runtime Library Exception, version 3.1, as published by
 * the Free Software Foundation.
 *
 * You should have received a copy of the GNU General Public License and a copy of
 * the GCC Runtime Library Exception along with this program; see the files
 * COPYING3 and COPYING.RUNTIME respectively.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <sux/util/Vector.hpp>
#include <sux/support/common.hpp>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <vector>
#include <iostream>

namespace morphishash {

    using namespace std;
#define expect(x, y) __builtin_expect(x, y)

/** Storage for Golomb-Rice codes of a RecSplit bucket.
 *
 * This class exists solely to implement RecSplit.
 * @tparam AT a type of memory allocation out of util::AllocType.
 */
    template<sux::util::AllocType AT = sux::util::AllocType::MALLOC>
    class RiceBitVector {

    public:
        class Builder {
            sux::util::Vector<uint64_t, AT> data;
            size_t bit_count = 0;

        public:
            Builder() : Builder(16) {}

            Builder(const size_t alloc_words) : data(alloc_words) {}

            Builder &operator=(Builder &&builder) {
                data = std::move(builder.data);
                bit_count = builder.bit_count;
                return *this;
            }

            void appendFixed128(const __uint128_t v, const int width) {
                if (width > 64) {
                    appendFixed(uint64_t(v >> 64), width - 64);
                    appendFixed(uint64_t(v), 64);
                } else {
                    appendFixed(v, width);
                }
            }

            void appendFixed(const uint64_t v, const int log2golomb) {
                if(log2golomb == 0) {
                    return;
                }
                const uint64_t lower_bits = v & ((uint64_t(1) << (log2golomb - 1) << 1) - 1); // weird shifting to support 64 bit
                int used_bits = bit_count & 63;

                data.resize((((bit_count + log2golomb + 7) / 8) + 7 + 7) / 8);

                uint64_t *append_ptr = &data + bit_count / 64;
                uint64_t cur_word = *append_ptr;

                cur_word |= lower_bits << used_bits;
                if (used_bits + log2golomb > 64) {
                    *(append_ptr++) = cur_word;
                    cur_word = lower_bits >> (64 - used_bits);
                    used_bits += log2golomb - 64;
                }
                *append_ptr = cur_word;
                bit_count += log2golomb;
            }

            void appendUnaryAll(const std::vector<uint32_t> &unary) {
                size_t bit_inc = 0;
                for (const auto &u: unary) {
                    bit_inc += u + 1;
                }

                data.resize((((bit_count + bit_inc + 7) / 8) + 7 + 7) / 8);

                for (const auto &u: unary) {
                    bit_count += u;
                    uint64_t *append_ptr = &data + bit_count / 64;
                    *append_ptr |= uint64_t(1) << (bit_count & 63);
                    ++bit_count;
                }
            }

            void appendRiceBitVector(const Builder &other) {
                data.resize((((bit_count + other.bit_count + 7) / 8) + 7 + 7) / 8);
                if (bit_count % 8 == 0) {
                    std::memcpy((uint8_t *) &data + bit_count / 8, (uint8_t *) &other.data, (other.bit_count + 7) / 8);
                } else {
                    const int used_bits = bit_count & 63;
                    uint64_t *append_ptr = &data + bit_count / 64;
                    uint64_t cur_word = *append_ptr;
                    for (size_t i = 0; i < (other.bit_count + 63) / 64; ++i) {
                        const uint64_t value = other.data[i];
                        cur_word |= value << used_bits;
                        *(append_ptr++) = cur_word;
                        cur_word = value >> (64 - used_bits);
                    }
                    if (cur_word != 0)
                        *append_ptr = cur_word;
                }
                bit_count += other.bit_count;
            }

            uint64_t getBits() { return bit_count; }

            RiceBitVector <AT> build() {
                data.trimToFit();
                return RiceBitVector(std::move(data));
            }
        };

    private:
        sux::util::Vector<uint64_t, AT> data;

        friend std::ostream &operator<<(std::ostream &os, const RiceBitVector<AT> &rbv) {
            os << rbv.data;
            return os;
        }

        friend std::istream &operator>>(std::istream &is, RiceBitVector<AT> &rbv) {
            is >> rbv.data;
            return is;
        }

    public:
        RiceBitVector() {}

        RiceBitVector(sux::util::Vector<uint64_t, AT> data) : data(std::move(data)) {}

        size_t getBits() const { return data.size() * sizeof(uint64_t) * 8; }

        class Reader {
            size_t curr_fixed_offset = 0;
            uint64_t curr_window_unary = 0;
            uint64_t *curr_ptr_unary;
            int valid_lower_bits_unary = 0;
            sux::util::Vector<uint64_t, AT> &data;

        public:
            Reader(sux::util::Vector<uint64_t, AT> &data) : data(data) {}

            uint64_t readNext(const int log2golomb) {
                uint64_t result = 0;

                if (curr_window_unary == 0) {
                    result += valid_lower_bits_unary;
                    curr_window_unary = *(curr_ptr_unary++);
                    valid_lower_bits_unary = 64;
                    while (expect(curr_window_unary == 0, 0)) {
                        result += 64;
                        curr_window_unary = *(curr_ptr_unary++);
                    }
                }

                const size_t pos = sux::rho(curr_window_unary);

                curr_window_unary >>= pos;
                curr_window_unary >>= 1;
                valid_lower_bits_unary -= pos + 1;

                result += pos;
                result <<= log2golomb;


                // In the worst case, curr_fixed_offset is 7.
                // Loading only 64 bit here would then fail for log2golomb>(64-7)=57
                __uint128_t fixed;
                memcpy(&fixed, (uint8_t *) &data + curr_fixed_offset / 8, 16);
                result |= (fixed >> curr_fixed_offset % 8) & ((uint64_t(1) << log2golomb) - 1);
                curr_fixed_offset += log2golomb;
                return result;
            }

            __uint128_t readFixed128(const int width) {
                if (width > 64) {
                    return (__uint128_t(readFixed(width - 64)) << 64) | readFixed(64);
                } else {
                    return readFixed(width);
                }
            }

            uint64_t readFixed(const int width) {
                // In the worst case, curr_fixed_offset is 7.
                // Loading only 64 bit here would then fail for width>(64-7)=57
                __uint128_t fixed;
                memcpy(&fixed, (uint8_t *) &data + curr_fixed_offset / 8, 16);
                uint64_t res = (fixed >> curr_fixed_offset % 8) & ((__uint128_t(1) << width) - 1);
                curr_fixed_offset += width;
                //std::cout<<width<<" WIDTH "<< res<<std::endl;
                return res;
            }

            void skipSubtree(const size_t nodes, const size_t fixed_len) {
                assert(nodes > 0);
                size_t missing = nodes, cnt;
                while ((cnt = sux::nu(curr_window_unary)) < missing) {
                    curr_window_unary = *(curr_ptr_unary++);
                    missing -= cnt;
                    valid_lower_bits_unary = 64;
                }
                cnt = sux::select64(curr_window_unary, missing - 1);
                curr_window_unary >>= cnt;
                curr_window_unary >>= 1;
                valid_lower_bits_unary -= cnt + 1;

                curr_fixed_offset += fixed_len;
            }

            void toFixedPos(size_t w, size_t index) {
                curr_fixed_offset=w*index;
            }

            void readReset(const size_t bit_pos, const size_t unary_offset) {
                // assert(bit_pos < bit_count);
                curr_fixed_offset = bit_pos;
                size_t unary_pos = bit_pos + unary_offset;
                curr_ptr_unary = &data + unary_pos / 64;
                curr_window_unary = *(curr_ptr_unary++) >> (unary_pos & 63);
                valid_lower_bits_unary = 64 - (unary_pos & 63);
            }
        };

        Reader reader() { return Reader(data); }
    };

} // namespace morphishash
