#pragma once
#include <vector>
#include <cassert>
#include <queue>
#include <bytehamster/util/Function.h>
#include <bytehamster/util/MurmurHash64.h>
#include <cstring>

namespace morphishash {
struct HashedKey {
    uint64_t mhc;

    HashedKey() {
        this->mhc = 0;
    }

    explicit HashedKey(uint64_t mhc) : mhc(mhc) {
    }

    explicit HashedKey(const std::string &element, uint32_t seed = 0) {
        uint64_t stringHash = bytehamster::util::MurmurHash64(element.data(), element.length());
        uint64_t modified = stringHash + seed;
        mhc = bytehamster::util::MurmurHash64(&modified, sizeof(uint64_t));
    }

    [[nodiscard]] inline uint64_t hash(int hashFunctionIndex, size_t range) const {
        return bytehamster::util::fastrange64(bytehamster::util::remix(mhc + hashFunctionIndex), range);
    }
};

/**
 * Tiny binary cuckoo hash table. Construction needs multiple tries before succeeding.
 */
class TinyBinaryCuckooHashTable {
    public:
        struct TableEntry {
            HashedKey hash;
            uint32_t candidateCellsXor = 0;
        };
        TableEntry *heap;
        TableEntry** cells;
        const size_t NMax;
    private:
        size_t seed = 0;
        size_t numEntries = 0;
    public:
        explicit TinyBinaryCuckooHashTable(size_t NMax) : NMax(NMax) {
            heap = new TableEntry[NMax];
            cells = new TableEntry*[NMax];
        }

        ~TinyBinaryCuckooHashTable() {
            delete[] heap;
            delete[] cells;
        }

        void clear() {
            numEntries = 0;
        }

        void prepare(HashedKey hash) {
            assert(numEntries < NMax);
            heap[numEntries].hash = hash;
            numEntries++;
        }

        bool testPassesFilter(size_t seed) {
            assert(numEntries <= 64);
            uint64_t used = 0;
            for (size_t i = 0; i < numEntries; i++) {
                CandidateCells hash = getCandidateCells(heap[i].hash, seed, numEntries);
                used |= (1ul << hash.cell1) | (1ul << hash.cell2);
            }
            return used == (1ul << numEntries) - 1;
        }

        bool construct(size_t seed_) {
            seed = seed_;

            /*
            // Filter based on unused cells
            if (!testPassesFilter(seed_)) {
                return false;
            }

            // Filter based on tree/pseudotree
            unionFind.clear();
            for (size_t i = 0; i < numEntries; i++) {
                Union64 hash = getCandidateCells(heap[i].hash, seed, M);
                if (!unionFind.unionIsStillPseudoforest(hash.halves.high, hash.halves.low)) {
                    return false;
                }
            }*/

            // Actual cuckoo hashing, we know that this will succeed at this point
            clearPlacement();
            for (size_t i = 0; i < numEntries; i++) {
                if (!insert(&heap[i])) {
                    return false;
                }
            }
            return true;
        }

        void clearPlacement() {
            memset(cells, 0, numEntries * sizeof(void*)); // Fill with nullpointers
        }

        [[nodiscard]] size_t size() const {
            return numEntries;
        }

        static inline size_t hashToCell(HashedKey key, size_t seed, size_t range, size_t hashFunctionIndex) {
            CandidateCells hash = getCandidateCells(key, seed, range);
            if (hashFunctionIndex == 0) {
                return hash.cell1;
            } else {
                return hash.cell2;
            }
        }

        typedef struct {
            uint32_t cell1;
            uint32_t cell2;
        } CandidateCells;

        template <size_t range>
        static inline CandidateCells getCandidateCells(const HashedKey key, size_t seed) {
            uint64_t remixed = bytehamster::util::remix(key.mhc + seed);
            const uint32_t hash1 = bytehamster::util::fastrange32<range / 2>(remixed);
            const uint32_t hash2 = bytehamster::util::fastrange32<(range + 1) / 2>(remixed >> 32) + range / 2;
            return {hash1, hash2};
        }

        static inline CandidateCells getCandidateCells(const HashedKey key, size_t seed, size_t range) {
            uint64_t remixed = bytehamster::util::remix(key.mhc + seed);
            const uint32_t hash1 = bytehamster::util::fastrange32(remixed, range / 2);
            const uint32_t hash2 = bytehamster::util::fastrange32(remixed >> 32, (range + 1) / 2) + range / 2;
            /* // Performance is slightly worse, but produces slightly better hash distribution (0.005 bpk better)
            const uint32_t hash1 = bytehamster::util::fastrange32(remixed, range);
            uint32_t hash2 = bytehamster::util::fastrange32(remixed >> 32, range);
            if (hash1 == hash2) {
                hash2 = bytehamster::util::fastrange32((remixed >> 32) ^ (remixed & 0xffffffff), range);
            }*/
            return {hash1, hash2};
        }

        bool insert(TableEntry *entry) {
            return insert(entry, getCandidateCells(entry->hash, seed, numEntries));
        }

        bool insert(TableEntry *entry, CandidateCells candidates) {
            TableEntry *origEntry = entry;
            assert(candidates.cell1 < numEntries && candidates.cell2 < numEntries);
            assert((unsigned long)(entry - heap) < numEntries * sizeof(TableEntry));
            entry->candidateCellsXor = candidates.cell1 ^ candidates.cell2;
            if (cells[candidates.cell1] == nullptr) {
                cells[candidates.cell1] = entry;
                return true;
            }
            if (cells[candidates.cell2] == nullptr) {
                cells[candidates.cell2] = entry;
                return true;
            }
            uint32_t currentCell = candidates.cell1;
            uint32_t origCell = currentCell;

            do {
                uint32_t alternativeCell = entry->candidateCellsXor ^ currentCell;
                std::swap(entry, cells[alternativeCell]);
                if (entry == nullptr) {
                    return true;
                }
                currentCell = alternativeCell;
            } while (currentCell != origCell || entry != origEntry);
            return false;
        }
};
} // Namespace morphishash
