#include <iostream>
#include <vector>
#include <cmath>
#include "../include/TinyBinaryCuckooHashTable.h"
#include "../include/UnionFind.h"
#include <bytehamster/util/XorShift64.h>
#include <chrono>
#include <set>
#include <unordered_set>
#include <MorphisHash-precompiled.h>
#include <PairingFunction.h>
#include <thread>

#include <iostream>
#include <sstream>
#include <mutex>

#define terr ThreadStream(std::cerr)
#define tout ThreadStream(std::cout)

class ThreadStream : public std::ostringstream
{
public:
    ThreadStream(std::ostream& os) : os_(os)
    {
        // copyfmt causes odd problems with lost output
        // probably some specific flag
//            copyfmt(os);
        // copy whatever properties are relevant
        imbue(os.getloc());
        precision(os.precision());
        width(os.width());
        setf(std::ios::fixed, std::ios::floatfield);
    }

    ~ThreadStream()
    {
        std::lock_guard<std::mutex> guard(_mutex_threadstream);
        os_ << this->str();
    }

private:
    static std::mutex _mutex_threadstream;
    std::ostream& os_;
};

std::mutex ThreadStream::_mutex_threadstream{};
static constexpr uint64_t rotate(size_t l, uint64_t val, uint32_t x) {
    return ((val << x) | (val >> (l - x))) & ((1ul << l) - 1);
}

double lower(size_t l) {
    double res = 1;
    for (size_t i = 1; i <= l; ++i) {
        res = res * double(l) / i;
    }
    return log2(res);
}

double geometricEntropy(double p) {
    if (p <= 0.0 || p > 1.0) {
        throw std::invalid_argument("The probability p must be in the range (0, 1].");
    }
    double q = 1.0 - p;
    return -(q * log2(q) / p) - log2(p);
}

void testRotationFitting(size_t l) {
    size_t numIterations = std::max(2ul, (size_t) 7e7 / (1 << l));
    size_t totalTries = 0;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<morphishash::HashedKey> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(std::to_string(i) + " " + std::to_string(iteration));
        }
        while (true) {
            uint64_t taken1 = 0;
            uint64_t taken2 = 0;
            for (size_t i = 0; i < l; i++) {
                size_t pos = keys[i].hash(totalTries, l);
                hfEvals++;
                if (keys[i].mhc % 2 == 0) {
                    taken1 |= (1 << pos);
                } else {
                    taken2 |= (1 << pos);
                }
            }
            bool canBeRotated = false;
            for (size_t r = 0; r < l; r++) {
                if ((taken1 | rotate(l, taken2, r)) == (1U << l) - 1U) {
                    canBeRotated = true;
                    break;
                }
            }
            if (canBeRotated) {
                break;
            }
            totalTries++;
        }
    }
    tout<< "RESULT"
              << " method=rotations"
              << " l=" << l
              << " hfEvals=" << (double) hfEvals / (double) numIterations
              << " tries=" << (double) totalTries / (double) numIterations
              << " iterations=" << numIterations
              << " spaceEstimate=" << log2((double) totalTries / (double) numIterations * l) / l
              << std::endl;

}

void testBruteForce(size_t l) {
    size_t numIterations = size_t(5e9 / pow(2.71, double(l)));
    if (numIterations < 4) {
        return;
    }
    size_t totalTries = 1;
    size_t hfEvals = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<morphishash::HashedKey> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(std::to_string(i) + " " + std::to_string(iteration));
        }
        while (true) {
            uint64_t taken = 0;
            size_t i = 0;
            for (; i < l; i++) {
                size_t pos = keys[i].hash(totalTries, l);
                hfEvals++;
                if (taken & (1 << pos)) {
                    break;
                }
                taken |= (1 << pos);
            }
            if (i == l) {
                break;
            }
            totalTries++;
        }
    }
    double p = double(numIterations) / double(totalTries);
    double space = (geometricEntropy(p)) / double(l);
    double space2 = (-log2(p)) / double(l);
    tout << "RESULT"
              << " method=bruteforce"
              << " l=" << l
              << " hfEvals=" << (double) hfEvals / (double) numIterations
              << " tries=" << (double) totalTries / (double) numIterations
              << " iterations=" << numIterations
              << " spaceEstimate1=" << space
              << " spaceEstimate2=" << space2
              << " spaceEstimate1o=" << (space * l - lower(l))
              << " spaceEstimate2o=" << (space2 * l - lower(l))
              << std::endl;


}


void testBipMorphisHash(size_t l, size_t w) {
    assert(l < 120);
    size_t numIterations = 2 + size_t(1e8 / pow(1.166, double(l)));
    if (numIterations < 4) {
        return;
    }
    double totalLargerPart = 0;
    double totalSeed = 0;
    for (size_t iteration = 0; iteration < numIterations; iteration++) {
        std::vector<uint64_t> keys;
        for (size_t i = 0; i < l; i++) {
            keys.emplace_back(bytehamster::util::MurmurHash64(std::to_string(i) + " " + std::to_string(iteration)));
        }
        // WARNING: To use this, switch to BasicSeedCandidateFinder in MorphisHash-precompiled.h!
        std::pair<uint64_t, __uint128_t> seed = morphishash::morphisHashconstruct(l, l - w, keys);
        auto [largerPart, smallerPart] = morphishash::unpairTriangular(seed.first);
        totalSeed += seed.first;
        totalLargerPart += largerPart;
    }
    double p = numIterations / totalSeed;
    double space = (geometricEntropy(p) + double(l) - double(w)) / double(l);
    double space2 = (-log2(p) + double(l) - double(w)) / double(l);
    tout << "RESULT"
              << " method=morphisHash"
              << " l=" << l
              << " w=" << w
              << " tries=" << (double) totalLargerPart / (double) numIterations
              << " iterations=" << numIterations
              << " spaceEstimate1=" << space
              << " spaceEstimate2=" << space2
              << " spaceEstimate1o=" << (space * l - lower(l))
              << " spaceEstimate2o=" << (space2 * l - lower(l))
              << std::endl;
}

int main() {
    auto runb = [&](uint64_t l) {
        testBruteForce(l);
    };
    auto runm = [&](uint64_t l, uint64_t  w) {
        testBipMorphisHash(l,w);
    };

    std::vector<std::thread> threads;
    for (size_t l = 10; l <= 80; l += 2) {
        threads.emplace_back(runb, l);
        for (int w = 0; w < 7; ++w) {
            threads.emplace_back(runm, l,w);
        }
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }


}
