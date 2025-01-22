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

static constexpr uint64_t rotate(size_t l, uint64_t val, uint32_t x) {
    return ((val << x) | (val >> (l - x))) & ((1ul << l) - 1);
}

double lower(size_t l) {
    double res=1;
    for (size_t i = 1; i <= l; ++i) {
        res=res*double(l)/i;
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
    size_t numIterations = std::max(2ul, (size_t) 7e7 / (1<<l));
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
                    taken1 |= (1<<pos);
                } else {
                    taken2 |= (1<<pos);
                }
            }
            bool canBeRotated = false;
            for (size_t r = 0; r < l; r++) {
                if ((taken1 | rotate(l, taken2, r)) == (1<<l) - 1) {
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
    std::cout<<"RESULT"
             <<" method=rotations"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate="<<log2((double)totalTries / (double)numIterations * l) / l
             <<std::endl;

}

void testBruteForce(size_t l) {
    size_t numIterations = std::max(2ul, (size_t) 1e8 / (1<<l));
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
                if (taken & (1<<pos)) {
                    break;
                }
                taken |= (1<<pos);
            }
            if (i == l) {
                break;
            }
            totalTries++;
        }
    }
    double p = double(numIterations)/double(totalTries);
    double space = (geometricEntropy(p)) / double(l);
    double space2 = (-log2(p)) / double(l);
    std::cout<<"RESULT"
             <<" method=bruteforce"
             <<" l="<<l
             <<" hfEvals="<<(double)hfEvals/(double)numIterations
             <<" tries="<<(double)totalTries / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate1="<<space
             <<" spaceEstimate2="<<space2
             <<" spaceEstimate1o="<<(space*l- lower(l))
             <<" spaceEstimate2o="<<(space2*l- lower(l))
             <<std::endl;


}


void testBipMorphisHash(size_t l, size_t w) {
    assert(l < 120);
    size_t numIterations = l <= 40 ? 500000 : 50000;
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
    double p = numIterations/totalSeed;
    double space = (geometricEntropy(p) + double(l) - double(w)) / double(l);
    double space2 = (-log2(p) + double(l) - double(w)) / double(l);
    std::cout<<"RESULT"
             <<"method = morphisHash"
             <<" l="<<l
             <<" w="<<w
             <<" tries="<<(double)totalLargerPart / (double)numIterations
             <<" iterations="<<numIterations
             <<" spaceEstimate1="<<space
             <<" spaceEstimate2="<<space2
             <<" spaceEstimate1o="<<(space*l- lower(l))
             <<" spaceEstimate2o="<<(space2*l- lower(l))
             <<std::endl;
}

int main() {
    for (size_t l = 10; l <= 80; l += 2) {
        if (l <= 25) {
            testRotationFitting(l);
        }
        if (l <= 22) {
            testBruteForce(l);
        }
        for (int w = 0; w < 7; ++w) {
            testBipMorphisHash(l, w);
        }
    }
}
