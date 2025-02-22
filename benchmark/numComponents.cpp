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



int main() {
    size_t threadscnt=128;
    size_t iters=10000;
    std::vector<std::vector<size_t>> allcounts;
    allcounts.resize(threadscnt);

    auto runb = [&](size_t id) {
        for (size_t l = 2; l <= 60; l+=2) {
            size_t counts = 0;
            for (size_t iter = 0; iter < iters; ++iter) {
                std::vector<uint64_t> keys;
                for (size_t i = 0; i < l; i++) {
                    keys.emplace_back(bytehamster::util::MurmurHash64(std::to_string(id) +" "+std::to_string(i) + " " + std::to_string(l)));
                }
                std::pair<uint64_t, __uint128_t> seed = morphishash::morphisHashconstruct(l, l, keys);
                counts+=seed.first;
            }
            allcounts[id].push_back(counts);
        }
    };

    std::vector<std::thread> threads;
    for (size_t l = 0; l < threadscnt; l += 1) {
        threads.emplace_back(runb, l);
    }

    for (auto &t: threads) {
        if (t.joinable()) t.join();
    }

    size_t i=0;
    for (size_t l = 2; l <= 60; l+=2) {
        size_t total=0;
        for (size_t j = 0; j < threadscnt; ++j) {
            total += allcounts[j][i];
        }
        std::cout<<"RESULT n="<<l<<" ec="<<(double(total)/double(iters*threadscnt))<<std::endl;

        i++;
    }

    return 0;
}
