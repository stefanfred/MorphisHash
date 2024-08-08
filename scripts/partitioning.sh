#!/bin/bash

hostname
for leafSize in $(seq 8 4 42); do
    params="--numObjects 10M --numQueries 10M -w 4 --leafSize $leafSize"
    ./Benchmark $params --shockhash2 --bucketSize 2000
done
