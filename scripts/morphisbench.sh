#!/bin/bash

hostname
for leafSize in $(seq 20 4 50); do
  for mw in $(seq 1 2 7); do
    params="--numObjects 100K --numQueries 1M --leafSize $leafSize -w $mw"
    ./Benchmark $params -2
    ./Benchmark $params -f
  done
done
