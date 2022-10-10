# GPURecSplit / SIMDRecSplit

Parallelization (Threads, SIMD, GPU) of the Minimal Perfect Hash Function
[RecSplit](https://arxiv.org/abs/1910.06416).
The code in this repository achieves a speedup of up to 333 on SIMD machines, and 1873 on GPUs when compared
to the original [RecSplit implementation](https://github.com/vigna/sux/blob/master/sux/function/RecSplit.hpp).
