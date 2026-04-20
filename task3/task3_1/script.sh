#!/bin/bash
for n in 1 2 4 7 8 16 20 40; do
    export export NTHREADS=$n
    echo $n потоков
    ./build/dgemv_thread
    ./build/dgemv_jthread
    ./build/dgemv_async
done