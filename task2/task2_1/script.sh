#!/bin/bash
for n in 40; do
    export OMP_NUM_THREADS=$n
    echo $n
    ./build/dgemv
done