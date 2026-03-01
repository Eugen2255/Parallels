#!/bin/bash
for n in 1 2 4 7 8 16 20 40; do
    export OMP_NUM_THREADS=$n
    ./build/dgemv
done