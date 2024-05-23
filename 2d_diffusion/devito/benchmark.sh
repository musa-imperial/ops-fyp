#!/bin/bash

max_thread_count=2
max_core_count=8

compiler="icx"

parallel="openmp"

echo "Benchmark started"
if [ "$compiler" = "icx" ]; then
    echo "compiler: icx"
    if [ "$parallel" = "openmp" ]; then
            echo -n " Runtime, Max, Thread count" > benchmark_output.txt
            for (( t=1; t<=$max_thread_count; t++ )); do
                export OMP_NUM_THREADS=$t
                sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) devitocodes/devito:icx-latest python 2d_diffusion.py
                # echo -n "$pi $t" >> benchmark_output.txt
                # echo -n >> benchmark_output.txt
                # echo "Runtime, Max"
                # echo "$pi"
                echo "Thread count: $t"
            done
        fi
elif [ "$compiler" = "gnu" ]; then
    echo "compiler: gnu"
    if [ "$parallel" = "openmp" ]; then
        echo -n " Runtime, Max, Thread count" > benchmark_output.txt
        for (( t=1; t<=$max_thread_count; t++ )); do
            export OMP_NUM_THREADS=$t
            sudo docker run --rm -it -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) devitocodes/devito:gcc-latest python 2d_diffusion.py
            # echo -n "$pi $t" >> benchmark_output.txt
            # echo -n >> benchmark_output.txt
            # echo "Runtime, Max"
            # echo "$pi"
            echo "Thread count: $t"
        done
    fi
fi


