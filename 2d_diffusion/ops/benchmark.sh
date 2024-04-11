#!/bin/bash

max_thread_count=48
max_core_count=4

compiler="gnu"

parallel="openmp"

export OPS_INSTALL_PATH=/home/musa/apps/OPS/ops

if [ "$compiler" = "icx" ]; then
    export OPS_COMPILER=icx
    export INTEL_INSTALL_PATH=/home/musa/intel/oneapi

    export PATH=$INTEL_INSTALL_PATH/compiler/latest/bin:$PATH
    export MPICH_CC=icx
    export MPICH_CXX=icpx

    
    source $INTEL_INSTALL_PATH/compiler/latest/env/vars.sh
    source $INTEL_INSTALL_PATH/mpi/latest/env/vars.sh

    export MPI_INSTALL_PATH=$INTEL_INSTALL_PATH/mpi/latest
    export MPICPP=$MPI_INSTALL_PATH/bin/mpicxx

    unset HDF5_INSTALL_PATH
    export HDF5_INSTALL_PATH=/home/musa/apps/asl/install/build_hdf5/intel
    export LD_LIBRARY_PATH=$HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH
    export PATH=$HDF5_INSTALL_PATH/bin:$PATH

elif [ "$compiler" = "gnu" ]; then
    export OPS_COMPILER=gnu
    export HDF5_INSTALL_PATH=$HOME/apps/hdf5
fi

./removescript.sh
make diffusion_seq

if [ "$parallel" = "openmp" ]; then
    echo -n " Runtime, Max, Thread count" > benchmark_output.txt
    for (( t=1; t<=$max_thread_count; t++ )); do
        export OMP_NUM_THREADS=$t
        pi=$(./diffusion_seq)
        echo -n "$pi, $t" >> benchmark_output.txt
        echo -n >> benchmark_output.txt
        echo "Runtime, Max"
        echo "$pi"
        echo "Thread count: $t"
    done
elif [ "$parallel" = "mpi" ]; then
    echo -n " Runtime, Max, Core count" > benchmark_output.txt
    for (( t=1; t<=$max_core_count; t++ )); do
        export OMP_NUM_THREADS=1
        pi=$(mpiexec -np $t ./diffusion_seq)
        for (( s=$t; s>=2; s-- )); do
            pi=$(sed '$d' <<< "$pi")
        done 
        echo -n "$pi, $t" >> benchmark_output.txt
        echo -n >> benchmark_output.txt
        echo "Runtime, Max"
        echo "$pi"
        echo "Core count: $t"
    done
fi