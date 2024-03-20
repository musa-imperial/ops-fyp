#!/bin/bash

./removescript.sh
./removetxt.sh

export OPS_COMPILER=icx
export OPS_INSTALL_PATH=/home/musa/apps/ops_installations/OPS/ops
export PATH=/home/musa/intel/oneapi/compiler/2024.0/bin:$PATH

source /home/musa/intel/oneapi/compiler/2024.0/env/vars.sh
source /home/musa/intel/oneapi/mpi/2021.11/env/vars.sh


export MPICH_CC=icx 
export MPICH_CXX=icpx

export MPI_INSTALL_PATH=/home/musa/intel/oneapi/mpi/2021.11
export MPICPP=$MPI_INSTALL_PATH/bin/mpicxx

unset HDF5_INSTALL_PATH
export HDF5_INSTALL_PATH=/home/musa/apps/asl/install/build_hdf5/intel
export LD_LIBRARY_PATH=$HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$HDF5_INSTALL_PATH/bin:$PATH

make diffusion_seq
#./diffusion_seq