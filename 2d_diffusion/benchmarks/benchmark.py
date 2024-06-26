import subprocess, os


max_thread_count=16

compiler = "icx"

if compiler == "icx":

        command = """
export OPS_COMPILER=icx
export OPS_INSTALL_PATH=/home/musa/apps/OPS/ops
export PATH=/home/musa/intel/oneapi/compiler/2024.0/bin:$PATH

# source /home/musa/intel/oneapi/compiler/2024.0/env/vars.sh
# source /home/musa/intel/oneapi/mpi/2021.11/env/vars.sh

export MPICH_CC=icx 
export MPICH_CXX=icpx

export MPI_INSTALL_PATH=/home/musa/intel/oneapi/mpi/2021.11
export MPICPP=$MPI_INSTALL_PATH/bin/mpicxx

unset HDF5_INSTALL_PATH
export HDF5_INSTALL_PATH=/home/musa/apps/asl/install/build_hdf5/intel
export LD_LIBRARY_PATH=$HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH
export PATH=$HDF5_INSTALL_PATH/bin:$PATH
"""
       

elif compiler == "gnu":
        command = """
export OPS_COMPILER=gnu
export OPS_INSTALL_PATH=/home/musa/apps/OPS/ops
export HDF5_INSTALL_PATH=$HOME/apps/hdf5
"""

subprocess.run('./removescript.sh', shell=True)

makeoutput = subprocess.run(command+'make diffusion_seq', shell=True, stdout=subprocess.PIPE)
#print(makeoutput)
with open('benchmark_output.txt', 'w') as f:
        f.write("Runetime, Thread count")  
        for t in range(1, max_thread_count+1):
                subprocess.run(f'export OMP_NUM_THREADS={t}', shell=True)
                #pi = subprocess.run(['./quick_run.sh'], shell=True,  stdout=subprocess.PIPE, text=True)
                pi = subprocess.run(command+'./diffusion_seq', shell=True,  stdout=subprocess.PIPE, text=True, env={'OMP_NUM_THREADS': f'{t}'})
                f.write(pi.stdout)
                f.write(f", {t}")
                print("Runtime, "+pi.stdout)
                print(f'Thread count: {t}')

