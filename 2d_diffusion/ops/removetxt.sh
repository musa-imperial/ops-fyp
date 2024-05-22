rm -rf CUDA
rm -rf HIP
rm -rf MPI_inline
rm -rf MPI_OpenMP
rm -rf OpenACC
rm -rf OpenCL
rm -rf OpenMP_offload
rm -rf SYCL
rm -rf .generated
rm -rf *_ops.cpp
rm -rf *_ops.c
rm -rf *.txt
rm -rf *.txt.*

#The grep -v says "only allow filenames that don't contain a dot"
#the xargs rm says "then pass the list of filenames to rm".
#ls | grep -v "\." | xargs rm #removes execuable
#ls | grep -v "\." | grep -v Makefile | xargs rm

