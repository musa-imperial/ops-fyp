rm -rf cuda
rm -rf hip
rm -rf mpi_inline
rm -rf mpi_openmp
rm -rf openacc
rm -rf opencl
rm -rf openmp_offload
rm -rf sycl
rm -rf .generated
rm -rf *_ops.cpp
rm -rf *.txt
rm -rf *.txt.*

#The grep -v says "only allow filenames that don't contain a dot"
#the xargs rm says "then pass the list of filenames to rm".
#ls | grep -v "\." | xargs rm #removes execuable
ls | grep -v "\." | grep -v Makefile | xargs rm

