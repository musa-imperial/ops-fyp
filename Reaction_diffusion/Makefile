#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
#include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
#include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hip
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5




HEADERS=reaction_diffusion_kernels.h

OPS_FILES=reaction_diffusion_parallel.cpp

OPS_GENERATED=reaction_diffusion_parallel_ops.cpp

OTHER_FILES=


APP=reaction_diffusion_parallel
MAIN_SRC=reaction_diffusion_parallel

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app