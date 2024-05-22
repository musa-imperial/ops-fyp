#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
#include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hip
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5




HEADERS=diffusion_kernels.h

OPS_FILES=diffusion.cpp

OPS_GENERATED=diffusion_ops.cpp

OTHER_FILES=


APP=diffusion
MAIN_SRC=diffusion

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app_legacy
