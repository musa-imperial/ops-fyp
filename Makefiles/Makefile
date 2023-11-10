#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
#include $(OPS_INSTALL_PATH)/../makefiles/Makefile.cuda
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hip
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5




HEADERS=cw_kernels.h

OPS_FILES=coursework_ops.cpp

OPS_GENERATED=coursework_ops_ops.cpp

OTHER_FILES=


APP=coursework_ops
MAIN_SRC=coursework_ops

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app
