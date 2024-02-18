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




HEADERS=2d_compressible_kernels.h derivatives.h fluxx.h global_decl.h \
		   global_ops_vars.h global_params.h

OPS_FILES=2d_compressible.cpp build_datasets.cpp

OPS_GENERATED=2d_compressible_ops.cpp build_datasets_ops.cpp

OTHER_FILES=


APP=2d_compressible
MAIN_SRC=2d_compressible

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app