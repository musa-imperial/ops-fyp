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




HEADERS=avg_kernel.h derivatives.h fluxx.h global_decl.h \
		   global_ops_vars.h global_params.h initial_conditions_kernel.h \
		   avg_kernel.h 


OPS_FILES=2d_compressible.cpp initial_conditions.cpp avg.cpp

OPS_GENERATED=2d_compressible_ops.cpp initial_conditions_ops.cpp \
				avg_ops.cpp

OTHER_FILES=build_datasets.cpp #initial_conditions.cpp avg.cpp


APP=2d_compressible
MAIN_SRC=2d_compressible

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app