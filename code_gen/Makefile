#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.common
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.mpi
include $(OPS_INSTALL_PATH)/../makefiles/Makefile.hdf5




HEADERS=reaction_diffusion_kernels.h

OPS_FILES=reaction_diffusion.cpp

OPS_GENERATED=reaction_diffusion_ops.cpp

OTHER_FILES=


APP=reaction_diffusion
MAIN_SRC=reaction_diffusion

include $(OPS_INSTALL_PATH)/../makefiles/Makefile.c_app_legacy