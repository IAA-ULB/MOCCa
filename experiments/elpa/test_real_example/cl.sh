#!/bin/bash

# # Following the manual "4 Compiling and linking against ELPA" p28 
# # FCFLAGS, CFLAGS, CXXFLAGS, and LDFLAGS.
# # 1. Extend the PKG_CONFIG_PATH environment variable to point to the subfolder lib/pkgconfig
# export PKG_CONFIG_PATH=/data/antwerpen/201/vsc20170/buildingElpa/elpa/_install/lib/pkgconfig
# # 2. Fetch the correct flags for Fortran (FCFLAGS)
# export FCFLAGS=$(pkg-config --variable=fcflags elpa)
# # 4. Fetch the correct linker flags (LDFLAGS)
# export LDFLAGS=$(pkg-config --libs elpa)
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/antwerpen/201/vsc20170/buildingElpa/elpa/_install-intel/lib

# echo '${FCFLAGS}=' ${FCFLAGS}
# echo "LDFLAGS="${LDFLAGS}

make && mpirun -np 4 ./exe