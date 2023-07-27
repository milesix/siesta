#!/bin/bash

# Set the paths and environment variables specific to your local system
mkdir -p siesta-max1-berkeley_install
INSTALL_DIR=$(realpath ../siesta-max1-berkeley_install)
export NETCDF_ROOT="/usr"
export NETCDF_INCFLAGS="-I/usr/include"

# ScaLAPACK
export SCALAPACK_ROOT=/usr
export SCALAPACK_LIB=${SCALAPACK_ROOT}/lib/x86_64-linux-gnu/
# export SCALAPACK_INC=${SCALAPACK_ROOT}/include

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCALAPACK_LIB

LDFLAGS = -L$(SCALAPACK_LIB) -lscalapack

export PEXSI_ROOT=$(realpath ../pexsi_install/v2.1_gnu)

HDIR=$(realpath .)
BDIR=$HDIR/build
# /usr/lib/x86_64-linux-gnu/libscalapack-mpich.so

# rm -rf ${BDIR}/*
# Configure CMake with the desired options

cmake --trace-expand -H${HDIR} -B${BDIR} -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DSCALAPACK_LIBRARY=${SCALAPACK_LIB}/libscalapack-mpich.so -DWITH_PEXSI=TRUE -DWITH_ELSI=OFF ..

cd ${BDIR}
make -j
make install
cd ..

