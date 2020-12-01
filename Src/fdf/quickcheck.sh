#!/bin/bash
#
# Copyright (C) 2017 Yann Pouillon <devops@materialsevolution.es>

# Note: this script is for maintainers and testers working with GCC

# Stop at first error and echo commands
set -ev

# Check that we are in the correct directory
test -s "configure.ac" -a -s "src/fdf.F90" || exit 0

# Init build parameters
export DBGFLAGS="-O0 -g3 -ggdb -Wall -Wextra -fbounds-check -fno-inline"

# Prepare source tree
./wipeout.sh
./autogen.sh

# Check default build
mkdir tmp-minimal
cd tmp-minimal
../configure \
  CC="gcc" CFLAGS="${DBGFLAGS}" FC="gfortran" FCFLAGS="${DBGFLAGS}"
sleep 3
make dist
make
make check
mkdir install-minimal
make install DESTDIR="${PWD}/install-minimal"
ls -lR install-minimal >install-minimal.log
cd ..

# Make distcheck
mkdir tmp-distcheck
cd tmp-distcheck
../configure \
  CC="gcc" CFLAGS="${DBGFLAGS}" FC="gfortran" FCFLAGS="${DBGFLAGS}"
sleep 3
make distcheck -j4
make distcleancheck

# Clean-up the mess
cd ..
rm -rf tmp-minimal tmp-distcheck
