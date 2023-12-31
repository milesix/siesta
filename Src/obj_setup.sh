#!/bin/sh
#
# ---
# Copyright (C) 1996-2014	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
##set -x
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
user_specified_dir=$(dirname $0)
testdir=$(dirname $srcdir)/Tests
#
destdir=$(pwd)
#
# Replicate the hierarchy of makefiles
#
(cd $srcdir;
  for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp $relpath/*akefile ${destdir}/$relpath
  done
)
# Replicate any .h files
# This is needed in some systems with broken include file import heuristics
# (e.g., CSCS blanc)
#
(cd $srcdir;
  for i in $(find . -name '*.h' ); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp -f $relpath/*.h ${destdir}/$relpath
  done
)
#
sed "s#VPATH=\.#VPATH=${srcdir}#g" ${srcdir}/Makefile > ${destdir}/Makefile

#
# Tests directory
# Create a list of files and use tar to process the list and copy the files
# to the destination directory
#
( cd ${testdir} ; cd .. ; find Tests  \
              -path *Reference -prune -o  \
              -path *Reference-xml -prune -o  \
              -path *work -prune      -o  \
              -path *.arch-ids  -prune -o -print \
              | tar -cf - --no-recursion -T- )   | ( cd ${destdir} ; tar xf -)
#
echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file or run configure as:"
echo "    ${user_specified_dir}/configure [configure_options]"
