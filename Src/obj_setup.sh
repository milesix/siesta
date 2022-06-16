#!/bin/sh
#
##set -x
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
objdir=$(
cd -P -- "$(pwd)" &&
pwd -P
)
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
testdir=$(dirname $srcdir)/Tests
utildir=$(dirname $srcdir)/Util
pseudodir=$(dirname $srcdir)/Pseudo
topdir=$(dirname $srcdir)
#
destdir=$(pwd)
#
# Copy build.mk and the checker snippet
#
cp -p ${srcdir}/../Config/mk-build/build.mk ${destdir}
cp -p ${srcdir}/../Config/mk-build/check_for_build_mk.mk ${destdir}
#
# Replicate the hierarchy of makefiles
#
(cd $srcdir;
  for i in $(find . -name \[mM\]akefile ) ; do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    filename=$(basename $i)
    sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
    sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${destdir}/$relpath/$filename
  done
)
(cd $utildir;
  for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile) ; do
    relpath=${i%/*}
    mkdir -p ${destdir}/Util/$relpath
    filename=$(basename $i)
    sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
    sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${destdir}/Util/$relpath/$filename
  done
)
(cd $pseudodir;
  for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile) ; do
    relpath=${i%/*}
    mkdir -p ${destdir}/Pseudo/$relpath
    filename=$(basename $i)
    sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
    sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${destdir}/Pseudo/$relpath/$filename
  done
)
# Replicate any .inc files
#
(cd $srcdir;
  for i in $(find . -name '*.inc' ); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp -fp $relpath/*.inc ${destdir}/$relpath
  done
)
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
# (deactivated for now)
echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file into the directory."
echo " *** These files are template arch.make files:"
echo " ***    gfortran.make (for gfortran compiler)"
echo " ***    intel.make (for intel compiler)"
echo " ***    DOCUMENTED-TEMPLATE.make (requires customization)"

