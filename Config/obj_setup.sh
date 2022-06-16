#!/bin/sh
#
##set -x
#
#
# Absolute path of the MAIN_OBJDIR (compilation directory)
#
objdir=$(
cd -P -- "$(pwd)" &&
pwd -P
)
#
# Get absolute path of this script, as that will be the Config directory to use
# as reference when copying files.
# 
configdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
topdir=$(dirname $configdir)

# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
srcdir=${topdir}/Src
testdir=${topdir}/Tests
utildir=${topdir}/Util
pseudodir=${topdir}/Pseudo
#
# Copy build.mk and the checker snippet
#
cp -p ${configdir}/mk-build/build.mk ${objdir}
cp -p ${configdir}/mk-build/check_for_build_mk.mk ${objdir}
#
# Replicate the hierarchy of makefiles
#
(cd $srcdir;
  for i in $(find . -name \[mM\]akefile ) ; do
    relpath=${i%/*}
    mkdir -p ${objdir}/Src/$relpath
    filename=$(basename $i)
    sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
    sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${objdir}/Src/$relpath/$filename
  done
)
(cd $utildir;
  for i in $(find . -name \[mM\]akefile ) ; do
    relpath=${i%/*}
    mkdir -p ${objdir}/Util/$relpath
    filename=$(basename $i)
    sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
    sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${objdir}/Util/$relpath/$filename
  done
)
(cd $pseudodir;
  for i in $(find . -name \[mM\]akefile ) ; do
    relpath=${i%/*}
    mkdir -p ${objdir}/Pseudo/$relpath
    filename=$(basename $i)
    sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
    sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${objdir}/Pseudo/$relpath/$filename
  done
)
# Replicate any .inc files
#
(cd $srcdir;
  for i in $(find . -name '*.inc' ); do
    relpath=${i%/*}
    mkdir -p ${objdir}/Src/$relpath
    cp -fp $relpath/*.inc ${objdir}/Src/$relpath
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
              | tar -cf - --no-recursion -T- )   | ( cd ${objdir} ; tar xf -)
#
echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file into the directory."
echo " *** Build Siesta in Src."

