#!/bin/sh
#
# This is a dummy script to redirect users to the new one
#
echo "***"
echo "***  Some details of the building procedure have changed. ***"
echo "***"
echo "***  Please use the new script: Config/obj_setup.sh       ***"
echo "***"
echo "***  and see the INSTALL file in the top directory        ***"
echo "***"

#
#  Put here the right code to emulate the other script, taking care
#  to use the right paths. This is for backwards compatibility.
#
# Absolute path of the MAIN_OBJDIR (compilation directory)
#
objdir=$(
cd -P -- "$(pwd)" &&
pwd -P
)
# 
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
topdir=$(dirname $srcdir)
#
configdir=${topdir}/Config
testdir=${topdir}/Tests
utildir=${topdir}/Util
pseudodir=${topdir}/Pseudo
#
# Copy build.mk and its subordinate files and the checker snippet
#
cp -p ${configdir}/mk-build/build.mk ${objdir}
cp -p ${configdir}/mk-build/extlibs.mk ${objdir}
cp -p ${configdir}/mk-build/check_for_build_mk.mk ${objdir}
#
# Replicate the hierarchy of makefiles
#
# First the top-level Makefile
#
sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" ${topdir}/Makefile | \
sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${objdir}/Makefile
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
              -path *.mk -prune      -o  \
              -path *.arch-ids  -prune -o -print \
              | tar -cf - --no-recursion -T- )   | ( cd ${objdir} ; tar xf -)
(cd ${testdir};
    for i in *.mk ; do
	sed "s#TOPDIR=\.#TOPDIR=${topdir}#g" $i | \
        sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${objdir}/Tests/$i
    done
)    
 
#
echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file into this directory."
echo " *** Build Siesta in Src here."

