#!/bin/sh
#
# Creates the 'serial' distribution of siesta. Needs python.
#
#  This script should be run on an "EXPORTED" directory:
#
#  cvs export { -r <appropriate tag>, -d <date> } [ -d directory ] siesta
#

dir=siesta-1.1-f1s

if [ -d $dir ] ; then
      echo "Directory exists"
      exit
fi

mkdir $dir
echo "Copying Pseudo Examples Docs Util..."
cp -rp Pseudo Examples Docs Util $dir
cp -rp README send-bug $dir
rm -rf $dir/Docs/Tech

mkdir $dir/Src
cd Src
echo "Copying fdf Sys NetCDF Tests and Libs to Src..."
cp -rp fdf Sys NetCDF Libs Tests ../$dir/Src
cp -rp Makefile ../$dir/Src

echo "De-MPI'ng .F and .F90 files..."
for i in *.F *.F90; do
  python ../dempi.py $i > ../$dir/Src/$i
done
echo "Copying serial files..."
for i in *.f *.f90 ; do
  cp -p $i ../$dir/Src
done
rm -rf $dir/Src/Include

echo ""
echo "**** Remember to include the .ps.gz User guide in 2up form"

