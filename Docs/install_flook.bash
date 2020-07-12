#!/bin/bash

# Installation script for zlib, hdf5, netcdf-c and netcdf-fortran
# with complete CDF-4 support (in serial).
# This installation script has been written by:
#  Nick R. Papior, 2016-2017.
#
# (Some comments by Alberto Garcia follow)  ----begin----
#
# The version string has been changed to 0.8.1
#
# To use with the Intel Fortran compiler, a manual step of creation
# of a setup.make file in the top level of the created flook-0.8.1 directory
# might be needed. A minimal example:
#
# CC = icc
# FC = ifort
# CFLAGS = -g
# FFLAGS = -g
# .f90.o:
# 	$(FC) -c $(FFLAGS) $(INC) $<
# .F90.o:
# 	$(FC) -c $(FFLAGS) $(INC) $<
# .c.o:
# 	$(CC) -c $(CFLAGS) $(INC) $<
#
# Make sure to type 'make clean' before typing 'make' again.
# You then need to install the libraries manually.
#
# ---- end ----
#
# The author takes no responsibility of damage done to your hardware or
# software. It is up to YOU that the script executes the correct commands.
#
# This script is released under the LGPL license.

# VERY BASIC installation script of required libraries
# for installing these packages:
#   flook-0.8.1
# If you want to change your compiler version you should define the
# global variables that are used for the configure scripts to grab the
# compiler, they should be CC and FC. Also if you want to compile with
# different flags you should export those variables; CFLAGS, FFLAGS.

# If you have downloaded other versions edit these version strings
f_v=0.8.1

# Install path, change accordingly
# You can change this variable to control the installation path
# If you want the installation path to be a "packages" folder in
# your home directory, change to this:
# ID=$HOME/packages
if [ -z $PREFIX ]; then
    ID=$(pwd)/build
else
    ID=$PREFIX
fi

echo "Installing libraries in folder: $ID"
mkdir -p $ID

# First we check that the user have downloaded the files
function file_exists {
    if [ ! -e $(pwd)/$1 ]; then
	echo "I could not find file $1..."
	echo "Please download the file and place it in this folder:"
	echo " $(pwd)"
	exit 1
    fi
}

# Check for function $?
function retval {
    local ret=$1
    local info="$2"
    shift 2
    if [ $ret -ne 0 ]; then
	echo "Error: $ret"
	echo "$info"
	exit 1
    fi
}

file_exists flook-${f_v}.tar.gz
unset file_exists

#################
# Install flook #
#################
[ -d $ID/flook/${f_v}/lib64 ] && flook_lib=lib64 || flook_lib=lib
if [ ! -d $ID/flook/${f_v}/$flook_lib ]; then
    tar xfz flook-${f_v}.tar.gz
    cd flook-${f_v}
    make liball
    retval $? "flook make liball"
    make install PREFIX=$ID/flook/${f_v}
    retval $? "flook make install"
    cd ../
    rm -rf flook-${f_v}
    echo "Completed installing flook"
    [ -d $ID/flook/${f_v}/lib64 ] && flook_lib=lib64 || flook_lib=lib
else
    echo "flook directory already found."
fi

##########################
# Completed installation #
##########################

echo ""
echo "##########################"
echo "# Completed installation #"
echo "#    of flook package    #"
echo "#  and its dependencies  #"
echo "##########################"
echo ""
echo ""

echo "Please add the following to the BOTTOM of your arch.make file"
echo ""
echo "INCFLAGS += -I$ID/flook/${f_v}/include"
echo "LDFLAGS += -L$ID/flook/${f_v}/$flook_lib -Wl,-rpath,$ID/flook/${f_v}/$flook_lib"
echo "LIBS += -lflookall -ldl"
echo "COMP_LIBS += libfdict.a"
echo "FPPFLAGS += -DSIESTA__FLOOK"
echo ""
