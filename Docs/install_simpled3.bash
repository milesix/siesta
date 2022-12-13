#!/bin/bash

# Installation script for simple-dftd3 library.
# Based on the installation script for psml.
#
# The author takes no responsibility of damage done to your hardware or
# software. It is up to YOU that the script executes the correct commands.
#
# This script is released under the LGPL license.

# VERY BASIC installation script of required libraries
# If you want to change your compiler version you should define the
# global variables that are used for the configure scripts to grab the
# compiler, they should be CC and FC. Also if you want to compile with
# different flags you should export those variables; CFLAGS, FFLAGS.

# If you have downloaded other versions edit these version strings
simpled3_v=0.6.0

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
# Decide whether everything is installed in 1 directory
_single_dir=1
_make_j=2

while [ $# -gt 0 ]; do
    opt=$1 ; shift
    case $opt in
	--prefix|-p)
	    ID=$1 ; shift
	    ;;
	--d3-version|-d3-v)
	    simpled3_v=$1 ; shift
	    ;;
	--single-directory)
	    _single_dir=1
	    ;;
	--separate-directory)
	    _single_dir=0
	    ;;
	--make-j)
	    _make_j=$1 ; shift
	    ;;
	--help|-h)
	    echo " $0 --help shows this message"
            echo " simple-d3 installation requires both CMake and Ninja."
	    echo ""
	    echo "These options are available:"
	    echo ""
	    echo "  --prefix|-p <>: specify the installation directory of the libraries"
	    echo "  --d3-version|-d3-v <>: specify the simple-d3 version (default: $simpled3_v)"
	    echo "  --single-directory : all libraries are installed in --prefix/{bin,lib,include} (default: YES)"
	    echo "  --separate-directory : all libraries are installed in --prefix/<package>/<version>/{bin,lib,include} (default: NO)"
	    echo "  --make-j <>: run make in parallel using <> number of cores (default: $_make_j)"
	    echo ""
	    echo "To customize compilers and flags please export these environment variables:"
	    echo "  CC"
	    echo "  FC"
	    echo "  CFLAGS"
	    echo "  FFLAGS"
	    echo ""
	    exit 0
	    ;;
    esac
done


echo "Installing libraries in folder: $ID"
echo ""
echo "Using these software version:"
echo "  simple-d3 : $simpled3_v"
echo ""
echo "If you want other versions, please try $0 --help and check the options"
echo "  Will wait 2 sec before proceeding"
echo ""
sleep 2

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

# Download a file, if able and the file does not exist
which wget > /dev/null
if [ $? -eq 0 ]; then
    # success we can download using wget
    function _dwn_file {
	wget -O $1 $2
    }
else
    function _dwn_file {
	curl -o $1 $2
    }
fi

# Use download function
#  $1 is name of file
#  $2 is URL
function download_file {
    if [ ! -e $(pwd)/$1 ] ; then
	# Try and download
	_dwn_file $1 $2
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

# Download files if they can
download_file s-dftd3-${simpled3_v}.tar.xz https://github.com/awvwgk/simple-dftd3/releases/download/v${simpled3_v}/s-dftd3-${simpled3_v}-source.tar.xz
file_exists s-dftd3-${simpled3_v}.tar.xz
unset file_exists

########################
# Install simple-dftd3 #
########################
if [ $_single_dir -eq 0 ]; then
    sd3_dir=$ID/simple-d3/$simpled3_v
else
    sd3_dir=$ID
fi
[ -d $sd3_dir/lib64 ] && sd3_lib=lib64 || sd3_lib=lib
if [ ! -e $sd3_dir/$sd3_lib/libs-dftd3.a ]; then
    rm -rf s-dftd3-${simpled3_v}
    tar xvf s-dftd3-${simpled3_v}.tar.xz
    d=s-dftd3-${simpled3_v}
    cd $d

    cmake -B _build -G Ninja -DWITH_OpenMP=0 -DCMAKE_INSTALL_PREFIX=${sd3_dir}
    retval $? "simple-d3 cmake setup"
    cmake --build _build/
    retval $? "simple-d3 cmake build"
    pushd _build && ctest && popd
    retval $? "simple-d3 tests"
    cmake --install _build/
    retval $? "simple-d3 make install"
    cd ../
    rm -rf $d
    echo "Completed installing simple-dftd3"
    [ -d $sd3_dir/lib64 ] && sd3_lib=lib64 || sd3_lib=lib
else
    echo "simple-d3 directory already found."
fi

##########################
# Completed installation #
##########################

echo ""
echo "##########################"
echo "# Completed installation #"
echo "#    of simple-dftd3     #"
echo "#  and its dependencies  #"
echo "##########################"
echo ""
echo ""

echo "Please add the following to the BOTTOM of your arch.make file:"
echo ""
echo "DFTD3_ROOT = $sd3_dir"
echo ""
