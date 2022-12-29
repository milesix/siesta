.. sectionauthor:: Yann Pouillon <y.pouillon@simuneatomistics.com>

Preparing the environment
=========================

The following actions are essential for the build of SIESTA to succeed. They
have to be performed at least once. Then, depending on the installation method
used, updating the available toolchains and libraries will be either automatic
or manual.


On Linux
--------

Debian-based distributions with GCC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On distributions like Debian, Mint or Ubuntu, the minimum requirements can be
installed system-wide with the *apt* package manager. To enable the default
GCC-based toolchain, it is as simple as typing::

   sudo apt install build-essential gfortran

The *build-essential* package makes sure that the C and C++ compilers from
GCC, as well as GNU Make, are always installed.

If you want to benefit from enhanced performance, you can also install a MPI
distribution. On Debian and its derivatives, the most widely used
implementation is OpenMPI::

   sudo apt install libopenmpi-dev openmpi-bin

We will assume in the following that you have installed OpenMPI. If this is
not the case, simply remove the `-mpi-` packages from the install arguments of
*apt*.

Although SIESTA 4.x comes with an embedded linear algebra implementation, it
is highly recommended to use a native one, e.g. a recent version from Netlib.
With *apt*, this reads::

   sudo apt install libblas-dev liblapack-dev libscalapack-mpi-dev

For platform-independent I/O with HDF5 and NetCDF, the following packages are
necessary::

   sudo apt install libhdf5-dev libhdf5-mpi-dev libnetcdf-dev libnetcdff-dev netcdf-bin

Please note that HDF5 is a very complex set of libraries, in particular when it
comes to configuring and using the correct ones. You can expect a few teething
problems if you are not yet familiar with it.

If you intend to download the source code of SIESTA from internet, we recommend
you to install *wget* or *curl*. Both work very well, choosing one over the
other is mostly a question of personal preference. If you are not sure, you can
install both::

   sudo apt install curl wget

In some cases, you may have to apply patches to some source file before compiling SIESTA. For this, you need difference-related utilities::

   sudo apt install diffutils patch

Although this is not something you will use daily, these tools are still quite
handy.


RPM-based distributions with GCC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On distributions like RedHat, CentOS or Fedora, the minimum requirements can be installed system-wide with the *yum* package manager.

On CentOS, the first step is to install the *epel-release* repository::

   sudo yum install epel-release

It will be useful to get up-to-date libraries. After that, you can install the *Development Tools* set of packages and *gfortran*::

   sudo yum groupinstall "Development Tools"
   sudo yum install gcc-gfortran

In order to benefit from a MPI distribution like OpenMPI, you have to install
the corresponding packages::

   sudo yum install openmpi-devel

Linear algebra libraries can then be installed with::

   sudo yum install blas-devel lapack-devel blacs-openmpi-devel scalapack-openmpi-devel

For platform-independent I/O with HDF5 and NetCDF, you have to install the
following packages::

   sudo yum install hdf5-devel hdf5-openmpi-devel netcdf-devel netcdf-fortran-devel

Please note that HDF5 is a very complex set of libraries, in particular when it
comes to configuring and using the correct ones. You can expect a few teething
problems if you are not yet familiar with it.

If you intend to download the source code of SIESTA from internet, we recommend
you to install *wget* or *curl*. Both work very well, choosing one over the
other is mostly a question of personal preference. If you are not sure, you can
install both::

   sudo yum install curl wget

The RedHat Development Tools come with utilities to patch source files when
necessary. There is thus no need to install any additional package.


Common issues
~~~~~~~~~~~~~

What happens if you cannot become root on the computer where you would like to
use SIESTA? When possible, you can point your system administrator to these
instructions. Otherwise, if you have enough disk space available, you can try
the approach proposed by EasyBuild, as explained in
:doc:`build-easybuild`.


On macOS
--------

Building and using Fortran programs on macOS is not the easiest of things,
since it departs strongly from the usage policy of Apple. However, if you are
ready for a few bumps on the road, in particular when updating the system, it
is possible to develop scientific software on a macOS-based computer.

Please note that the recommended solution to set the build environment varies with time, depending on the
issues encountered on the current releases of the operating systems. In 2021,
it seems that Homebrew is the most satisfactory solution. In any case, there
are steps that you have to perform beforehand, no matter how you install the
rest of the packages.


Prerequisites
~~~~~~~~~~~~~

Before doing anything, you have to check that you have enough free disk space on your computer. You will need 30 to 40 Gb free for everything to work. If you use a Macbook Air, this can be problematic, since you might have to extend the available storage space, which has some limitations (XXX put link here).

Before developing anything, you have to install and activate XCode. Please be aware that it requires a lot of disk space to be installed and work properly. There are several ways to do it. If you follow the instructions from `FreeCodeCamp`_, you'll be ready to go.


Installing packages with Homebrew
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

XXX Help wanted!


Common issues
~~~~~~~~~~~~~

The most frequent issue encountered with macOS is getting failures due to inconsistent permission settings. This can be remedied with system utilities like `Onyx`_, which are beyond the scope of the current document.


On Windows
----------

Coming soon ...


.. _FreeCodeCamp: https://www.freecodecamp.org/news/how-to-download-and-install-xcode/
.. _Onyx: https://titanium-software.fr/en/onyx.html
