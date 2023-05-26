# Installation of Siesta with CMake

Siesta requires CMake >= 3.17, and (if used) the ninja (>=1.10) backend.
Both cmake and ninja can be installed easily in most systems, in
particular with conda or pip.

The CMake approach facilitates the handling of the external libraries
that must be installed before Siesta can be compiled. Depending on the
needs and experience of the user, several modes of operation are
available.

Siesta CMake configurations have heavily borrowed ideas from the
DFTB+ and SIRIUS codes, and we fully acknowledge their contributions
and efforts in streamlining CMake infrastructure.


## Quick and go

The most basic compilation of Siesta can be done simply by:

```shell
cmake -S. -B_build -DCMAKE_INSTALL_PREFIX=/path/to/installation
cmake --build _build -j 4
cmake --install _build
```
If all required dependencies are found this will succeed, otherwise
follow below instructions.


## Compilation flags

Compilation flags are generally managed through the environment
variables (NOT CMake variables).

- `FC` for specifying the fortran compiler
- `FFLAGS` for specifying the compilation flags

An invocation might be:
```shell
FC=gfortran FFLAGS='-O3 -march=native' cmake ...
```

Alternatively, the flags can be supplied on the command line
```shell
cmake -DFortran_FLAGS=-Os -DC_FLAGS=-Os
```
This enables fine tuning of the compiler flags.

> Customarily, CMake uses the `CMAKE_<LANG>_FLAGS`.
> These may still be used, but the shorter, simpler flags
> allows less typing and faster proto-typing.

Siesta's infrastructure also allows the usage of toolchain files.
This can either be set in 2 different ways:
```shell
# Use default toolchain files located in Config/cmake/toolchains
cmake ... -DSIESTA_TOOLCHAIN=gnu
# or a full path (for local edited files)
cmake ... -DSIESTA_TOOLCHAIN=/path/to/toolchain/file/gnu.cmake

# Direct usage of the toolchain file
cmake ... -C Config/cmake/toolchains/gnu.cmake
# or equivalently
cmake ... -DCMAKE_TOOLCHAIN_FILE=Config/cmake/toolchains/gnu.cmake
```
When using `SIESTA_TOOLCHAIN` one can use multiple toolchains.
This can be valuable for overwriting or adding variables from various
toolchains. Mainly useful for developers.
```shell
cmake -DSIESTA_TOOLCHAIN=gnu;local ...
```
to use `./Config/cmake/toolchains/gnu.cmake` and `./local.cmake`.


These toolchain files may be used to default your variables and caching of the flags.

Currently the default toolchain will be decided with:
- GNU compilers will use the `Config/cmake/toolchains/gnu.cmake`
  toolchain file.
- Intel (and the newer Intel LLVM backend) compilers will use the
  `Config/cmake/toolchains/intel.cmake` toolchain file.
- Otherwise a _generic_ toolchain file will be used, which uses the
  default CMake variables.


To gain complete control of the compiler flags (without adding the toolchain
ones) you will have to select the `none` toolchain and set the flags.
```shell
cmake -DSIESTA_TOOLCHAIN=none -DFortran_FLAGS="-Os -Dasheusatoehu"
```

A custom toolchain may contain any setting of variables. They can
be thought of as an `arch.make` file with default parameters.
Parameters that exists in a toolchain file can be overwriting on
the command-line with `cmake -D<VAR>=<VALUE>` for temporary
changing its value.


### Build type

CMake compilation infrastructure utilizes a build-type to determine the
flags used.

These build-types are primarily used for experienced users, the default
build type (`Release`) should be sufficient for most (if not all users).

A specific build-type can be enabled with:
```shell
cmake -DCMAKE_BUILD_TYPE=Debug
```

Currently the default Siesta toolchain files allows these different
build types:

- `Release`: the default and recommended build type, it uses a high optimization
  level without sacrifycing accuracy.
- `Debug`: used for debugging Siesta, or if there are runs that shows problems
  this build-type may be useful.
  *Bug reports* should use this build
- `Check`: used for debug + checking code execution runs, primarily
  useful for developers; equally good for bug-reports.
- `RelWithDebInfo`: a release mode with debug mode.
- `MinSizeRel`: optimizes the executables for minimum size (`-Os`)

One can specify different compiler flags for different build types
to more easily switch between them, for instance:
```shell
cmake -DFortran_FLAGS=-Os -DFortran_FLAGS_DEBUG=-g -DCMAKE_BUILD_TYPE=Debug
```
will use the `Fortran_FLAGS_DEBUG` flags while omitting the `Fortran_FLAGS`.
This allows toolchain files to be self-contained and contain multiple
user-configurations.

The currently supported build-types in the shipped toolchain files are:
- `Fortran_FLAGS`
- `Fortran_FLAGS_RELEASE`
- `Fortran_FLAGS_DEBUG`
- `Fortran_FLAGS_CHECK`
- `Fortran_FLAGS_RELWITHDEBINFO`
- `Fortran_FLAGS_MINSIZEREL`


### Developers

Developers are suggested to create custom toolchain files with the appropriate
compiler flags and linker flags to sustain a quick and easy turn-around for
the compilation procedure.

> Bash scripts are notorious for omitting quotation marks when passing
> variables to CMake.
> For instance a small script like this will fail due to the quotation
> marks being disconnected when passed as arguments to the `cmake` executable
> ```shell
> opts="-DFortran_FLAGS='-Os -g'
> cmake $opts
> ```
> Full control is easier to gain by using custom toolchain files.


## Building in parallel (recommended!)

To build in parallel simply add these flags:
```shell
cmake ... -j 4
```
to build using 4 processes.



## Options

Siesta provides a set of options that controls the capabilities
or some intricate feature of Siesta. The generic Siesta executable
should be sufficient for most, but some may need different details.

- `WITH_GRID_SP=OFF|ON` use single-precision grid operations (`ON`).
  Can greatly reduce the memory requirements for large mesh-cutoffs
  and/or large unit-cells. At the expense of some precision.
  The default is to use double precision `-DWITH_GRID_SP=OFF`



## Dependencies

Siesta heavily relies on numerous dependencies, some are required
while some are optional.

To ease the installation several of the packages are shipped in the
Siesta source tree. These can be checked out by doing:
```shell
git submodule update --init --recursive
```
to fetch all of them.  
If users do not have internet access on the compiling machine one must
send the sources by other means. To aid this procedure one may use
the `stage_submodules.sh` script to gather all sources for later uploading.

Ensure that the required packages are present in these environment variables:

- `CMAKE_PREFIX_PATH` variable
- `PKG_CONFIG_PATH` variable

For instance:
```shell
CMAKE_PREFIX_PATH=/path/libxc/share/cmake:/path/libgridxc/share/cmake
PKG_CONFIG_PATH=/path/libxc/lib/pkgconfig:/path/libgridxc/pkgconfig
cmake ...
```
which ensures that CMake can search in the appropriate directories.
Alternatively one can put `CMAKE_PREFIX_PATH` as a CMake variable:
```shell
cmake ... -DCMAKE_PREFIX_PATH=/path/libxc/share/cmake;/path/libgridxc/share/cmake
```
Note the different delimiters, `:` (Unix OS) vs. `;` (CMake list separator).


Here they are listed together with their options:


#### BLAS (required)

- `BLAS_LIBRARY=<name of library>|NONE` specifies the library name
  for linking. If `NONE` BLAS is implicitly linked through other
  libraries/flags or the compiler itself (e.g. Cray or
  for instance in OpenBLAS LAPACK can be implicitly contained and
  the BLAS library is not needed).
- `BLAS_LIBRARY_DIR=<path to library>` place where to find the library
  `BLAS_LIBRARY`
- `BLAS_LINKER_FLAG` flags to use when linking

Example:
```shell
cmake ... -DBLAS_LIBRARY=blis \
          -DBLAS_LIBRARY_DIR=/opt/blis/lib
```


#### LAPACK (required)

- `LAPACK_LIBRARY=<name of library>|NONE` specifies the library name
  for linking. If `NONE` LAPACK is implicitly linked through other
  libraries/flags or the compiler itself (e.g. Cray).
- `LAPACK_LIBRARY_DIR=<path to library>` place where to find the library
  `LAPACK_LIBRARY`
- `LAPACK_LINKER_FLAG` flags to use when linking

Example:
```shell
cmake ... -DLAPACK_LIBRARY=openblas \
          -DLAPACK_LIBRARY_DIR=/opt/openblas/lib \
	  -DBLAS_LIBRARY=NONE
```


#### ScaLAPACK (required for MPI support)

- `SCALAPACK_LIBRARY=<name of library>|NONE` specifies the library name
  for linking. If `NONE` ScaLAPACK is implicitly linked through other
  libraries/flags or the compiler itself.
- `SCALAPACK_LIBRARY_DIR=<path to library>` place where to find the library
  `SCALAPACK_LIBRARY`
- `SCALAPACK_LINKER_FLAG` flags to use when linking

Example:
```shell
cmake ... -DSCALAPACK_LIBRARY="-lmkl=cluster" \
	  -DBLAS_LIBRARY=NONE -DLAPACK_LIBRARY=NONE
```


#### MPI (highly recommended)

- `WITH_MPI=ON|OFF` to enable, disable support respectively.

MPI will defaulted to be `ON` when an MPI compiler is found.
Using MPI will forcefully require ScaLAPACK (see above).



#### OpenMP

Enable threading support using OpenMP.

Far from all of Siesta is done using OpenMP.
Both TranSiesta/TBtrans may benefict when running large systems
which may result in performance gains plus memory reductions.

- `WITH_OPENMP=OFF|ON` to disable, enable support respectively.
  By default it is `OFF`.

Users are recommended to test whether it makes sense for them to
utilize the threading support.

Be aware of `OMP_NUM_THREADS` and `OMP_PROC_BIND` variables
which may highly influence the performance gains.



#### xmlf90 (required)

Contained in the External/xmlf90 folder. Can be pre-installed,
installed from custom source directory, fetched at compile time.

- `XMLF90_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `XMLF90_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `XMLF90_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `XMLF90_FIND_METHOD=source`
- `XMLF90_GIT_TAG` when `XMLF90_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `XMLF90_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.


#### libfdf (required)

Contained in the External/libfdf folder. Can be pre-installed,
installed from custom source directory, fetched at compile time.

- `LIBFDF_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `LIBFDF_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `LIBFDF_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `LIBFDF_FIND_METHOD=source`
- `LIBFDF_GIT_TAG` when `LIBFDF_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `LIBFDF_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.


#### libpsml (required)

libpsml enables the reading of pseudopotential files in the `PSML` file
format. It thus enables re-use of PSML files from www.pseudo-dojo.org,
amongst others.

- `LIBPSML_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `LIBPSML_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `LIBPSML_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `LIBPSML_FIND_METHOD=source`
- `LIBPSML_GIT_TAG` when `LIBPSML_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `LIBPSML_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.


#### libgridxc (required)

libgridxc enables the calculation of the XC functionals on the grid
where density is calculated. It can leverage the libxc library, which
is also highly recommended.

libgridxc depends on the `WITH_GRID_SP` flag which controls the
precision of the grid operations.


<to be filled>



#### libxc (highly recommended)

libxc is an XC functional library implementing a very large variety of functionals.

libgridxc can leverage libxc and use the functionals from there.

- `WITH_LIBXC=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if the library can be found.
- `LIBXC_ROOT=<path to installation>` where libxc has been installed.
- `LIBXC_Fortran_INTERFACE=f03;f90` to search for a specific interface,
  defaults to both, but prefers the f03 interface. To only search for
  f90, do `-DLIBXC_Fortran_INTERFACE=f90`.



#### simple-DFTD3 (optional)

Add support for DFTD3 dispersion corrections as suggested by Grimme et.al.

- `WITH_DFTD3=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if the `./External/DFTD3/` directory contains
  directories with the appropriate sources.
- `S-DFTD3_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `S-DFTD3_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `S-DFTD3_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `S-DFTD3_FIND_METHOD=source`
- `S-DFTD3_GIT_TAG` when `S-DFTD3_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `S-DFTD3_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.


##### mctc, mstore, test-drive, toml-f

These packages are dependencies of the `simple-DFTD3` library.
Generally one need not change these unless one changes the `S-DFTD3_*`
flags in which case dependencies may require manual changes.

Here are flags for each of these sub-dependencies.


- `MCTC-LIB_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `MCTC-LIB_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `MCTC-LIB_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `MCTC-LIB_FIND_METHOD=source`
- `MCTC-LIB_GIT_TAG` when `MCTC-LIB_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `MCTC-LIB_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.




- `MSTORE_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `MSTORE_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `MSTORE_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `MSTORE_FIND_METHOD=source`
- `MSTORE_GIT_TAG` when `MSTORE_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `MSTORE_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.




- `TEST-DRIVE_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `TEST-DRIVE_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `TEST-DRIVE_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `TEST-DRIVE_FIND_METHOD=source`
- `TEST-DRIVE_GIT_TAG` when `TEST-DRIVE_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `TEST-DRIVE_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.




- `TOML-F_FIND_METHOD=cmake/pkgconf/fetch/source`
  a CMake list of multiple ways to check for library existance,
  `cmake;fetch` will first search using CMake `find_package`, if
  that fails it will fetch it from the `TOML-F_GIT_REPOSITORY`
  variables
  `cmake` and `pkgconf` are generically implemented using
  package finders shipped with CMake. Ensure `CMAKE_PREFIX_PATH`
  and `PKG_CONFIG_PATH` are (prepended)/appended with the
  directories that should be searched.
  This value defaults to `SIESTA_FIND_METHOD`.
  `SIESTA_FIND_METHOD` defaults to `cmake;pkgconf;source;fetch`.
- `TOML-F_SOURCE_DIR` should point to a directory where the
  sources are present, ither manually cloned on unpacked from
  a release archive.
  Applicable when `TOML-F_FIND_METHOD=source`
- `TOML-F_GIT_TAG` when `TOML-F_FIND_METHOD=fetch` this
  revision of the source will be checked out.
- `TOML-F_GIT_REPOSITORY` is the URL of the Git repository
  when cloning the sources.
  Is defaulted to the original development site, may be
  useful for testing clones with fixes/changes or.





#### NetCDF (highly recommended)

Enable writing NetCDF files for faster (parallel) IO and
also for easier post-processing utilities.


- `WITH_NETCDF=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if NetCDF can be found (i.e. specifying `NetCDF_PATH`
  should be enough).
- `NetCDF_ROOT|PATH=<path to installation>` to specificy the location
  of the NetCDF installation. Generally should be enough with this
  flag.
- `NetCDF_INCLUDE_DIR` to manually specify include directories
  for modules etc.


In conjunction with NetCDF there are supporter libraries shipped
with Siesta which are required for TBtrans.

- `WITH_NCDF=ON|OFF` enable NCDF support (a wrapper around NetCDF
  that makes it easier to work with).
  This is automatically detected.
  The default is sufficient.

- `WITH_NCDF_PARALLEL=ON|OFF` allow parallel IO through NCDF library.
  This is automatically detected.
  The default is sufficient.



#### ELPA (recommended)

See `Config/cmake/Modules/FindCustomElpa.cmake` for details on
how to link against ELPA.



#### FFTW

The FFTW library is only used in the Util/STM/ol-stm utility.

- `WITH_FFTW=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if found.



#### FLOOK

A library for interacting with the internal Siesta variables
on the fly and/or create custom molecular dynamics trajectories.
It exposes the Lua language for scripting capabilities.

- `WITH_FLOOK=ON|OFF` to enable, disable support respectively.
  Defaults to `ON` if found.



## Tests

CMake integrates a testing framework.

Not all tests are present, this is a work-in-progress.
```shell
ctest <options>
```
If the required external libraries have been compiled as part of
the current CMake invocation, installation tests for them will
also be executed.





_Experimental_) SPACK packages are available in Config/spack_package_defs

   After setting up spack and defining compilers, etc, a user can simply install
   a new repo with the 'siesta-project' namespace:

     spack repo add /path/to/spack_package_defs

   and issue commands such as:

     spack install xmlf90
     spack info siesta
     spack spec siesta +mpi +netcdf +libxc +elpa
     spack install siesta -mpi build_type=Debug

   Note that the spack builtin repo *might have* other Siesta-related
   recipes prepared in the past by other members of the community.
   By making sure that the 'siesta-project' repo is listed first in the
   spack repository chain, those can be avoided. Check:

     spack repo list
     
   At this point, the spack recipes pull sources from development
   branches of Siesta, without a well-defined source id. This is
   a temporary situation during the final stages of development of
   the CMake framework and the update of the dependency libraries.
 
