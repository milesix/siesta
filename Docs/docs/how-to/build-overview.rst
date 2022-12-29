.. sectionauthor:: Yann Pouillon <y.pouillon@simuneatomistics.com>

Overview of building SIESTA
===========================

Overall process
---------------

There are many ways to build SIESTA. What they have in common is the steps to
follow to get there:

#. Check the minimum requirements (once per computer).
#. Prepare the build environment (once or twice a year).
#. Build SIESTA and its utilities.
#. Test the built executables.
#. Install the executables and possible libraries.
#. Run the executables in the correct environment.

The content of each step will of course depend a lot on your hardware,
operating system, the installed software dependencies, the selected build
framework/methodology, as well as the purpose for which the executables will
be used. This will even influence your choice of the version of SIESTA to
build.

Here are the options available to you when building SIESTA:

- Hardware: `ARM`_, `Intel`_, `PowerPC`_, ...
- Operating system: `Linux`_, `macOS`_, `Windows`_.
- MPI distributions: `Intel MPI`_, `MPICH`_, `OpenMPI`_, ...
- Linear algebra libraries: `Netlib`_, `MKL`_, ...
- Platform-independent I/O: `HDF5`_ + `NetCDF`_.
- Methodologies: manual, `Docker`_, `EasyBuild`_, `ESL Bundle`_, 
  `Singularity`_, `Spack`_, ...
- Desired enhancements: `OpenMP`_, `GPU`_, linear scaling, ...
- Will the executables used system-wide or just by you?
- Do you have administrator privileges on your computer?

Depending on the ingredients you select, the simplicity or difficulty of the
build may vary a lot. In particular, if the above names do not mean anything
to you, we recommend that you have a look at the corresponding websites before
starting.


.. _ARM: https://www.arm.com/
.. _Intel: https://www.intel.com/
.. _PowerPC: https://en.wikipedia.org/wiki/Power_ISA

.. _Linux: https://www.linuxfoundation.org/
.. _macOS: https://www.apple.com/macos/
.. _Windows: https://www.microsoft.com/en-us/windows

.. _`Intel MPI`: https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/mpi-library.html
.. _MPICH: https://www.mpich.org/
.. _OpenMPI: https://www.open-mpi.org/

.. _Netlib: https://www.netlib.org/
.. _MKL: https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5
.. _NetCDF: https://www.unidata.ucar.edu/software/netcdf/

.. _Docker: https://www.docker.com/
.. _EasyBuild: https://easybuild.io/
.. _`ESL Bundle`: https://esl.cecam.org/bundle/
.. _Singularity: https://sylabs.io/singularity/
.. _Spack: https://spack.io/

.. _OpenMP: https://www.openmp.org/
.. _GPU: https://en.wikipedia.org/wiki/Graphics_processing_unit
