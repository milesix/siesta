README for FDF
==============

FDF stands for Flexible Data Format and is used in High-Performance Computing
for Electronic Structure Calculations.

LibFDF is the official implementation of the FDF Specifications.

The LibFDF package, designed by Alberto Garcia and Jose Soler, was
reimplemented by Raul de la Cruz (Barcelona Supercomputer Center) in 2007
to make it memory-side instead of file-oriented.

The former fully MPI-enabled library has now been turned into a
single-node smaller library with the basic functionality. MPI
operation is delegated to the client code, which just has
to coordinate the broadcast of the (serialized) internal data
structure containing the fdf data. (See doc/README_MPI.md)

We request that the copyright notices in the code be kept if the package
is eventually used as part of another program.

