This MR implements three changes:

## Re-design of the matel tables, including a new approach to their initialization

(By Rogeli Grima and Alberto Garcia)

Previously, the matel tables were initialized as each element was
first needed.  In particular, most of the matrix elements were
computed first in the nlefsm routine, which unfortunately had a large
degree of redundancy when running in parallel.

The matrix elements are now pre-computed in parallel, followed by a
globalization of the data among all the MPI processes.The matel tables
are initialized as needed: first the main tables, and the rest, if
necessary, at the point just before they are needed (currently calls
to 'optical', 'born_charge', 'ksv', or 'wannier').

The data structures have been redesigned. The MATEL type tries to be a
replica of the data in the old module, but reducing the dimensionality
associated with the various operations ('S', 'T', 'U', ...). Eachh
operation is now handled by as many instances ("tables" with prefix
tab_) of the MATEL type as are necessary. The MATEL type has
associated a series of procedures:

     * INIT (IOPER .. ): Creates and initializes an IOPER type matrix

       
     * GET_MATEL (ig1, ig2, ...): Extract the element associated with the
                                  functions tagged with ig1 and ig2 in the
				  matel registry.
				  
    
The rest of MATEL procedures are used to build parts of the matrix,
communicate those parts between the various processes or to index the
elements of said matrix.

A new wrapper routine matel_mod::new_matel holds the logic to dispatch
the right evaluators, while maintaining the old calling sequence
and code.

The code now runs much faster, but it can use more memory to hold the
tables, as data vectors which are proportional to existing ones (which
appear due to the intrinsic symmetries of some matrix elements) are
only removed if detected by a given MPI rank. A further step of
collective filtering will be needed to fix this.

Also needed: Clarify the conditions to initialize the "optical",
"polarization", and "wannier" tables, if the analysis stage is
requested by means other than the fdf file (e.g., a Lua script).


##  Consolidate parameters in spher_harm.f and avoid allocations
    
    Define MAXL=10 at the top level. Make Y(:) and dydr(:,:) static
    arrays dimensioned to the maximum.


##  Exploit more opportunities for parallelization over atoms

Some routines (kinefsm, dnafsm, etc) perform loops over atoms. These
loops can be distributed among different MPI ranks, and a final "reduce"
operation performed at the end.

A new module m_mpi_inplace provides interfaces for the "in place"
version of mpi_all_reduce, which was implicitly used by Rogeli Grima
in his work. The calls might have been re-worked, but the new module
will be useful in the future to eliminate the need for extra variables
and for extra copy operations.
    
As the 'in-place' option is not supported by the 'classic' MPI
interface contained in mpi_siesta, the new module does not use the
latter. As side effects, the 'comm' argument is now mandatory, and no
'MPI' timing information is provided.

## Things to do in this branch before merging:

* Add code to show the contents of the ylmk and matel tables.
* Add code to (re)-initialize the tables (currently the data
  is not released at the end of the program).
* Add a 'name' field to the tables.
* Adhere to the coding guidelines in Docs/developer/...
* Make sure that possible Lua-driven cases are treated correctly in
  the initialization.
* Perform a collective filtering of redundant data.
* Maybe eliminate the vna special case in the evaluation routine in
  matel_registry (factor of Y_00=1/sqrt(4*pi)).
* Maybe use "object-oriented" techniques in matel_registry to handle
  radfuncs and trial_orbitals.
* Add FORD documentation for the new suite of modules.
* Add "unit tests"

