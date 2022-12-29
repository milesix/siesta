:sequential_nav: next

..  _tutorial-basic-first-encounter:

A first encounter with Siesta
=============================

..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.

In this exercise we will get a first acquaintance with ``SIESTA`` by
studying two simple molecules, CH4 and CH3. We will cover quite a lot
of features and concepts, without worrying too much about issues of
convergence.

Basic execution. Input and output files
---------------------------------------


.. hint::
   Enter the directory 'CH4'

You will find an input file named `ch4.fdf` along with the files
`C.psf` and `H.psf` containing the information about the
pseudopotentials.  The ch4.fdf file sets the value of several
parameters which specify both the system we want to study and the
accuracy of the calculation. We find first the inputs that specify the
system::

 SystemName          CH4 molecule
 SystemLabel         ch4

 %block ChemicalSpeciesLabel
  1  6 C   # Species index, atomic number, species label
  2  1 H   # Species index, atomic number, species label
 %endblock ChemicalSpeciesLabel

 AtomicCoordinatesFormat  Ang

 %block AtomicCoordinatesAndAtomicSpecies
   0.000      0.000      0.000   1
   1.219     -0.284     -0.377   2 
  -0.284      1.219     -0.377   2 
  -0.140     -0.140      1.219   2 
  -0.833     -0.833     -0.503   2 
 %endblock AtomicCoordinatesAndAtomicSpecies

Pay special attention to the block ChemicalSpeciesLabel. In
this block you assign an index and a label to each atomic species. The
label will allow to recognize the files containing the information
about the pseudopotential and the basis set (when provided).
 
Check the input of the coordinates (they are just some guess
coordinates, not the optimized equilibrium ones).

The file ch4.fdf contains also the most important parameters to take into
account to perform a molecular calculation. Namely:
 
* Those defining the size and localization of the basis set. Here the
  number of orbitals per atom is defined by the parameter
  ``PAO.BasisSize``, and we have set it to select a minimal basis (SZ)
  for quick, cheap calculation, looking for qualitative results,
  rather than quantitative results.


* The parameter ``MeshCutoff``, controlling the fineness of the real-space
  grid used to compute the integrals for the matrix elements of the
  Hamiltonian.

* Those that control the self-consistent cycle.
 
While the parameters specifying the system are mandatory, all other
parameters have some default values and, in principle, it is not
necessary to explicitly include them in the input file. However, it is
important to note that the default values do not always guarantee a
converged calculation.

Run the program::

   siesta < ch4.fdf > ch4.out           # traditional way

   or simply:
   
   siesta ch4.fdf > ch4.out            
   
Take a look at `ch4.out` once the program has finished. 
The file contains a human-readable account of the doings of SIESTA for
this calculation:

* A header with the version and compilation options, and
  a copy (in the first mode above) or a mention (in the second mode)
  of the input file. 

* Details about the pseudopotential read and the basis set generated
  for each species.

* A summary of the values used in the calculation for the
  most important parameters.

* A log of the self-consistent-field (SCF) cycle

* A summary of the total energy decomposition, forces, and stress

.. hint::
   You can play with the values of a few parameters and check their
   effect on the output results:

   * The basis-set size: Set PAO.Basis-size to any of: DZ, SZP, DZP,
     TZP. More on this topic in :ref:`this tutorial<tutorial-basic-basis-sets>`. 

   * The fineness of the real-space grid: Use a unreasonable low value
     for the the parameter MeshCutoff (may be 10-30 Ry) and check the
     resulting total energy and forces (you can also find the forces in
     the file `ch4.FA`). Try to determine the minimum value of the
     MeshCutoff parameter that gives an energy converged to 0.1 eV.
     More on this topic in :ref:`this tutorial<tutorial-basic-grid-convergence>`. 

Periodic boundary conditions
----------------------------

You might have wondered about the appearance of this block in the
input file::

  #Unit cell for the calculation
  LatticeConstant 15 Ang
  %block LatticeVectors 
  1.000 0.000 0.000
  0.000 1.000 0.000
  0.000 0.000 1.000
  %endblock LatticeVectors

which does not seem to make sense for a 'molecule' calculation. In
fact, SIESTA uses periodic boundary conditions (PBC), and this means
in this case that we are doing a calculation for an infinite
collection of regularly spaced molecules. If we want to simulate an
isolated molecule it is important to have enough distance between the
molecule and its neighboring images. At the very minimum, there should
not be overlap between the orbitals on different image molecules. This
can be actually automatically checked by the program, so the block is
not strictly necessary. However, in the general case it might be
important to have more control over the separation. (This is quite
important for molecules with a dipole, for example, which will have a
long-range interaction with their images in PBC.)

.. hint:: You can play with the size of the lattice parameters to go
   from 'interacting' molecules to effectively isolated ones. Look at
   the variation in the total energy as a function of the cell size,
   to see how the interaction between molecules decreases with
   increasing distance between images. For this non-polar molecule,
   the interaction should be very small. (But see the case of the water
   molecule `here
   <https://personales.unican.es/junqueraj/JavierJunquera_files/Metodos/Basic/H2O/Exercise-H2O.pdf>`_.
 
DFT functional
--------------
 
Up to now we have been implicitly using LDA for our
calculations. However, it is also possible to use other functionals,
such as those of GGA type. Edit the ch4.fdf file to include this block::

  #Density functional (Notice that Xc.authors and XC.functional
  #are both needed and must be consistent)
  XC.functional GGA
  XC.authors  PBE

Note the use of '#' to mark comments, and, once again, the fact that
Siesta uses defaults (in this case an LDA functional) for certain
parameters if they are not specified in the input.

Run the program again and look for possible lines with 'WARNING' or
'ERROR' in them. You will see
that there is a warning. The code does not like that you are using a
GGA functional with a pseudopotential generated using LDA, as this is not
consistent!. Fortunately, we have produced also the pseudopotentials
using GGA for you. They are in the files C.gga.psf and
H.gga.psf. You can modify the input file again to use these files by
simply changing the 'species' strings::

  %block ChemicalSpeciesLabel
  1  6 C.gga   # Species index, atomic number, species label
  2  1 H.gga   # Species index, atomic number, species label
  %endblock ChemicalSpeciesLabel

Run the program again and check whether the warning disappears from
the output.
 
Structural optimization
-----------------------

Now add to `ch4.fdf` the following lines::

  #Geometrical optimization
  MD.TypeOfRun CG
  MD.NumCGsteps 50
  MD.MaxCGDispl         0.1 Bohr
  MD.MaxForceTol        0.04d0 eV/Ang

to instruct Siesta to perform a structural optimization using the
conjugate gradient algorithm (you can check the manual to understand
the meaning of the lines added). If you run the program again you will
notice that the output file contains several new sections, each
corresponding to a different structure, in a series that should
converge to an optimal configuration with zero forces. There are
defaults for the tolerance in convergence. (We will cover relaxation
in more detail in :ref:`this tutorial<tutorial-basic-structure-optimization>`. 

.. hint::
   Relax the structure for various basis set sizes (SZ, DZ,
   DZP) and check the differences on geometry and total energy.

.. note::
   The file ch4.ANI contains all the structures generated during the
   relaxation in XYZ format. It can be processed by various graphical
   tools (**more refs**).

Spin polarization in the CH3 molecule
-------------------------------------

.. hint::
   Enter the directory 'CH3'

Now we are going to perform calculations for the molecule CH3. If you
look at the input file `ch3.fdf`` you should realize that we are
requesting, within the LDA, the optimized geometry of the molecule,
using an automatically generated unit-cell.

However, this molecule contains an unpaired electron.  Therefore the
system should show some spin polarization. We can request 
a spin-polarized calculation by including the line::

  Spin-polarized T

If you compare the results of this calculation with those of the
previous one you will see that there is extra output regarding the
total spin moment during the scf cycle. The final energy should be
lower than for the calculation without spin.

.. note::
   To obtain spin polarization we need to break the symmetry
   between the up and down spins. If spin symmetry is somehow imposed
   or assumed in the initial configuration the results will not be
   spin-polarized. You can check this by adding the block::

     %block DM.InitSpin 
      1 0.0
      2 0.0
      3 0.0
      4 0.0
     %endblock DM.InitSpin

   which will set the initial spins to zero on all atoms (check the
   manual to see the meaning of the block contents). When not using
   the block, the built-in Siesta heuristics prepared an initial
   density-matrix with a spin imbalance.

   
 
Plotting densities
------------------

You might  have noticed the lines::


  SaveRho .true.                 # -- Output of charge density
  %block LocalDensityOfStates    # -- LDOS ('charge' for an energy window)
  -6.00  -3.00 eV
  %endblock LocalDensityOfStates

In response to these lines, SIESTA produced two extra files:

* ``ch3.RHO``: contains the values of the self-consistent electronic
  density on the real-space mesh.

* ``ch3.LDOS``: contains the 'charge density' associated only to the HOMO
  of the molecule. We had to specify an energy window (-6 to -3 eV) in
  which we know that there is only this state. (We can get this
  information by looking at the ch3.EIG file).

.. hint::
   You can modify the block ``LocalDensityOfStates`` to plot
   the 'density' associated with different molecular orbitals lying in
   different energy windows.

More details on how to visualize the charge density and other
quantities represented on the real-space grid are given in :ref:`this
how-to <how-to-analysis-tools>`. For this tutorial we will
use Xcrysden.  Execute::

 rho2xsf < rho2xsf.inp

to generate a 'ch3.XSF' file that contains both the total charge
density and the LDOS information.

Then::

 xcrysden --xsf ch3.XSF

will open the Xcrysden window. You need to go into the 'Grid' section
and set the options to select the data-set (ch3.RHO or ch3.LDOS). For
the charge density, you can select your preferred combination of the 'up' and 'down'
charge densities. If you give them factors of '+1' and '-1' you will get the
spin-density, showing in essence the 'unpaired electron'.

 



