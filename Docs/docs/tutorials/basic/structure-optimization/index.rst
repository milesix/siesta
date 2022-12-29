:sequential_nav: next

..  _tutorial-basic-structure-optimization:

Structural optimization using forces and stresses
=================================================

..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.

This tutorial provides an overview of the methods available in Siesta
for structural optimization, which are quite varied in their nature,
scope, and sophistication.

.. note::
   For background, you can refer to these slide presentations from
   past Siesta Schools:

   * Presentation by `Emilio Artacho (2014) <https://siesta-project.org/SIESTA_MATERIAL/Docs/Tutorials/tlv14/slides15.pdf>`_
   * Presentation by `Marivi Fdez-Serra (2017) <https://siesta-project.org/SIESTA_MATERIAL/Docs/Tutorials/max-2017-bsc/Talks/MVS-ForcesRelaxation.pdf>`_ 

Most DFT codes use the Hellman-Feynman theorem to compute forces and
stresses directly, without having to perform several energy
calculations and use finite differences. This affords tremendous gains
in efficiency. With forces and stresses at hand given the structure,
one can implement algorithms that will search for the optimal structure with
zero forces and stresses.

.. note::
   This is a good time to review the issues related to the use
   of a :ref:`real-space grid<tutorial-basic-grid-convergence>`, and how to
   get appropriately accurate forces and stresses.
   
We have already seen geometry relaxation in action in the context of
the CH4 molecule in :ref:`tutorial-basic-first-encounter`. We will
re-take molecules further below, when we discuss constraints, but now
we turn to examples with crystals.

.. hint::
   Enter directory SiH
   
We will use the (slightly artificial) example of a cubic box of
64-atoms of Si with an extra H atom placed in the middle of a Si-Si
bond. That is not quite optimal, and the H atom 'pushes out' its
neighbors to lower the energy and minimize the forces.

.. (this is a very special example because, due to symmetry, the forces
   on the H atom are always zero!).

For the purposes of running this example quickly, we are using a SZ basis set and a
very low cutoff. You might want to move beyond this later, but there
should not be qualitative changes in the results.

We look first at the options for relaxation::

  MD.TypeOfRun         CG
  MD.NumCGsteps        50
  MD.MaxForceTol       0.05 eV/Ang

CG stands for "conjugate gradients", and is one of the standard
optimizers within Siesta.

For this example the program will run for around 15 *geometry
iterations* (that is, steps in which the structure is modified in
response to the forces, until the cycle converges). The tolerance for
convergence is given by the last line above (note the explicit units),
and is slighty above the default in Siesta (0.04 eV/Ang).

You might want to see the evolution of the total energy of the system
during the run. You can scroll through the file or, as a shortcut, use
the command (assuming you have redirected the output to file *OUT*::

  grep enth OUT

There is a gain of more than 5 eV in the relaxation.

To check the enviroment of the H atom before and after the relaxation,
you can look in the *.BONDS* and *.BONDS_FINAL* files::

  tail *.BONDS        # distances to the neighbors of H
                      # at the start

  tail *.BONDS_FINAL  # distances to the neighbors of H
                      # at the end

You can also use a visualizer program that can read the *.xyz* file
produced by the option ``write-coor-xmol T`` in the fdf file.

It is possible to use a variety of optimization engines within
Siesta. Above we used ``MD.TypeOfRun CG`` to select the
conjugate-gradients engine, but one could use also ``MD.TypeOfRun
Broyden``, or ``MD.TypeOfRun FIRE`` to employ other algorithms, which
might be more efficient than `CG`. (Test this by replacing the option
and starting the relaxation again -- you should see, for example, that
the `Broyden` optimizer needs fewer steps (about half!) to reach convergence.)

These optimization engines are coded in Siesta itself. It is possible
to use the :ref:`Lua scripting engine<tutorial-lua-engine>` embedded in
Siesta to access other algorithms. In particular, a version of the
very efficient LBFGS algorithm can be used for structural
optimization.

Variable-cell optimization
--------------------------

In the example above the atoms are moved in response to the forces,
but the lattice vectors stay fixed. If you have scrolled through the
output file you might have seen lines mentioning "stress". For
example::

  Stress tensor Voigt[x,y,z,yz,xz,xy] (kbar):      -70.56      -70.56      -70.56        0.37        0.37        0.37

A non-zero stress means that the lattice vectors are not optimal, and
the energy can be lowered by optimizing them. The relaxation algorithm uses
now two sets of variables: atomic coordinates and lattice vectors, and
moves them in response to the forces and the stresess, until these are
below a set tolerance.

Variable-cell optimization is non-trivial to implement or run
properly, as the forces and the stresses are not only dependent on the
coordinates and lattice vectors, respectively. Both sets are
coupled. Besides, the dimensions (energy/length and energy/volume,
respectively for forces and stresses) and the 'spring constants' that
determine the changes in them in response to distortions are of course
different.

.. hint::
   Enter directory *varcell_cg*

The example we will use is an 8-atom Si cell (the 'conventional cubic
cell'). Note the options::

  MD.TypeOfRun          CG
  MD.NumCGSteps         100
  MD.VariableCell       T   
  MD.MaxForceTol        0.1 eV/Ang
  MD.MaxStressTol       0.1 GPa
  MD.TargetPressure     0.0 GPa

As expected, the tolerances have different physical dimensions. The
last line tells the program to aim for a zero (diagonal) stress (in
practice, "atmospheric pressure"), but it can be changed to obtain the
optimum geometry of the system under applied pressure (hydrostatic in
this case, but an arbitrary tensor target can be also specified for
non-hydrostatic conditions -- see the manual).

The example is set-up to start with a structure determined by::

  LatticeConstant     5.535 Ang
  %block LatticeVectors
   1.150  0.200  0.000
   0.000  1.050  0.000
  -0.100  0.000  0.900
  %endblock LatticeVectors

We see that the first two cell vectors are too large, so the system is
under tensile stress along the x and y directions. Conversely, the
third cell vector is too small (compressive stress). In addition,
there are off-diagonal components of the strain, leading to 'shear'
components of the stress. Note the signs and sizes of these initial
stresses, and watch how they move towards zero. (You can do this by
``grep oigt OUT``, where *OUT* is the output file you have chosen.)
Also, check the evolution of the energy (``grep enth OUT``)

There are a number of things to note:

* The stresses do not decrease monotonically to zero.
* The energy does not decrease monotonically either.
* The final lattice vectors look "funny" when one looks at their
  cartesian components (for example, in file *si8.STRUCT_OUT*), but
  their modules are roughly the same and the angles between them are
  very close to 90 degrees (grep for `modules` and `angles` in the
  output file). The cell might have rotated during the process of
  relaxation, but it is basically cubic at the end.
* The atomic coordinates (in fractional form) are very close to their
  initial values, which are the standard sites in the Si diamond
  structure (in the conventional cell). This is to be expected, since
  these are high-symmetry positions.

Some of these might be related to, or made worse by, the ridiculously low
mesh-cutoff chosen (30 Ry) and by the small basis set (try improving
these, but note the increased cpu time).

Relaxation with "quenched" molecular dynamics
---------------------------------------------

There is an alternative relaxation method that uses a physically
motivated scheme, rather than a purely mathematical search for a 'zero
forces and stresses' configuration. Imagine that we perform a
molecular dynamics simulation in which, rather than relaxing, we move
the atoms, and the cell vectors, according to the (classical)
equations of motion, using the forces and stresses. For more
information about this, see :ref:`this
tutorial<tutorial-molecular-dynamics>`. The trick in this case is
that, every time an atom senses a force 'opposite' its velocity (in
the sense that their scalar product is negative), the velocity is set
to zero. This roughly corresponds to the idea: "since the atom seems to
be moving away from its equilibrium point, we rather stop it". The
same can be done with strains and stresses in the case of variable
cell.

.. hint::
   Enter the directory *varcell_md*

Look now at the "relaxation section" of the *si8.fdf* file::

 MD.TypeOfRun          ParrinelloRahman
 MD.InitialTimeStep    1
 MD.FinalTimeStep      200
 MD.LengthTimeStep           3.0 fs
 MD.ParrinelloRahmanMass  10.0 Ry*fs**2
 MD.Quench             T
 
The ``ParrinelloRahman`` scheme is a combined atoms+cell
microcanonical scheme (see MD tutorial). We allow it to run for up to
200 steps, with a time-step of three femtoseconds. Note also the
appearance of a "mass" with the dimensions of energy*time^2: this is
used to homogeneize the dynamics of the system, which has to deal, as
we indicated earlier, with fundamentally different sets of variables.
The final line request the "quenching".

If you now run the example, you will notice that it converges quite
nicely, with monotonic decrease in the energy and cleaner evolution of
the stresses (even if not monotonic). This method is always quite
robust, and in this particular case of variable-cell, particularly
efficient. It can generally be counted on to bring systems closer to
the optimal structure, and the final relaxation can be done with a
faster method.

.. note::
   We can even start with some extra kinetic energy, in the
   form of some 'starting temperature', to wiggle things around and
   free the system from any undesired local minimum.

Constraints
-----------

Our next system is a very simple model of the H-terminated (100)
as-cut surface of Si. It is a periodically-repeated slab, with three
layers of Si atoms (six atoms in total), and four H atoms. We want to
know how the top-most section relaxes, while maintaining the
"bulk-like" bottom layers fixed to simulate the connection to the bulk
below.

.. hint::
   Enter the *si100_constrained* directory

Note the block::

  %block GeometryConstraints
        position from 1 to 4
  %endblock GeometryConstraints

which requests that the positions of the first four atoms (those at
the bottom of the slab) are kept fixed.

.. the model surface is not very clear. We should get a new one.

The system converges (with the Broyden algorithm) in about 15
steps. During the process, forces on all atoms are computed, but only
those on non-constrained atoms enter into the check for convergence
(the maximum absolute value of any component of these forces is what
``grep constrained OUT`` prints).

There are many options for the specification of relaxation constraints
in Siesta, including features such as fixing a group of atoms
('molecule') to move rigidly together, constraining cell vector sizes
or angles, etc. You should check the manual to get more information.

In the next section we will explore another way to express
constraints, by using reduced coordinates.

Z-matrix constraints
....................

.. hint::
   Enter the *h2o_zmatrix* directory

Sometimes the important structural degrees of freedom that we want to
optimize are not easily represented in cartesian coordinates. In
molecules, for example, bond lengths, angles, and torsion angles are
much more relevant than particular values of the cartesian
coordinates. For many years, chemists have used more appropriate ways
to represent molecular structure, and the `Zmatrix
<https://en.wikipedia.org/wiki/Z-matrix_(chemistry)>`_ is one of
them. Consider this block in file *h2oZ.fdf*::

  %block Zmatrix
  molecule_cartesian
    1 0 0 0   0.0 0.0 0.0 0 0 0
    2 1 0 0   HO1 90.0 37.743919 1 0 0
    2 1 2 0   HO2 HOH 90.0 1 1 0
  variables
      HO1 1.0
      HO2 1.0
      HOH 106.0
  %endblock Zmatrix

which represents the structure of a water molecule by giving the
cartesian coordinates of the oxygen atom (placed at the origin), but
using essentially bond-lengths (HO1 and HO2), and the HOH angle for
the H atoms, instead of coordinates. (There are other pieces of data
that can be explained by looking in the manual.).

Furthermore, these symbols represent *variables*, which can vary
during a relaxation. We have set initially, for example, the
bond-length to the bball-park value of 1 Ang, and the HOH angle to 106
degrees.  If you run the example you will see how this variables are
changed until relaxation within the tolerances is achieved. The
tolerances themselves have a new, more appropriate form::

  ZM.ForceTolLength 0.04 eV/Ang
  ZM.ForceTolAngle 0.0001 eV/deg

You might want to play further with this example in several
directions:

* Fix the HOH bond angle. For this,
  introduce a new section *constants:* in the block (see the manual)
  and place it there. 

* Compare the results to experimental data. Maybe you need to use a
  GGA functional for better results. Try to get (PSML)
  pseudopotentials from Pseudo-Dojo for this.

* You can use an extendend Zmatrix format to study molecules near
  surfaces. See the manual for an example. 
 







   















   
   

   
  
  
   

   
