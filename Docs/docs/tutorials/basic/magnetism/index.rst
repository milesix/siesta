:sequential_nav: next

..  _tutorial-basic-magnetism:

Magnetism
=========

..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.

:Author: Daniel Sanchez-Portal and Andrei Postnikov (2010)

Initialization of different magnetic orders in a magnetic insulator
(MnO), and analysis of the results.

.. As it is, this exercise is a bit "hard", because it demands some
   skill in setting up the AF structures, and in analyzing the
   output. There are hints in Answers for now
      
The aims of the exercise:

* To learn to set up a calculation in a bulk (crystal) system composed
  of "magnetic" (i.e., `3d`) and "non-magnetic" (`s-p`) atoms, whereby
  spin moments of "magnetic" atoms may be arranged in different types
  of magnetic orderings, and hence allow different (metastable)
  magnetic solutions.

* To analyze the results of converged calculations, demonstrating the type
  of magnetic ordering in each of the cases considered.

.. note::
   3d oxydes are notorious examples of systems in which Coulomb correlation
   effects play an important role, so that the "conventional" DFT treatment
   is in some senses misledaing: the insulating band gaps are underestimated,
   and the placement of main features of the band structure not consistent
   with spectroscopic experiments. A large part of these inconsistencies
   can be fixed by applying the "DFT+U" formalism, as discussed in
   :ref:`another tutorial<tutorial-dft+u>`. 

.. note::
   Another complication is that the crystal structure of real
   3d oxydes, basically very close to the B1 (NaCl) type, undergoes
   slight but noticeable distortions which help to stabilize one or
   another magnetic phases.  However, in the hystorical context as
   well as for didactic purposes, conventional DFT calculations, done
   in a nominally cubic lattice, play an important role.

In MnO, the system under study in the present exercise, the origin
of antiferromagnetism can be understood; the correct antiferromagnetic
phase has the lowest total energy and a band gap, following already
from the results of conventional GGA calculation.

Neglecting very small distortions which occur in different magnetic
phases, MnO has a B1 (NaCl) structure. The latter is described in
the input file as follows::

  LatticeConstant     4.43 Ang
  %block LatticeVectors
  0.00     0.50      0.50
  0.50     0.00      0.50
  0.50     0.50      0.00
  %endblock LatticeVectors

  AtomicCoordinatesFormat ScaledCartesian
  %block AtomicCoordinatesAndAtomicSpecies
  0.000   0.000   0.000  1   # Mn
  0.500   0.500   0.500  2   # O
  %endblock AtomicCoordinatesAndAtomicSpecies

The lattice constant is set at the experimental value, kept constant
throughout the present exercise. 

In order to run a spin-polarized calculation, we set::

  SpinPolarized           .true

and moreover initialize non-zero (in fact, maximal) spin on Mn site::

  %block DM.InitSpin       # Initial magnetic order (on Mn only)
  1   +
  %endblock DM.InitSpin

The rest of the input file sets calculation parameters (GGA, k-mesh,
mesh cutoff, ...) and notably provides printing out Mulliken
populations, saving charge (and magnetic) density, and writing down
the density of states::

  WriteMullikenPop  1
  SaveRho           
  %block ProjectedDensityOfStates
  -25.0  10.0  0.1   700   eV
  %endblock ProjectedDensityOfStates

This output part will be identical for calculations with different
magnetic orderings.

Tasks:

1. When the first calculation is done, check the local (on Mn and O sites)
and total (per unit cell) magnetic moments. Plot the spin-resolved
density of states (total and that for Mn and O sites).

2. Prepare the input file(s) for antiferromagnetc structures of MnO
and run the calculations. Note that an antiferromagnetic ordering 
doubles the unit cell: it contains now two Mn atoms (with opposite
spins) and, correspondingly, two O atoms. Don't forget to introduce
changes in the lattice vectors, and atomic coordinates.

We'll consider two different antiferromagnetic structures:

* AF1 has, what is called, the [001] ordering, that is, all Mn atoms
  in a given [001] plane are equivalent (has the same spin
  orientation), and between the consecutive [001] planes the spin
  orioentations alternate.

* AF2 has the [111] ordering, that is, equivalent spin orientation is
  throughout a given [111] plane, and alternates between such planes.

Prepare the ``%block LatticeVectors``
for AF1 and AF2 cases. Check that the unit cell volume (mixed product
of three lattice vectors) is twice that of the FM case (which equals :math:`a_0^3/4`).
Provide coordinates of (four) atoms, initialize spins in opposite sense
at two Mn atoms, and run the calculations.

3. Check Mulliken populations; plot total (from corresponding DOS files)
and (optionally) atom-resolved densities of states. 
Check that the AF structures have pronounced band gaps.

4. Using the calculated RHO, visualize the spin density in each
magnetic structure. You'll see that the spin density has nearly perfect
spherical distribution in the vicinity of Mn sites; why?

5. For the AF2 structure, you can find in the occupied part
of the Mn3d-related density of states two groups of peaks, related to
t2g and eg states of Mn. Try to vizualize them separately, selecting
for each of these groups a corresponding energy interval and calculating
LDOS. Why the separation into t2g and eg states is not so neat
in the AF1 structure? 
   



