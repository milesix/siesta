:sequential_nav: next

..  _tutorial-basic-first-crystals:

First crystals
==============

..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.

In this exercise we will move beyond molecules to treat crystals, that
is, periodic solids, also known as `bulk systems`. 


.. hint::
   Enter the directory 'MgO'


In bulk systems the electronic states can be characterized by a
continuous quantum number k, the so-called Bloch vector or crystalline
quasi-momentum of the electrons.  The space of allowed Bloch vector
(the Brillouin zone) is usually sampled using a finite grid. The most
popular scheme to generate this integration mesh is the Monkhorst-Pack
algorithm. This is also used in SIESTA.  The integration grid in the
Brillouin zone is specified using the block kgrid_Monkhorst_Pack:


.. code-block:: bash

  %block kgrid_Monkhorst_Pack
     6  0  0  0.5
     0  6  0  0.5
     0  0  6  0.5
  %endblock kgrid_Monkhorst_Pack

where the 3x3 matrix of integer numbers defines the mesh and the last
number (an optional float) in each line defines the displacement,
which might be useful to reduce the number of (symmetry irreducible)
k-points used in the calculation (for this, 0.5 is a good default),
or, for some applications, to exclude the Gamma point from
consideration. If not specified, the displacement will be zero in each
direction. The program will warn the user if the displacements chosen
are not optimal.

Explore the `Mg.fdf` file in search for the Monkhorst Pack block.

Run the simulation, which will take a few seconds:

.. code-block:: bash

  siesta  MgO.fdf > MgO.out

We will now analyse the band-structure for MgO and look at the density of
states.

To plot the DOS you will use the utility program 'Eig2DOS' (see
:ref:`this how-to<how-to-dos-pdos>` for background and installation
notes, if needed). This program reads the file MgO.EIG (always
produced by Siesta), which contains all the eigenvalues for each k-point used
to sample the BZ, and the file MgO.KP, which contains the k-point
sampling information.  Useful options to the program (type 'Eig2DOS
-h' for a full list) are the broadening for each state in eV (a value
of the order of the default (0.2 eV) is usually reasonable), the
number of energy points where the DOS will be calculated (200 by
default) and the Emin and Emax of the energy window where the DOS will
be calculated (the default is to compute the DOS in the whole range of
energies available in the EIG file).

For example::

   Eig2DOS -k MgO.KP MgO.EIG > dos 

will compute the DOS in the whole range, using a grid of 200 points, a
broadening ("smearing") of 0.2 eV, and using the k-point info from
MgO.KP. (The k-point file information is important when different
k-points have different weights).

Plot the dos using gnuplot, type ``gnuplot`` and enter the commands:

.. code-block:: gnuplot

  plot "dos" with lines

The result is a bit confusing, but it is just because the range of energies
is very large and the number of grid points not so big.
Open the `MgO.out` file and look for the Fermi energy.
Now generate the dos only for a range of energies of about 30 eV around the Fermi
energy (use the ``-E`` and ``-e`` options of Eig2DOS). The result will be
more familiar.

The file MgO.fdf will also produce a file `MgO.bands` containing the
band structure along the several high-symmetry lines in the Brillouin zone (BZ).
This file is produced because in input the
the block BandLines was included::

  BandLinesScale       pi/a
  %block BandLines
  1   1.5   1.5   0.0   K             # Begin at K
  38  0.0   0.0   0.0   \Gamma        # 38 points from K to Gamma
  36  0.0   2.0   0.0   X             # 36 points from Gamma to X
  18  1.0   2.0   0.0   W             # 18 points from X to W
  26  1.0   1.0   1.0   L             # 26 points from W to L
  31  0.0   0.0   0.0   \Gamma        # 31 points from L to Gamma
  %endblock BandLines

Let's understand these few lines:

* The **BandLinesScale** entry specifies the scale of the k-vectors given in the BandLines block.
  The possible options are "pi/a" and "ReciprocalLatticeVectors".

* Each line of the BandLines block represents the end-point of a segment to be spanned in the
  reciprocal space. The comments on the code block above help to familiarize with the format
  of these k-space segments.

To learn more about the special high-symmetry points and lines of the
BZ for this case, you can visit the `Bilbao Crystallographic Server
<https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-kv-list?gnum=225&fig=fm3qmf>`_
(click on 'Optimized listing of k-vector types...'; then on 'Brillouin
Zone').

.. (Note that the final Gamma point in the BandLines sequence is actually
   an image along (1,1,1) of the Gamma point at the origin.)

The bands are stored in the file `MgO.bands`.
To plot the band structure you need to use the utility program ``gnubands``:

.. code-block:: bash
 
   gnubands < MgO.bands > bands.dat

and you can plot `bands.dat` with ``gnuplot``, using the `bands.gplot`
file provided (which has extra information to place and label the
'ticks' for each symmetry point in the BZ):

.. code-block:: bash

  gnuplot -persist bands.gplot


.. hint::
   Enter the directory 'Al'


For metals such as Al there are electronic bands that are not completely filled
and therefore for an accurate description of the total energy, forces
and all properties of the materials it is necessary to use a better
sampling in reciprocal space (Bloch vectors) than for insulators. In
the input file Al_bulk.fdf a 4x4x4 grid is used. This might be
insufficient for a good description of aluminium.

.. Use eig2bxsf to get the Fermi surface?
   
.. You should explore the convergence of total energy, the lattice
   parameter, and density of states respect to the fineness of the
   k-sampling.

Run the input as it is.

.. code-block:: bash

  siesta < Al_bulk.fdf > Al.out

And create the DOS

.. code-block:: bash

   Eig2DOS -k Al.KP Al.EIG > dos 

If you plot the DOS with gnuplot the result does not look ok.

Create a new folder

.. code-block:: bash

  mkdir MORE_KP
  cp Al_bulk.fdf Al.psf MORE_KP
  cd MORE_KP

Change in the `Al_bulk.fdf` file the mesh to [14 14 14] and run siesta.

.. code-block:: bash

  siesta < Al_bulk.fdf > Al.out

Then:

.. code-block:: bash

  Eig2DOS -k Al.KP Al.EIG > dos

Plotting the DOS, we know recognize the "free-electron-like" curve we see in
textbooks.
For coarse samplings, instead, the DOS was not at all
like the "free-electron-like" curve since too few points were considered.

Another point to notice is that the DOS seems to have an upper limit
(DOS is zero after a certain energy). This is due to the limited number
of orbitals in the basis that, in turn, limits the number of produced bands.
A richer basis will increase the number of bands (see tutorial 
:ref:`on basis<tutorial-basic-basis-sets>`.)

Try to compare also the bands among the case with few kpoints and more kpoints.

A detailed tutorial on kpoints convergence is :ref:`here<tutorial-basic-kpoint-convergence>`.

.. A SZ basis set is specified in the file Al_bulk.fdf.  It might be
   quite interesting to see how the band structure changes when more
   complete basis sets are used (DZ,DZP). You might defer this for the
   :ref:`tutorial on Basis Sets <tutorial-basic-basis-sets>`.
