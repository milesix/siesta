:sequential_nav: next

..  _tutorial-basic-vibrational-properties-benzene:

Modes of vibration of the benzene molecule
==========================================

:Author: Javier Junquera and Andrei Postnikov
	 
..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.

The goal of this exercise is to compute the vibrational frequency of
a molecule (benzene).

Before running the calculation to compute the vibrational frequencies,
the first step is to relax the geometrical structure of the system under
study. 

So, to start with, we run a conjugate gradient minimization to relax
the atomic positions.  The input file has been prepared for you in the
file benzene.relax.fdf.  See how the structure of the benzene molecule
has been introduced in the **Z-matrix format** (especification of internal
variables, such as distances, angles, and torsional angles).  This
allows the minimization including some constraints in the symmetry in
a trivial way. We just allow the C-C and C-H distances to relax::

  %block Zmatrix
  molecule
  2 0 0 0  xm1 ym1  zm1   0 0 0
  2 1 0 0  CC  90.0 60.0  0 0 0
  2 2 1 0  CC  CCC  90.0  0 0 0
  2 3 2 1  CC  CCC  0.0   0 0 0
  2 4 3 2  CC  CCC  0.0   0 0 0
  2 5 4 3  CC  CCC  0.0   0 0 0
  1 1 2 3  CH  CCH  180.0 0 0 0
  1 2 1 7  CH  CCH  0.0   0 0 0
  1 3 2 8  CH  CCH  0.0   0 0 0
  1 4 3 9  CH  CCH  0.0   0 0 0
  1 5 4 10 CH  CCH  0.0   0 0 0
  1 6 5 11 CH  CCH  0.0   0 0 0
  constants
    ym1 5.00
    zm1 0.00
    CCC 120.0
    CCH 120.0
  variables
    CC 1.390
    CH 1.090
  constraints
    xm1 CC -1.0 3.903229
  %endblock Zmatrix

We run Siesta to carry out the relaxation::
  
    siesta  benzene.relax.fdf > benzene.relax.out
  

Computing the force constants in real space
--------------------------------------------

*  We have prepared an input file, called benzene.ifc.fdf
   to run Siesta and compute the interatomic force
   constants in real space:

   * We have copied the relaxed coordinates and unit cell from
     the *benzene.XV* generated after the relaxation to the 
     `AtomicCoordinatesAndAtomicSpecies` block. 

   * We have included the atomic masses after the coordinates 
     of each atom. This will be useful for the next step with the
     *vibra* code
  
   Here are the relevant sections::
     
     LatticeConstant     1.0 Bohr
     %block LatticeVectors
     20.932528150       0.000000000       0.000000000
     0.000000000      19.551203193       0.000000000
     0.000000000       0.000000000      10.714661844
     %endblock LatticeVectors

     AtomicCoordinatesFormat NotScaledCartesianBohr
     %block AtomicCoordinatesAndAtomicSpecies
     4.738724869       9.448634389       0.000000000     2    12.0107
     6.057380810      11.732613477      -0.000000000     2    12.0107
     8.694692693      11.732613477      -0.000000000     2    12.0107
     10.013348634       9.448634389      -0.000000000     2    12.0107
     8.694692693       7.164655301      -0.000000000     2    12.0107
     6.057380810       7.164655301      -0.000000000     2    12.0107
     2.647979028       9.448634389      -0.000000000     1     1.00794
     5.012007889      13.543252488      -0.000000000     1     1.00794
     9.740065613      13.543252488      -0.000000000     1     1.00794
     12.104094475       9.448634389      -0.000000000     1     1.00794
     9.740065613       5.354016289      -0.000000000     1     1.00794
     5.012007889       5.354016289      -0.000000000     1     1.00794
     %endblock AtomicCoordinatesAndAtomicSpecies

To compute the interatomic force constant in real space, we have
to run Siesta::

     siesta < benzene.ifc.fdf > benzene.ifc.out 

The interatomic force constant matrix in real space is stored
in a file with name of the form  *SystemLabel.FC*.
 
.. Again, the explanation of the different entries of this file can
   be found in the theoretical lectures.


Computing the dynamical matrix at the Gamma point, and (phonon) modes
---------------------------------------------------------------------

Once the interatomic force constants in real space have been computed,
a discrete Fourier transform is performed to compute the dynamical
matrix in reciprocal space.  Then, the dynamical matrix is
diagonalized and its eigenfrequencies and eigenvectors are computed.
This is done using the vibra code.

In the case of a molecule, only the Gamma point is relevant.  It is
specified in the same way as to compute the electronic band structure,
in the same file benzene.ifc.fdf::

     Eigenvectors    .true.        # Compute both phonon eigenvalues and eigenvectors
     BandLinesScale  pi/a
     %block BandLines
     1   0.0   0.0   0.0   \Gamma  # Only the Gamma point (enough for a molecule)
     %endblock BandLines

To compute the vibrational frequencies::

  Your_siesta_directory/Util/Vibra/Src/vibra < benzene.ifc.fdf > vibra.out

The output of this code is:

* *SystemLabel.bands*: with the different mode frequencies (in cm^-1).
  They are stored in the same way as the electronic band structure.

* *SystemLabel.vectors*: with the eigenmodes at Gamma
  (the format is self-explained).

.. note::
   Some of the modes might have negative frequencies. How could that
   be?
   
How to visualize the normal modes
---------------------------------

After getting the .vectors file (calculated by vibra) and the .XV file
(computed in Siesta), run the vib2xsf program.

You have to answer a few question on the fly, regarding the name of
the files where the .vectors are stored, the units to be used to
introduce the lattice vectors (Bohrs or Angstroms), the zero of
coordinates, the unit cell lattice vectors, the first mode to
visualize, the last mode to visualize, the amplitude of the modes to
be visualized, and the number of steps in the movie.

You can play a little bit, but to save time we have prepared
all the answers in the file vib2xsf.dat for you.   Just run::

      vib2xsf < vib2xsf.dat

This will produce two files per mode:

* .XSF file: contains a static structures (as in .XV), 
  with arrors added to each atom to indicate displacement pattern.

* .AXSF file: contains the animation of a phonon, for a (user-chosen) 
  amplitude and number of steps.

They can be visualized using XCRYSDEN::

      xcrysden

      Select "File"
      Open Structure
      Open AXSF (Animation XCrySDen Structure File)

The same can be done to visualize the XSF file, but just choosing::

      Select "File"
      Open Structure
      Open XSF file (XCrySDen Structure File)

.. note::
   It might be interesting to analyze those modes with negative or small
   frequencies.
   


     
  







   

   
    
	       
   
   

   
  
  
   

   
