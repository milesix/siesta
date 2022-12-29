:sequential_nav: next

..  _tutorial-basic-vibrational-properties-si-bulk:

Phonon dispersion of bulk Si
============================

:Author: Javier Junquera and Andrei Postnikov

..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.


The goal of this exercise is to compute the phonon band structure of
bulk Si in the diamond structure.

.. note::
   Before computing the vibrational structure, it is always necessary
   to relax the geometrical structure of the system under
   study, so that the reference configuration truly corresponds to a
   minimum of the energy, with zero forces.

   In this case, the forces on the Si atoms are always zero by
   symmetry, and the only structural parameter left is the lattice
   constant. We can run phonon calculations for any (reasonable) value
   of the lattice constant. In fact, useful information about
   average anharmonicity (the "Gruneisen parameter") can be obtained
   by monitoring the dependency on volume of the phonon frequencies.

   For this exercise we will just use a lattice constant obtained from
   a previous optimization with the same basis set and parameters
   (the value is slightly larger than the experimental one).


Building the supercell to compute the force-constant matrix in real space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  As explained in the slide presentation
   `<https://drive.google.com/file/d/1SRhSOWwFVPmS3vr4REJSvqgXf2Tic8fl/view?usp=sharing>`_,
   if we know the force constant matrix in real space we can compute
   the dynamical matrix for every q-point in reciprocal space. For a
   bulk system, we need a way to record the forces felt by atoms in
   neighboring unit cells when we displace an atom in the unit
   cell. Hence the need for a *supercell*.

*  Run the fcbuild program to generate the supercell that will be used
   to compute the interatomic force constant matrices in real space.
   From the time being, we will generate a supercell replicating the
   unit cell three times (named -1, 0, and 1) along each cartesian direction
   (edit the Si.fcbuild.fdf and see the lines::

     #
     # Options to generate the supercell
     #

     SuperCell_1    1     # number of shells in which the unit cell is
     #   repeated in the direction of the first lattice vector.
     SuperCell_2    1     # Idem for the second lattice vector.
     SuperCell_3    1     # Idem for the third  lattice vector.

To generate the supercell run::

        fcbuild < Si.fcbuild.fdf
 
This code dumps the information of the Supercell in an output file,
called FC.fdf, that contains

   - The structural data of the supercell, including

        * The number of atoms. 
        * The lattice constant. 
        * The lattice vectors. 
        * The atomic coordinates and the atomic species of all the atoms

   - The variables required to compute the interatomic force constants.
     in real space.

        * Atoms that will be displaced (those belonging to the home unit cell).
        * The amount by which the atoms will be displaced.

For a final analysis of the results, the supercell should contain
enough atoms so that all non-negleable elements of the force constant
matrix are computed. The range in real space in which the force
constant matrix decays to zero varies widely from system to system!!.


Computing the force-constant matrix in real space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


*  We prepare an input file, called Si.ifc.fdf,
   to run Siesta and compute the interatomic force
   constant in real space. 

*  Many variables are taken directly from the file where the supercell is
   described (FC.fdf)

*  To compute the interatomic force constant in real space, we have
   to run Siesta

          siesta < Si.ifc.fdf > Si.ifc.111.out 

*  The interatomic force constant matrix in real space are stored
   in a file called SystemLabel.FC
 
..   Again, the explanation of the different entries of this file can
     be found in the theoretical lectures.


Computing the q-dependent dynamical matrix, and phonon modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


*  Once the interatomic force constants in real space have been computed,
   a discrete Fourier transform is performed to compute the dynamical matrix
   in reciprocal space. 
   Then, the dynamical matrix is diagonalized and its eigenfrequencies and
   eigenvectors are computed.
   This is done using the vibra code.
 
   The k-points are defined in the same way as to compute the electronic
   band structure, in the same file used to define the supercell
   
       vibra < Si.fcbuild.fdf

*   The output of this code is:

    SystemLabel.bands: with the different mode frequencies (in cm^-1).
    They are stored in the same way as the electronic band structure.

    SystemLabel.vectors: with the eigenmodes for each k-points 
    (the format is self-explained).

*   To plot the phonon band structure, proceed in the same way as 
    to plot the electronic band structure (using the gnubands.x, etc).
   
       gnubands  < Si.bands > Si.phonon-bands.111.dat

       gnuplot
       gnuplot> plot "Si.phonon-bands.111.dat" using 1:2 with lines

*   To produce a postscript file of this figure to be included in a document

gnuplot> set terminal postscript
gnuplot> set output "Si.phonon-bands.111.ps" 
gnuplot> replot


Checks of the convergence with the supercell size
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
   It might take too long to generate the FC files for the 222
   and 333 cases. We provide them in the FILES directory. You will
   need to adapt the instructions below to this case (rename or copy
   files as needed).

*  One should always check the convergence of the computed phonon 
   band structure with respect the size of the supercell, to be sure
   that all the relevant interatomic force constant matrix elements
   are included.

   (Note: the simulations for larger cells require more than hour of CPU
   time to generate the force constant matrix. You can either repeat
   the procedure explained below or directly take the force constant
   matrix prepared for you, direct output of the proposed simulations.
   The name of the output files are Si.222.FC and Si.333.FC respectively).

*  To do this:

   - First, we save all the input and output files used upto now
     in order to be overwritten::

       $ cp Si.fcbuild.fdf Si.fcbuild.111.fdf
       $ mv FC.fdf FC.111.fdf
       $ mv Si.FC Si.111.FC
       $ mv Si.vectors Si.111.vectors
       $ mv Si.bands Si.111.bands

   - Edit the file Si.fcbuild.fdf and increase the size of the supercell,
     adding up to 5 periodic repetitions of the unit cell in each direction
     (named -2, -1, 0, 1, 2) ::

       #
       # Options to generate the supercell
       #

       SuperCell_1    2     # number of shells in which the unit cell is
       #   repeated in the direction of the first lattice vector.
       SuperCell_2    2     # Idem for the second lattice vector.
       SuperCell_3    2     # Idem for the third  lattice vector.

    - Repeat the previous procedure for SuperCell_1,2,3 = 2::

	fcbuild < Si.fcbuild.fdf
	siesta < Si.ifc.fdf > Si.ifc.222.out 
	vibra < Si.fcbuild.fdf
	gnubands < Si.bands > Si.phonon-bands.222.dat
	gnuplot
	gnuplot> plot "Si.phonon-bands.222.dat" using 1:2 with lines

	$ cp Si.fcbuild.fdf Si.fcbuild.222.fdf
	$ mv FC.fdf FC.222.fdf
	$ mv Si.FC Si.222.FC
	$ mv Si.vectors Si.222.vectors
	$ mv Si.bands Si.222.bands

    - Repeat the previous procedure for SuperCell_1,2,3 = 3::

	fcbuild < Si.fcbuild.fdf
	siesta < Si.ifc.fdf > Si.ifc.333.out
	vibra < Si.fcbuild.fdf
	gnubands < Si.bands > Si.phonon-bands.333.dat
	gnuplot
	gnuplot> plot "Si.phonon-bands.333.dat" using 1:2 with lines

	$ cp Si.fcbuild.fdf Si.fcbuild.333.fdf
	$ mv FC.fdf FC.333.fdf
	$ mv Si.FC Si.333.FC
	$ mv Si.vectors Si.333.vectors
	$ mv Si.bands Si.333.bands

      (for this, you might have to edit the vibra.h file in the
      Vibra/Src directory, change the values of ::

          parameter (maxx = 3)
          parameter (maxy = 3)
          parameter (maxz = 3)

      and recompile the code typing ``make``).

*  To compare the results obtained with the three superlattices::

     $ gnuplot
     gnuplot> plot "Si.phonon-bands.111.dat" using 1:2 with lines,
                    "Si.phonon-bands.222.dat" using 1:2 w l,
		    "Si.phonon-bands.333.dat" u 1:2 with lines

*   You can produce postscript files as indicated above





     
  







   

   
    
	       
   
   

   
  
  
   

   
