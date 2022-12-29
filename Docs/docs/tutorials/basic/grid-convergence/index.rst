:sequential_nav: next

..  _tutorial-basic-grid-convergence:

The real-space grid
===================

(More material to be added)

In this tutorial we will look more closely at the real-space grid
and how to monitor its adequacy for a given system.

Generalities
------------

The real-space grid is the home of densities, potentials, and other
real-space magnitudes, and is also the scaffolding on which to compute
some contributions to matrix elements in Siesta. It is specified quite
simply. For example::

  Mesh-Cutoff  100 Ry

Internally, the program will translate this "energy yardstick" into a
mesh of points, information about which will appear in the output file
in a way similar to this::

  InitMesh: MESH = 18 x 18 x 30 = 9720
  InitMesh: Mesh cutoff (required, used) =   100.000   101.039 Ry

giving the number of mesh points along each of the lattice
vectors. Note that the actual cutoff used is (in this and in most
cases) higher than requested. This is due to geometrical issues, as
well as to the fact that the grid algorithms used internally by Siesta (in
particular the implementation of the fast fourier transform) need
'magic numbers' (multiples of 2, 3, or 5).


Convergence of properties with mesh-cutoff
------------------------------------------

We will initially use the Methane (CH4) molecule (directory *ch4*)

Check the input file *ch4.fdf*. The system is a methane molecule, inside a cubic cell of 10 Ang
on a side. Look for the the mesh cutoff chosen.
 
Run the example.
 
* Can you calculate what is the distance between mesh points?
 
* Can you estimate how many mesh points fit inside the (first zeta) 2s orbital of C? 
  
Now do a series of calculations increasing the value of the MeshCutoff
parameter, in 10 Ry steps. From the output files, extract:
 
* The actual MeshCutoff used (which as we mentioned above is usually different from the one
  required in the input file)
* The total energy
* The total force (that is, the sum of forces over all the atoms). It
  will *not* be zero!
* The CPU  time. You might want to turn on the option ``use-tree-timer
  T`` for a hierarchical display of the different sections.
 
Plot the total energy and total force versus the requested meshcutoff
and the used meshcutoff. Try to decide what would be a reasonable
value for a real calculation in this system.

Note the cpu-time implications. In particular, compare the time taken
by the diagonalization ("compute_dm" in the timings) and the setup of
the hamiltonian ("setup_H"), for various basis cardinalities (SZ,
DZP...) and values of mesh-cutoff. You might want to check these
timings for other systems as you explore them in these tutorials.

In this case, since we have comparatively a lot of vacuum, and very
few orbitals, the cost of the grid operations is substantially higher
than that of the diagonalization step. 
 

The egg-box effect
------------------

.. note::
   For background, you can see `this slide presentation
   <https://personales.unican.es/junqueraj/JavierJunquera_files/Metodos/Convergence/Eggbox-MgO/Exercise-eggbox-MgO.pdf>`_
   by Javier Junquera.

   The tutorial will offer similar worked examples, but the files and
   description are yet to be uploaded.

..  focus on the egg-box effect and two methods to
    alleviate it: increase mesh-cutoff and use of 'grid-cell-sampling'.

..  Optionally: suggest using Lua to automate the displacements.
  
  
   

   
