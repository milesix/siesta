.. _reference_macroave:

Macroave User Guide
===================

   **MACROAVE 1.2.1**

   *November 18, 2003*

   Javier Junquera

   *Departamento de Ciencias de la Tierra y Fı́sica de la Materia
   Condensada
   Universidad de Cantabria, Santander, E-39005, Spain*

   **javier.junquera@unican.es**

   Pablo Ordejón

   *Institut de Ciència de Materials de Barcelona - CSIC, Campus de la
   U.A.B., 08193 Bellaterra, Barcelona, Spain*

.. _sec:intro:

Introduction
------------

The Macroave program implements the *macroscopic average technique*,
introduced by A. Baldereschi and coworkers (A. Baldereschi, S. Baroni,
and R. Resta, Phys. Rev. Lett. **61**, 734 (1988) ). This is an
extremely powerful method that relates microscopic quantities, typical
outputs of first-principles codes, with macroscopic magnitudes, needed
to perform electrostatic analysis. Within this methodology, we will be
able of washing out all the wiggles of the rapidly-varying functions of
position (resembling the underlying atomic structure) of the microscopic
quantities, blowing up only the macroscopic features.

It is a basic tool to calculate some important magnitudes in surface or
interface-related problems, such as:

-  | Band offsets and Work functions:
   | L. Colombo, R. Resta and S. Baroni, Phys Rev B **44**, 5572 (1991)

-  | Effective charges:
   | R. Martin and K. Kunc, Phys Rev B **24**, 2081 (1981)

-  | High-frequency dielectric constants:
   | F. Bernardini and V. Fiorentini, Phys. Rev. B **58**, 15292 (1998)

Macroave reads the magnitude, :math:`f \left( \vec{r} \right)`, whose
macroscopic average will be calculated (typically, charge densities or
potentials) at the points of a three-dimensional uniform real space
grid, as it is dumped into output files by standard first-principle
codes. Then it performs the macroscopic average in a two step process:

-  | First: a planar average of :math:`f \left( \vec{r} \right)` on
     planes parallel to the interface.
   | (To establish the notation, we will call the plane parallel to the
     surface or the interface the :math:`(x,y)` plane, whereas the
     perpendicular direction will be referred to as the :math:`z` axis).

   .. math::

      \overline{f} \left( z \right) =
                   \frac{1}{S} \int_{S}
                   f \left( \vec{r} \right) dx dy
                   \label{eq:planar}

   where :math:`S` is the surface of the unit cell perpendicular to the
   given direction.

-  Second: a final convolution of :math:`\overline{f} \left( z \right)`
   with filter functions. We choose step functions, :math:`\Theta`, of
   length :math:`l`.

   .. math::

      \omega_{l} \left( z \right) 
                       = \frac{1}{l} \Theta\left( \frac{l}{2} - |z| \right)
                      \label{eq:step}

   .. math::

      \overline{ \overline{f}} \left( z \right) =
                             \int dz' \int dz'' \omega_{l_{1}} \left( z-z' \right)
                             \omega_{l_{2}} \left( z'-z'' \right)
                             \overline{f} \left( z'' \right)
                     \label{eq:macro}

Currently, Macroave can handle directly the microscopic information
provided by Siesta and Abinit, but it should be easily adapted to any
other first-principle code.

Coded by J. Junquera and P. Ordejón, April 1999

Adapted for Abinit by J. Junquera, October 2002

This is a short description of the compilation procedures and of the
datafile format for the Macroave code. This version is a very
preliminary release of the code. Please report problems, bugs and
suggestions to javier.junquera@ulg.ac.be

Compilation
------------

Everything is automated within the Siesta distribution. Just go to
Macroave/Src and type ’make’. Optionally, compile the ’permute’ program
if needed.

Running the program
-------------------

As it was mentionned in the introduction (Section  `1 <#sec:intro>`__),
Macroave needs as input the microscopic magnitude,
:math:`f \left( \vec{r} \right)`, whose macroscopic average we want to
calculate. :math:`f \left( \vec{r} \right)` will be, typically, a charge
density or a given potential (electrostatic, exchange-correlation only,
total,...). This information, that will be supplied by a
first-principles electronic-simulation code, is usually stored at the
points of a three-dimensional real-space grid.

Obviously, the first thing we must do is to run the
electronic-simulation program for the system we are interested in,
setting up the variables that instruct to write the corresponding
magnitude. At the current time, Macroave is able to digest directly the
output files supplied by Siesta and Abinit. The relevant input variables
*in these first-principles codes* are:

-  Siesta

   -  SaveRho

   -  SaveDeltaRho

   -  SaveElectrostaticPotential

   -  SaveTotalPotential

   -  SaveIonicCharge

   -  SaveTotalCharge

   -  LocalDensityOfStates

-  Abinit

   -  prtpot

   -  prtvha

   -  prtvhxc

   -  prtvxc

   -  prtden

We refer the reader to the User’s Guide of Siesta or Abinit to learn
more about these different options.

Once the simulation is finished, and the relevant output files written,
then move to the directory where the job was run (let’s call it
``\sim/rundir``)

``cd \sim/rundir``

Edit the macroave’s input file (called *macroave.in*) and set up the
right values for the different variables. This file will be fully
explained in section  `4 <#section:input>`__.

**NOTE**: If you have chosen ``x`` as the direction perpendicular to the
slab, instead of ``z``, you need to run the ``Src/permute`` program to
permute the axes.

Execute macroave

``$ \sim/Macroave/Src/macroave``

The output is dumped in files which will be described in Section
 `5 <#section:output>`__.

.. _section:input:

Input data file
---------------

Apart from the information taken from the electronic-simulation code,
Macroave requires only an input data file, named *macroave.in*.

This input file has eigth lines:

**first line**
   (*string*): Name of the first-principles code used to generate the
   microscopic magnitude, :math:`f \left( \vec{r} \right)`. At present,
   it only accepts two options:

   -  Siesta

   -  Abinit

**second line**
   (*string*): Microscopic magnitude whose macroscopic average will be
   calculated:

   -  Potential

   -  Charge

**third line**
   (*string*): Name of the file (output of the first-principles code)
   where the magnitude :math:`f \left( \vec{r} \right)` is stored. In
   the case of Siesta, only the **SystemLabel** is required (see Siesta
   User’s Guide).

**fourth line**
   (*integer*): Number of convolutions with step functions required to
   perform the macroscopic average. It can take only two different
   values:

   -  1 (for surface-related problems).

   -  2 (for interface-related problems).

**fifth line**
   (*real*): Length of the first step function used to perform the
   macroscopic average (see Eq.  `[eq:step] <#eq:step>`__)

   *Units:* bohrs

**sixth line**
   (*real*): Length of the second step function used to perform the
   macroscopic average (see Eq.  `[eq:step] <#eq:step>`__)

   *Units:* bohrs

   *Use:* Only use if the number of convolutions is equal to 2.

**seventh line**
   (*integer*): Electronic charge of the sistem

   *Units:* electrons

   *Use:* Only use if we are computing the macroscopic average of charge
   densities.

**eigth line**
   (*string*): Kind of interpolation to get
   :math:`f \left( \vec{r} \right)` at a fine FFT grid, starting from
   the grid used in the first-principles code.

   At the current time, it only accepts two different values

   -  Spline

   -  Linear

.. _section:output:

Output files
-------------

Two output files are produced, containing the information about the
planar (see Eq.  `[eq:planar] <#eq:planar>`__) and the macroscopic
average (see Eq.  `[eq:macro] <#eq:macro>`__) of
:math:`f \left( \vec{r} \right)`.

Contains, in two colums, values of :math:`z` and the profile of the
planar or macroscopic average.

The name of these output files is the same as the one introduced in the
third line of the input, plus an extension:

.PAV
   for the planar average.

   *Units:*

   -  electrons/bohr\ :math:`^3` if :math:`f \left( \vec{r} \right)` is
      a charge density.

   -  eV if :math:`f \left( \vec{r} \right)` is a potential. density.

.MAV
   for the macroscopic average.

   *Units:*

   -  electrons/bohr\ :math:`^3` if :math:`f \left( \vec{r} \right)` is
      a charge density.

   -  eV if :math:`f \left( \vec{r} \right)` is a potential. density.

Examples
---------

In directory ``\sim/Macroave/Examples`` you will find some examples of
input files.

Known bugs and errors
---------------------

-  The code only works for orthorrombic unit cells.

-  Spin polarization not implemented yet. The planar average, and the
   corresponding macroscopic average are only implemented for the first
   component of array RHO.
