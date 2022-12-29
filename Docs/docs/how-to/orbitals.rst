:sequential_nav: next

..  _how-to-orbital-visualization:

Visualizing orbital shape
=========================


Using the ioncat/ionplot tools to process .ion files
----------------------------------------------------

You can look at the shape of the orbitals by plotting the
contents of the .ion files produced by Siesta. These files are not
easily readable, but the 'ioncat' program can extract the relevant 
pieces of information, and the "ionplot" script can drive 'ioncat' to
plot the desired graphs. For example::

        ionplot -o 1 O

will plot the orbital with number "1" in the O.ion file.::

        ioncat -i O

will print the numbers of the representative orbitals of each nlz
shell (i.e., disregarding the 'm' quantum number, which does not
affect the radial part).::

        ioncat -o 1 O 

will output the data for the first orbital in O.ion.

Typing ``ioncat -h`` produces a display of the full set of options for
the program::

  Usage: ioncat [options] Species_Label
  Options:

  -s          : Show header information
  -i          : Print indexes of unique orbitals
  -j          : Print indexes of unique KB projectors

  -o ORBINDEX : Generate table for orbital ORBINDEX
  -k KBINDEX  : Generate table for KB proj KBINDEX
  -v          : Generate table for Vna potential
  -c          : Generate table for pseudocore charge
  -l          : Generate table for Vlocal charge
  
  -Z          : Zoom in near rc for table generation
  -O rlog     : Zoom near zero (up to 10^rlog)
  
  -h          : Print this help message

The script ``ionplot`` accepts the options (such as -o)
that generate data and drives also gnuplot to plot the information.

By combining the two one could plot all the orbitals::

 for i in $(ioncat -i O); do ionplot -o $i; done

Note that each orbital will appear in a different window.

More examples::

  sh ionplot.sh -o 3 H.ion

  for i in $(ioncat -i H.ion) ; do
    ionplot.sh -Z -o $i H.ion
  done

.. note::
   The ionplot script is not installed by default. Simply copy this
   file (name it *ionplot*) to your *$HOME/.local/bin* directory, or anywhere in your *PATH*:

   .. code:: shell

	     #!/bin/sh
	     #
	     string="$@"

	     ioncat $@  > .tmp_ioncat

	     cat > .plot.g << EOF
	     set title "$string"
	     plot ".tmp_ioncat" using 1:2 w l title "f"
	     replot ".tmp_ioncat" using 1:3 w l title "grad f"
	     EOF
	     #
	     gnuplot -persist .plot.g

	     
Using sisl to process .ion.nc files
-----------------------------------

When Siesta is compiled with netcdf support, information about basis
orbitals, KB projectors, etc, is also written to .ion.nc files. These
netcdf files can be read by `sisl <http://zerothi.github.io/sisl>`_.

In particular, plots that include orbital information as one of the
panels can be obtained, for example, by executing the ``sdata`` tool
of the sisl distribution::

   sdata Al.ion.nc --plot

.. note::
   You might need to install sisl and/or activate a specific
   python environment on your machine to access its functionality.
   In the School, type ``workon`` to see the environments, and then
   ``workon school`` or the appropriate environment.
   


