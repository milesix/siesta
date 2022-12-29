.. sectionauthor:: Yann Pouillon <y.pouillon@simuneatomistics.com>

Installing and running executables
==================================

In this section, we will assume that you have successfully built SIESTA, as
well as the utilities you need. If this not the case, yet, you can do it by
following the instructions in one of the related sections:

- :doc:`build-manually` (all SIESTA versions)
- :doc:`build-esl-bundle` (SIESTA 5.x and Max-x versions)
- :doc:`build-easybuild` (SIESTA 5.x and Max-x versions)

As a starting point, we will suppose that you have just built the *siesta*,
*gnubands* and *denchar* executables from SIESTA 4.1.5 and that they are still
in your build directory, located at *$HOME/src/siesta-4.1.5/Obj/*.


Preparing an installation directory
-----------------------------------

One thing is to build software, another is to use it reliably. For the latter,
we have to decide where we want to store what we have built. Here we will
provide an example of good practice that we invite you to tune to best suit
your own needs.

In this example, we have chosen to install all files related to SIESTA in a
``$HOME/siesta`` directory. Inside, we first have to create a basic
structure, so that we can build and install several versions of SIESTA at the
same time. Starting with our example, we would type::

    mkdir -p $HOME/siesta/4.1.5/bin

As you can see, we are organising files by version (*4.1.5*) and by type
(*bin*). Feel free to adapt this structure to your preferences. The most
important aspect is to make a choice and stay consistent with it over time.

.. note::

   Documenting your choices is a very good practice. A simple way to do it is
   to create a README file in the *siesta/* directory and summarise there how
   you install files in its subdirectories.

For the moment, even though we are only installing executables, we will store
them in a *bin/* subdirectory, so that we don't have to change our directory
structure when we install other kinds of files, e.g. libraries or data files.


Copying the files to the installation location
----------------------------------------------

In our example, we have built executables from the
*$HOME/src/siesta-4.1.5/Obj/* build directory. To make them available in a
permanent way and free the space for new builds, we just have to copy the
executables to our installation directory. In this case, we type::

    cp $HOME/src/siesta-4.1.5/Obj/siesta $HOME/siesta/4.1.5/bin
    cp $HOME/src/siesta-4.1.5/Obj/Util/Bands/gnubands $HOME/siesta/4.1.5/bin
    cp $HOME/src/siesta-4.1.5/Obj/Util/Denchar/Src/denchar $HOME/siesta/4.1.5/bin

From now on, we can safely run ``make clean`` in the build directory without
loosing the fruit of our efforts.

Saving the *arch.make* that was used to build teh executables is usually a
good idea if we want to keep track of how we have built them. In this case, we
can create a *conf/* subdirectory where we will store them. In the process, it might also be useful to rename the arch.make to reflect our configuration
choices explicitly. If we have built a SIESTA executable with GCC, OpenMPI and
the Netlib linear algebra libraries, we can rename the file accordingly::

    mkdir -p $HOME/siesta/4.1.5/conf
    cp $HOME/src/siesta-4.1.5/Obj/arch.make \
        $HOME/siesta/4.1.5/conf/gcc-openmpi-netlib.make

Now, if we want to build different flavours of the same SIESTA version using
different *arch.make* files, we could also rename the SIESTA executable
following the same conventions::

    mv $HOME/siesta/4.1.5/bin/siesta \
        $HOME/siesta/4.1.5/bin/siesta-gcc-openmpi-netlib

In this case, it would be advisable to update the *$HOME/siesta/README* file
to keep track of these conventions.


Accessing the executables
-------------------------

One important step before being able to run the executables reliably is to
update the environment variables. Two of them are of special importance:

- *LD_LIBRARY_PATH* for the libraries.
- *PATH* for the executables.

A good practice is to group these settings in a shell script that will be
*sourced* whenever you want to run a specific version of SIESTA.

The first thing to decide is where to store the script. The most common used
options are:

1. Store the script in the version-specific installation subdirectory, in
   which case a single name can be used for all the instances of the script.
   In our example, this would be *$HOME/siesta/4.1.5/siesta-vars.sh*. If, in
   the future, we install SIESTA 5.0.0, the corresponding script would be
   *$HOME/siesta/5.0.0/siesta-vars.sh*.
2. Store all the scripts in the same directory, e.g. *$HOME/siesta/env/*, in
   which case they would have to be named differently. In the above example,
   this would be *$HOME/siesta/env/siesta-4.1.5-vars.sh* and
   *$HOME/siesta/env/siesta-5.0.0-vars.sh*.

No solution is better than the other. What is important is that you feel at
ease with it.

In the following, we will illustrate option 1. Adapting the instructions to
option 2 is relatively easy. For SIESTA 4.1.5, the
*$HOME/siesta/4.1.5/siesta-vars.sh* script will look like the following:

.. code-block:: shell

   #!/bin/sh

   PATH="$HOME/siesta/4.1.5/bin:$PATH"
   export PATH

This is a minimalistic content. If your are familiar with shell scripting, you
can of course refine on this.

From SIESTA 5.x on, as well as for the MaX-* releases, you have to install
dependencies before compiling SIESTA. In this case, the recommended
installation directory for e.g. SIESTA 5.0.0 would be
*$HOME/siesta/5.0.0/lib/* and the environment script would look like the
following:

.. code-block:: shell

   #!/bin/sh

   LD_LIBRARY_PATH="$HOME/siesta/5.0.0/lib:$LD_LIBRARY_PATH"
   export LD_LIBRARY_PATH

   PATH="$HOME/siesta/5.0.0/bin:$PATH"
   export PATH

Once the environment is correctly prepared, the SIESTA executables can be run
in all kinds of circumstances in a reliable way.


Running the executables
-----------------------

Whenever you want to run a SIESTA executable you have installed as recommended
here, you will have to **source the corresponding script once per terminal
session**. For instance, if you want to run SIESTA 4.1.5, you will open a new
terminal and type the following::

    source $HOME/siesta/4.1.5/siesta-vars.sh

From that moment on, you can run SIESTA and its utilities by typing their
names. If you want to understand better what is happening to the environment,
you can start a terminal session and type ``which siesta``, then source the
script and type ``which siesta`` again. Your terminal will look like this::

    $ which siesta
    $ source $HOME/siesta/4.1.5/siesta-vars.sh
    $ which siesta
    $HOME/siesta/4.1.5/bin/siesta

The ``which`` command tells you which executable is actually run when you type
the corresponding command. The first time, the SIESTA command is not available
in the environment, hence ``which`` returns nothing. However, after sourcing
the script, the environment has been updated and the command is available.

.. note::

   It is always a good idea to check commands with ``which`` when you install
   software, no matter the method you use to install it.

These instructions are also valid when you run SIESTA through a batch
scheduler on a HPC cluster. In this case, your batch job would look like the
following:

.. code-block:: shell

   #!/bin/bash
   #PBS -q dft
   #PBS -l nodes=1:ppn=12
   #PBS -V
   #PBS -N "Test"

   # Set environment
   source $HOME/siesta/4.1.5/siesta-vars.sh
   which siesta

   # Run SIESTA
   cd my_job_dir
   mpirun -np 12 siesta <my_input.fdf >my_output.log

One last important piece of advice: it is essential to avoid mixing the
environments for different versions of SIESTA within the same shell session,
as the results are highly impredictable and will likely end up in random
crashes or garbled SIESTA output. **If you want to run different versions of
SIESTA, always prefer running them in different terminal sessions.**

You're now done. Enjoy your shiny new SIESTA installation!
