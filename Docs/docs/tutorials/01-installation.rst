.. _local_installation:

Setting up the local working environment for the tutorial exercises
===================================================================

Execution environment
---------------------

For the 2021 Siesta School we will be using the Siesta Mobile, a
virtual machine based on the `Quantum Mobile
<https://quantum-mobile.readthedocs.io>`_ sponsored by the Marvel and
MaX projects.

The virtual machine contains an Ubuntu operating system, plus
school-related executables, python packages, and data files. Users
just need to download the VM image and run it using the VirtualBox
framework, which is available for most operating systems.

For more details about the setup and a link to download the VM image,
please see `this link
<https://drive.google.com/drive/folders/14V50YRuJfW1jxdWkQzZPnTx0TIa10ftX>`_.

.. note::
   Using the Siesta Mobile is strongly recommended, since:

   * It provides a cutting-edge version of Siesta with a number of
     features that are not yet in the latest officially released version
     (4.1.5)::

       TD-DFT
       Spin-orbit coupling without the 'on-site approximation'
       Support for PSML pseudopotentials

   * Contains the most relevant visualizers and analysis tools that
     are needed for the tutorials

   * There is no need to deal with possibly complex installation
     issues, allowing full focus on learning how to use Siesta.

   However, if using the Siesta Mobile is not possible for technical
   reasons, there are (less optimal)
   :ref:`alternatives<how-to-alt-setup>`.
	

Working files for the tutorials
-------------------------------

They are currently embedded in the source for the document repository
you are reading, for consistency. One can obtain the whole
distribution and generate an independent tree of working files::

     git clone   https://gitlab.com/garalb/siesta-docs.git
     cd siesta-docs
     cd work-files
     sh link.sh

During the school, new working files might be added, or typos or
inconsistencies fixed. To get updated versions, one should do::

     cd siesta-docs     # To be inside the git-controlled domain
     git fetch
     git pull

Now you are ready to update the working-files hierarchy in directory
``work-files`` for the actual
work, but notice that the next step will potentially over-write files
that you might have changed. For example, if you were working in the
*tutorials/basic/first-encounter/CH4* directory and had edited ch4.fdf
to try some new feature, your changes will be lost. This only happens
for files with the same name as the `official` files, so it is good
practice to rename files in which you want to try new things that you
might want to keep.  In any case, you have also the option of copying
your working directory to some other place for safekeeping.

After this step, you can do::

     cd work-files
     sh link.sh

You should read the docs online  `here <https://docs.siesta-project.org/projects/siesta>`_, but could
also regenerate them locally if desired using the above repo. (Check
for availability of sphinx et al)

Things to add or download/compile later
---------------------------------------

Some items are not included in the Siesta Mobile for technical
reasons:

* Extra programs

  The programs to process the PDOS/PDOS.xml files can be obtained and
  compiled very easily. See :ref:`this how-to<how-to-dos-pdos>`.

* The FLOS library

  This is needed for the simulation code implemented in Lua
  (e.g. NEB, new relaxation algorithms, etc). It is easily installed
  following the instructions in :ref:`this how-to<how-to-flos>`.

* Other tutorial-specific material (such as extra packages and
  potentially large data files).

When needed, each tutorial will provide information about any
additional material needed.
