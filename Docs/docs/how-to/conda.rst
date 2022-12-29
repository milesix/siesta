:sequential_nav: next

..  _building_with_conda:

Building Siesta with conda
==========================

`Conda <https://coda.io>`_ is a package management system that runs on
multiple operating systems.

We have prepared recipes to install Siesta with it. Setup is
simple. You need at least to install the miniconda framework (provide
more details) and then do::

  conda install -c conda-forge siesta

to install a serial Siesta executable (current version is 4.1.5).

To install a parallel version, you need to use, depending on your
choice of parallel environment::

  conda install -c conda-forge "siesta=*=*openmpi*"
  conda install -c conda-forge "siesta=*=*mpich*"


  

  








