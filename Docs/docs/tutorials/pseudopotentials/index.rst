:sequential_nav: next

..  _tutorial-pseudopotentials:

Generation and testing of pseudopotentials
==========================================

..  sidebar:: **Have you set up the local environment?**

    If not, :ref:`do that now <local_installation>` before proceeding.

This tutorial covers a subject that used to be
very important: one needed to carefully control the quality of the
pseudopotentials used in a calculation. Nowadays, with the
availability of databases of curated pseudopotentials, such as the
`Pseudo Dojo <https://www.pseudo-dojo.org>`_, the need to look under the
hood is lessened (but has not disappeared completely!).

This module uses the ATOM program and its suite of tutorials.
The ATOM package itself can be installed simply by executing::

  git clone https://gitlab.com/garalb/atom-test
  cd atom-test
  ln -sf arch.make.QuantumMobile arch.make # this step might vary
  make

The documentation and tutorials can be accesed `here <https://docs.siesta-project.org/projects/atom>`_.
There is a lot to explore.


   
   
