:sequential_nav: next

..  _how-to-dos-pdos:

DOS and projected DOS
=========================


DOS from EIG file using Eig2DOS
-------------------------------

Given the set of eigenvalues :math:`\epsilon_n({\bf k})`, the DOS can
be estimated as:

.. math::

   DOS(\omega) = \sum_{n{\bf k}} { \delta(\epsilon_n({\bf k}) - \omega)}
  
replacing the delta function with some suitable gaussian broadening.

The advantage of this method is that is very simple, and the EIG file
containing the BZ eigenvalues is always available after the end of the
scf cycle.

Since the algorithm is very simple minded (just sum over energy points
with some smearing, the DOS might be spiky (coarse k-sampling, small
broadening) or too broad. It is possible to increase the number of
k-points in the sampling of the BZ, but then a new scf cycle is
needed.

.. note::
   If not available, Eig2DOS can be built by doing::

     git clone https://gitlab.com/siesta-project/analysis-tools/eig2dos.git
     cd eig2dos; make
     cp -p Eig2DOS $HOME/.local/bin/

   (The last statement might need to be adapted, depending on the
   system)

Projected DOS from PDOS and PDOS.xml files
------------------------------------------

When a ``Projected-Density-of-States`` block is used in Siesta, such
as::

  %block Projected-density-of-states
  -26.00 4.00 0.200 500 eV
  %endblock Projected-density-of-states

Siesta will compute a full decomposition of the DOS over all orbitals,
in the energy range provided (above: -26.00 to 4.00 eV), using a given
broadening (0.2 eV above), and a given number of energy points in the
range (500 in the above example). The output is placed, for historical
reasons, in two files: `SystemLabel.PDOS` and
`SystemLabel.PDOS.xml`. Both contain the same (XML) data, but the
first one is formatted so that the data items appear in separate
lines. The second is more compact. We have kept the original .PDOS
file since it can be processed by the ``fmpdos`` program by Andrei
Postnikov, which does not really parse the XML, and depends on the
special format.

The program ``pdosxml`` uses an XML parser to process either file, and
is therefore more robust. However, it currently has a drawback: the
specification of which orbitals to take into account for the
projection of the DOS is done *in the code*, so the program must be
recompiled for each use. Not very convenient, but in practice quite
fast and providing full control of the ingredients of the projection.

.. note::
   If not available, *fmpdos* and *pdosxml* can be downloaded and
   built by doing::

     git clone https://gitlab.com/siesta-project/analysis-tools/pdos-xml.git
     cd pdos-xml/pdosxml ; make XMLF90_ROOT=/usr/local
     cd ../fmpdos; make
     cp -p fmpdos $HOME/.local/bin/

   The last statement might need to be adapted, depending on the
   system. Note that we do not install *pdosxml*, as it needs to be
   recompiled for each use.
   
The PDOS file can also be processed by `sisl
<http://zerothi.github.io/sisl>`_, using ???.


Note that if a different energy range, smearing, or number of points
is desired, a new Siesta run must be used. And since the PDOS
information is computed internally from the coefficients of the wave
functions, the program needs the Hamiltonian, which cannot currently
be read from a file. So even if a converged density-matrix file is
used to restart the calculation, some substantial computing is still
involved to get the new projected DOS data.

One can specify a denser BZ sampling for the PDOS calculation by
using a special block::

  %block PDOS.kgrid_Monkhorst_Pack
   8  0  0  0.5
   0  8  0  0.5
   0  0  8  0.5
  %endblock PDOS.kgrid_Monkhorst_Pack


Projected DOS from stored wavefunctions
---------------------------------------

An alternative way to process the projected DOS is to use directly the
wavefunctions associated to a full sampling of the BZ. These can be
produced by Siesta if the option::

  COOP.write T

is used. The name of the option is related to the framework for the
analysis of the :ref:`crystal-overlap populations<tutorial-coop-cohp>`
(COOP/COHP), which includes the processing of the pDOS as a
by-product.

This framework allows the interactive control of all the parameters
(energy range, broadening, orbitals involved, etc), at the expense of
having potentially large wavefunction files around.



