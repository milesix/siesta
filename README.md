SIESTA is a program for efficient electronic structure calculations
and ab initio molecular dynamics simulations of molecules and
solids in the framework of Density-Functional Theory (DFT).
SIESTA's efficiency stems from the use of a basis set of
strictly-localized atomic orbitals. A very important feature of the
code is that its accuracy and cost can be tuned in a wide range, from
quick exploratory calculations to highly accurate simulations matching
the quality of other approaches, such as plane-wave methods.

The main web page for the project is at [icmab.es/siesta](https://icmab.es/siesta).

Further information:

* Siesta development is multi-pronged, with stable, beta, and various
other branches with the latest features. A guide to the different
versions can be found
[in the project wiki](https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions).

* TranSiesta is now part of the executable; see the documentation for details.

* The LaTeX source for the manual is in Docs/siesta.tex. Assuming you have
a working Tex/LaTeX installation, you can type `make final` to generate a pdf file.
Alternatively, manuals in pdf can be found in the Documentation section of the main web page.

* Tutorial material can be found also in the Documentation section of the main web page.

For bug reports, and other code suggestions, please follow the guidelines
in the file Docs/REPORTING_BUGS

### MaX branch
This is a special branch that contains all the features currently
implemented within the [MaX Center of
Excellence](http://www.max-centre.eu) EU H2020 project. It is
addressed to HPC users that want to try the new optimizations and
enhancements, and give feedback on them. It currently contains, in
addition to the features in `master`:

    * An interface to the [ELSI](http://www.elsi-interchange.org) solver interface library 
    * Support for PSML pseudopotentials
    * Extra optimizations to the parallel operation
    * A more modular architecture, using libraries and other modules developed within MaX.

Note that the libPSML and libGridXC libraries are now
mandatory. Please see the 'Compilation' section of the manual for
installation instructions and extra building advice.

The ELSI library (and ELPA) are recommended. Please see 000_INSTALL for installation
instructions, and README_ELSI for important information.




