SIESTA is a program for efficient electronic structure calculations
and ab initio molecular dynamics simulations of molecules and
solids in the framework of Density-Functional Theory (DFT).
SIESTA's efficiency stems from the use of a basis set of
strictly-localized atomic orbitals. A very important feature of the
code is that its accuracy and cost can be tuned in a wide range, from
quick exploratory calculations to highly accurate simulations matching
the quality of other approaches, such as plane-wave methods.

The main web page for the project is at [siesta-project.org](https://siesta-project.org).

Online documentation, tutorials, and how-to's are being developed in [docs.siesta-project.org](https://docs.siesta-project.org).

Further information:

* Siesta development is multi-pronged, with stable, beta, and various
other branches with the latest features. A guide to the different
versions can be found
[in the project wiki](https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions).

* The LaTeX source for the manual is in Docs/siesta.tex. Assuming you have
a working Tex/LaTeX installation, you can type `make final` to generate a pdf file.

For bug reports, and other code suggestions, please follow the guidelines
in the file Docs/REPORTING_BUGS

### max1-berkeley branch

This is a special branch that contains, in addition to the features in `master`:

    * An interface to the [ELSI](http://www.elsi-interchange.org) solver interface library 
    * Extra optimizations to the parallel operation (including asynchronous grid comms)
    * Updated native PEXSI interface (version 2.X)
    
This is work-in-progress for development, so buckle up!



