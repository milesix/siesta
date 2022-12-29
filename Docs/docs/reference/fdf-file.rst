:sequential_nav: next

..  _reference-fdf-file:

The Siesta FDF input file
=========================

The main input for Siesta is a text file in FDF format, which contains
a description of the system to be simulated, a (maybe implicit) script
of actions to be taken, and a specification of the quality parameters
of the simulation.

The FDF format allows data to be given in any order, or to be omitted in favor of
default values. Here we offer a glimpse of it through the following rules:

*  The basic syntax is a ’data label’ followed by its value. Values that are
   not specified in the datafile are assigned a default value.

*  Labels are case insensitive, and the three characters *- _ .*
   (hyphen, underscore, dot)  in a data label
   are ignored. Thus, LatticeConstant, lattice-constant and LAttice_ConstaNT represent the
   same label.

*  All text following the *#* character is taken as comment.

*  Logical values can be specified as T, true, .true., yes, F, false,
   .false., no. Blank is also equivalent to *true*.

*  Character strings should **not** be in apostrophes.

*  Real values which represent a physical magnitude **must be followed by
   its units** (see :ref:`this section<fdf_units>`).
   
*  In some cases it is important to include a decimal point in a real
   number to distinguish it from an integer, in order to prevent
   ambiguities when mixing the types on the same input line.

*  Complex data structures are called blocks and are placed between
   *%block label* and a *%endblock label*.

*  You may ‘include’ other FDF files and redirect the search for a
   particular data label to another file. If a data label appears more
   than once, its first appearance is used.

*  If the same label is specified several times, **the first one takes
   precedence**, since the FDF parser stops looking for a label when
   it finds it.

*  If a label is misspelled it will not be recognized (there is no
   internal list of “accepted” tags in the program). You can check the
   actual value used by siesta by looking for the label in the output
   *fdf.log* file.

These are some examples::

              SystemName      Water molecule  # This is a comment
              SystemLabel     h2o
              Spin polarized
              SaveRho
              NumberOfAtoms         64
              LatticeConstant       5.42 Ang
              %block LatticeVectors
                       1.000  0.000  0.000
                       0.000  1.000  0.000
                       0.000  0.000  1.000
              %endblock LatticeVectors
              KgridCutoff < BZ_sampling.fdf

              # Reading the coordinates from a file
              %block AtomicCoordinatesAndAtomicSpecies < coordinates.data

              # Even reading more FDF information from somewhere else
              %include mydefaults.fdf

The file *fdf-XXXX.log* (with XXXX standing for some process-dependent
numbers) contains all the parameters used by the program in a given run,
both those specified in the input FDF file and those taken by default.
They are written in fdf format, so that you may reuse them as input
directly. Input data blocks are copied to the fdf.log file only if you
specify the *dump* option for them.

.. _fdf_units:

Physical unit names recognized by FDF
-------------------------------------

.. container:: center

   ========= ========= ============
   Magnitude Unit name MKS value
   ========= ========= ============
   mass      kg        1.E0
   mass      g         1.E-3
   mass      amu       1.66054E-27
   length    m         1.E0
   length    cm        1.E-2
   length    nm        1.E-9
   length    Ang       1.E-10
   length    Bohr      0.529177E-10
   time      s         1.E0
   time      fs        1.E-15
   time      ps        1.E-12
   time      ns        1.E-9
   time      mins      60.E0
   time      hours     3.6E3
   time      days      8.64E4
   energy    J         1.E0
   energy    erg       1.E-7
   energy    eV        1.60219E-19
   energy    meV       1.60219E-22
   energy    Ry        2.17991E-18
   energy    mRy       2.17991E-21
   energy    Hartree   4.35982E-18
   energy    Ha        4.35982E-18
   energy    K         1.38066E-23
   energy    kcal/mol  6.94780E-21
   energy    mHartree  4.35982E-21
   energy    mHa       4.35982E-21
   energy    kJ/mol    1.6606E-21
   energy    Hz        6.6262E-34
   energy    THz       6.6262E-22
   energy    cm-1      1.986E-23
   energy    cm**-1    1.986E-23
   energy    cm        1.986E-23
   force     N         1.E0
   force     eV/Ang    1.60219E-9
   force     Ry/Bohr   4.11943E-8
   ========= ========= ============

   ========= ========== ============
   Magnitude Unit name  MKS value
   ========= ========== ============
   pressure  Pa         1.E0
   pressure  MPa        1.E6
   pressure  GPa        1.E9
   pressure  atm        1.01325E5
   pressure  bar        1.E5
   pressure  Kbar       1.E8
   pressure  Mbar       1.E11
   pressure  Ry/Bohr**3 1.47108E13
   pressure  eV/Ang**3  1.60219E11
   charge    C          1.E0
   charge    e          1.602177E-19
   dipole    C*m        1.E0
   dipole    D          3.33564E-30
   dipole    debye      3.33564E-30
   dipole    e*Bohr     8.47835E-30
   dipole    e*Ang      1.602177E-29
   MomInert  Kg*m**2    1.E0
   MomInert  Ry*fs**2   2.17991E-48
   Efield    V/m        1.E0
   Efield    V/nm       1.E9
   Efield    V/Ang      1.E10
   Efield    V/Bohr     1.8897268E10
   Efield    Ry/Bohr/e  2.5711273E11
   Efield    Har/Bohr/e 5.1422546E11
   Efield    Ha/Bohr/e  5.1422546E11
   angle     deg        1.d0
   angle     rad        5.72957795E1
   torque    eV/deg     1.E0
   torque    eV/rad     1.745533E-2
   torque    Ry/deg     13.6058E0
   torque    Ry/rad     0.237466E0
   torque    meV/deg    1.E-3
   torque    meV/rad    1.745533E-5
   torque    mRy/deg    13.6058E-3
   torque    mRy/rad    0.237466E-3
   ========= ========== ============





