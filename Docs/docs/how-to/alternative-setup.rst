:sequential_nav: next

..  _how-to-alt-setup:

(Less optimal) Alternatives to the SiestaMobile
===============================================

You might not be able to use the Siesta Mobile for technical reasons,
or maybe you are a knowledgeable user that is fine without it. In that
case, to follow and profit from the tutorials, you are going to need:

* A compiled version of Siesta (See the `Guide to Siesta Versions
  <https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions>`_).

  The tutorials can be followed with the 'psml' version, or
  (without PSML support), with the 'master' version.  The released
  4.1.X versions are also missing TD-DFT and full-SOC support.

  The 4.1.5 version can be installed with
  :ref:`conda<building_with_conda>`.
       
  To compile any of the other versions, you can get inspiration from
  :ref:`this tutorial<tutorial-deployment>` and from other material in
  the Siesta documentation tree.

  Note that you also need most of the utility programs in the Util
  directory.
  
* Visualizers

  * gnuplot
  * XcrysDen
  * ...

* Extra analysis tools

  * `sisl <http://zerothi.github.io/sisl>`_. (It can be installed with
    conda)

  * Andrei Postnikov's 'visualization' and 'lattice-dynamics'
    packages, which can be obtained from::

      https://gitlab.com/siesta-project/analysis-tools/visualization
      https://gitlab.com/siesta-project/analysis-tools/lattice-dynamics


* Extra optional packages for some tutorials (as listed therein)

* The AiiDA and aiida-siesta packages if you plan to try the AiiDA
  tutorials. 




  

  








