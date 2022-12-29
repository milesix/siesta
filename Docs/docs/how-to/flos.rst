:sequential_nav: next

..  _how-to-flos:

Installing the FLOS library for Lua
===================================

.. note::
   For background on FLOS, please see the `flos documentation site
   <https://flos.readthedocs.io/>`_
   
To take advantage of the simulation functionality implemented with Lua
scripts (e.g., NEB, new structural relaxation algorithms, etc), the
FLOS library needs to be installed and the LUA_PATH appropriately
configured.

The procedure is really simple. First, choose a place for the FLOS
library, and change your working directory to that location. Then::
  
  git clone https://github.com/siesta-project/flos.git
  cd flos
  git submodule init
  git submodule update

That's it for the installation, but now we need to tell Lua where to
find the library. For this, you need to put a line in your .bashrc
file, or issue the command for every work session in which you intend
to use Lua functionality::

  export LUA_PATH="${FLOS_PLACE}/flos/?.lua;${FLOS_PLACE}/flos/?/init.lua;$LUA_PATH;;"

where you need to substitute `${FLOS_PLACE}` by the right directory
where you invoked the installation command. (For example, if you chose
your home directory, just replace `${FLOS_PLACE}` by `$HOME` in the
above export command).


  

  








