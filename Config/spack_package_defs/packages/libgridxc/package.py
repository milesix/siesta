from spack import *


class Libgridxc(CMakePackage):
    """A library for exchange-correlation calculations on radial and
       solid-state grids"""

    homepage = "https://gitlab.com/siesta-project/libraries/libgridxc"

    git = 'https://gitlab.com/siesta-project/libraries/libgridxc.git'

    version('master', branch='master')

    variant('mpi', default=False, description='Use MPI.')
    variant('libxc', default=False, description='Use libxc')
    

    depends_on('cmake@3.14.0:', type='build')

    depends_on('mpi', when='+mpi')
    depends_on('libxc@4:', when='+libxc')
    

    def cmake_args(self):
       args = [
            self.define_from_variant('WITH_MPI', 'mpi'),
            self.define_from_variant('WITH_LIBXC', 'libxc'),
       ]

       return args
