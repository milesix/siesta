from spack import *


class Siesta(CMakePackage):
    """An efficient DFT simulation code (recipe from Siesta-Project)"""

    homepage = "https://siesta-project.org/siesta"

    git = 'https://gitlab.com/garalb/siesta.git'

    # Only a development version for now
    version('master-spack', branch='master-spack')

    variant('mpi', default=False, description='Use MPI')
    variant('netcdf', default=False, description='Use NetCDF')
    variant('libxc', default=False, description='Use libxc')
    variant('elpa', default=False, description='Use ELPA library (native interface)')
    variant('fftw', default=True, description='Use FFTW library (needed only for STM/ol-stm)')
    

    depends_on('cmake@3.14.0:', type='build')

    # generator = 'Ninja'
    # depends_on('ninja', type='build')
    
    depends_on('lapack')
    depends_on('xmlf90')
    depends_on('libpsml')
    depends_on('mpi', when='+mpi')
    depends_on('scalapack', when='+mpi')
    depends_on('netcdf-fortran', when='+netcdf')
    depends_on('libxc@4:', when='+libxc')
    depends_on('libgridxc+libxc', when='+libxc')
    depends_on('libgridxc~libxc', when='-libxc')
    depends_on('libgridxc+mpi', when='+mpi')
    depends_on('libgridxc~mpi', when='-mpi')
    depends_on('elpa', when='+elpa')
    depends_on('fftw@3.3.0:', when='+fftw')
    

    def cmake_args(self):
       args = [
            self.define_from_variant('WITH_MPI', 'mpi'),
            self.define_from_variant('WITH_LIBXC', 'libxc'),
            self.define_from_variant('WITH_NETCDF', 'netcdf'),
            self.define_from_variant('WITH_ELPA', 'elpa'),
            self.define_from_variant('WITH_FFTW', 'fftw'),
       ]

       return args

