from spack import *


class Xmlf90(CMakePackage):
    """XML package for Fortran."""

    homepage = "https://gitlab.com/siesta-project/libraries/xmlf90"

    git = 'https://gitlab.com/siesta-project/libraries/xmlf90.git'

    version('master', branch='cmake')
    
    depends_on('cmake@3.14.0:', type='build')

    def cmake_args(self):
       args = []

       return args
