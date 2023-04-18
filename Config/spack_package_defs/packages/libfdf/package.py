from spack import *


class Libfdf(CMakePackage):
    """Flexible Data Format library."""

    homepage = "https://gitlab.com/siesta-project/libraries/libfdf"

    git = 'https://gitlab.com/siesta-project/libraries/libfdf.git'

    version('master', branch='master')
    
    depends_on('cmake@3.14.0:', type='build')

    def cmake_args(self):
       args = []

       return args
