from spack import *


class Libfdf(CMakePackage):
    """Flexible Data Format library."""

    homepage = "https://gitlab.com/siesta-project/libraries/libfdf"

    git = 'https://gitlab.com/siesta-project/libraries/libfdf.git'

    #version("0.5.0", tag="0.5.0")
    version("0.5.0", commit="d7e7178eadf5")   # same as above, but trusted
    
    depends_on('cmake@3.17.0:', type='build')

    def cmake_args(self):
       args = []

       return args
