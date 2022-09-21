from spack import *


class Libpsml(CMakePackage):
    """A library to process PSeudopotential Markup Language files"""

    homepage = "https://gitlab.com/siesta-project/libraries/libpsml"

    git = 'https://gitlab.com/siesta-project/libraries/libpsml.git'

    version('master', branch='cmake')
    
    depends_on('cmake@3.14.0:', type='build')

    depends_on('xmlf90')
    
    def cmake_args(self):
       args = []

       return args
