# Run this with:
#  python setup.py install --install-lib=.
from RDBuild import *
from distutils.core import setup,Extension

setup(name="SPtrTestModule",version="1.0",
  ext_modules=[Extension("SPtrTestModule",["module.cpp"],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_link_args=linkArgs,
                         extra_compile_args=compileArgs)])


