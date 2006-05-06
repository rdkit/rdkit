# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension

#libraries = ["RDBoost",'RDGeneral'] + libraries


setup(name="rdBase",version="1.0",
  ext_modules=[Extension("rdBase",
                         ["RDBase.cpp",
                          ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_link_args=linkArgs,
                         extra_compile_args=compileArgs)])

